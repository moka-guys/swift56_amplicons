import math
import MySQLdb
import sqlalchemy
import requests
import xml.etree.ElementTree as ET

class Genes(object):
	'''
	Takes a lits of genes. 
	Methods to:
	Pull out all RefSeq transcripts and details from UCSC. 
	Map cDNA exon coordinates to genomic coordinates.
	Get canonical transcript for each gene from UCSC
	Get the most reported transcript in ClinVar
	Results stored in nested dictionary of format:
	gene_dict{gene: {canonical: NM_xxxxx, clinvar_most_reported: NM_xxxxx, transcripts: {NM_xxxxx: {cds_start: xxxx, ....etc.}}} ...}
	'''
	def __init__(self, gene_list):
		self.engine = None
		self.connection = None
		self.gene_transcript_dict = {}
		for gene in gene_list:
			self.gene_transcript_dict[gene] = {'transcripts': {}}

	def connect_db(self):
		# Opens connection to UCSC
		self.engine = sqlalchemy.create_engine('mysql+mysqldb://genome:@genome-mysql.soe.ucsc.edu/hg19')
		self.connection = self.engine.connect()

	def disconnect_db(self):
		# CLoses connection to UCSC
		self.connection.close()	

	def get_UCSC_refGene(self):
		#Get all refSeq transcript rows from UCSC for each gene
		self.connect_db()
		all_transcripts = self.connection.execute("SELECT * FROM refGene WHERE refGene.name2 IN ('" + "','".join(self.gene_transcript_dict.keys()) + "');")
		self.disconnect_db()
		# Convert results to dictionary
		for tx in all_transcripts:
			t_dict = dict(tx)
			t_dict['exonStarts'] = map(long, t_dict['exonStarts'].rstrip(',').split(',')) #Convert comma separated string to list of numbers
			t_dict['exonEnds'] = map(long, t_dict['exonEnds'].rstrip(',').split(',')) #Convert comma separated string to list of numbers
			t_dict['exonFrames'] = map(int, t_dict['exonFrames'].rstrip(',').split(',')) #Convert comma separated string to list of numbers
			gene = t_dict['name2']
			self.gene_transcript_dict[gene]['transcripts'][t_dict['name']] = t_dict			

	def get_cds(self):
		# Trims the 5' and 3' UTRs from the exons to leave the cDNA exon coordinates.
		# Converts exon starts to 1 based
		for gene in self.gene_transcript_dict:
			for transcript, tx_details in self.gene_transcript_dict[gene]['transcripts'].items():				
				all_exons = zip(tx_details['exonStarts'], tx_details['exonEnds']) # These exons include UTR
				cds_start = tx_details['cdsStart']
				cds_end = tx_details['cdsEnd']	
				cds_exons = []
				for exon_start, exon_end in all_exons:
					#Exons lying completely outside cds will not meet any of the following conditions so will be skipped
					#if coding sequence starts within the exon, change start site to cds start site
					if exon_start <= cds_start and exon_end >= cds_start:
						cds_exons.append((cds_start+1, exon_end)) #Increase start by 1 so that it's 1 based instead of 0 based
					#if all of exon within coding sequence, keep exon as it is
					elif exon_start > cds_start and exon_end < cds_end:
						cds_exons.append((exon_start+1, exon_end)) #Increase start by 1 so that it's 1 based instead of 0 based
					#if coding sequence ends within the exon, change end site to cds end site
					elif exon_start <= cds_end and exon_end >= cds_end:
						cds_exons.append((exon_start+1, cds_end)) #Increase start by 1 so that it's 1 based instead of 0 based
				tx_details['cds_exons'] = cds_exons
		self.map_g_cds()

	def map_g_cds(self):
		##Map the cDNA coordinates to the genomic coordinates [((g_start, g_end), (c_start, c_end)), ((g_start, g_end), (c_start, c_end))...]
		for gene in self.gene_transcript_dict:
			for transcript, tx_details in self.gene_transcript_dict[gene]['transcripts'].items():	
				g2cDNA = []
				cDNA_start = 0
				cDNA_end = 0
				if tx_details['strand'] == '+':
					for g_start, g_end in tx_details['cds_exons']:
						cDNA_start = cDNA_end + 1 #cDNA start of one exon is 1 higher than cDNA end of the previous exon
						cDNA_end = cDNA_start + (g_end - g_start) #cDNA end of one exon is the cDNA start + difference between genomic start and end
						g2cDNA.append(((g_start, g_end), (cDNA_start, cDNA_end)))
				# If negative strand, loop through exons in reverse order...
				elif tx_details['strand'] == '-':
					for g_start, g_end in reversed(tx_details['cds_exons']):
						cDNA_start = cDNA_end + 1 #cDNA start of one exon is 1 higher than cDNA end of the previous exon
						cDNA_end = cDNA_start + (g_end - g_start) #cDNA end of one exon is the cDNA start + difference between genomic start and end
						g2cDNA.append(((g_start, g_end), (cDNA_end, cDNA_start)))
					g2cDNA.reverse()
				tx_details['g2cDNA'] = g2cDNA

	def get_canonical(self):
		# Retrieve the canonical transcript for each gene from UCSC
		self.connect_db()
		canonical_transcripts = self.connection.execute("SELECT kgXref.geneSymbol, kgXref.refseq FROM knownCanonical JOIN kgXref ON knownCanonical.transcript = kgXref.kgID WHERE kgXref.geneSymbol IN ('" + "','".join(self.gene_transcript_dict.keys()) + "');")		
		self.disconnect_db()
		for tx in canonical_transcripts:
			self.gene_transcript_dict[tx['geneSymbol']]['canonical'] = tx['refseq']

	def get_clinVar(self):
		# Get the most reported transcript from ClinVar for each gene
		# This can take a few minutes if gene list is big
		for gene in self.gene_transcript_dict:
			payload = {'db': 'clinvar',
 			'term': gene + '[gene]',
			'retmax': 500,
			'retmode': 'json'}
			# First request pulls out JSON containing the IDs for each ClinVar variant record for a gene (max 500 entries)
			r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params = payload)
			j = r.json()
			ids = j['esearchresult']['idlist']
			# Second request pulls out an XML containing the variant reports in the list of IDs retrieved in first request
			# (This endpoint doesn't support JSON)	
			payload = {'db': 'clinvar',
					   'rettype': 'variation',
					   'id': ','.join(ids)}
			r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params = payload)
			root = ET.fromstring(r.content)
			transcripts = {}
			# loop through XML and pull out transcript from variant name.
			# Keep a count of how many times a transcript is used in dictionary
			for variationReport in root.iter('VariationReport'):
				if variationReport.attrib['VariationName'].startswith('NM_'):
					transcript = variationReport.attrib['VariationName'].split('.')[0]
					transcripts[transcript] = transcripts.setdefault(transcript, 0) + 1 #Increase count by 1
			# Find transcript with the highest count
			counts_transcripts = [(count, tx) for tx, count in transcripts.items()]
			most_reported_tx = max(counts_transcripts)[1]
			self.gene_transcript_dict[gene]['clinvar_most_reported'] = most_reported_tx



class Amplicon(object):

	def __init__(self, gene, chr, start, end):
		self.gene = gene
		self.chr = chr
		self.start = start
		self.end = end
		self.cds_aa_coords = {}

	def get_cDNA_aa(self, gene_transcript_dict):
		#Takes the gene_transcript_dict from Genes object and uses it to calculate the cDNA and aa start and end positions for the amplicon
		for transcript, tx_details in gene_transcript_dict[self.gene]['transcripts'].items():
			#Call appropriate method (dependent on +ve or -ve strand) to generate dictionary containing cDNA and aa starts, ends etc	
			if tx_details['strand'] == "+":
				t_dict = self.get_cDNA_aa_positive(tx_details)
			elif tx_details['strand'] == "-":
				t_dict = self.get_cDNA_aa_negative(tx_details)
			self.cds_aa_coords[transcript] = t_dict
		return self.cds_aa_coords

	def get_cDNA_aa_positive(self, tx_details):
		##Convert the genomic start and end to cDNA start and end
		amp_cDNA_start = None
		amp_cDNA_end = None
		amp_cDNA_start_diff = None
		amp_cDNA_end_diff = None
		aa_start = None
		aa_start_diff = None
		aa_end = None
		aa_end_diff = None		
		amp_g_start = self.start
		amp_g_end = self.end
		g2cDNA = tx_details['g2cDNA']
		#If the amplicon lies completely outside the transcript, set cDNA start and end to 0
		if ((amp_g_start < g2cDNA[0][0][0] and amp_g_end < g2cDNA[0][0][0])
			or (amp_g_start > g2cDNA[-1][-2][1] and amp_g_end > g2cDNA[-1][-2][1])):
			amp_cDNA_start = 0
			amp_cDNA_end = 0
			amp_cDNA_start_diff = 0
			amp_cDNA_end_diff = 0
			aa_start = 0
			aa_start_diff = 0
			aa_end = 0
			aa_end_diff = 0
		#If amp_cDNA_start is still None after above step, loop through exons and find exon (or intron) containing start position
		exon = 0
		while amp_cDNA_start == None:
			g_exon, c_exon = g2cDNA[exon]
			#If the genomic start lies upstream of the exon (i.e. in an intron), set cDNA start position to cDNA start of the exon
			if amp_g_start < g_exon[0]:
				amp_cDNA_start = c_exon[0] #First cDNA base included in the transcript
				amp_cDNA_start_diff = g_exon[0] - amp_g_start #Calculate number of bases upstream of the start cDNA base
			#Else if the genomic start lies within the intron, calculate the cDNA position and set it as amp_cDNA_start
			elif amp_g_start >= g_exon[0] and amp_g_start <= g_exon[1]:
				amp_cDNA_start = c_exon[0] + (amp_g_start - g_exon[0])
				amp_cDNA_start_diff = 0 #cDNA is not upstream exon, so diff is 0
			exon += 1 #Move to next exon
		#To get amino acid start pos add 2 to cDNA start, divide by 3 and round up
		#Need to add 2 because the start of the codon won't be divisible by 3, only the end 
		aa_start = int(math.ceil((amp_cDNA_start+2)/3.0))
		if (amp_cDNA_start+2) % 3:
			aa_start_diff =  3 - ((amp_cDNA_start+2) % 3) ##If the amplicon starts partway through a codon, this number states number of bases from start of amplicon to first codon
		else:	
			aa_start_diff = 0

		#Next find the cDNA end position. Needs work backwards from the last exon for the transcript.
		exon = len(g2cDNA) - 1 #index of the last exon in list, subtract 1 to convert from 1 based to 0 based
		while amp_cDNA_end == None:
			g_exon, c_exon = g2cDNA[exon]
			#If the end is downstream of the exon...	
			if amp_g_end > g_exon[1]: 
				amp_cDNA_end = c_exon[1] #Set the cDNA base position to the last base in the exon
				amp_cDNA_end_diff = amp_g_end - g_exon[1] #Calculate number of bases downstream of the last cDNA base
			#Else if the genomic end lies within the exon, calculate cDNA position and set it as amp_cDNA_end
			elif amp_g_end >= g_exon[0] and amp_g_end <= g_exon[1]: 
				amp_cDNA_end = c_exon[0] + (amp_g_end - g_exon[0])
				amp_cDNA_end_diff = 0 #cDNA end is not downstream of the exon, so diff is 0
			exon -= 1 #Move to next exon (reverse order)
		#Get the last amino acid fully included in the amplicon
		aa_end = int(math.floor((amp_cDNA_end)/3.0))
		aa_end_diff = amp_cDNA_end % 3 #If the amplicon ends partway through a codon, this number states number of bases from last codon to end of amplicon
		t_dict = {}
		t_dict['cDNA_start'] = amp_cDNA_start
		t_dict['amp_cDNA_start_diff'] = amp_cDNA_start_diff
		t_dict['cDNA_end'] = amp_cDNA_end
		t_dict['amp_cDNA_end_diff'] = amp_cDNA_end_diff
		t_dict['aa_start'] = aa_start
		t_dict['amp_aa_start_diff'] = aa_start_diff
		t_dict['aa_end'] = aa_end
		t_dict['amp_aa_end_diff'] = aa_end_diff
		#t_dict['transcript'] = transcript
		return t_dict

	def get_cDNA_aa_negative(self, tx_details):
		##Convert the genomic start and end to cDNA start and end
		amp_cDNA_start = None
		amp_cDNA_end = None
		amp_cDNA_start_diff = None
		amp_cDNA_end_diff = None
		aa_start = None
		aa_start_diff = None
		aa_end = None
		aa_end_diff = None		
		amp_g_start = self.start
		amp_g_end = self.end
		g2cDNA = tx_details['g2cDNA']
		#If the amplicon lies completely outside the transcript, set cDNA start and end to 0
		if ((amp_g_start < g2cDNA[0][0][0] and amp_g_end < g2cDNA[0][0][0])
			or (amp_g_start > g2cDNA[-1][-2][1] and amp_g_end > g2cDNA[-1][-2][1])):
			amp_cDNA_start = 0
			amp_cDNA_end = 0
			amp_cDNA_start_diff = 0
			amp_cDNA_end_diff = 0
			aa_start = 0
			aa_start_diff = 0
			aa_end = 0
			aa_end_diff = 0

		#If amp_cDNA_start is still None after above step, find the cDNA start position. 
		#Needs work backwards from the last exon for the transcript because it is on reverse strand.
		exon = len(g2cDNA) - 1 #index of the last exon in list, subtract 1 to convert from 1 based to 0 based
		while amp_cDNA_start == None:
			g_exon, c_exon = g2cDNA[exon]
			#If the start is downstream of the exon in an intron...
			if amp_g_end > g_exon[1]: 
				amp_cDNA_start = c_exon[1] #Set the cDNA base position to the last base in the exon
				amp_cDNA_start_diff = amp_g_end - g_exon[1] #Calculate number of bases downstream of the last cDNA base
			#Else if the genomic end lies within the exon, calculate cDNA position and set it as amp_cDNA_start
			elif amp_g_end >= g_exon[0] and amp_g_end <= g_exon[1]: 
				amp_cDNA_start = c_exon[1] + (g_exon[1] - amp_g_end)
				amp_cDNA_start_diff = 0 #cDNA end is not downstream of the exon, so diff is 0
			exon -= 1 #Move to next exon (reverse order)
		#To get amino acid start pos add 2 to cDNA start, divide by 3 and round up
		#Need to add 2 because the start of the codon won't be divisible by 3, only the end 
		aa_start = int(math.ceil((amp_cDNA_start+2)/3.0))
		if (amp_cDNA_start+2) % 3:
			aa_start_diff =  3 - ((amp_cDNA_start+2) % 3) ##If the amplicon starts partway through a codon, this number states number of bases from start of amplicon to first codon
		else:	
			aa_start_diff = 0


		#If amp_cDNA_end is still None after the first step, loop through exons and find exon (or intron) containing end position
		exon = 0
		while amp_cDNA_end == None:
			g_exon, c_exon = g2cDNA[exon]
			#If the genomic start lies upstream of the exon in an intron, set cDNA end position to cDNA end of the exon
			if amp_g_start < g_exon[0]:
				amp_cDNA_end = c_exon[0] #First cDNA base included in the transcript
				amp_cDNA_end_diff = g_exon[0] - amp_g_start #Calculate number of bases upstream of the start cDNA base
			#Else if the genomic start lies within the exon, calculate the cDNA position and set it as amp_cDNA_end
			elif amp_g_start >= g_exon[0] and amp_g_start <= g_exon[1]:
				amp_cDNA_end = c_exon[1] + (g_exon[1] - amp_g_start)
				amp_cDNA_end_diff = 0 #cDNA is not upstream exon, so diff is 0
			exon += 1 #Move to next exon
		#Get the last amino acid fully included in the amplicon
		aa_end = int(math.floor((amp_cDNA_end)/3.0))
		aa_end_diff = amp_cDNA_end % 3 #If the amplicon ends partway through a codon, this number states number of bases from last codon to end of amplicon

		t_dict = {}
		t_dict['cDNA_start'] = amp_cDNA_start
		t_dict['amp_cDNA_start_diff'] = amp_cDNA_start_diff
		t_dict['cDNA_end'] = amp_cDNA_end
		t_dict['amp_cDNA_end_diff'] = amp_cDNA_end_diff
		t_dict['aa_start'] = aa_start
		t_dict['amp_aa_start_diff'] = aa_start_diff
		t_dict['aa_end'] = aa_end
		t_dict['amp_aa_end_diff'] = aa_end_diff
		return t_dict