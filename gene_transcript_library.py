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
			# The name2 field if the gene name
			gene = t_dict['name2']
			# The name field is the transcript accession (NM...)
			# Store the dictionary created above within the nested dictionary 
			self.gene_transcript_dict[gene]['transcripts'][t_dict['name']] = t_dict			

	def get_cds(self):
		# Trims the 5' and 3' UTRs from the exons to leave the cDNA exon coordinates.
		# Converts exon starts to 1 based
		for gene in self.gene_transcript_dict:
			for transcript, tx_details in self.gene_transcript_dict[gene]['transcripts'].items():	
				# transcript is the NM accession, tx_details is the dictionary containing data retrieved from UCSC refGene
				# Zip the exon starts and ends together to create list of tuples, where each tuple is the start and end of an exon
				# i.e. [(exon1_start, exon1_end), (exon2_start, exon2_end)...]
				# At this stage the UTRs are still included 			
				all_exons = zip(tx_details['exonStarts'], tx_details['exonEnds']) # These exons include UTR
				# Capture cdsStart and cdsEnd (i.e. start and end of transcript without UTR)
				cds_start = tx_details['cdsStart']
				cds_end = tx_details['cdsEnd']
				# Empty list which will be used to create list in same format as above but without the UTRs	
				cds_exons = []
				# Loop through all_exons list and use cdsStart and cdsEnd to trim the UTRs
				for exon_start, exon_end in all_exons:
					#Exons lying completely outside cds (i.e. whole exon is in UTR) will not meet any of the following conditions so will be skipped
					# if the entire coding sequence is contained within this exon, change the start and end of the exon to cds start and end sites
					if exon_start < cds_start and exon_end > cds_end:
						cds_exons.append((cds_start + 1, cds_end)) # Increase start by 1 so that it's 1 based instead of 0 based
					#else if coding sequence (cdsStart) starts within the exon, change the start site of the exon to cds start site
					elif exon_start <= cds_start and exon_end >= cds_start:
						cds_exons.append((cds_start+1, exon_end)) #Increase start by 1 so that it's 1 based instead of 0 based
					#else if all of exon within coding sequence, keep exon as it is
					elif exon_start >= cds_start and exon_end <= cds_end:
						cds_exons.append((exon_start+1, exon_end)) #Increase start by 1 so that it's 1 based instead of 0 based
					#else if coding sequence ends within the exon, change end site to cds end site
					elif exon_start <= cds_end and exon_end >= cds_end:
						cds_exons.append((exon_start+1, cds_end)) #Increase start by 1 so that it's 1 based instead of 0 based
				# Store the cds_exons list in the dictionary
				tx_details['cds_exons'] = cds_exons
		# Call the map_g_cds() method 
		self.map_g_cds()

	def map_g_cds(self):
		##Map the cDNA coordinates of the exons onto the genomic coordinates and produce a list in following format:
		# [((g_start, g_end), (c_start, c_end)), ((g_start, g_end), (c_start, c_end))...]
		# The list will be in order of increasing genomic coordinates with respect to forward strand
		# Therefore for transcripts on positive strand, the genomic and cDNA coordinates will increase from start to end of list
		# Whereas for transcripts on negative strand, the genomic coordinates will increase but the cDNA coordinates will decrease from start to end of list
		for gene in self.gene_transcript_dict:
			for transcript, tx_details in self.gene_transcript_dict[gene]['transcripts'].items():
				# transcript is the NM number
				# tx_details is the dictionary containing data for that transcript	
				g2cDNA = [] # This will hold the list of mapped genomic and cDNA coordinates described above
				cDNA_start = 0
				cDNA_end = 0
				# If the transcript is on positive strand...
				if tx_details['strand'] == '+':
					for g_start, g_end in tx_details['cds_exons']:
						cDNA_start = cDNA_end + 1 #cDNA start of one exon is 1 higher than cDNA end of the previous exon
						cDNA_end = cDNA_start + (g_end - g_start) #cDNA end of one exon is the cDNA start + difference between genomic start and end of the exon
						g2cDNA.append(((g_start, g_end), (cDNA_start, cDNA_end))) # Add to the list
				# If transcript on negative strand, the first exon will be at the end of the cds_exons list
				# Therefore need to loop through list in reverse order...
				elif tx_details['strand'] == '-':
					for g_start, g_end in reversed(tx_details['cds_exons']):
						cDNA_start = cDNA_end + 1 #cDNA start of one exon is 1 higher than cDNA end of the previous exon
						cDNA_end = cDNA_start + (g_end - g_start) #cDNA end of one exon is the cDNA start + difference between genomic start and end
						# Unlike positive strand transcripts, need to switch the cDNA_end and cDNA_start around, so that g_start maps to cDNA end and vice versa 
						g2cDNA.append(((g_start, g_end), (cDNA_end, cDNA_start)))
					# Because we worked backwards, the g2cDNA list will be in reverse genomic coordinate order. 
					# Therefore need to reverse it so that it is in	order of increasing genomic coordinates
					g2cDNA.reverse()
				# Store the 'g2cDNA' list in the dictionary for that transcript
				tx_details['g2cDNA'] = g2cDNA

	def get_canonical(self):
		# Retrieve the canonical transcript for each gene from UCSC knownCanonical table (https://www.ensembl.org/Help/Glossary?id=346)
		self.connect_db() 
		canonical_transcripts = self.connection.execute("SELECT kgXref.geneSymbol, kgXref.refseq FROM knownCanonical JOIN kgXref ON knownCanonical.transcript = kgXref.kgID WHERE kgXref.geneSymbol IN ('" + "','".join(self.gene_transcript_dict.keys()) + "');")		
		self.disconnect_db()
		# Add the canonical transcript for each gene into the dictionary
		for tx in canonical_transcripts:
			self.gene_transcript_dict[tx['geneSymbol']]['canonical'] = tx['refseq']

	def get_clinVar(self):
		# Get the most reported transcript from the last 500 records in ClinVar for each gene
		# This can take a few minutes if gene list is big
		for gene in self.gene_transcript_dict:
			payload = {'db': 'clinvar',
 			'term': gene + '[gene]',
			'retmax': 500,
			'retmode': 'json'}
			# First request pulls out JSON containing the IDs for each ClinVar variant record for a gene (max 500 entries, most recent first)
			r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params = payload)
			j = r.json()
			ids = j['esearchresult']['idlist']
			# Second request pulls out an XML containing the variant reports in the list of IDs retrieved in first request
			# (This endpoint doesn't support JSON)	
			payload = {'db': 'clinvar',
					   'rettype': 'variation',
					   'id': ','.join(ids)}
			r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params = payload)
			# read the XML content into an element tree object so that it can be parsed
			root = ET.fromstring(r.content)
			transcripts = {}
			# loop through XML and pull out transcript from variant name.
			# Keep a count of how many times a transcript is used in dictionary
			for variationReport in root.iter('VariationReport'):
				# Pull out the 'NM' accession from the transcript
				if variationReport.attrib['VariationName'].startswith('NM_'):
					# Remove the version number of the transcript by splitting on '.' and taking first element
					transcript = variationReport.attrib['VariationName'].split('.')[0]
					# Increase the count for that transcript in the dictionary by 1
					# Set default sets the count to 0 if the transcript is a new entry in the dictionary
					transcripts[transcript] = transcripts.setdefault(transcript, 0) + 1
			# Find transcript for this gene with the highest count
			# Build a list in format [(count, transcript), (count, transcript)...]
			counts_transcripts = [(count, tx) for tx, count in transcripts.items()]
			# max will order the list by the first element of each tuple (count) and return the tuple with highest value (i.e. highest count)
			# Capture the transcript accession from that tuple and add to dictionary
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
		# Loop through each transcript for the gene
		#Takes the gene_transcript_dict from Genes object and uses it to calculate the cDNA and aa start and end positions for the amplicon
		# tx_details is a dictionary containing the data for that transcript
		for transcript, tx_details in gene_transcript_dict[self.gene]['transcripts'].items():
			#Call appropriate method (dependent on +ve or -ve strand) to generate dictionary containing cDNA and aa starts, ends etc
			# tx_details is a dictionary containing the information generated from the Genes object for the transcript 
			if tx_details['strand'] == "+":
				t_dict = self.get_cDNA_aa_positive(tx_details)
			elif tx_details['strand'] == "-":
				t_dict = self.get_cDNA_aa_negative(tx_details)
			self.cds_aa_coords[transcript] = t_dict
		# return the cds_aa_coords dictionary
		return self.cds_aa_coords

	def get_cDNA_aa_positive(self, tx_details):
		##Convert the genomic start and end to cDNA start and end for POSITIVE strand transcripts
		# Initialise variables with None values
		amp_cDNA_start = None
		amp_cDNA_end = None
		amp_cDNA_start_diff = None
		amp_cDNA_end_diff = None
		aa_start_full = None
		aa_starts_in = None
		aa_end_full = None
		aa_ends_in = None
		# Store the amplicon start and end genomic coordinates
		amp_g_start = self.start
		amp_g_end = self.end
		# g2cDNA contains mappings from genomic to cds coordinates for this transcript in format [((exon1_g_start, exon1_g_end), (exon1_cDNA_end, exon1_cDNA_start)), ((exon2_g_start, exon2_g_end), (exon2_cDNA_end, exon2_cDNA_start)), ...etc]
		g2cDNA = tx_details['g2cDNA']
		#If there is no overlap between the transcript and amplicon, set cDNA start and end to 0
		# If the amplicon genomic start and end coords are both LESS than the FIRST exon START coord, 
		# OR the amplicon genomic start and end coords are both GREATER than the LAST exon END coord
		# then there is no overlap between the amplicon and transcript, so set all values to 0 
		if ((amp_g_start < g2cDNA[0][0][0] and amp_g_end < g2cDNA[0][0][0])
			or (amp_g_start > g2cDNA[-1][-2][1] and amp_g_end > g2cDNA[-1][-2][1])):
			amp_cDNA_start = 0
			amp_cDNA_end = 0
			amp_cDNA_start_diff = 0
			amp_cDNA_end_diff = 0
			aa_start_full = 0
			aa_starts_in = 0
			aa_end_full = 0
			aa_ends_in = 0
		# If the amplicon does overlap with the transcript, the value for amp_cDNA_start will still be None (rather than 0).
		# If this is the case, loop through exons and find exon (or intron) containing start position
		exon = 0 # Counter for the exon number
		# Loop until amp_cDNA_start is set to a value
		while amp_cDNA_start == None:
			# Capture the exon start and end genomic coordinates tuple from g2cDNA and store in g_exon e.g. (exon1_g_start, exon1_g_end)
			# Capture the exon start and end cDNA coordinates tuple from g2cDNA and store in c_exon e.g. (exon1_cDNA_start, exon1_cDNA_end)
			# Use the exon counter as an index, will increase each time the loop executes
			g_exon, c_exon = g2cDNA[exon]
			# If the amplicon start lies upstream of the exon start (e.g. in an intron), the first coding base in the amplicon is the first base of the current exon
			# Therefore set cDNA start position of the amplicon to the cDNA start position of the exon
			if amp_g_start < g_exon[0]:
				amp_cDNA_start = c_exon[0] #First cDNA base included in the transcript
				# Next work out how many bases upstream of the exon the amplicon actually starts so this can be included in output
				# This is the difference between the genomic start of the amplicon and the genomic start of the exon 
				amp_cDNA_start_diff = g_exon[0] - amp_g_start #Calculate number of bases upstream of the start cDNA base
			# Else if the genomic start lies within the exon, calculate the cDNA position and set it as amp_cDNA_start
			# If the genomic start of the amplicon is >= the genomic start of the exon AND the genomic start of the amplicon is <= the genomic end of the exon, the amplicon starts within the exon
			elif amp_g_start >= g_exon[0] and amp_g_start <= g_exon[1]:
				# The cDNA start for the amplicon is equal to the cDNA start of the exon + the difference between genomic starts of the amplicon and exon
				amp_cDNA_start = c_exon[0] + (amp_g_start - g_exon[0])
				# The amplicon starts within the exon, so there are no bases included that are upstream of the exon. Therefore set amp_cDNA_start_diff to 0
				amp_cDNA_start_diff = 0 #cDNA is not upstream exon, so diff is 0
			# Increase the exon count by 1. If neither of the above rules were met, amp_cDNA_start will still be None, so the while loop will continue, moving to next exon
			exon += 1 #Move to next exon
		
		# To get the amino acid start pos add 2 to cDNA start and divide by 3 e.g. if the amplicon started at cDNA position 10. Would need to add 2 (= 12), then divide by 3 to get amino acid start (4)
		# If the amplicon starts partway through a codon, the result of dividing by 3 will be a float
		# To get the first full codon included in the amplicon, this number should be rounded up. 
		# To get the codon that the amplicon starts within, this number should be rounded down. 
		# e.g if the transcript started at cDNA position 5 (so partway through codon 2). Would add 2 (=7), / 3 (=2.3)
		# Round 2.3 up to get the first full codon included in the amplicon (3)
		# Round 2.3 down to get the codon that the amplicon actually starts in (2)
		aa_start_full = int(math.ceil((amp_cDNA_start+2)/3.0)) #math.ceil rounds up
		aa_starts_in = int(math.floor((amp_cDNA_start+2)/3.0)) #math.floor round down

		#Next find the cDNA end position. Works backwards from the last exon of the transcript.
		# Set exon counter to the last exon in list
		exon = len(g2cDNA) - 1 #index of the last exon in list, subtract 1 to convert from 1 based to 0 based
		# If the amplicon does overlap with the transcript, the value for amp_cDNA_end will still be None (rather than 0).
		# If this is the case, loop through exons and find exon (or intron) containing end position
		while amp_cDNA_end == None:
			# Capture the exon start and end genomic coordinates tuple from g2cDNA and store in g_exon e.g. (exon1_g_start, exon1_g_end)
			# Capture the exon start and end cDNA coordinates tuple from g2cDNA and store in c_exon e.g. (exon1_cDNA_start, exon1_cDNA_end)
			# Use the exon counter as an index, will decrease each time the loop executes (because we're working backwards from last exon)			
			g_exon, c_exon = g2cDNA[exon]
			# If the amplicon end lies downstream of the exon end (e.g. in an intron), the last coding base in the amplicon is the last base of the current exon
			# Therefore set cDNA end position of the amplicon to the cDNA end position of the exon
			if amp_g_end > g_exon[1]:
				amp_cDNA_end = c_exon[1] #Set the cDNA base position to the last base in the exon
				# Next work out how many bases downstream of the exon the amplicon actually ends so this can be included in output
				# This is the difference between the genomic end of the amplicon and the genomic end of the exon 
				amp_cDNA_end_diff = amp_g_end - g_exon[1] #Calculate number of bases downstream of the last cDNA base
			#Else if the genomic end lies within the exon, calculate cDNA position and set it as amp_cDNA_end
			# If the genomic end of the amplicon is >= the genomic start of the exon AND the genomic end of the amplicon is <= the genomic end of the exon, the amplicon ends within the exon
			elif amp_g_end >= g_exon[0] and amp_g_end <= g_exon[1]: 
				# The cDNA end for the amplicon is equal to the cDNA start of the exon + the difference between genomic end of the amplicon and the genomic start of the exon
				amp_cDNA_end = c_exon[0] + (amp_g_end - g_exon[0])
				# The amplicon ends within the exon, so there are no bases included that are downstream of the exon. Therefore set amp_cDNA_end_diff to 0
				amp_cDNA_end_diff = 0 #cDNA end is not downstream of the exon, so diff is 0
			exon -= 1 #Move to next exon (reverse order)
		# Get the last amino acid/codon fully included in the amplicon
		# Here can just divide the cDNA end position by 3. e.g. if the the cDNA end position was 9, the codon end woruld be 9/3 (=3)
		# If the amplicon ends partway through a codon, the result of dividing by 3 will be a float
		# To get the last full codon included in the amplicon, this number should be rounded down. 
		# To get the codon that the amplicon ends within, this number should be rounded up. 
		# e.g if the transcript ended at cDNA position 10 (so partway through codon 4). Would / 3 (= 3.3).
		# Round 3.3 down to get the last full codon included in the amplicon (3)
		# Round 3.3 up to get the codon that the amplicon actually ends in (2)
		aa_end_full = int(math.floor((amp_cDNA_end)/3.0))
		aa_ends_in = int(math.ceil((amp_cDNA_end)/3.0))

		# Add the calculated values to a dictionary and return them
		t_dict = {}
		t_dict['cDNA_start'] = amp_cDNA_start
		t_dict['amp_cDNA_start_diff'] = amp_cDNA_start_diff
		t_dict['cDNA_end'] = amp_cDNA_end
		t_dict['amp_cDNA_end_diff'] = amp_cDNA_end_diff
		t_dict['aa_start_full'] = aa_start_full
		t_dict['aa_starts_in'] = aa_starts_in
		t_dict['aa_end_full'] = aa_end_full
		t_dict['aa_ends_in'] = aa_ends_in
		return t_dict

	def get_cDNA_aa_negative(self, tx_details):
		##Convert the genomic start and end to cDNA start and end for NEGATIVE strand transcripts
		# Initialise variables with None values
		amp_cDNA_start = None
		amp_cDNA_end = None
		amp_cDNA_start_diff = None
		amp_cDNA_end_diff = None
		aa_start_full = None
		aa_starts_in = None
		aa_end_full = None
		aa_ends_in = None	
		# Store the amplicon start and end genomic coordinates	
		amp_g_start = self.start
		amp_g_end = self.end
		# g2cDNA contains mappings from genomic to cds coordinates for this transcript in format [((exonN_g_start, exonN_g_end), (exonN_cDNA_end, exonN_cDNA_start)), ((exonN-1_g_start, exonN-1_g_end), (exonN-1_cDNA_end, exonN-1_cDNA_start)), ...etc]
		# Note that g2cDNA is ordered by increaseing genomic coordinates. This means that the FIRST exon in the transcript is the LAST element of the list.
		# The genomic coordinates will INCREASE, whilst the cDNA coordinates will DECREASE as you iterate through the list.
		g2cDNA = tx_details['g2cDNA']
		#If there is no overlap between the transcript and amplicon, set cDNA start and end to 0
		# If the amplicon genomic start and end coords are both LESS than the FIRST exon START coord, 
		# OR the amplicon genomic start and end coords are both GREATER than the LAST exon END coord
		# then there is no overlap between the amplicon and transcript, so set all values to 0
		if ((amp_g_start < g2cDNA[0][0][0] and amp_g_end < g2cDNA[0][0][0])
			or (amp_g_start > g2cDNA[-1][-2][1] and amp_g_end > g2cDNA[-1][-2][1])):
			amp_cDNA_start = 0
			amp_cDNA_end = 0
			amp_cDNA_start_diff = 0
			amp_cDNA_end_diff = 0
			aa_start_full = 0
			aa_starts_in = 0
			aa_end_full = 0
			aa_ends_in = 0
		# If the amplicon does overlap with the transcript, the value for amp_cDNA_start will still be None (rather than 0).
		# If this is the case, loop through exons and find exon (or intron) containing start position
		# Needs to work backwards because transcript is on reverse strand, so FIRST exon is LAST element in list.
		exon = len(g2cDNA) - 1 #index of the last exon in list, subtract 1 to convert from 1 based to 0 based
		while amp_cDNA_start == None:
			g_exon, c_exon = g2cDNA[exon]
			# If amplicon genomic end is greater than the genomic start of the first exon of transcript
			# the start of the amplicon is therefore upstream (on reverse strand) of the cDNA start 
			if amp_g_end > g_exon[1]:
				# Set the amplicon cDNA start position to the first base of the first exon  (reverse strand so this is actually the last cDNA coordinate in the list)
				amp_cDNA_start = c_exon[1]
				# Next work out how many bases upstream (on rev strand) of the exon the amplicon actually starts so this can be included in output
				# This is the difference between the genomic end of the amplicon and the genomic start of the exon
				amp_cDNA_start_diff = amp_g_end - g_exon[1]
			#Else if the genomic end lies within the exon, calculate cDNA position and set it as amp_cDNA_start
			# If the genomic end of the amplicon is >= the genomic end of the exon AND the genomic end of the amplicon is <= the genomic start of the exon, the amplicon starts within the exon
			elif amp_g_end >= g_exon[0] and amp_g_end <= g_exon[1]:
				# The cDNA start for the amplicon is equal to the cDNA start of the exon + the difference between genomic start of the exon and the genomic end of the amplicon
				amp_cDNA_start = c_exon[1] + (g_exon[1] - amp_g_end)
				# The amplicon starts within the exon, so there are no bases included that are upstream of the exon. Therefore set amp_cDNA_start_diff to 0
				amp_cDNA_start_diff = 0 #cDNA end is not downstream of the exon, so diff is 0
			exon -= 1 #Move to next exon (reverse strand so have to work backwards from end of list)
		# To get the amino acid start pos add 2 to cDNA start and divide by 3 e.g. if the amplicon started at cDNA position 10. Would need to add 2 (= 12), then divide by 3 to get amino acid start (4)
		# If the amplicon starts partway through a codon, the result of dividing by 3 will be a float
		# To get the first full codon included in the amplicon, this number should be rounded up. 
		# To get the codon that the amplicon starts within, this number should be rounded down. 
		# e.g if the transcript started at cDNA position 5 (so partway through codon 2). Would add 2 (=7), / 3 (=2.3)
		# Round 2.3 up to get the first full codon included in the amplicon (3)
		# Round 2.3 down to get the codon that the amplicon actually starts in (2)
		aa_start_full = int(math.ceil((amp_cDNA_start+2)/3.0)) #math.ceil rounds up
		aa_starts_in = int(math.floor((amp_cDNA_start+2)/3.0)) #math.floor round down

		#If amp_cDNA_end is still None after the first step, loop through exons and find cDNA end position
		exon = 0
		while amp_cDNA_end == None:
			# Capture the exon start and end genomic coordinates tuple from g2cDNA and store in g_exon e.g. (exon1_g_start, exon1_g_end)
			# Capture the exon start and end cDNA coordinates tuple from g2cDNA and store in c_exon e.g. (exon1_cDNA_start, exon1_cDNA_end)
			# Use the exon counter as an index, will increase each time the loop executes			
			g_exon, c_exon = g2cDNA[exon]
			# If amplicon genomic start is less than the genomic end of the last exon of transcript
			# the end of the amplicon is therefore downstream (on reverse strand) of the cDNA end 
			if amp_g_start < g_exon[0]:
				# Set the amplicon cDNA end position to the last base of the last exon (reverse strand so this is actually the first cDNA coordinate in the list)
				amp_cDNA_end = c_exon[0] #First cDNA base included in the transcript
				# Next work out how many bases downstream (on rev strand) of the exon the amplicon actually ends so this can be included in output
				# This is the difference between the genomic end of the exon and the genomic end of the amplicon
				amp_cDNA_end_diff = g_exon[0] - amp_g_start
			#Else if the genomic amplicon start lies within the exon, calculate cDNA position and set it as amp_cDNA_end
			# If the genomic start of the amplicon is >= the genomic end of the exon AND the genomic start of the amplicon is <= the genomic start of the exon, the amplicon ends within the exon
			elif amp_g_start >= g_exon[0] and amp_g_start <= g_exon[1]:
				# The cDNA end for the amplicon is equal to the cDNA start of the exon + the difference between genomic starts of the amplicon and exon
				amp_cDNA_end = c_exon[1] + (g_exon[1] - amp_g_start)
				# The amplicon starts within the exon, so there are no bases included that are upstream of the exon. Therefore set amp_cDNA_start_diff to 0
				amp_cDNA_end_diff = 0 #cDNA is not upstream exon, so diff is 0
			exon += 1 #Move to next exon
		# Get the last amino acid/codon fully included in the amplicon
		# Here can just divide the cDNA end position by 3. e.g. if the the cDNA end position was 9, the codon end woruld be 9/3 (=3)
		# If the amplicon ends partway through a codon, the result of dividing by 3 will be a float
		# To get the last full codon included in the amplicon, this number should be rounded down. 
		# To get the codon that the amplicon ends within, this number should be rounded up. 
		# e.g if the transcript ended at cDNA position 10 (so partway through codon 4). Would / 3 (= 3.3).
		# Round 3.3 down to get the last full codon included in the amplicon (3)
		# Round 3.3 up to get the codon that the amplicon actually ends in (2)
		aa_end_full = int(math.floor((amp_cDNA_end)/3.0))
		aa_ends_in = int(math.ceil((amp_cDNA_end)/3.0))

		# Add the calculated values to a dictionary and return them
		t_dict = {}
		t_dict['cDNA_start'] = amp_cDNA_start
		t_dict['amp_cDNA_start_diff'] = amp_cDNA_start_diff
		t_dict['cDNA_end'] = amp_cDNA_end
		t_dict['amp_cDNA_end_diff'] = amp_cDNA_end_diff
		t_dict['aa_start_full'] = aa_start_full
		t_dict['aa_starts_in'] = aa_starts_in
		t_dict['aa_end_full'] = aa_end_full
		t_dict['aa_ends_in'] = aa_ends_in
		return t_dict
