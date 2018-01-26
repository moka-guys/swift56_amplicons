from gene_transcript_library import Genes, Amplicon
# input file
i_file = '/home/andy/Dropbox/STP/Oncology/swift_56_amplicons.txt'
# output file
o_file = '/home/andy/Dropbox/STP/Oncology/SWIFT56_output.txt'

# Create list of genes from input file 
all_genes = set([])
with open(i_file, 'r') as i:
	for line in i:
		if not line.startswith('#'):
			all_genes.add(line.split('\t')[1].rstrip())

# Create Genes object
g = Genes(all_genes)
# Get data for refseq transcripts from UCSC 
g.get_UCSC_refGene()
# Find the transcript that is reported most frequently in ClinVar
g.get_clinVar()
# get exon coordinates and map genomic posititons to cDNA positions for each transcript. (This info is stored in the object's gene_transcript_dict dictionary)
g.get_cds()

# Open the input and output files
with open(i_file, 'r') as i, open(o_file, 'w') as o:
    # Write the header row to output file
	o.write("\t".join([
        '#Swift_Amplicon_Name', 
        'Chr', 
        'Start', 
        'End', 
        'Gene', 
        'Transcript', 
        'Orientation', 
        'cDNA_start', 
        'upstream_bases', 
        'cDNA_end', 
        'downstream_bases', 
        'starts_in_codon',
        'first_full_codon',
        'ends_in_codon',
        'last_full_codon'])
        + '\n'
        )
	# Loop through the amplicons in the input file
	for line in i:
		# If the line isn't a header row
		if not line.startswith('#'):
			# Capture the fields stored in the input file
			gene = line.split('\t')[1].rstrip()
			# First field is [chr]:[start]-[stop]
			# Need to split out into chr, start and stop
			chrom_start_end = line.split('\t')[0].rstrip()
			chrom, start_end = chrom_start_end.split(':') # SPlit on colon to separate chrom from start-stop
			# Split start-stop on - to separate, and convert from strings to integers
			start, end = map(int, start_end.split('-'))
			swift_amplicon_name = line.split('\t')[2].rstrip()
			# Create and amplicon object by passing the gene, chrom, start and end
			a = Amplicon(gene, chrom, start, end)
			# Capture the most reported transcript (NM...) for that gene from the Genes object gene_transcript_dict dictionary
			clinvar_most_reported = g.gene_transcript_dict[gene]['clinvar_most_reported']
			# Capture the dictionary containing transcript level data from gene_transcript_dict
			transcript_data = g.gene_transcript_dict[gene]['transcripts'][clinvar_most_reported]
			# Call the get_cDNA_aa() method of the Amplicon object which calculates the cDNA and aa positions of the amplicon in each transcript for that gene and returns a nested dictionary
			# Capture the amplicon level data for the most reported transcript in clinvar from the dictionary
			amplicon_data = a.get_cDNA_aa(g.gene_transcript_dict)[clinvar_most_reported]
			# Write out the fields to output file
			o.write("\t".join(map(str, [
				swift_amplicon_name,
				chrom,
				start,
				end,
				gene,
				clinvar_most_reported,
				transcript_data['strand'],
				amplicon_data['cDNA_start'],
				amplicon_data['amp_cDNA_start_diff'],
				amplicon_data['cDNA_end'],
				amplicon_data['amp_cDNA_end_diff'],
				amplicon_data['aa_starts_in'],
				amplicon_data['aa_start_full'],
				amplicon_data['aa_ends_in'],
				amplicon_data['aa_end_full']
			]))
			+ '\n'
			)
