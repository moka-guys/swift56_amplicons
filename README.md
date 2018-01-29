# swift56_amplicons

The oncology team have provided a file containing the chromosome, genomic start, genomic stop and gene symbol for the 264 amplicons in the Swift56 oncology panel. They requested that this be annotated with cDNA and amino acid start and stop coordinates.

The code contained in the `Genes` class of `gene_transcript_library.py` finds all RefSeq transcripts for a given list of gene symbols. It also identifies the canonical transcript and the transcript most reported in ClinVar for each gene.

The code contained in the `Amplicon` class of `gene_transcript_library.py` takes a chromosome, start, stop and gene for each amplicon from the input file, and uses the data contained in the Genes object to calculate the cDNA and amino acid coordinates for each transcript for the gene.

There were a few differences between the 'canonical' and 'most reported in ClinVar' transcripts. The oncology team reviewed both lists of transcripts and decided to go with the 'most reported in ClinVar' transcripts.

`swift_56_get_results.py` parses the input file (`swift_56_amplicons.txt`) and uses the `Genes` and `Amplicon` objects from the `gene_transcript_library` module to find the cDNA and amino acid coordinates for each amplicon. This data is then output to `SWIFT56_output.txt` (included in this repo).

NOTES:
The script was tested using the swift5 amplicons that had previously been calculated manually. See 180124_swift56_amplicons folder on MokaNAS2.

Explanation of some of the (less obvious) output fields:
* upstream_bases = bases upstream of the cDNA start site
* downstream_bases = bases downstream of the cDNA stop site
* starts_in_codon = the codon the amplicon starts in
* first_full_codon = the first codon that is fully covered in the amplicon. e.g. if the amplicon starts 2 bases into codon 10, the start_in_codon would be 10, but the first_full_codon would be 11
* ends_in_codon & last_full_codon = equivalent of above

