# swift56_amplicons
This repo contains code to find start and stop cDNA and codon coordinates for the 264 amplicons in the Swift56 oncology panel.

`swift_56_get_results.py` parses the input file (`swift_56_amplicons.txt`) containing genomic coordinates and Gene names for each amplicon. It then uses the Genes and Amplicon objects in `gene_transcript_library.py` to find the required information and outputs `SWIFT56_output.txt`.

NOTES:
The script was tested using the swift5 amplicons that had previously been calculated manually. See 180124_swift56_amplicons folder on MokaNAS2.

Explanation of some of the output fields:
* upstream_bases = bases upstream of the cDNA start site
* downstream_bases = bases downstream of the cDNA stop site
* starts_in_codon = the codon the amplicon starts in
* first_full_codon = the first codon that is fully covered in the amplicon. e.g. if the amplicon starts 2 bases into codon 10, the start_in_codon would be 10, but the first_full_codon would be 11
* ends_in_codon & last_full_codon = equivalent of above








