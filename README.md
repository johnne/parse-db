# parse-db
This repository contains scripts to process information in sequence and metabolic pathway datases, as well as scripts to parse output from bioinformatic tools.

## Cluster of Orthologous Groups (COG) database
The COG database contains several orthologous groups of genes which are categorized into one or more functional categories (currently 25 in total). The [process_cog_db.py](https://github.com/johnne/parse-db/blob/master/cog/process_cog_db.py) script will simply read from two files on the [COG ftp server](https://ftp.ncbi.nih.gov/pub/COG/COG2014/data/) (fun2003-2014.tab and cognames2003-2014.tab) and produce a tab-delimited file with each COG listed with its one-letter code and functional category name. **Note that for COGs classified into more than one functional category the COG is printed on multiple lines, one per functional category.**

## HMMER
The [HMMER](http://hmmer.org/) suite uses hidden markov models of protein alignments to find homologs for sequences. The [parse_hmmout.py](https://github.com/johnne/parse-db/blob/master/hmmer/parse_hmmout.py) script parses output from hmmsearch and hmmscan and prints best scoring hits for each query sequence. The user can set cutoffs using an e-value or by supplying "trusted score cutoffs" (see [the example](https://github.com/johnne/parse-db/blob/master/hmmer/example/example.trusted_cutoffs.tab)) for each HMM. This is analogous to running hmmsearch with --cut_tc which uses the trusted cutoff listed for HMMs in databases such as PFAM and TIGRFAM. By default the script prints the best scoring hit that passed the threshold per query but if there are several hits on the same query sequence that pass the threshold they can all be stored. The allowed overlap (as a fraction of the shorter HMM alignment) between multiple hits on the same query can be set by the user.

A simple run that stores only the highest scoring hit with an e-value < 1e-5 per query sequence:

`python parse_hmmout.py -i example/example.hmmout > parsed.tab`

Using trusted cutoffs instead of e-value:

`python parse_hmmout.py -i example/example.hmmout -t example/example.trusted_cutoffs.tab > parsed.tab`

Allowing multiple hits to the same query sequence (no overlap allowed):

`python parse_hmmout.py -i example/example.hmmout --overlap 0 > parsed.tab`

Allowing multiple hits to the same query sequence with 10% overlap on the shortest alignment allowed:

`python parse_hmmout.py -i example/example.hmmout --overlap 0.1 > parsed.tab`

## Metacyc
The [Metacyc](http://metacyc.org/) database contains information on enzymes and metabolic pathways. Flatfiles can be downloaded after [applying](http://metacyc.org/download.shtml) for a (free) license. After that, download the MetaCyc tarball and run [process_metacyc_db.py](https://github.com/johnne/parse-db/blob/master/metacyc/process_metacyc_db.py) from the unzipped data/ directory, or provide paths to the classes.dat, compounds.dat, pathways.dat and reactions.dat files to create a tab-delimited file linking enzyme EC numbers to pathways and higher order metabolic categories.
