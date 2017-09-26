RNAv

GENERAL INFORMATION
-------------------
RNAv is a program that searches for RNA secondary structure variations based on the notion of a structure graph to specify the consensus structure of an RNA family.

This program is implemented in C++. It has been compiled and tested on several systems, including Desktop Linux computers, a Linux cluster, and a SUN workstation running SunOS 5.1.

REFERENCE
---------
1.	Zhibin Huang, Russell L. Malmberg, Mohammad Mohebbi, Liming Cai. 
	RNAv: Non-coding RNA Secondary Structure Variation Search via Graph Homomorphism, 
	In Proceedings of Computational Systems Bioinformatics Conference (CSB 2010), 
	August, 2010. Vol. 9, p. 56-69.

2.	Yingfeng Wang, Zhibin Huang, Yong Wu, Russell L. Malmber and Liming Cai, 
	RNATOPS-W: a web server for RNA structure searches of genomes, 
	Bioinformatics, vol. 25, no. 8, pp. 1080-1081, 2009.

3.	Huang, Z., Wu, Y., Robertson, J., Feng, L., Malmberg, R., and Cai, L. 
	Fast and accurate search for non-coding RNA pseudoknot structures in genomes, 
	Bioinformatics, 2008 24(20):2281-2287.

Copyright (C) 2010 University of Georgia. All Rights Reserved.

WHAT IS CONTAINED IN THIS PACKAGE
---------------------------------
The software package of RNAv and one example are contained in this package.

INSTALLATION
------------
Compile RNAv.
Type "make" to build an executable program called "rnav".

USAGE
-----
1. How to run this program
Usage: rnav <-tf training_file> <-gf genome_file> [other options]

2. Command line parameter list:
[-tf]			training_file
[-gf]			genome_file
[-hmmfilter]	search with automatic-hmmfilter
[-pcnt]			pseudocount_value [DEFAULT=0.001]
[-k]			k_value [DEFAULT=10]
[-th]			threshold_value_for_candidate_hit [DEFAULT=0.0]
[-no]			number of overlapping nts between two stems when overlap is allowed [DEFAULT=2]
[-mc]			to merge candidates in preprocessing
[-ms]			when taking the merge strategy, take the candidate with the (m)ax score 
				or the one with the (s)hortest length (m|s) [DEFAULT=s]
[-ns]			number of shift positions allowed in the merge strategy [DEFAULT=0]
[-ni]			number of insertions allowed in a null loop [DEFAULT=3]
[-pcv]			weight for prior base pairing matrix [DEFAULT=0.4]
[-pf]			prior_file_name [DEFAULT='./base_pair_prior.txt']
[-pv]			pcoeff value [DEFAULT=2.0]
[-st]			split threshold [DEFAULT=6]
[-js]			stepsize in jump strategy [DEFAULT=1]
[-jt]			score_filtering threshold in jump strategy [DEFAULT=0.0]
[-r]			to also search reverse complement strand
[-ps]			to print structure alignment info
[-d]			to print out the debug info
[-time]			to print out the time info
[-pscore]		to print score info for stem and loop
[-sopt]			score option [0 - stem only; 1 - stem + positive-loop; 2 - stem + loop; [DEFAULT=0]
[-swt]			weight factor for the stem score [DEFAULT=1.0]
[-lwt]			weight factor for the loop score [DEFAULT=1.0]
[-ept]			max number of empty stem [DEFAULT=0]
[-dfilter]		to print out the debug info for the filter based search
[-ddpinput]		to print out the debug info of the input in dp search
[-ddptree]		to print out the debug info of the tree in dp search
[-ddpwin]		to print out the debug info of the window in dp search
[-ddptd]		to print out the debug info of the td in dp search
[-ddpsearch]		to print out the debug info of the dp search

3. Command line parameter explanation:
[-tf]	training_file
File of RNA training data in the pasta format to produce the RNA structure model.
Two examples (one pk-free and one pk) are included in the directory of "example".

[-gf]	genome_file
File to be searched for the structure modeled with the training data.
Current version requires that genome_file only contain ACGT nucleotide characters.
See Section of "PREPARATION OF THE GENOME FILE".

[-hmmfilter]	search with automatic-hmmfilter
This option is to do the automatically selected hmm filter-based fast search,
If not specified this option, then the program will do the whole-structure search.

[-pcnt] pseudocount_value [DEFAULT=0.001]
The value of pseudocount.
Default value is 0.001.

[-k]	k_value [DEFAULT=10]
The number of candidates for each stem in the structure to be searched in genome_file.
Default value is 10.

[-th]	threshold_value_for_candidate_hit [DEFAULT=0.0]
the threshold value which is used to determine whether the current structure is hit or not.
Default value is 0.0.

[-no]	number of overlapping nts between two stems when overlap is allowed [DEFAULT=2]
The number of overlapping nucleotides allowed between adjacent candidates.
Default value is 2.

[-mc]	to merge candidates in preprocessing
If you use this option, the strategy of merging similar candidates will be taken.
By default, it is used when doing the structure search.

[-ms]	when taking merge strategy, take the candidate with (m)ax score or one with the (s)hortest length (m|s) [DEFAULT=s]
When using the -mc option, merging candidates, then you need to choose which merging strategy you will choose.
Choosing the representive candidate with max score or shortest length.
m: max score;  s: shortest length.
Default value is s.

[-ns]	number of shift positions allowed in the merge strategy [DEFAULT=0]
When merging candidates, the candidates within the number of shift positions will be merged together.
Default value is 0.

[-ni]	number of insertions allowed in null loop [DEFAULT=3]
For the null loop model, the number of inserted nucleotides will be allowed. 
Default value is 3.

[-pcv]	prior weight value [DEFAULT=0.4]
The weight for prior base pairing frequency matrix. 
Default value is 0.4.

[-pf]	prior_file [DEFAULT='./base_pair_prior.txt']
The file containing prior pair frequency matrix (5x5). 
base_pair_prior.txt.
(Note: User can put this file any places as long as he needs to specify the absolute path of this file or 
he can use his own matrix file)

[-pv]	pcoeff value [DEFAULT=2.0]
The weight of distance penalty. 
Default value is 2.0.

[-st]	split threshold [DEFAULT=6]
The standard deviation value to group the training sequences.  
Default value is 6.

[-js]	stepsize in jump strategy [DEFAULT=1]
This program supports the skip-and-jump strategy.
The scanning window is shifted by the 'stepsize' nucleotides to speed up search.
2 or 3 for the stepsize is suggested. 

[-jt]	score_filtering threshold in jump strategy [DEFAULT=0.0]
Default value is 0.0.

[-r]	to search the complementary strand 
This option can be used to search the complementary genome strand.

[-ps]	to print out info of the folded structure and structure alignment.
This option can give you more info of the candidate hit the program finds.
Definitions of folded structure and structure alignment are given in "EXPLANATION OF THE RESULTS AND OUTPUT" in this file.

[-time]		to print out the time info
This option can give you the time info about your running RNATOPS.

[-pscore]	to print score info for stem and loop
This option can give you the score info for every stem and loop.

[-sopt]	the strategy used to compute the score, the first one is totally based on stem score					(0 - stem only) [DEFAULT]
												the second one is based on stem score and positive loop score;	(1 - stem + positive-loop)
												the third one is based on stem score and loop score;			(2 - stem + loop) 
[-swt]			weight factor for the stem score [DEFAULT=1.0]
[-lwt]			weight factor for the loop score [DEFAULT=1.0]
[-ept]			max number of empty stem [DEFAULT=0]
This option is to set the max number of empty stem in the structure you want to search.

[-dfilter]	to print out the debug info for the filter based search
[-ddpinput]	to print out the debug info of the input in dp search
[-ddptree]	to print out the debug info of the tree in dp search
[-ddpwin]	to print out the debug info of the window in dp search
[-ddptd]	to print out the debug info of the td in dp search
[-ddpsearch]	to print out the debug info of the dp search
[-d]		to print out the debug info

If you want to debug the program, you can use -d option.
We suggest you not debugging this program until you have a good understanding of the code 
because it will produce a large volume of data especially for the tree decomposition based dynamic programming.

4.  How to use RNAv to do the search.
Generally user can do two kinds of search: automatic hmm-filter-based search and whole-structure search in a stream line or direct whole-structure search. 

The data directory contains one sample using the "-hmmfilter" option and the script for running it is
./rnav -tf train_000.txt -gf test_000.fasta -hmmfilter -r -th 0 -k 10 -pscore -no 0 -ept 5 -sopt 1 -o train_000_result_debug.txt >train_000_result.txt

EXPLANATION OF THE RESULTS AND OUTPUT
-------------------------------------
The following information of every found structure candidate (of a score above the threshold) is generated as a result: 
	search direction;
	total alignment score, begin and end positions in the genome; 
	folded structure of the candidate; 
	structural alignment with the model;

The following is an example of running result from test_000.fasta.

------------------------------------------------------------
- Whole Structure Result                                   -
-                                                          -
- Profile file : train_000.txt                             -
- Profile length : 119                                     -
-                                                          -
- Filter type : HMM                                        -
- Filter info : positions from 0      to 118               -
-                                                          -
- Genome file : test_000.fasta                             -
- Number of sequences : 1                                  -
- Total length of sequences : 2074                         -
-                                                          -
- Search parameters setting:                               -
- Pseudocount = 0.001                                      -
- Num of stem candidates = 10                              -
- Score threshold for hits = 0                             -
- Num of nt overlap between stems = 0                      -
- Candidate representatives only = Yes                     -
- Shortest candidate representatives = Yes                 -
- IShiftNumMergeCand = No                                  -
- Nts allowed in null loops = 3                            -
- Pcoeff = 2                                               -
- Search with jump = Yes                                   -
- Search step size = 1                                     -
- Search reversed complement sequence = Yes                -
- Max num of empty stem = 5                                -
------------------------------------------------------------

Whole structure search hit 1
----------------------------
000
Plus search result
Hit Positions: 1000-1072
HIT Alignment score = 68.2662
Num of empty stem   = 0

Folded structure

         1 DDDDDDD..AAAA.........aaaa.BBBBB.......bbbbb.....CCCCC...... 61
         1 1111111..1111.........1111.11111.......11111.....11111...... 61
      1001 GCCGCCGTAGCTCAGCCCGGGAGAGCGCCCGGCTGAAGACCGGGTTGTCCGGGGTTCAAG 1061

        61 .cccccddddddd 74
        61 .111111111111 74
      1061 TCCCCGCGGCGGC 1074

Structure alignment

         1 DDDDDDD..AAAA...--....aaaa.BBBBB.......bbbbb...-.CCCCC...... 61
         1 1111111..1111.........1111.11111.......11111.....11111...... 61
      1001 GCCGCCGUAGCUCAGCCCGGGAGAGCGCCCGGCUGAAGACCGGGUUGUCCGGGGUUCAAG 1061

        61 .cccccddddddd 74
        61 .111111111111 74
      1061 UCCCCGCGGCGGC 1074

Whole structure search hit 2
----------------------------
000
Minus search result
Hit Positions: 986-1046
HIT Alignment score = 24.6886
Num of empty stem   = 1

Folded structure

         1 ^^^^^^^^AAAA......aaaaBBBBB.......bbbbb.....CCCCC.......cccc 61
         1 ^^^^^^^^1111......111111111.......11111.....11111.......1111 61
       987 GAGTTAAAGGGTTATGCCGCCGCGGGGACTTGAACCCCGGACAACCCGGTCTTCAGCCGG 1047

        61 c 62
        61 1 62
      1047 G 1048

Structure alignment

         1 ^^^^^^^^AAAA.......aaaaBBBBB.......bbbbb...-.CCCCC.......ccc 61
         1 ^^^^^^^^1111.......111111111.......11111.....11111.......111 61
       987 GAGTTAAAGGGUUAUGC-CGCCGCGGGGACUUGAACCCCGGACAACCCGGUCUUCAGCCG 1047

        61 cc 63
        61 11 63
      1047 GG 1049

------------------------------------------------------------
- Searched done : with RNATOPS V2.0                        -
- By RNA-Informatics @ UGA                                 -
- Total no of hits: 2                                      -
- Total time used 6.38889e-05     hours                    -
- Time: 19:14:05 EDT 2010-10-31                            -
------------------------------------------------------------

Note:
The folded structure is the predicted secondary structure of the found RNA candidate 
in the genome, where a base pair is annotated by a pair of letter. (See referenced pasta paper). 

The structure alignment shows the optimal alignment between the model and the 
RNA sequence found on the genome. It contains additional information of matches,
deletions, and insertions computed by the RNAv program.

CONTACT INFORMATION
-------------------
For suggestions, questions, requests and bug report:
please contact Liming Cai at cai@cs.uga.edu.
