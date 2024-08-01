# local_align
PLAZZOTTA GIOVANNI, MAT. 232312, ALGORITHMS FOR BIOINFORMATICS (QCB), 2023/2024, UNITN
PERSONAL IMPLEMENTATION OF THE SMITH - WATERMAN PAIRWISE LOCAL ALIGNMENT ALGORITHM.

The algorithm relies on the construction of two different matrices: the 'scoring matrix' is populated with scores according to the strategy developed by Smith and Waterman, 
while the 'origin matrix' is populated with pointers (both in the form of explicit coordinates and with an arrow symbol that can be horizontal '-', vertical '|' or diagonal '\\') 
that point to the nearby cell in the scoring matrix that allowed to achieve the (highest) score in the current cell of the scoring matrix. In the current implementation, the 
user can provide the two sequences to locally align, and a series of parameters (seed, match_score, mismatch_score, gap_opening_penalty, gap_extension_penalty, cost_model, 
num_results, verbose). All parameters are optional and have default values, except for the sequences, that are compulsory and must be fed as ungapped strings of ACTG characters. 
The seed is used for reproducible results (if multiple 'arrows' are available from the same position in the scoring matix, only one arrow is chosen, randomly). The parameter 
num_results allows to set the desired number of alignments to print in output, by default ordered in decreasing value for the final alignment score. Please provide all numerical 
parameters as positive integers (the code will take care of handling these values in the correct way). Please provide all parameters in the format -parameter=value.

Required packages:
- random
- argparse

Please read the instruction before use, with:
./plazzotta_giovanni_localalign.py -h 
