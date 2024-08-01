# local_align
PLAZZOTTA GIOVANNI, MAT. 232312, ALGORITHMS FOR BIOINFORMATICS (QCB), 2023/2024, UNITN
PERSONAL IMPLEMENTATION OF THE SMITH - WATERMAN PAIRWISE LOCAL ALIGNMENT ALGORITHM.

The algorithm relies on the construction of two different matrices: the 'scoring matrix' is populated with scores according to the strategy developed by Smith and Waterman, 
while the 'origin matrix' is populated with pointers (both in the form of explicit coordinates and with an arrow symbol that can be horizontal '-', vertical '|' or diagonal '\\') 
that point to the nearby cell in the scoring matrix that allowed to achieve the (highest) score in the current cell of the scoring matrix. In the current implementation, the 
user can provide the two sequences to locally align, and a series of parameters (seed, match_score, mismatch_score, gap_opening_penalty, gap_extension_penalty, cost_model, 
num_results, verbose). All parameters are optional and have default values, except for the sequences, which are compulsory and must be fed as ungapped strings of ACTG characters. 
The seed is used for reproducible results (if multiple 'arrows' are available from the same position in the scoring matrix, only one arrow is chosen, randomly). The parameter 
num_results allows to set the desired number of alignments to print in output, by default ordered in decreasing value for the final alignment score. Please provide all numerical 
parameters as positive integers (the code will take care of handling these values in the correct way). Please provide all parameters in the format -parameter=value.

The local aligner is provided in two versions. In the first (default) one, termed 'original', the backtracking procedure is repeated only num_results times, starting only from the
highest scoring positions in the scoring matrix (in decreasing order of score). The alignments produced this way are directly given in output and are naturally ordered according to
their alignment score. The original version thus allows to quickly recover and output the desired number of local alignments, starting from the highest scoring ones. The second version,
termed 'alternative', is more versatile but computationally expensive. Instead of producing and giving in output only the num_results best alignments ordered by score, the code runs the 
backtracking procedure from all cells in the scoring matrix, always respecting the order instructed by the alignment score (so always starting the very first backtracking from the cell 
carrying the highest score, and then performing the backtracking from the cell carrying the second highest score, if this has not been visited yet and if this contains at least one arrow,
and so on...). Each time, the backtracking procedure is performed exactly as described in the first version of the code. Generating all of these alignments allows the user to sort them on 
the basis of multiple statistics, that may be different from the simple alignment score: for example, the user may be interested in getting in output a certain number of alignments, but 
sorted according to the length of the alignment, or according to the number of gaps, and not according to the alignment score. To recover the alignment with the greatest length or with the
highest number of gaps, it is not wise to limit the search to the num_results alignments having the highest alignment score: it is necessary to look through the entire matrix. 

Required packages:
- random
- argparse

Please read the instruction before use, with:
./plazzotta_giovanni_localalign.py -h 
