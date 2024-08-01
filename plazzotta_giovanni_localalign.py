#!/usr/bin/env python3

# PLAZZOTTA GIOVANNI, MAT. 232312, ALGORITHMS FOR BIOINFORMATICS (QCB), 2023/2024, UNITN
# PERSONAL IMPLEMENTATION OF THE SMITH - WATERMAN PAIRWISE LOCAL ALIGNMENT ALGORITHM.

# The algorithm relies on the construction of two different matrices: the 'scoring matrix' is populated with scores according to the strategy developed by Smith and Waterman, 
# while the 'origin matrix' is populated with pointers (both in the form of explicit coordinates and with an arrow symbol that can be horizontal '-', vertical '|' or diagonal '\\') 
# that point to the nearby cell in the scoring matrix that allowed to achieve the (highest) score in the current cell of the scoring matrix. In the current implementation, the 
# user can provide the two sequences to locally align, and a series of parameters (seed, match_score, mismatch_score, gap_opening_penalty, gap_extension_penalty, cost_model, 
# num_results, verbose). All parameters are optional and have default values, except for the sequences, that are compulsory and must be fed as ungapped strings of ACTG characters. 
# The seed is used for reproducible results (if multiple 'arrows' are available from the same position in the scoring matix, only one arrow is chosen, randomly). The parameter 
# num_results allows to set the desired number of alignments to print in output, by default ordered in decreasing value for the final alignment score. Please provide all numerical 
# parameters as positive integers (the code will take care of handling these values in the correct way). Please provide all parameters in the format -parameter=value.

# Importing required packages
import random
import argparse

# Welcome message
print('PLAZZOTTA GIOVANNI, MAT. 232312, ALGORITHMS FOR BIOINFORMATICS (QCB), 2023/2024, UNITN')
print('PERSONAL IMPLEMENTATION OF THE SMITH - WATERMAN PAIRWISE LOCAL ALIGNMENT ALGORITHM')
print('PLEASE USE -H TO VISUALIZE A DESCRIPTION AND DETAILED INSTRUCTIONS')
print(r"""
 ___       ________  ________  ________  ___               ________  ___       ___  ________  ________   _____ ______   _______   ________   _________   
|\  \     |\   __  \|\   ____\|\   __  \|\  \             |\   __  \|\  \     |\  \|\   ____\|\   ___  \|\   _ \  _   \|\  ___ \ |\   ___  \|\___   ___\ 
\ \  \    \ \  \|\  \ \  \___|\ \  \|\  \ \  \            \ \  \|\  \ \  \    \ \  \ \  \___|\ \  \\ \  \ \  \\\__\ \  \ \   __/|\ \  \\ \  \|___ \  \_| 
 \ \  \    \ \  \\\  \ \  \    \ \   __  \ \  \            \ \   __  \ \  \    \ \  \ \  \  __\ \  \\ \  \ \  \\|__| \  \ \  \_|/_\ \  \\ \  \   \ \  \  
  \ \  \____\ \  \\\  \ \  \____\ \  \ \  \ \  \____        \ \  \ \  \ \  \____\ \  \ \  \|\  \ \  \\ \  \ \  \    \ \  \ \  \_|\ \ \  \\ \  \   \ \  \ 
   \ \_______\ \_______\ \_______\ \__\ \__\ \_______\       \ \__\ \__\ \_______\ \__\ \_______\ \__\\ \__\ \__\    \ \__\ \_______\ \__\\ \__\   \ \__\
    \|_______|\|_______|\|_______|\|__|\|__|\|_______|        \|__|\|__|\|_______|\|__|\|_______|\|__| \|__|\|__|     \|__|\|_______|\|__| \|__|    \|__|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
""")

# Argparse setup
my_parser = argparse.ArgumentParser(description= ("""The algorithm relies on the construction of two different matrices: the 'scoring matrix' is populated with 
                                                  scores according to the strategy developed by Smith and Waterman, while the 'origin matrix' is populated with 
                                                  pointers (both in the form of explicit coordinates and with an arrow symbol that can be horizontal '-', vertical
                                                 '|' or diagonal '\\') that point to the nearby cell in the scoring matrix that allowed to achieve the (highest) 
                                                  score in the current cell of the scoring matrix. In the current implementation, the user can provide the two 
                                                  sequences to locally align, and a series of parameters (seed, match_score, mismatch_score, gap_opening_penalty, 
                                                  gap_extension_penalty, cost_model, num_results, verbose). All parameters are optional and have default values, 
                                                  except for the sequences, that are compulsory and must be fed as ungapped strings of ACTG characters. The seed 
                                                  is used for reproducible results (if multiple 'arrows' are available from the same position in the scoring matix, 
                                                  only one arrow is chosen, randomly). The parameter num_results allows to set the desired number of alignments 
                                                  to print in output, ordered in decreasing value for the final alignment score. Please provide all numerical 
                                                  parameters as positive integers (the code will take care of using these values in the correct way). Please 
                                                  provide all parameters in the format -parameter=value."""))
my_parser.add_argument("sequence1", type = str, help = "The first sequence to be aligned (ungapped string of ACTG characters). Upper of lower case.")           # Horizontal sequence
my_parser.add_argument("sequence2", type = str, help = "The first sequence to be aligned (ungapped string of ACTG characters). Upper or lower case.")           # Vertical sequence
my_parser.add_argument("-seed", type = int, default = 1, help = "Seed for reproducible results (please assign an integer value). Default is 1.")
my_parser.add_argument("-match_score", type = int, default = 3, help = "Score to assign for a match (please assign a positive value). Default is 3.")
my_parser.add_argument("-mismatch_score", type = int, default = 3, help = "Penalty to assign for a mismatch (please assign a positive value). Default is 3.")
my_parser.add_argument("-gap_opening_penalty", type = int, default = 2, help = "Penalty to assign every time a gapping section is opened (please assign a positive value). Default is 2.")
my_parser.add_argument("-gap_extension_penalty", type = int, default = 2, help = "Penalty to assign every time an exsisting gapping section is elongated (please assign a positive value). Default is 2.")
my_parser.add_argument("-cost_model", type = str, default = 'linear', help = "Cost model for populating the scoring matrix: can be 'linear' or 'affine'. Default is 'linear'.")
my_parser.add_argument("-num_results", type = int, default = 1, help = "Number of results (alignments) to display, in decreasing order of alignment score. Default is 1.")
my_parser.add_argument("-verbose", action = "store_true", help = "Print the scoring matrix and the origin matrix.")
args = my_parser.parse_args()

# Parameters and inputs: seq1, seq2, seed, match score, mismatch score, gap_opening_penalty and gap_extension_penalty parameters are considered. The latter can be combined using two 
# different cost models: a liner cost model and an affine cost model. The algorithm is able to output not only the best alignment in terms of alignment score, but also the second best 
# and the third best and so on, depending on the value of the num_results parameter. If verbose is selected, the scoring matrix and the origin matrix are printed. 
for character in args.sequence1:
    if character not in ['A', 'C', 'T', 'G', 'a', 'c', 't', 'g']:              # Checking for correct nature of the input sequences (must be ungapped DNA strings)
        raise Exception("Invalid character in sequence 1")
    
for character in args.sequence2:
    if character not in ['A', 'C', 'T', 'G', 'a', 'c', 't', 'g']:
        raise Exception("Invalid character in sequence 2")

seq1 = args.sequence1.upper()                                                  # Reassigning parameters for interal use
seq2 = args.sequence2.upper()
chosen_seed = args.seed
match_score = args.match_score
mismatch_score = args.mismatch_score
gap_opening_penalty = args.gap_opening_penalty                              
gap_extension_penalty = args.gap_extension_penalty
cost_model = args.cost_model                                         
num_results = args.num_results
verbose = args.verbose

# Parameters and inputs check
print('\n')
print('Input sequences provided:')
print(seq1)
print(seq2)
print('\n')
print('Parameters selected (if not provided, default values are used):')
print('Seed:', chosen_seed)
print('Match score:', match_score)
print('Mismatch score:', mismatch_score)
print('Gap opening penalty:', gap_opening_penalty)
print('Gap extension penalty:', gap_extension_penalty)
print('Cost model:', cost_model)
print('Number of results to output:', num_results)
print('Verbose output:', verbose)
print('\n')

import random

# For debugging from Python
# seq1 = 'ATGGGTCTGAAGATAGTCTTATATTATTTTAATTGTGGGTTTTCTTTATATTATTGTTCTGATAATCGTTTGTTTGAAGGTCAAGTAGATACTATATATA'                                              
# seq2 = 'ATGCGAAGTAACTGGCGTATGGTTGCCGCGTTTGTCAGTCACCACCACCCTTGGCCAAGGGTCGTGTACGGTGCATACCT'
# chosen_seed = 1
# match_score = 3
# mismatch_score = 3
# gap_opening_penalty = 2                             
# gap_extension_penalty = 2
# cost_model = 'linear'                                        
# num_results = 5
# verbose = True

# Setting the seed for reproducible results
random.seed(chosen_seed)

# Scoring matrix and origin matrix initializer
m = len(seq1)
n = len(seq2)
scoring_matrix = [[0]*(m+1) for i in range(n+1)]       # Contains the scores
origin_matrix = [[0]*(m+1) for i in range(n+1)]        # Contains, for each score in scoring matrix, the source of that score (as indexes) and the directionality of the movement 

# Gap penality calculator function: as mentioned above, two cost models are considered: linear and affine
def gamma(i, k):
    if cost_model == 'linear':
        return (i - k)*gap_extension_penalty
    if cost_model == 'affine':
        return (gap_opening_penalty + (i - k - 1)*gap_extension_penalty)

# Scoring matrix and origin matrix builder
for i in range(1, n+1):                                                      # Range starts form 1 to exclude the first row and the first column, since their value is initialized at zero
    for j in range(1, m+1):
        if seq1[j-1] == seq2[i-1]:
            s = match_score
        else:
            s = - mismatch_score
        upper_values_list = []
        leftward_values_list = []
        for k in range(0, i):                                                # Last element is always excluded by python, so instead of stopping at i-1 we stop at i
            upper_values_list.append(scoring_matrix[k][j] - gamma(i, k))
        max_upper = max(upper_values_list)
        for k in range(0, j):
            leftward_values_list.append(scoring_matrix[i][k] - gamma(j, k))
        max_leftward = max(leftward_values_list)
        cases = [scoring_matrix[i-1][j-1] + s,                               # The value in each cell of the scoring matrix will be the maximum among these four values
                 max_upper, 
                 max_leftward,
                 0]
        scoring_matrix[i][j] = max(cases)                           # Scoring matrix populator
        origin_matrix[i][j] = []
        if cases[0] == max(cases) and cases[0] >= 0:                # Origin matrix populator: in this implementation, it allows for multiple arrows to be put in the same cell
            origin_matrix[i][j].append((i-1, j-1, "\\"))            # Each cell in the origin matrix contains not only the coordinates of origin, but also an 'arrow' representing the        
        if cases[1] == max(cases) and cases[1] >= 0:                # movement: a diagonal arrow '\\' suggest that the maximum value in the current cell was achieved from the cell in the
            origin_matrix[i][j].append((i-1, j, "|"))               # upper left corner; an horizontal arrow '-' suggests the same but form the right, and a vertical arrow '|' suggests
        if cases[2] == max(cases) and cases[2] >= 0:                # the same but from above.
            origin_matrix[i][j].append((i, j-1, "-"))               
        if cases[0] < 0 and cases[1] < 0 and cases[2] < 0:          # If the first three options are all negative, the algorithm puts 0 in the current cell of the scoring matrix;
            origin_matrix[i][j] = '#'                               # if this is the case however, no arrow is stored in the corresponding cell of the origin matrix. Given that
                                                                    # these zeroes have no origin in terms of directions, they are stored as '#' in the origin matrix.

# If verbose is True, the script prints the scoring matrix and the origin matrix
if verbose:
    print('Scoring matrix:') 
    for row in scoring_matrix:
        print(row)
    print('\n')
    print('Origin matrix:')
    for row in origin_matrix:
        print(row)
    print('\n')

############################################################## ORIGINAL VERSION ######################################################################################
# Backtracking: the code below flattens and sorts the scoring matrix in order for the cells carrying the highest score to be in the first positions. Then, for each cell in the sorted 
# and flattened scoring matrix, if the cell contains at least one arrow and if the cell has not been visited or used as a startig position by a previous backtracking instance (this is to
# avoid generating subalignments of already produced alignments), a backtracking instance is initiated. This procedure is repeated until the script has generated a number num_results of 
# alignments. The backtracking is performed by following iteratively the directions encoded in the origin matrix, until a cell with no arrows is reached. The symbols '-', '|' and '\' allow
# to quickly reconstruct the alignment: we will append a base for both sequences if the arrow is diagonal '\', while if the arrow is horitzontal '-' or vertical '|' we will append a base only
# in one sequence, and we will append a gap in the other. All cells visited during the first backtracking instance are added to the visited_positions list, in order for the subsequent backtracking
# instances (needed if the user asks for more than one result) not to start from these cells. In this implementation, it is possible for a certain cell in the origin matrix to carry more than 
# one arrow: if this is the case, when backtracking one arrow is chosen randomly. The 'alignments' object is a list of num_results lists: each sublist contains the aligned first sequence and 
# the aligned second sequence. 
flattened_scoring_matrix = [(scoring_matrix[i][j], i, j) for i in range(len(scoring_matrix)) for j in range(len(scoring_matrix[0]))]      # touple(alignment score, row_index, column_index)
flattened_scoring_matrix.sort()
inverted_sorted_flattened_scoring_matrix = flattened_scoring_matrix[::-1]
visited_positions = []                                                                          # Positions in the matrix that have already been used
alignments = []
for element in inverted_sorted_flattened_scoring_matrix:
    # We continue only if the current cell has not been involved in an alignment yet, if there is at least an arrow in the current cell and if the number of results generated so far 
    # is smaller than the required num_results.
    if [element[1], element[2]] not in visited_positions and origin_matrix[element[1]][element[2]] not in [0, '#'] and len(alignments) < num_results:   # we repeat this process only num_results times
        alignment_score = element[0]                                                                                  
        starting_position = [element[1], element[2]]        
        empty = False                                                                           # If there is at least one arrow in the current cell of the origin matrix, empty is set to False                                                                                                           
        current_coordinates = starting_position                                          
        seq1_aligned = []
        seq2_aligned = []
        while empty == False:
            visited_positions.append(current_coordinates) 
            num_arrows = len(origin_matrix[current_coordinates[0]][current_coordinates[1]])     # Records the number of arrows available from the current position in the origin matrix
            random_index = random.randint(0, num_arrows-1)                                      # Chooses randomly one direction (arrow) among the ones available
            chosen_arrow_directionality = origin_matrix[current_coordinates[0]][current_coordinates[1]][random_index][2]
            if chosen_arrow_directionality == '\\':
                seq1_aligned.append(seq1[current_coordinates[1]-1])    # Depending on the directionality of the arrow available from the current position, we append only bases or also gaps
                seq2_aligned.append(seq2[current_coordinates[0]-1])
            if chosen_arrow_directionality == '|':
                seq1_aligned.append('-')
                seq2_aligned.append(seq2[current_coordinates[0]-1])
            if chosen_arrow_directionality == '-':
                seq1_aligned.append(seq1[current_coordinates[1]-1])
                seq2_aligned.append('-')
            next_coordinates = list(origin_matrix[current_coordinates[0]][current_coordinates[1]][random_index][0:2])     
            if origin_matrix[next_coordinates[0]][next_coordinates[1]] in [0, '#']:     # If the origin matrix contains no arrows in the next position, we preventively set empty = False
                empty = True
                visited_positions.append(next_coordinates)
            else:
                current_coordinates = next_coordinates
        seq1_aligned.reverse()
        seq2_aligned.reverse()
        alignments.append([seq1_aligned, seq2_aligned, alignment_score, starting_position])

# If the user asks for a number of results higher than the number of alignments in the 'alignments' list, a warning message is printed, and all the avilable alignments 
# in the 'alignments' list are printed.
if num_results > len(alignments):
    print('Warning: the desired number of results is higher than the number of cells in the scoring matrix from which an alignment instance can be started. Showing', len(alignments), 'results.')
    print('\n')

# Alignment packer: for each alignment in the 'alignments' list, the alignment packer builds a list containing the aligned sequences (as strings), the 'alignment bars' (used to print 
# the alignment in a more visually appealing way), the total alignment score, the number of matches, mismatches, the gaps, the length of the alignment and its starting position in the
# scoring matrix. The list is then appended to the alignment_pack list. Remember that by default the aligments in the 'alignments' list are already sorted by the alignment score, by
# construction. 
alignment_pack = []
for alignment in alignments:
    first_seq = ''
    connection_bars = ''
    second_seq = ''
    num_matches = 0
    tot_gap_num = 0
    num_mismatches = 0
    alignment_length = len(alignment[0])                             # The two aligned sequences have the same length, so we just pick one
    alignment_score = alignment[2]
    position = alignment[3]
    for i in range(len(alignment[0])):                               # The two aligned sequences have the same length, so we just pick one
        first_seq = first_seq + alignment[0][i]
        second_seq = second_seq + alignment[1][i]
        if alignment[0][i] == alignment[1][i]:
            connection_bars = connection_bars + '|'
            num_matches = num_matches + 1
        else:
            connection_bars = connection_bars + ' '
        if alignment[0][i] == '-' or alignment[1][i] == '-':
            tot_gap_num = tot_gap_num + 1
        if alignment[0][i] != alignment[1][i] and (alignment[0][i] != '-' and alignment[1][i] != '-'):
            num_mismatches = num_mismatches + 1
    alignment_pack.append([first_seq, second_seq, connection_bars, alignment_score, alignment_length, num_matches, num_mismatches, tot_gap_num, position])

# Results printer: prints the alignments contained in alignment_pack in a visually appealing way, providing also statistics on each alignment such as the alignment score, 
# length, number of matches, of gaps and of mismatches. Remember that by default the aligments in the 'alignments' list are already sorted by the alignment score, by construction. 
print('RESULTS')
print('\n')
counter = 1
for entry in alignment_pack:
    print('########################################')
    print('ALIGNMENT NUMBER', counter)
    print('Starting position on the scoring matrix:', entry[8])
    counter = counter + 1
    print('\n')
    print(entry[0])
    print(entry[2])
    print(entry[1])
    print('\n')
    print('Alignment score:', entry[3])     
    print('Alignment length:', entry[4])                                                                
    print('Number of matches:', entry[5])
    print('Number of mismatches:', entry[6])
    print('Total number of gaps:', entry[7])
    print('########################################')
    print('\n')

# ############################################################### ALTERNATIVE VERSION ######################################################################################
# # In this alternative (more versatile but computationally expensive) version, instead of producing and giving in output only the num_results best alignments ordered by score, 
# # the code runs the backtracking procedure from all cells in the scoring matrix, always respecting the order instructed by the score (so we always start the very first backtracking
# # from the cell having the highest score, and then we perform the backtracking from the cell carrying the second highest score, if this has not been visited yet and if this contains 
# # at least one arrow, and so on...). Each time, the backtracking procedure is performed exactly as described in the first version of the code. Generating all of these alignments allows
# # the user to sort them on the basis of multiple statistics, that may be different from the simple alignment score: for example, the user may be interested in getting in output a certain
# # number of alignments, but sorted according to the length of the alignment, or according to the number of gaps, and not according to the aligment score. To recover the alignment with the 
# # greatest length or with the highest number of gaps, we cannot limit our search to the num_results aligment having the highest alignment score: we need to look through the entire matrix. 

# # Backtracking 
# flattened_scoring_matrix = [(scoring_matrix[i][j], i, j) for i in range(len(scoring_matrix)) for j in range(len(scoring_matrix[0]))]      # touple(alignment score, row_index, column_index)
# flattened_scoring_matrix.sort()
# inverted_sorted_flattened_scoring_matrix = flattened_scoring_matrix[::-1]
# visited_positions = []                                                                          # positions in the matrix that have already been used
# alignments = []
# for element in inverted_sorted_flattened_scoring_matrix:
#     if [element[1], element[2]] not in visited_positions and origin_matrix[element[1]][element[2]] not in [0, '#']:   # we continue only if the current cell has not been involved in an alignment
#         alignment_score = element[0]                                                                                  # and only if there is at least an arrow in the current cell.
#         starting_position = [element[1], element[2]]        
#         empty = False                                                                           # If there is at least one arrow in the current cell of the origin matrix, empty is set to False                                                                                                           
#         current_coordinates = starting_position                                          
#         seq1_aligned = []
#         seq2_aligned = []
#         while empty == False:
#             visited_positions.append(current_coordinates) 
#             num_arrows = len(origin_matrix[current_coordinates[0]][current_coordinates[1]])     # Records the number of arrows available from the current position in the origin matrix
#             random_index = random.randint(0, num_arrows-1)                                      # Chooses randomly one direction (arrow) among the ones available
#             chosen_arrow_directionality = origin_matrix[current_coordinates[0]][current_coordinates[1]][random_index][2]
#             if chosen_arrow_directionality == '\\':
#                 seq1_aligned.append(seq1[current_coordinates[1]-1])    # Depending on the directionality of the arrow available from the current position, we append only bases or also gaps
#                 seq2_aligned.append(seq2[current_coordinates[0]-1])
#             if chosen_arrow_directionality == '|':
#                 seq1_aligned.append('-')
#                 seq2_aligned.append(seq2[current_coordinates[0]-1])
#             if chosen_arrow_directionality == '-':
#                 seq1_aligned.append(seq1[current_coordinates[1]-1])
#                 seq2_aligned.append('-')
#             next_coordinates = list(origin_matrix[current_coordinates[0]][current_coordinates[1]][random_index][0:2])     
#             if origin_matrix[next_coordinates[0]][next_coordinates[1]] in [0, '#']:     # If the origin matrix contains no arrows in the next position, we preventively set empty = False
#                 empty = True
#                 visited_positions.append(next_coordinates)
#             else:
#                 current_coordinates = next_coordinates
#         seq1_aligned.reverse()
#         seq2_aligned.reverse()
#         alignments.append([seq1_aligned, seq2_aligned, alignment_score, starting_position])

# # If the user asks for a number of results higher than the number of alignments in the 'alignments' list, a warning message is printed, and all the avilable alignments 
# # in the 'alignments' list are printed.
# if num_results > len(alignments):
#     print('Warning: the desired number of results is higher than the number of cells in the scoring matrix from which an alignment instance can be started. Showing', len(alignments), 'results.')
#     print('\n')

# # Alignment packer: for each alignment in the 'alignments' list, the alignment packer builds a list containing the aligned sequences (as strings), the 'alignment bars' (used to print 
# # the alignment in a more visually appealing way), the total alignment score, the number of matches, mismatches, the gaps, the length of the alignment and its startign position in the
# # scoring matrix. The list is then appended to the alignment_pack list. 
# alignment_pack = []
# for alignment in alignments:
#     first_seq = ''
#     connection_bars = ''
#     second_seq = ''
#     num_matches = 0
#     tot_gap_num = 0
#     num_mismatches = 0
#     alignment_length = len(alignment[0])                             # The two aligned sequences have the same length, so we just pick one
#     alignment_score = alignment[2]
#     position = alignment[3]
#     for i in range(len(alignment[0])):                               # The two aligned sequences have the same length, so we just pick one
#         first_seq = first_seq + alignment[0][i]
#         second_seq = second_seq + alignment[1][i]
#         if alignment[0][i] == alignment[1][i]:
#             connection_bars = connection_bars + '|'
#             num_matches = num_matches + 1
#         else:
#             connection_bars = connection_bars + ' '
#         if alignment[0][i] == '-' or alignment[1][i] == '-':
#             tot_gap_num = tot_gap_num + 1
#         if alignment[0][i] != alignment[1][i] and (alignment[0][i] != '-' and alignment[1][i] != '-'):
#             num_mismatches = num_mismatches + 1
#     alignment_pack.append([first_seq, second_seq, connection_bars, alignment_score, alignment_length, num_matches, num_mismatches, tot_gap_num, position])

# # Alignment sorter: up to now, the alignments stored in alignment_pack are ALREADY sorted for the alignment score, since the first cells to be subjected to the backtracking were the ones
# # carrying the highest scores. The code below can be quickly modified to sort the alignments in alignment_pack according to some variable (for example, length of the alignment, or 
# # number of gaps ecc...), if this is required.
# alignments_sorted_by_score = sorted(alignment_pack, key=lambda x: x[3], reverse=True)  # The index 4 is the index of the statistics 'alignment_length', so if we want to sort for length, we switch 3 with 4 
#                                                                                        # (at the moment aligment_pack is already sorted for the score, so this is redundant, but may be used to sort according to something else)
# desired_number_alignments_sorted_by_score = alignments_sorted_by_score[0:num_results]      

# # Results printer: prints the alignments contained in desired_number_alignments_sorted_by_len in a visually appealing way, providing also statistics on each alignment such as the alignment
# # score, length, number of matches, of gaps and of mismatches.
# print('RESULTS (alternative version)')
# print('\n')
# counter = 1
# for entry in desired_number_alignments_sorted_by_score:
#     print('########################################')
#     print('ALIGNMENT NUMBER', counter)
#     print('Starting position on the scoring matrix:', entry[8])
#     counter = counter + 1
#     print('\n')
#     print(entry[0])
#     print(entry[2])
#     print(entry[1])
#     print('\n')
#     print('Alignment score:', entry[3])     
#     print('Alignment length:', entry[4])                                                                
#     print('Number of matches:', entry[5])
#     print('Number of mismatches:', entry[6])
#     print('Total number of gaps:', entry[7])
#     print('########################################')
#     print('\n')
