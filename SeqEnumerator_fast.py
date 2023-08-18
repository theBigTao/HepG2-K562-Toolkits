'''
A universal genome sequence tailor tool. 
Capable of enumerating all variants combinations within a user-specified window.

Written by Tao Wang, Department of Genetics, School of Medicine, Stanford University, CA, USA.

Contributors: 
    Jay Luo, Johns Hopkins University, Baltimore, MD, USA 
    Shannon White, Department of Genetics, School of Medicine, Stanford University, CA, USA.
'''

import os
import argparse
import gzip
from Bio import SeqIO

import pandas as pd
from tqdm import tqdm
import copy

import time

BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
data_all_report_submit = dict()

CHROM_SEQ_dict = dict()

# alt is a list of lists [[alt sequence, sequence len], [alt sequence, sequence len]...]
input_df = pd.DataFrame(columns=['chrom', 'position', 'ref', 'ref_len', 'alt_lst', 'vcf_row_num'])

def argument_parser():
    '''
    parse the input parameters
    '''

    parser = argparse.ArgumentParser(description="SeqTailor DNA sequence extraction for genomic variants (neighbourhood) in VCF format")
    parser.add_argument("-g", "--genome", help="genome sequence filename (FASTA / FASTA.GZ)")
    parser.add_argument("-c", "--coordinate", choices=['0', '1'], default='0', help="coordinate indexing")
    parser.add_argument("-s", "--strand", choices=['BOTH', 'FORWARD', 'REVERSE'], default='BOTH', help="strand")
    parser.add_argument("-wd", "--window_down", type=int, default=25, help="window size downstream")
    parser.add_argument("-wu", "--window_up", type=int, default=25, help="window size upstream")
    parser.add_argument("-ws", "--window_size", type=int, default=2114, help="total window size")
    parser.add_argument("-wp", "--window_prediction", type=int, default=1000, help="the prediction result window size")
    parser.add_argument("-q", "--seq_type", choices=['BOTH', 'WT', 'MT'], default='BOTH', help="output sequence type")
    parser.add_argument("-i", "--input", help="input filename (VCF)")
    parser.add_argument("-o", "--output", help="output filename. default: suffix (.DNA.fasta)")
    parser.add_argument("-r", "--report", help="report filename. default: suffix (.report.txt)")
    parser.add_argument("-f", "--tsv_file", help="tsv format for each sub-sequence, default: suffix (.tsv)")
    parser.add_argument('--is_debug', default=False, type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
            help="Debug (true) or production (false, default to production)")
    parser.add_argument('--independent', default=True, type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
            help="Enumerate variant combinations (True, default to not enumerate all the combinations but to" 
            "treat each line of input independently.)")
    parser.add_argument('--ref_as_variant', default=False, type=lambda x: (str(x).lower() in ['true', '1', 'yes']),
            help="Use reference sequence as one of the variants in enumeration (false, default to not use)")
    
    return parser

def intialize_chrom_set_and_dataframe(FILE_INPUT):
    '''
    Initialize chrom set and the internal dataframe for the entire input
    dataset. This is to facilitate the forward and backward walking in the
    dataset by lines, which is very inconvenient to do with file I/O.
    
    Parameters:
        FILE_INPUT: the handler to the input file.

    Return values:
        input_df: the dataframe for the input file, after eliminating invalid entries.
    '''

    assert FILE_INPUT != None, "File must have been opened and ready for reading."

    line_num = 0

    vcf_row_num = 0
    try:
        FILE_INPUT.seek(0)
        for eachline in FILE_INPUT:
            vcf_row_num += 1
            if eachline.strip('\n'):
                eachline = eachline.strip('\n')
                field = eachline.split('\t')

                chrom = field[0]
                pos = int(field[1])
                ref = field[3]
                alt = field[4]
                
                each_variant = chrom + '_' + str(pos) + '_' + ref + '_' + alt
                
                # check each line of input, and recored the comment, the action, and the submitted data
                # dict may be faster than set: O(1) vs O(1) average or O(n) worst case.
                # reference: https://www.geeksforgeeks.org/internal-working-of-set-in-python/
                #if each_variant not in data_variant_set:
                #    data_variant_set.add(each_variant)
                if ref == '.' or ref == '-':
                    continue

                elif ref != CHROM_SEQ_dict[chrom][pos: pos+len(ref)].upper():
                    continue

                elif ref == alt:
                    continue

                elif ',' in alt:
                    alt_column = []
                    alt_list = alt.split(',')
                    each_variant_new = ''
                    for each_alt in alt_list:
                        each_variant_temp = chrom + '_' + str(pos) + '_' + ref + '_' + each_alt
                        each_variant_new += each_variant_temp + ','
                        alt_column.append([each_alt, len(each_alt)])
                    each_variant_new = each_variant_new[0:-1]
                    
                    input_df.loc[line_num] = [chrom, pos, ref, len(ref), alt_column, vcf_row_num ] 
                    line_num += 1

                else:
                    input_df.loc[line_num] = [chrom, pos, ref, len(ref), [[alt, len(alt)]], vcf_row_num] 
                    line_num += 1

    except IOError:
        assert False, "Error happens in file I/O."
    
    except Exception:
        pass

    return input_df

def compute_windows(WINDOW_SIZE, focus_len):
    '''
    WINDOW_SIZE: window size in bps
    focus_len: length in bps for reference sequence or alt alleles

    Return
        window_up: window size for all upstream variants, excluding the center upstream variants
        window_down: window size for downstream variants, excluding the center downstream variants
    '''
    # Compute upstream and downstream window size
    window_up = WINDOW_SIZE//2 - focus_len//2
    if WINDOW_SIZE % 2 == 0:
        if focus_len % 2 == 0:
            window_down = window_up
        else:
            window_down = window_up - 1
    else:
        if focus_len % 2 == 0:
            window_down = window_up + 1
        else:
            window_down = window_up
    
    return window_up, window_down
        
def initialize_variants(focus_row_num, focus_alts_lst, focus_ref_len, WINDOW_SIZE):
    '''
    **Goal**: Give this focus variant at the focus position, generate all valid combinations
    
    **Analysis**
    Search upstream variants in [the start of current chrom, focus_row_num): overlap test and range test
    Each sub-sequence can be represented as combinations of variants and coodinates in the reference genome.
    Each sub-sequence could have different lengths at the same step, some may complete earlier than others.
    Each sub-sequence should have its one flag to determine whether it is complete or not.
    At each step, only the incomplete sub-sequence will be considered to increment.
    
    For upstream variants: the overlap test is done between the left-hand position (the start-point) of 
    the current sub-sequence and the right-hand position (the end-point) of the current variant's reference 
    sequence. The overlap test should not be beween the start-point of the subsequence and each variant's end position,
    because the variant biologically replaces the reference sequence.

    For downstream variants: the overlap test is done between the right-hand position (the end-point) of 
    the current sub-sequence and the left-hand position (the start-point) of the current variant's reference 
    sequence. The overlap test should not be beween the end-point of the current subsequence and each variant's 
    start position, because the variant biologically replaces the reference sequence.
    
    **Design**
    Each sub-sequence should have a flag to mark it as complete or incomplete, current length, start-point
    coordinate, end-point coordinate, (variant start-point coordinate, variant sequence, length),
    (reference start-point coordinate, reference end-point coordinate, length).

    '''
    
    upstream_variants = []
    downstream_variants = []
    focus_variant = input_df.loc[focus_row_num]
    focus_position = focus_variant['position']

    identifier = 1
    for focus_alt in focus_alts_lst:
        focus_sequence = focus_alt[0]
        focus_len = focus_alt[1]

        # Assumption: 
        # 1. WINDOW_PREDICTION is smaller than WINDOW_SIZE
        # 2. WINDOW_PREDICTION is greater than the max variant length.
        # 3. WINDOW_PREDICTION and WINDOW_SIZE are both even.
        # window_up, window_down, chrom are shared by all upstream or downstream variants.
        # if (window_up - bp_length) > 0, the variant is incomplete; otherwise, it is complete.
        # the number of combinations can be calculated as len(upstream_variants) or len(downstream_variants).
        window_up, window_down = compute_windows(WINDOW_SIZE, focus_len)
        
        sub_sequence_upstream = focus_sequence[0:focus_len//2]
        sub_sequence_downstream = focus_sequence[focus_len//2:focus_len]
        upstream_variants.append({'upstream': True,
                'bp_length': focus_len//2,
                'window_size': window_up,
                'identifier': identifier,
                'num_sequences': 1,
                'complete': False,
                'sequence': [{'is_ref': False, 'is_center': True,
                    'row_num': focus_row_num,
                    'start': focus_position, 
                    'end': focus_position+focus_len//2 - 1,
                    'length': focus_len//2,
                    'sub_sequence': sub_sequence_upstream}]
                })
        
        downstream_variants.append({'upstream': False, 
                'bp_length': focus_len - focus_len//2,
                'window_size': window_down,
                'identifier': identifier,
                'num_sequences': 1,
                'complete': False,
                'sequence': [{'is_ref': False, 'is_center': True, 
                    'row_num': focus_row_num,
                    'start': focus_position + focus_len//2, 
                    'end': focus_position+focus_ref_len - 1,
                    'length': focus_len - focus_len//2,
                    'sub_sequence': sub_sequence_downstream}] 
                })

        identifier += 1

    return upstream_variants, downstream_variants

def increment_upstream_variant(upstream_variant, left_over, current_row_num, 
                                current_right_bound, increment_sequence, is_ref=False):
    '''
    Given an increment_sequence, update the current upstream variant.
    upstream_variant: the current upstream variant.
    left_over: the left over window size to be determined
    current_row_num: current row in input_df
    current_right_bound: start-point of the current upstream variant in reference genome coordinate
    increment_sequence: the sequence to prepend to the current upstream variant

    Return:
        upstream_variant: the updated upstream variant
    '''

    increment_sequence = increment_sequence[len(increment_sequence) - left_over:]
    
    increment_variant = { 'is_ref': is_ref, 'is_center': False, 'row_num': current_row_num,
        'start': current_right_bound - left_over, 'end': current_right_bound - 1,
        'length': left_over,
        'sub_sequence': increment_sequence}
    
    upstream_variant['sequence'] = [increment_variant] + upstream_variant['sequence']
    upstream_variant['bp_length'] += left_over
    upstream_variant['window_size'] -= left_over
    upstream_variant['num_sequences'] += 1
    upstream_variant['complete'] = True

    return upstream_variant
                
def increment_downstream_variant(downstream_variant, left_over, current_row_num, 
                                current_left_bound, increment_sequence, is_ref=False):
    '''
    Given an increment_sequence, update the current downstream variant.
    downstream_variant: the current downstream variant.
    left_over: the left over window size to be determined
    current_row_num: current row in input_df
    current_left_bound: end-point of the current downstream variant in reference genome coordinate
    increment_sequence: the sequence to prepend to the current downstream variant

    Return:
        downstream_variant: the updated downstream variant
    '''

    increment_sequence = increment_sequence[0:left_over]
    
    increment_variant = { 'is_ref': is_ref, 'is_center': False, 'row_num': current_row_num,
        'start': current_left_bound + 1, 'end': current_left_bound + left_over,
        'length': left_over,
        'sub_sequence': increment_sequence}
    
    downstream_variant['sequence'] = downstream_variant['sequence'] + [increment_variant]
    downstream_variant['bp_length'] += left_over
    downstream_variant['window_size'] -= left_over
    downstream_variant['num_sequences'] += 1
    downstream_variant['complete'] = True

    return downstream_variant
                
def cross_chrom_update_upstream_variants(upstream_variants, chrom_sequence, current_row_num, base=1):
    '''
    when the current chromesome is different than the focus chromesome,
    all incomplete upstream variants shall be completed with reference sequence.
    If the reference is insufficently long, use 'N' to extend it.

    upstream_variants: list of all multiplicative upstream variants
    chrome_sequence: the focus chromesome
    current_row_num: the current row number in input_df
    base: 0-base (0) or 1-base (1)

    Return:
        updated: whether any upstream variants have been updated.
        updated_upstream_variants: the updated upstream variants. Note that even the
        variants that are not updated also get added.
    '''

    updated_upstream_variants = []
    updated = False
    for upstream_variant in upstream_variants:
        if not upstream_variant['complete']: # incomplete
            # window_up: the left_over upstream window size
            window_up = upstream_variant['window_size']
            left_over = window_up
            
            # current_right_bound: the ref start of the current upstream variant. It is also the exclusive
            # right boundary of the new variant.
            current_right_bound = upstream_variant['sequence'][0]['start']
            increment_sequence = ''
            if (current_right_bound - left_over) >= base: # the chromesome has sufficient nucleutides
                increment_sequence = str(chrom_sequence[(current_right_bound - left_over):current_right_bound])
            else:
                # Any nucleutides: 'N'
                any_nucleutides = 'N'*(left_over+1 - current_right_bound)
                increment_sequence = any_nucleutides + str(chrom_sequence[1:current_right_bound])
            
            increment_variant = { 'is_ref': True, 'is_center': False, 'row_num': current_row_num,
                        'start': current_right_bound - left_over, 'end': current_right_bound - 1,
                        'length': left_over,
                        'sub_sequence': increment_sequence}
            
            upstream_variant['sequence'] = [increment_variant] + upstream_variant['sequence']
            upstream_variant['bp_length'] += left_over
            upstream_variant['window_size'] -= left_over
            upstream_variant['num_sequences'] += 1
            upstream_variant['complete'] = True
            updated = True

        updated_upstream_variants.append(upstream_variant)

    return updated, updated_upstream_variants

def cross_chrom_update_downstream_variants(downstream_variants, chrom_sequence, current_row_num):
    '''
    when the current chromesome is different than the focus chromesome,
    all incomplete downstream variants shall be completed with reference sequence.
    If the reference is insufficently long, use 'N' to extend it.

    downstream_variants: list of all multiplicative downstream variants
    window_down: the downstream window size, excluding the downstream center sequence
    chrome_sequence: the focus chromesome
    current_row_num: the current row number in input_df

    Return:
        updated: whether any downstream variants have been updated.
        updated_downstream_variants: the updated downstream variants. Note that even the
        variants that are not updated also get added.
    '''

    updated_downstream_variants = []
    updated = False
    chrom_length = len(chrom_sequence)
    for downstream_variant in downstream_variants:
        if not downstream_variant['complete']: # incomplete
            left_over = downstream_variant['window_size']
            current_left_bound = downstream_variant['sequence'][-1]['end']
            
            increment_sequence = ''
            if (chrom_length - current_left_bound) > left_over: # the chromesome has sufficient nucleutides
                increment_sequence = str(chrom_sequence[current_left_bound+1: (current_left_bound + left_over + 1)])
            else:
                # Any nucleutides: 'N'
                any_nucleutides = 'N'*(current_left_bound + 1 + left_over - chrom_length)
                increment_sequence = str(chrom_sequence[current_left_bound+1:]) + any_nucleutides 
            
            increment_variant = { 'is_ref': True, 'is_center': False, 'row_num': current_row_num,
                        'start': current_left_bound + 1, 'end': current_left_bound + left_over,
                        'length': left_over,
                        'sub_sequence': increment_sequence}
            
            downstream_variant['sequence'] = downstream_variant['sequence'] + [increment_variant]
            downstream_variant['bp_length'] += left_over
            downstream_variant['window_size'] -= left_over
            downstream_variant['num_sequences'] += 1
            downstream_variant['complete'] = True
            updated = True

        updated_downstream_variants.append(downstream_variant)

    return updated, updated_downstream_variants

def enumerate_upstream_variants(focus_position, focus_chrom, current_row_num, upstream_variants, chrom_sequence,
                                    ref_as_variant, base=1):
    '''
    Description:
        enumerate all combinations of upstream variants.

    Parameters:
        focus_position: the start-position of the focus variant
        focus_chrom: the chromesome of the focus variant
        current_row_num: the row number (0-based) of the current variant being considered to multiply
        upstream_variants: the initialized upstream variants, only containing the center upstream variants
        chrom_sequence: the entire sequence of the focus chromesome
        ref_as_variant: whether the reference sequence is considered as one of the multiplicative variants.
        base: 0-base (0) or 1-base (1)

    Return:
        upstream_variants: all the enumerated upstream variants
    '''

    # upstream:
    # for one step, the right-hand boundary is the same for all variants.
    current_right_bound = focus_position
    to_update = False
    
    while current_row_num >= 0:
        current_row = input_df.loc[current_row_num]
        current_chrom = current_row['chrom']

        # because the processing is line-by-line, the current_chrom is the same for all upstream variants
        if current_chrom != focus_chrom:
            # current upstream variant is from a different chromesome other than the focus
            # for all incomplete variants, complete their current sequence
            updated, updated_upstream_variants = cross_chrom_update_upstream_variants(upstream_variants, 
                                                                             chrom_sequence, current_row_num,
                                                                             base)
            if updated:
                upstream_variants = updated_upstream_variants
            
            # now upstream_variants are complete
            break;
        else:
            # prepare all the potential variants in a list
            current_variants_list = current_row['alt_lst']
            current_ref_len = current_row['ref_len']
            if ref_as_variant:
                current_ref = current_row['ref']
                current_variants_list = [[current_ref, current_ref_len]] + current_variants_list

            # increment all incomplete variants with alt alleles 
            updated_upstream_variants = []
            for upstream_variant in upstream_variants:
                left_over = upstream_variant['window_size']
                if left_over > 0:
                    to_update = True
                    current_right_bound = upstream_variant['sequence'][0]['start']

                    current_position = current_row['position']
                    current_ref_end = current_position + current_ref_len - 1

                    # do overlap and range test
                    # if current ref_end overlaps with current_right_bound, all the variants should be completed
                    # as if the current variant is from a different chromesome.
                    if current_ref_end < current_right_bound:   # no overlap
                        # calculate the gap between the right end of current ref and current_right_bound
                        # this is a sequence that will be included before any variants in current row would
                        # be added to the window.
                        # (current_ref_end, current_right_bound) = [current_ref_end+1, current_right_bound)
                        # gap_between_current_ref_end_and_right_bound = current_right_bound - current_ref_end-1
                        gap_sequence = str(chrom_sequence[(current_ref_end+1):current_right_bound])
                        increment_sequence = gap_sequence # the 2nd half of the increment_sequence

                        # the gap between the current row and the upstream_variant is larger than leftover of the window size
                        if left_over <= len(increment_sequence):
                            upstream_variant = increment_upstream_variant(upstream_variant, left_over, current_row_num, 
                                                            current_right_bound, increment_sequence, is_ref=True)

                            updated_upstream_variants.append(upstream_variant)
                            continue
                        else:
                            # add gap_sequence to the current upstream variant
                            increment_variant = { 'is_ref': True, 'is_center': False, 'row_num': current_row_num,
                                'start': current_ref_end+1, 'end': current_right_bound - 1,
                                'length': len(increment_sequence),
                                'sub_sequence': increment_sequence}
                            
                            upstream_variant['sequence'] = [increment_variant] + upstream_variant['sequence']
                            upstream_variant['bp_length'] = upstream_variant['bp_length'] + (current_right_bound - current_ref_end - 1)
                            upstream_variant['window_size'] -= len(increment_sequence)
                            upstream_variant['num_sequences'] += 1
                            upstream_variant['complete'] = False

                            # multiplication between the current variants and upstream variant
                            left_over = upstream_variant['window_size']
                            current_right_bound = current_ref_end + 1

                            for current_variant in current_variants_list:
                                upstream_variant_new_copy = copy.deepcopy(upstream_variant)
                                current_variant_sequence = current_variant[0]
                                increment_sequence = current_variant_sequence
                                current_variant_length = current_variant[1]

                                if left_over <= current_variant_length:
                                    upstream_variant_new_copy = increment_upstream_variant(upstream_variant_new_copy, left_over, current_row_num,
                                            current_right_bound, increment_sequence)

                                    updated_upstream_variants.append(upstream_variant_new_copy)
                                    continue
                                else:
                                    increment_variant = { 'is_ref': False, 'is_center': False, 'row_num': current_row_num,
                                        'start': current_position, 'end': current_right_bound - 1,
                                        'length': current_variant_length,
                                        'sub_sequence': increment_sequence}
                                    
                                    upstream_variant_new_copy['sequence'] = [increment_variant] + upstream_variant_new_copy['sequence']
                                    upstream_variant_new_copy['bp_length'] += current_variant_length
                                    upstream_variant_new_copy['window_size'] -= current_variant_length
                                    upstream_variant_new_copy['num_sequences'] += 1
                                    upstream_variant_new_copy['complete'] = False
                                    updated_upstream_variants.append(upstream_variant_new_copy)
                    else:
                        # overlap, the current variant should be skipped
                        print(f'Warning: the current variant at input line {current_row_num+1} should be skipped.\n')
                        continue
                else:
                    # even if the upstream variant is not updated, it should be kept.
                    updated_upstream_variants.append(upstream_variant)
        
            # update upstream variants with the new multiplicative combinations
            upstream_variants = updated_upstream_variants

        # all incomplete variants have been multiplied with the current row.
        if to_update:
            current_row_num -= 1
        else:
            break # break out while loop for upstream variants
    else:
        # current_row_num is -1, the focus variant row is 0, time to
        # complete all upstream variants
        updated, updated_upstream_variants = cross_chrom_update_upstream_variants(upstream_variants, 
                                                                        chrom_sequence, current_row_num+1, base)
        if updated:
            upstream_variants = updated_upstream_variants

    return upstream_variants

def enumerate_downstream_variants(focus_position, focus_chrom, current_row_num,
                                        downstream_variants, chrom_sequence, ref_as_variant):
    '''
    Description:
        enumerate all combinations of downstream variants.
        Note: the two function enumerate_upstream_variants and enumerate_downstream_variants have a duality relationship.
        up should be down, start is end, + is -, etc. It is very difficult to share the same code copy without confusion,
        so they end up with being two copies, even if they are similar.

    Parameters:
        focus_position: the start-position of the focus variant
        focus_chrom: the chromesome of the focus variant
        current_row_num: the row number (0-based) of the current variant being considered to multiply
        downstream_variants: the initialized downstream variants, only containing the center downstream variants
        chrom_sequence: the entire sequence of the focus chromesome
        ref_as_variant: whether the reference sequence is considered as one of the multiplicative variants.

    Return:
        downstream_variants: all the enumerated downstream variants
    '''
    
    # Search downstream variants in (focus_row_num, the end of current chrom]
    last_row_num = input_df.shape[0]
    while current_row_num < last_row_num:
        current_row = input_df.loc[current_row_num]
        current_chrom = current_row['chrom']
        to_update = False
        if current_chrom != focus_chrom:
            # current downstream variant is from a different chromesome other than the focus
            # for all incomplete variants, complete their current sequence
            updated, updated_downstream_variants = cross_chrom_update_downstream_variants(downstream_variants, 
                                                                            chrom_sequence, current_row_num)
            if updated:
                downstream_variants = updated_downstream_variants
            else:
                # now downstream_variants are complete
                break;
        else:
            # prepare all the potential variants in a list
            current_variants_list = current_row['alt_lst']
            current_ref_len = current_row['ref_len']
            if ref_as_variant:
                current_ref = current_row['ref']
                current_variants_list = [[current_ref, current_ref_len]] + current_variants_list
            
            # increment all incomplete variants with alt alleles 
            updated_downstream_variants = []
            for downstream_variant in downstream_variants:
                left_over = downstream_variant['window_size']
                if left_over > 0:
                    to_update = True
                    current_left_bound = downstream_variant['sequence'][-1]['end']

                    current_position = current_row['position']
                    current_ref_end = current_position + current_ref_len - 1
                    # do overlap and range test
                    # if current position overlaps with current_left_bound, all the variants should be completed
                    # as if the current variant is from a different chromesome.
                    if current_position > current_left_bound:   # no overlap
                        # calculate the gap between the start of current ref and current_left_bound,
                        # this is a sequence that will be included before any variants in current row would
                        # be added to the downstream window.
                        # (current_left_bound, current_position) = [current_left_bound+1, current_position)
                        # gap_between_current_ref_start_and_left_bound = current_position - current_left_bound-1
                        gap_sequence = str(chrom_sequence[(current_left_bound+1):current_position])
                        increment_sequence = gap_sequence # the 1st half of the increment_sequence
                        
                        # the gap between the current row and the downstream_variant is larger than leftover of the window size
                        if left_over <= len(increment_sequence):
                            downstream_variant = increment_downstream_variant(downstream_variant, left_over, current_row_num, 
                                                            current_left_bound, increment_sequence, is_ref=True)

                            updated_downstream_variants.append(downstream_variant)
                            continue
                        else:
                            # add gap_sequence to the current downstream variant
                            increment_variant = { 'is_ref': True, 'is_center': False, 'row_num': current_row_num,
                                'start': current_left_bound+1, 'end': current_position - 1,
                                'length': len(increment_sequence),
                                'sub_sequence': increment_sequence}
                            
                            downstream_variant['sequence'] = downstream_variant['sequence'] + [increment_variant]
                            downstream_variant['bp_length'] = downstream_variant['bp_length'] + len(increment_sequence)
                            downstream_variant['window_size'] -= len(increment_sequence)
                            downstream_variant['num_sequences'] += 1
                            downstream_variant['complete'] = False

                            # multiplication between the current variants and downstream variant
                            left_over = downstream_variant['window_size']
                            current_left_bound = current_position - 1

                            for current_variant in current_variants_list:
                                downstream_variant_new_copy = copy.deepcopy(downstream_variant)
                                current_variant_sequence = current_variant[0]
                                increment_sequence = current_variant_sequence
                                current_variant_length = current_variant[1]

                                # if current variant is longer than left_over
                                if left_over <= current_variant_length:
                                    downstream_variant_new_copy = increment_downstream_variant(downstream_variant_new_copy, left_over, current_row_num,
                                            current_left_bound, increment_sequence)

                                    updated_downstream_variants.append(downstream_variant_new_copy)
                                    continue
                                else:
                                    increment_variant = { 'is_ref': False, 'is_center': False, 'row_num': current_row_num,
                                        'start': current_position, 'end': current_ref_end,
                                        'length': current_variant_length,
                                        'sub_sequence': increment_sequence}
                                    
                                    downstream_variant_new_copy['sequence'] = downstream_variant_new_copy['sequence'] + [increment_variant]
                                    downstream_variant_new_copy['bp_length'] += current_variant_length
                                    downstream_variant_new_copy['window_size'] -= len(increment_sequence)
                                    downstream_variant_new_copy['num_sequences'] += 1
                                    downstream_variant_new_copy['complete'] = False
                                    updated_downstream_variants.append(downstream_variant_new_copy)
                    else:
                        # overlap, the current variant should be skipped
                        print(f'Warning: the current variant at input line {current_row_num+1} should be skipped.\n')
                        continue
                else:
                    # even if the downstream variant is not updated, it should be kept.
                    updated_downstream_variants.append(downstream_variant)

            # update downstream variants with the new multiplicative combinations
            downstream_variants = updated_downstream_variants

        # all incomplete variants have been multiplied with the current row.
        if to_update:
            current_row_num += 1
        else:
            break # break out while loop for downstream variants
    else:
        # current_row_num is -1, the focus variant row is 0, time to
        # complete all downstream variants
        updated, updated_downstream_variants = cross_chrom_update_downstream_variants(downstream_variants, 
                                                                        chrom_sequence, current_row_num+1)
        if updated:
            downstream_variants = updated_downstream_variants
    
    return downstream_variants

def enumerate_upstream_downstream_combinations(upstream_variants, downstream_variants, WINDOW_SIZE):
    '''
    Given a list of upstream variants and a list of downstream variants, enumerate all the combinations
    to generate a list of complete variants.
    
    # Enumerate all the combinations between all upstream and downstream variants
    # Design:
    # 1. an upstream variant can only pair up with a downstream variant that shares the same identfier
    # 2. the combinations could be generated on the fly, but that will makes the code less module-oriented.
    #    It will make code hard to understand and maintain, because the assumptions are implicit and control
    #    flows are mixed together. We decide to modularize the code in a way that all the combinations are
    #    generated first in one function, and pass down to the next function to parse the combinations for
    #    final results.

    upstream_variants: the list of upstream variants
    downstream_variants: the list of downstream variants

    Return:
        complete_variants: the list of all valid complete variants
    '''
    complete_variants = []

    for upstream_variant in upstream_variants:
        for downstream_variant in downstream_variants:
            if upstream_variant['identifier'] == downstream_variant['identifier']:
                #print('a valid combination')

                complete_variant = copy.deepcopy(upstream_variant)
                complete_variant['upstream'] = False
                complete_variant['bp_length'] += downstream_variant['bp_length']
                assert complete_variant['bp_length'] == WINDOW_SIZE, f"bp_length must be {WINDOW_SIZE}"

                complete_variant['num_sequences'] += downstream_variant['num_sequences']
                complete_variant['sequence'] = complete_variant['sequence'] + downstream_variant['sequence']

                complete_variants.append(complete_variant)


    return complete_variants


def reverse_sequence(sequence_input):
    '''
    # reverse the input sequence into complementary strand in reverse order
    sequence_input: a string of the sequence

    Return:
        sequence_reversed: complementary strand in reverse order
    '''

    sequence_length = len(sequence_input)
    sequence_reversed = ''
    for i in range(0, sequence_length):
        sequence_reversed += BASE_PAIRING[sequence_input[sequence_length - i - 1]]
    return sequence_reversed

def write_sequence_to_file(meta_info, sequence, STRAND, FILE_OUTPUT):
    '''
    write out the meta info and the complete sequence for the focus variant
    
    meta_info: meta information about the sequence, chromesome, variant type,
               start coordinate, length and row number in the original input file.
               meta_info is only for forward strand. For reverse strand, it is
               still meaningful but has to be carefully interpreted
    sequence:  the complete sequence for the focus variant
    STRAND: FORWARD/REVERSE/BOTH strand
    FILE_OUTPUT: the file handler for output sequences

    Return:
        None
    '''
                
    ###
    # process output sequences
    # wrap the sequences to 80nt per line
    ###
    if STRAND == 'BOTH' or STRAND == "FORWARD":
        # meta info
        FILE_OUTPUT.write(meta_info + '\n')

        # the forward sequence
        output_seq_forward = sequence
        output_seq_forward_wrap = ''
        for i in range(1, len(output_seq_forward)+1):
            if i % 80 == 0:
                output_seq_forward_wrap += output_seq_forward[i-1] + '\n'
            else:
                output_seq_forward_wrap += output_seq_forward[i-1]

        FILE_OUTPUT.write(output_seq_forward_wrap + '\n\n')
        
    if STRAND == 'BOTH' or STRAND == "REVERSE":
        output_seq_reverse = reverse_sequence(sequence.upper())

        output_seq_reverse_wrap = ''
        for i in range(1, len(output_seq_reverse)+1):
            if i % 80 == 0:
                output_seq_reverse_wrap += output_seq_reverse[i-1] + '\n'
            else:
                output_seq_reverse_wrap += output_seq_reverse[i-1]
        
        FILE_OUTPUT.write(output_seq_reverse_wrap + '\n\n')


def write_sequence_statistics(FILE_OUTPUT, num_upstream_variants, num_downstream_variants, total_num_variants):
    '''
    write the descriptive statistics for the sequence

    '''
    statistics = "# variants: " + str(num_upstream_variants) + " (upstream)\t" + str(num_downstream_variants) \
                    + " (downstream)\t" + str(total_num_variants) + " (total)"

    FILE_OUTPUT.write(statistics + '\n')

def write_complete_variants_to_file(FILE_OUTPUT, complete_variants, WINDOW_SIZE, STRAND, focus_chrom, focus_ref, focus_ref_len):
    '''
    Analysis: the output forms could be either tailored for our purpose or be consistent with original
    SeqTailor. The output could be in two files, one for all valid variants, one for invalid variants
    or inputs. The output should come with both the sequence and some meta data for the sequence.

    Design: the meta data lines begin with >, and are meant to give a quick overview of the sequence
            compisition, which part is composed of which type of variant (SNP/Insertion/Deletion/Reference, etc).

    FILE_OUTPUT: the file handler for output sequences
    complete_variants: the list of valid complete sequences (TODO: may include invalid ones later)
    WINDOW_SIZE: the window size in bp for each variant
    STRAND: FORWARD/REVERSE/BOTH strand. The sequence information in complete_variants is only for
            forward strand. To get reverse strand, apply reverse complement operations.
    focus_chrom: the chromesome of the current variants
    focus_ref: the focus reference sequence
    focus_ref_len: the length in bp for focus ref (the center)

    Return:
        None
    '''

    for complete_variant in complete_variants:
        assert complete_variant['bp_length'] == WINDOW_SIZE
        assert complete_variant['complete']
    
        # the center of a complete_variant is counted twice
        total_num_variants = complete_variant['num_sequences'] - 1
        identifier = complete_variant['identifier']
        
        meta_info = '>' + focus_chrom + '_' + focus_ref + '_' + str(focus_ref_len)
        total_num_variants = 0
        num_upstream_variants = 0
        num_downstream_variants = 0

        # collect the complete sequence out of the sub-sequences
        upstream = True
        complete_sequence = ""
        for sequence in complete_variant['sequence']:
            complete_sequence += sequence['sub_sequence']

            if not sequence['is_center']:
                if upstream:
                    # upstream variant
                    num_upstream_variants += 1
                else:
                    num_downstream_variants += 1
            else:
                upstream = False

            # sequence type
            if sequence['is_ref']:
                meta_info += "_REF"

            else:
                # check the type of the current variant
                current_variant = input_df.loc[(sequence['row_num'])]
                current_ref_len = current_variant['ref_len']

                if sequence['length'] == 1 and current_ref_len == 1:
                    meta_info += "_SNP"
                elif sequence['length'] > current_ref_len:
                    meta_info += "_INS"
                elif sequence['length'] < current_ref_len:
                    meta_info += "_DEL"
                else:
                    # substitution
                    meta_info += "_SUB"
                
                if sequence['is_center']:
                    if sequence['sub_sequence'] == '':
                        meta_info += "_center_" + "\'\'"
                    else:
                        meta_info += "_center_" + sequence['sub_sequence']
                else:
                    meta_info += '_' + sequence['sub_sequence']

            # row number, start and length
            meta_info += "_" + str(sequence['row_num']+1) # row number in 1-based input file
            meta_info += "_" + str(sequence['start']) 
            meta_info += "_" + str(sequence['length'])

        # write out the meta info and the complete sequence
        write_sequence_statistics(FILE_OUTPUT, num_upstream_variants, num_downstream_variants, total_num_variants)
        write_sequence_to_file(meta_info, complete_sequence, STRAND, FILE_OUTPUT)

def write_complete_variants_to_tsv(TSV_FILE_OUTPUT, complete_variants, 
                                    WINDOW_SIZE, WINDOW_PREDICTION, STRAND, 
                                    focus_chrom, focus_vcf_row_num, focus_ref, focus_ref_len, Window_reference):
    '''
    Analysis: the tsv file is a better way to represent the information because it is easier to process with
              pandas, e.g., filtering according to specific conditions.

    Design: the header is self-explanatory and specifies each column. 
            Variant_type is one of SNP/Insertion/Deletion/Reference;
            Multiple upstream (downstream) neigbors are separated by ","
            each variant is coded as varaint_sequence + "_" + start_position

    TSV_FILE_OUTPUT: the tsv file handler for output sequences
    complete_variants: the list of valid complete sequences (TODO: may include invalid ones later)
    WINDOW_SIZE: the window size in bp for each variant
    WINDOW_PREDICTION: the prediction result window size in bps
    STRAND: FORWARD/REVERSE/BOTH strand. The sequence information in complete_variants is only for
            forward strand. To get reverse strand, apply reverse complement operations.
    focus_chrom: the chromesome of the current variants
    focus_vcf_row_num: row num of the current variant in vcf file, used to identify corresponding PPM info.
    focus_ref: the focus reference sequence
    focus_ref_len: the length in bp for focus ref (the center)
    ref_sequence: the WINDOW_SIZE long reference sequence (placing the focus reference at the center)

    Return:
        None
    '''
    Row_num_tsv = -1
    Position_tsv = -1
    
    Chrom_tsv = focus_chrom
    Ref_tsv = focus_ref
    vcf_row_num_tsv = focus_vcf_row_num

    # the prediction window is [flank_left_bound, flank_right_bound), length: WINDOW_PREDICTION
    flank_size = (WINDOW_SIZE - WINDOW_PREDICTION) // 2
    flank_left_bound = flank_size # [0, flank_left_bound), length: flank_size
    flank_right_bound = WINDOW_SIZE - flank_size # [flank_right_bound, WINDOW_SIZE), length: flank_size
    Window_reference_prediction = Window_reference[flank_left_bound:flank_right_bound]

    for complete_variant in complete_variants:
        assert complete_variant['bp_length'] == WINDOW_SIZE
        assert complete_variant['complete']

        # focus variant-specific meta data
        Window_sequence_tsv = ''
        Focus_variant_tsv = ''
        Variant_type_tsv = ''
        Upstream_neighbors_tsv_lst = []
        Downstream_neighbors_tsv_lst = []

        Window_sequence_prediction_tsv = ''
        Upstream_neighbors_prediction_tsv_lst = []
        Downstream_neighbors_prediction_tsv_lst = []

        first_left = True
        upstream = True
        for sequence in complete_variant['sequence']:
            
            Window_sequence_tsv += sequence['sub_sequence']
            current_seq_length = len(Window_sequence_tsv)
            
            if sequence['is_center']:
                Focus_variant_tsv += sequence['sub_sequence']
                
                if first_left:
                    Row_num_tsv = sequence['row_num'] + 1
                    Position_tsv = sequence['start']
                    first_left = False
                    upstream = False
            else:
                if not sequence['is_ref']:
                    if upstream:
                        Upstream_neighbors_tsv_lst.append(sequence['sub_sequence'] + '_' + str(sequence['start']))
                        if current_seq_length > flank_size:
                            Upstream_neighbors_prediction_tsv_lst.append(sequence['sub_sequence'] + '_' + str(sequence['start']))
                    else:
                        Downstream_neighbors_tsv_lst.append(sequence['sub_sequence'] + '_' + str(sequence['row_num'] + 1) + '_' + str(sequence['start']))
                        if current_seq_length <= flank_right_bound:
                            Downstream_neighbors_prediction_tsv_lst.append(sequence['sub_sequence'] + '_' + str(sequence['row_num'] + 1) + '_' + str(sequence['start']))

        
        Upstream_neighbors_tsv = ",".join(Upstream_neighbors_tsv_lst)
        Downstream_neighbors_tsv = ",".join(Downstream_neighbors_tsv_lst)

        # check the type of the focus variant
        focus_variant_length = len(Focus_variant_tsv)
        if focus_variant_length == 1 and focus_ref_len == 1:
            Variant_type_tsv = "SNP"
        elif focus_variant_length > focus_ref_len:
            Variant_type_tsv = "INS"
        elif focus_variant_length < focus_ref_len:
            Variant_type_tsv = "DEL"
        else:
            # substitution
            Variant_type_tsv = "SUB"

        Window_sequence_prediction_tsv = Window_sequence_tsv[flank_left_bound:flank_right_bound]
        Upstream_neighbors_prediction_tsv = ",".join(Upstream_neighbors_prediction_tsv_lst)
        Downstream_neighbors_prediction_tsv = ",".join(Downstream_neighbors_prediction_tsv_lst)

        # write out one complete variant
        TSV_FILE_OUTPUT.write(f'{Row_num_tsv}\t{Chrom_tsv}\t{Position_tsv}\t{Ref_tsv}\t{Focus_variant_tsv}'
                                f'\t{Variant_type_tsv}\t{Window_reference}\t{Window_sequence_tsv}'
                                f'\t{Window_reference_prediction}\t{Window_sequence_prediction_tsv}'
                                f'\t{Upstream_neighbors_tsv}\t{Downstream_neighbors_tsv}'
                                f'\t{Upstream_neighbors_prediction_tsv}\t{Downstream_neighbors_prediction_tsv}'
                                f'\t{vcf_row_num_tsv}\n')

def get_ref_sequence_for_upstream(focus_position, focus_ref, chrom_sequence, WINDOW_SIZE):
    '''
    Get the pure reference sequence in the upstream window.

    focus_position: the focus ref start coordinate (1-based)
    focus_ref: the focus reference sequence
    chrom_sequence: the focus chromesome
    WINDOW_SIZE: window size in bps 
    Return:
        ref_sequence_upstream: the ref sequence
    '''
    ref_sequence_upstream = ""
    window_up, window_down = compute_windows(WINDOW_SIZE, len(focus_ref))

    left_over = window_up
    if left_over > 0:
        current_right_bound = focus_position
        if current_right_bound > left_over:
            ref_sequence_upstream = str(chrom_sequence[current_right_bound - left_over:current_right_bound])
        else:
            ref_sequence_upstream = 'N'*(left_over-current_right_bound+1) + str(chrom_sequence[1:current_right_bound])
    
    assert len(ref_sequence_upstream) == left_over, "ref_sequence_upstream should have left_over bps."
    
    return ref_sequence_upstream
    
def get_ref_sequence_for_downstream(focus_position, focus_ref, chrom_sequence, WINDOW_SIZE):
    '''
    Get the pure reference sequence in the downstream window.

    focus_position: the focus ref start coordinate (1-based)
    focus_ref: the focus reference sequence
    chrom_sequence: the focus chromesome
    WINDOW_SIZE: window size in bps 
    Return:
        ref_sequence_downstream: the ref sequence
    '''
    ref_sequence_downstream = ""
    chrom_length = len(chrom_sequence)
    window_up, window_down = compute_windows(WINDOW_SIZE, len(focus_ref))

    left_over = window_down
    if left_over > 0:
        current_left_bound = focus_position + len(focus_ref) - 1
        if (chrom_length - current_left_bound) > left_over:
            ref_sequence_downstream = str(chrom_sequence[current_left_bound+1:current_left_bound + left_over + 1])
        else:
            ref_sequence_downstream = str(chrom_sequence[current_left_bound+1:chrom_length]) \
                                        + 'N'*(left_over+current_left_bound+1-chrom_length) 
    
    assert len(ref_sequence_downstream) == left_over, "ref_sequence_downstream should have left_over bps."
    
    return ref_sequence_downstream


def main():
    parser = argument_parser()

    # paring the arguments
    args = parser.parse_args()
    GENOME = args.genome
    COORDINATE = args.coordinate
    STRAND = args.strand
    WINDOW_DOWN = args.window_down
    WINDOW_UP = args.window_up
    WINDOW_SIZE = args.window_size
    WINDOW_PREDICTION = args.window_prediction
    SEQUENCE_TYPE = args.seq_type
    INPUT = args.input
    OUTPUT = args.output
    REPORT = args.report
    TSV_OUTPUT = args.tsv_file
    # this flag controls whether to consider reference as one variant.
    ref_as_variant = args.ref_as_variant
    independent = args.independent

    # generate the filenames for output and report, when these are not provided by the user
    if not OUTPUT:
        OUTPUT = INPUT + '.DNA.fa'
    if not REPORT:
        REPORT = INPUT + '.report.txt'

    if not TSV_OUTPUT:
        TSV_OUTPUT = INPUT + '.tsv'
    
    # when the window size (downstream or upstream) is negative, skip it
    if WINDOW_DOWN < 0 or WINDOW_UP < 0:
        FILE_OUTPUT = open(OUTPUT, 'w')
        TSV_FILE_OUTPUT = open(TSV, 'w')
        FILE_REPORT = open(REPORT, 'w')
        FILE_OUTPUT.write('NEGATIVE WINDOW SIZE\n')
        FILE_REPORT.write('NEGATIVE WINDOW SIZE\n')
        TSV_FILE_OUTPUT.write('NEGATIVE WINDOW SIZE\n')

        return

    else:
        # when any provied window size is greater than 5,000bp, trim it to 5,000bp
        if WINDOW_DOWN > 5000:
            WINDOW_DOWN = 5000
        if WINDOW_UP > 5000:
            WINDOW_UP = 5000
        
        # provide all defined parameters and filenames in the output file
        SeqTailor_PARAMETER = '# DNA sequence extraction for genomic variants (neighbourhood) in VCF format\n' \
                            + '# GENOME: ' + GENOME + '\n' \
                            + '# COORDINATE: ' + COORDINATE + '\n' \
                            + '# STRAND: ' + STRAND + '\n' \
                            + '# WINDOW_DOWN: ' + str(WINDOW_DOWN) + ' bp\n' \
                            + '# WINDOW_UP: ' + str(WINDOW_UP) + ' bp\n' \
                            + '# SEQUENCE_TYPE: ' + SEQUENCE_TYPE + '\n' \
                            + '# INPUT FILE: ' + INPUT + '\n' \
                            + '# OUTPUT FILE: ' + OUTPUT + '\n' \
                            + '# REPORT FILE: ' + REPORT + '\n\n'
        
        SeqTailor_TSV_HEADERS = "Row_num\tChrom\tPosition\tFocus_ref\tFocus_varaint\tVariant_type" + \
                "\tWindow_reference\tWindow_sequence\tWindow_reference_prediction" + \
                "\tWindow_sequence_prediction\tUpstream_neighors\tDownstream_neighbors" + \
                "\tUpstream_neighors_prediction\tDownstream_neighbors_prediction\tvcf_row_num_tsv\n"

        try:
            ###
            # read input files, generate output files, identify chromosomes, and load sequences
            ###
            if GENOME.endswith('fasta') or GENOME.endswith('fa'):
                FILE_GENOME = SeqIO.parse(open(GENOME), 'fasta')
            elif GENOME.endswith('gz'):
                GENOME_gz = gzip.open(GENOME, 'r')
                FILE_GENOME = SeqIO.parse(GENOME_gz, 'fasta')

            if INPUT.endswith('vcf'):
                FILE_INPUT = open(INPUT, 'r')
            elif INPUT.endswith('vcf.gz'):
                FILE_INPUT = gzip.open(INPUT, 'r')
            
            FILE_OUTPUT = open(OUTPUT, 'w')
            FILE_OUTPUT.write(SeqTailor_PARAMETER)
            FILE_REPORT = open(REPORT, 'w')
            FILE_REPORT.write(SeqTailor_PARAMETER)
            TSV_FILE_OUTPUT = open(TSV_OUTPUT, 'w')
            TSV_FILE_OUTPUT.write(SeqTailor_TSV_HEADERS)

            # initial scan of the chromosomes of the given variants
            CHROM_set = set()
            for eachline in FILE_INPUT:
                if not eachline.startswith('#') and eachline.strip('\n'):
                    field = eachline.strip('\n').split('\t')
                    chrom = field[0]
                    CHROM_set.add(chrom)

            # read the genome sequence, and load the sequences of the given chromosomes
            for each_genome in FILE_GENOME:
                chromosome, sequence = each_genome.id, each_genome.seq
                if chromosome in CHROM_set:
                    if COORDINATE == '0':
                        CHROM_SEQ_dict[chromosome] = sequence
                    elif COORDINATE == '1':
                        CHROM_SEQ_dict[chromosome] = 'N' + sequence
            
            # initial scan of the chromosomes of the given variants
            start = time.time()
            intialize_chrom_set_and_dataframe(FILE_INPUT)
            done = time.time()
            elapsed = done - start
            print(f'intialize_chrom_set_and_dataframe: {elapsed} seconds')
        
        # if any IOerror, give an error message
        except IOError:
            FILE_OUTPUT = open(OUTPUT, 'w')
            FILE_REPORT = open(REPORT, 'w')
            TSV_FILE_OUTPUT = open(TSV_OUTPUT, 'w')
            FILE_OUTPUT.write('FILE OPEN ERROR\n')
            TSV_FILE_OUTPUT.write('FILE OPEN ERROR\n')
            FILE_REPORT.write('FILE OPEN ERROR\n')
        
        # for each row in input_df, do a multiplicative enumeration for its variants and neighboring variants.
        # neighboring variants are those which are within WINDOW_SIZE window. 
        # The center of each focus variant in the focus row is the center of the window.
        last_row_num = input_df.shape[0]
        rows_processed = 0
        with tqdm(total=last_row_num) as pbar:
            for focus_row_num in range(0, last_row_num):
                # Focus the variant (focus)
                focus_variant = input_df.loc[focus_row_num]
                focus_chrom = focus_variant['chrom']
                focus_position = focus_variant['position']
                focus_ref = focus_variant['ref']
                focus_ref_len = focus_variant['ref_len']
                focus_alts_lst = focus_variant['alt_lst']
                focus_vcf_row_num = focus_variant['vcf_row_num']
                chrom_sequence = CHROM_SEQ_dict[focus_chrom]
    
                # whether to use reference genome as one of the variants in the enumeration
                if ref_as_variant:
                    focus_alts_lst =  [[focus_ref, focus_ref_len]] + focus_alts_lst
    
                # find the center of the focus variants and initialize upstream variants and downstream variants
                # according to the focus variants.
                upstream_variants, downstream_variants = initialize_variants(focus_row_num, focus_alts_lst, focus_ref_len, WINDOW_SIZE)

                # ref_sequence: the WINDOW_SIZE long reference sequence (placing the focus reference at the center)
                ref_sequence_upstream = get_ref_sequence_for_upstream(focus_position, focus_ref, chrom_sequence, WINDOW_SIZE) 
                ref_sequence_downstream = get_ref_sequence_for_downstream(focus_position, focus_ref, chrom_sequence, WINDOW_SIZE) 
                ref_sequence = ref_sequence_upstream + focus_ref + ref_sequence_downstream
    
                # from now on, upstream and downstream variants are calculated separately.
                # enumerate all upstream variants
                if independent:
                    # fake the current row is the one row before the first row, so that enumerate_upstream_variants
                    # will just complete the all upstream variants with ref sequence.
                    current_row_num = -1
                else:
                    current_row_num = focus_row_num - 1
                upstream_variants = enumerate_upstream_variants(focus_position, focus_chrom, current_row_num, 
                                                                upstream_variants, chrom_sequence, ref_as_variant,
                                                                int(COORDINATE))

                # enumerate all downstream variants
                if independent:
                    # fake the current row is the last row, so that
                    # enumerate_downstream_variants will just complete the all downstream variants with ref sequence. 
                    current_row_num = input_df.shape[0]
                else:
                    current_row_num = focus_row_num + 1
                downstream_variants = enumerate_downstream_variants(focus_position, focus_chrom, current_row_num, 
                                                                downstream_variants, chrom_sequence, ref_as_variant)
                
                # Generate the combinations
                complete_variants = enumerate_upstream_downstream_combinations(upstream_variants, downstream_variants, WINDOW_SIZE)

                # Generate the report like SeqTailor
                write_complete_variants_to_file(FILE_OUTPUT, complete_variants, WINDOW_SIZE, STRAND, focus_chrom, focus_ref, focus_ref_len)

                # Write out the csv file
                write_complete_variants_to_tsv(TSV_FILE_OUTPUT, complete_variants, WINDOW_SIZE, WINDOW_PREDICTION,
                                                STRAND, focus_chrom, focus_vcf_row_num, focus_ref, focus_ref_len, ref_sequence)

                # progress report
                rows_processed += 1
                if rows_processed % 100 == 0:
                    pbar.update(100)


        FILE_OUTPUT.close()
        TSV_FILE_OUTPUT.close()

if __name__ == '__main__':
    main()
