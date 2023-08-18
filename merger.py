'''
A tool to merge overlapping variants into one.
- powerful than bedtools merge in overlapping variant merging
- merge all adjacent overlapping variants no matter how different they are into an uniform format

Written by Tao Wang, Department of Genetics, School of Medicine, Stanford University, CA, USA.

Contributors:
    Shannon White, Department of Genetics, School of Medicine, Stanford University, CA, USA.
'''

import argparse
import os
from Bio import SeqIO

import pandas as pd
from tqdm import tqdm
import random

BASE_PAIRING = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
CHROM_SEQ_dict = dict()

columns=['chrom', 'start_pos', 'end_pos', 'annotation', 'ref', 'alt_lst']

def non_negative_int(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is not a non-negative integer" % value)
    return ivalue

def argument_parser():
    parser = argparse.ArgumentParser(description="Merge adjacent overlapping variants (.bed, 0-based) into one variant (.bed, 0-based).")
    parser.add_argument("-g", "--genome", help="genome sequence filename (FASTA / FASTA.GZ)")
    parser.add_argument("-i", "--input", help="input filename (.bed)")
    parser.add_argument("-o", "--output", help="output filename (.bed)")
    parser.add_argument("-v", "--variant", action='store_true', default=False, 
                        help="target individual variant or uniform reference sequene")
    parser.add_argument("-r", "--random", action='store_true', default=False, help="randomize alternative sequence instead of simple deletion")
    parser.add_argument("-s", "--scramble", action='store_true', default=False, 
                        help="scramble alternative sequence instead of randomized w/ [A, T, C, G] for each base")
    parser.add_argument("-d", "--deletion_only", action='store_true', default=False, 
                        help="whether there is only DEL variant in input bed file")
    parser.add_argument("-n", "--num", type=non_negative_int, default=5, 
                        help='number of rounds to generate randomized sequences')

    return parser

def detect_overlapping(input_df, focus_row_num, focus_start_pos, focus_end_pos, focus_chrom):
    '''
    forward scan variants to identify those overlapping with the focus variant.
    Note that the scanning happens in a cascading/transitive fashion: a variant
    is added to the overlapping set and if a forth-coming variant overlaps with
    the overlapping set, then it is added to the overlapping set as well.
    This process is repeated until no new variant is overlapping with the focus
    overlapping set.
    '''
    left_pos = focus_start_pos # the input is already sorted with sort -k1,1V -k2,2n -k3,3n
                               # so that any forth-coming variant's start pos is >= focus'
    # initialization
    right_pos = focus_end_pos
    if focus_row_num + 1 == input_df.shape[0]:
        return left_pos, right_pos, focus_row_num

    exclude_last_row = True


    for variant_row_num in range(focus_row_num+1, input_df.shape[0]):

        variant = input_df.loc[variant_row_num]

        # overlapping
        if focus_chrom == variant['chrom']:
            assert variant['start_pos'] >= left_pos, "Error: Out of order variant."
            if variant['start_pos'] < right_pos:
                right_pos = max(right_pos, variant['end_pos'])
                exclude_last_row = False
            else:
                exclude_last_row = True
                break
        else:
            exclude_last_row = True
            break

    if exclude_last_row:
        variant_row_num = variant_row_num - 1

    return left_pos, right_pos, variant_row_num

def shuffle_with_no_reptition(ref, NUM):
    '''
    ref: the reference sequence to shuffle
    NUM: number of new sequences to generate
    '''

    alt_list = []
    char_list = list(ref)
    for i in range(0, NUM):
        tries = 100
        while tries > 0: # at most try 100 times
            random.shuffle(char_list)
            t = "".join(char_list)

            if t != ref and t not in alt_list:
                alt_list.append(t)
                break;

            tries -= 1

    return alt_list

def randomize_choice_with_no_reptition(ref, NUM):
    '''
    ref: the reference sequence
    NUM: number of sequences to produce
    '''
   
    alt_list = []
    nbases = len(ref)

    for i in range(0, NUM):
        tries = 100
        while tries > 0: # at most try 100 times
            char_list = random.choices(['A', 'T', 'C', 'G'], k=nbases)
            t = "".join(char_list)
            if t != ref and t not in alt_list:
                alt_list.append(t)
                break;

            tries -= 1

    return alt_list


def transform_overlapping_variants_with_randomized_unified_ref(left_pos, right_pos, focus_row_num, variant_end_row, 
                                                                chromsome, input_df, NUM, SCRAMBLE):
    '''
    Transform all overlapping variants into a uniform variant. 
    After merging, instead of using original individual variants within the resultant window,
    use scrambled NUM references as the alternative variants, per the discussion between Shannon and Tao
    on 05/01/2023: https://docs.google.com/document/d/1ubLjP6WGLonEaiOj8WcyfCOJgmBYI5EQgJJmvJvPco4/edit
    '''

    # step 2. Extract reference sequence
    ref_sequence = str(chromsome[left_pos:right_pos])

    # step 3. Transform
    uniform_variant = dict()

    uniform_variant['chrom'] = input_df.loc[focus_row_num, "chrom"]
    uniform_variant['start_pos'] = left_pos
    uniform_variant['end_pos'] = right_pos
    uniform_variant['annotation'] = input_df.loc[focus_row_num, "annotation"]
    uniform_variant['ref'] = ref_sequence

    random.seed(42)
    # since this function should be called on DEL variant only, there is only one variant
    # as "" indicating the deletion of the reference sequence in variant['ref'].
    if SCRAMBLE:
        alt_lst = shuffle_with_no_reptition(uniform_variant['ref'], NUM)
    else:
        alt_lst = randomize_choice_with_no_reptition(uniform_variant['ref'], NUM)

    new_alt_lst = ",".join(alt_lst)
    uniform_variant['alt_lst'] = new_alt_lst

    return uniform_variant

def transform_overlapping_variants_with_randomized_variant_ref(left_pos, right_pos, focus_row_num, variant_end_row, 
                                                        chromsome, input_df, NUM, SCRAMBLE):
    '''
    Transform all overlapping variants into a uniform variant, but scramble the variant's ref sequences NUM times.
    For all variants w/ the same ref, there is only NUM scrambled or randomized ref sequences.
    '''

    # step 2. Extract reference sequence
    ref_sequence = str(chromsome[left_pos:right_pos])

    # step 3. Transform
    uniform_variant = dict()

    uniform_variant['chrom'] = input_df.loc[focus_row_num, "chrom"]
    uniform_variant['start_pos'] = left_pos
    uniform_variant['end_pos'] = right_pos
    uniform_variant['annotation'] = input_df.loc[focus_row_num, "annotation"]
    uniform_variant['ref'] = ref_sequence

    random.seed(42)
    for row in range(focus_row_num, variant_end_row+1):
        variant = input_df.loc[row]
        prefix = ref_sequence[0 : variant['start_pos'] - left_pos]
        postfix = ref_sequence[variant['end_pos'] - left_pos : right_pos - left_pos]

        # since this function should be called on DEL variant only, there is only one variant
        # as "" indicating the deletion of the reference sequence in variant['ref'].
        # new alt: prefix + shuffle(variant['ref']) + postfix
        if SCRAMBLE:
            alt_lst = shuffle_with_no_reptition(variant['ref'], NUM)
        else:
            alt_lst = randomize_choice_with_no_reptition(variant['ref'], NUM)

        new_alt_lst = ",".join([prefix + alt + postfix for alt in alt_lst])
        if 'alt_lst' in uniform_variant:
            uniform_variant['alt_lst'] = uniform_variant['alt_lst'] + "," + new_alt_lst
        else:
            uniform_variant['alt_lst'] = new_alt_lst

    return uniform_variant

def transform_overlapping_variants_with_unified_ref_deletion(left_pos, right_pos, focus_row_num, variant_end_row, chromsome, input_df):
    '''
    Transform all overlapping variants into a uniform variant.
    Instead of deleting individual variant reference, delete the resultant unified reference.
    '''

    # step 2. Extract reference sequence
    ref_sequence = str(chromsome[left_pos:right_pos])

    # step 3. Transform
    uniform_variant = dict()

    uniform_variant['chrom'] = input_df.loc[focus_row_num, "chrom"]
    uniform_variant['start_pos'] = left_pos
    uniform_variant['end_pos'] = right_pos
    uniform_variant['annotation'] = input_df.loc[focus_row_num, "annotation"]
    uniform_variant['ref'] = ref_sequence
    uniform_variant['alt_lst'] = ''

    return uniform_variant

def transform_overlapping_variants_with_variant_ref_deletion(left_pos, right_pos, focus_row_num, variant_end_row, chromsome, input_df):
    '''
    Transform all overlapping variants into a uniform variant.
    Delete the reference for all variants with the same reference.
    '''

    # step 2. Extract reference sequence
    ref_sequence = str(chromsome[left_pos:right_pos])

    # step 3. Transform
    uniform_variant = dict()

    uniform_variant['chrom'] = input_df.loc[focus_row_num, "chrom"]
    uniform_variant['start_pos'] = left_pos
    uniform_variant['end_pos'] = right_pos
    uniform_variant['annotation'] = input_df.loc[focus_row_num, "annotation"]
    uniform_variant['ref'] = ref_sequence
    
    for row in range(focus_row_num, variant_end_row+1):
        variant = input_df.loc[row]
        prefix = ref_sequence[0 : variant['start_pos'] - left_pos]
        postfix = ref_sequence[variant['end_pos'] - left_pos : right_pos - left_pos]
        
        new_alt_lst = prefix + postfix # w/ variant reference deletion for all variants w/ the same reference
        if 'alt_lst' in uniform_variant:
            uniform_variant['alt_lst'] = uniform_variant['alt_lst'] + "," + new_alt_lst
        else:
            uniform_variant['alt_lst'] = new_alt_lst


    return uniform_variant

def transform_overlapping_variants(left_pos, right_pos, focus_row_num, variant_end_row, chromsome, input_df):
    '''
    Transform all overlapping variants into a uniform variant w/o removal/scrambling/randomization of variants
    or reference sequences.
    '''

    # step 2. Extract reference sequence
    ref_sequence = str(chromsome[left_pos:right_pos])

    # step 3. Transform
    uniform_variant = dict()

    uniform_variant['chrom'] = input_df.loc[focus_row_num, "chrom"]
    uniform_variant['start_pos'] = left_pos
    uniform_variant['end_pos'] = right_pos
    uniform_variant['annotation'] = input_df.loc[focus_row_num, "annotation"]
    uniform_variant['ref'] = ref_sequence

    for row in range(focus_row_num, variant_end_row+1):
        variant = input_df.loc[row]
        prefix = ref_sequence[0 : variant['start_pos'] - left_pos]
        postfix = ref_sequence[variant['end_pos'] - left_pos : right_pos - left_pos]

        # new alt: prefix + alt + postfix
        alt_lst = variant['alt_lst'].split(',')
        new_alt_lst = ",".join([prefix + alt + postfix for alt in variant['alt_lst'].split(',')])
        if 'alt_lst' in uniform_variant:
            uniform_variant['alt_lst'] = uniform_variant['alt_lst'] + "," + new_alt_lst
        else:
            uniform_variant['alt_lst'] = new_alt_lst

    return uniform_variant

def output_variant(FILE_OUTPUT, uniform_variant):
    '''
    output uniform variant into .bed file w/o header
    '''

    FILE_OUTPUT.write(f"{uniform_variant['chrom']}\t"
                      f"{uniform_variant['start_pos']}\t"
                      f"{uniform_variant['end_pos']}\t"
                      f"{uniform_variant['annotation']}\t"
                      f"{uniform_variant['ref']}\t"
                      f"{uniform_variant['alt_lst']}\n")


def main():
    parser = argument_parser()

    # paring the arguments
    args = parser.parse_args()
    GENOME = args.genome
    INPUT = args.input
    OUTPUT = args.output
    VARIANT = args.variant
    RANDOM=args.random
    SCRAMBLE=args.scramble
    DELETION_ONLY = args.deletion_only
    NUM=args.num

    if not OUTPUT:
        OUTPUT = INPUT.replace(".bed", "") + "_merged" + ".bed"

    if GENOME.endswith('fasta') or GENOME.endswith('fa'):
        FILE_GENOME = SeqIO.parse(open(GENOME), 'fasta')

    FILE_OUTPUT = open(OUTPUT, 'w')

    # read input into dataframe
    # https://cmdlinetips.com/2018/09/how-to-change-data-type-for-one-or-more-columns-in-pandas-dataframe/
    input_df = pd.read_csv(INPUT, sep="\t", header=None, names=columns, dtype={'alt_lst': str})
    input_df.alt_lst = input_df.alt_lst.fillna('') # https://stackoverflow.com/questions/26837998/pandas-replace-nan-with-blank-empty-string
    #input_df.columns = columns
    input_df['chrom'].unique().tolist()

    # initial scan of the chromosomes of the given variants
    CHROM_set = set()
    for chrom in input_df['chrom'].unique().tolist():
        CHROM_set.add(chrom)

    # read the genome sequence, and load the sequences of the given chromosomes
    for each_genome in FILE_GENOME:
        chromosome, sequence = each_genome.id, each_genome.seq
        if chromosome in CHROM_set:
            CHROM_SEQ_dict[chromosome] = sequence

    last_row_num = input_df.shape[0]
    rows_processed = 0
    variant_end_row = -1
    random.seed(42)
    with tqdm(total=last_row_num) as pbar:
        focus_row_num=0
        while focus_row_num < last_row_num:
            focus_variant = input_df.loc[focus_row_num]
            focus_chrom = focus_variant['chrom']
            focus_start_pos = focus_variant['start_pos']
            focus_end_pos = focus_variant['end_pos']
            focus_ref = focus_variant['ref']
            focus_alts_lst = focus_variant['alt_lst']
            chrom_sequence = CHROM_SEQ_dict[focus_chrom]

            # step 1. Figure out left-most coordinate and right-most coordinate
            left_pos, right_pos, variant_end_row = detect_overlapping(input_df, focus_row_num, focus_start_pos, focus_end_pos, focus_chrom)

            # step 2. transform to a uniform variant
            if focus_row_num < variant_end_row : # overlapping variants detected
                if DELETION_ONLY: # DEL only
                    if RANDOM: # scramble or randomize
                        if VARIANT: # target individual variant
                            uniform_variant = transform_overlapping_variants_with_randomized_variant_ref(left_pos, right_pos, focus_row_num, 
                                                                                              variant_end_row, chrom_sequence, input_df,
                                                                                              NUM, SCRAMBLE)
                        else:
                            uniform_variant = transform_overlapping_variants_with_randomized_unified_ref(left_pos, right_pos, focus_row_num,
                                                                                                    variant_end_row, chrom_sequence,
                                                                                                    input_df, NUM, SCRAMBLE)
                    else: # deletion
                        if VARIANT: # target individual variant
                            uniform_variant = transform_overlapping_variants_with_variant_ref_deletion(left_pos, right_pos, 
                                                                                                    focus_row_num, variant_end_row, 
                                                                                                    chrom_sequence, input_df)
                        else:
                            uniform_variant = transform_overlapping_variants_with_unified_ref_deletion(left_pos, right_pos, 
                                                                                                focus_row_num, variant_end_row, 
                                                                                                chrom_sequence, input_df)
                else: # DEL mixed with INS/SNP
                    uniform_variant = transform_overlapping_variants(left_pos, right_pos, focus_row_num, variant_end_row,
                                                                     chrom_sequence, input_df)

                rows_processed += variant_end_row - focus_row_num + 1
                focus_row_num = variant_end_row + 1
                
                # step 3. write out the uniform variant
                if uniform_variant['ref'] in uniform_variant['alt_lst'].split(","):
                    continue

                output_variant(FILE_OUTPUT, uniform_variant)
            else: # non-overlapping cases
                rows_processed += variant_end_row - focus_row_num + 1
                focus_row_num += 1
                focus_variant_copy = focus_variant.copy()
                
                # step 3. write out the uniform variant
                if DELETION_ONLY: # DEL only
                    if RANDOM: # scramble or randomize
                        alt_lst = []
                        
                        if SCRAMBLE: # scramble
                            alt_lst = shuffle_with_no_reptition(focus_ref, NUM) 
                        else: # randomize
                            alt_lst = randomize_choice_with_no_reptition(focus_ref, NUM)

                        if focus_ref in alt_lst:
                            continue

                        focus_variant_copy['alt_lst'] = ",".join(alt_lst)

                output_variant(FILE_OUTPUT, focus_variant_copy)

            # update the progress
            if rows_processed % 1000 == 0:
                pbar.update(1)
                rows_processed = 0

    # cleanup
    FILE_OUTPUT.close()


if __name__ == '__main__':
    main()
