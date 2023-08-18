'''
Model wrapper for basepairmodel

It takes a sequence file (in the format defined by
our universal genome sequence tailor tool) and makes prediction for TF-binding signal profile
for each input sequence.

Written by Tao Wang, Department of Genetics, School of Medicine, Stanford University, CA, USA.

Contributors: 
    Vivekanandan Ramalingam, Department of Computer Science, Stanford University, CA, USA.
'''

from tensorflow import keras
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import CustomObjectScope

import argparse
import time
import datetime
import json
import os
import pandas as pd
import sys
import pysam
from scipy.special import logsumexp
import imageio

from plotnine import *

import matplotlib as mpl
from tqdm import tqdm
import pandas as pd
import numpy as np
import random
from math import ceil

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


from tensorflow.keras import Model
import logging

import tensorflow as tf
import tensorflow.keras.backend as kb

import tensorflow_probability as tfp


def mse_loss_function(y_log_true, y_log_pred):
    # logcounts mse loss without sample weights
    mse_loss = tf.keras.losses.mean_squared_error(
            y_log_true, y_log_pred)
    return mse_loss

def poisson_loss_function(y_log_true, y_log_pred):
    # we can use the Possion PMF from TensorFlow as well
    # dist = tf.contrib.distributions
    # return -tf.reduce_mean(dist.Poisson(y_pred).log_pmf(y_true))

    # last term can be avoided since it doesn't depend on y_pred
    # however keeping it gives a nice lower bound to zero

    y_true = tf.math.exp(y_log_true) # -1? 
    y_pred = tf.math.exp(y_log_pred)
    y_true = tf.cast(y_true,tf.float32)
    y_pred = tf.cast(y_pred,tf.float32)
    loss = y_pred - y_true*tf.math.log(y_pred+1e-8) + tf.math.lgamma(y_true+1.0)

    return loss

def multinomial_nll(true_counts, logits):
    """Compute the multinomial negative log-likelihood
    Args:
        true_counts: observed count values
        logits: predicted logits values
    """
    counts_per_example = tf.reduce_sum(true_counts, axis=-1)
    dist = tfp.distributions.Multinomial(total_count=counts_per_example,
            logits=logits)
    return (-tf.reduce_sum(dist.log_prob(true_counts)) / 
            tf.cast(tf.shape(true_counts)[0], dtype=tf.float32))

    from tensorflow import keras
from tensorflow.keras import Model

class CustomModel(Model):

    def __init__(self, num_tasks, num_output_tracks, tracks_for_each_task, output_profile_len, loss_weights,counts_loss, orig_multi_loss=False, **kwargs):

        # call the base class with inputs and outputs
        super(CustomModel, self).__init__(**kwargs)

        # number of tasks
        self.num_tasks = num_tasks

        # number of output tracks used for original multinomial loss
        self.num_output_tracks = num_output_tracks

        self.orig_multi_loss = orig_multi_loss

        # number of tracks for each task
        self.tracks_for_each_task = tracks_for_each_task

        # output profile length
        self.output_profile_len = output_profile_len

        # weights for the profile mnll and logcounts losses
        self.loss_weights = loss_weights

        # logcounts loss funtion
        self.counts_loss = counts_loss

        # object to track overall mean loss per epoch
        self.loss_tracker = keras.metrics.Mean(name="loss")


    def _get_loss(self, x, y, sample_weights, training=True):
        # boolean mask for sample weights != 0                
        y_pred = self(x, training=training)  # Forward pass        

        def _count_loss_function(_y_log_true,_y_log_pred, count_loss_fn):
            # count_loss_fn: either poisson_loss_function or mse_loss_function
            total_count_loss = 0
            track_count_cuml = 0
            if self.orig_multi_loss:
                for i in range(self.num_output_tracks):
                    y_log_true = _y_log_true[:,i:(i+1)][:,-1]
                    y_log_pred = _y_log_pred[:,i:(i+1)][:,-1]
                    loss = count_loss_fn(y_log_true, y_log_pred)
                    total_count_loss += loss               
            else:
                for i in range(self.num_tasks):
                    num_of_tracks = self.tracks_for_each_task[i]
                    y_log_true = tf.reduce_logsumexp(_y_log_true[:,track_count_cuml:(track_count_cuml+num_of_tracks)],axis=1)
                    y_log_pred = _y_log_pred[:,i:(i+1)][:,-1]
                    loss = count_loss_fn(y_log_true, y_log_pred)
                    total_count_loss += loss
                    track_count_cuml += num_of_tracks
            return total_count_loss

        if self.counts_loss == "MSE":
            total_counts_loss = _count_loss_function(y['logcounts_predictions'],y_pred[1], mse_loss_function)

        elif self.counts_loss == "POISSON":        
            total_counts_loss = _count_loss_function(y['logcounts_predictions'],y_pred[1], poisson_loss_function)            

        else:
            raise Exception("Sorry, unknown loss funtion")


        # for mnll loss we mask out samples with weight == 0.0        
        boolean_mask = tf.math.greater_equal(sample_weights, 1.0)

        _y = tf.boolean_mask(y['profile_predictions'], boolean_mask)
        _y_pred = tf.boolean_mask(y_pred[0], boolean_mask)

        def _zero_constant():
            return kb.constant(0)

        def _multinomial_nll(_y,_y_pred):
            total_mnll_loss = 0
            track_count_cuml = 0

            if self.orig_multi_loss:
                for i in range(self.num_output_tracks):
                    loss = multinomial_nll(_y[..., i], _y_pred[..., i])
                    total_mnll_loss += loss

            else:
                for i in range(self.num_tasks):
                    num_of_tracks = self.tracks_for_each_task[i]
                    _y_reshape = tf.reshape(\
                            _y[:,:,track_count_cuml:(track_count_cuml+num_of_tracks)],\
                            [-1,(num_of_tracks)*(self.output_profile_len)]\
                            )
                    _y_pred_reshape = tf.reshape(\
                            _y_pred[:,:,track_count_cuml:(track_count_cuml+num_of_tracks)],\
                            [-1,(num_of_tracks)*(self.output_profile_len)]\
                            )

                    loss = multinomial_nll(_y_reshape, _y_pred_reshape)
                    track_count_cuml = track_count_cuml+num_of_tracks
                    total_mnll_loss += loss
            return total_mnll_loss

        # edge case where no example with non-zero weight             
        total_mnll_loss = tf.cond(tf.equal(tf.size(_y), 0), 
                _zero_constant,
                lambda:  _multinomial_nll(_y,_y_pred))

        if self.counts_loss == "MSE":
            loss =  (self.loss_weights[0] * total_mnll_loss) + \
                    (self.loss_weights[1] * total_counts_loss)   
        elif self.counts_loss == "POISSON":

            loss =  total_mnll_loss + total_counts_loss            
        else:
            raise Exception("Sorry, unknown loss funtion")

        return loss, total_mnll_loss, total_counts_loss

    def train_step(self, data):
        x, y, sample_weights = data       

        with tf.GradientTape() as tape:
            loss, total_mnll_loss, total_counts_loss = \
                    self._get_loss(x, y, sample_weights)
        # Compute gradients
        trainable_vars = self.trainable_variables
        gradients = tape.gradient(loss, trainable_vars)

        # Update weights
        self.optimizer.apply_gradients(zip(gradients, trainable_vars))

        # Compute our own metrics
        self.loss_tracker.update_state(loss)
        return {"loss": self.loss_tracker.result(),
                "batch_loss": loss,
                "profile_predictions_loss": total_mnll_loss, 
                "logcounts_predictions_loss": total_counts_loss}

    @property
    def metrics(self):
        # We list our `Metric` objects here so that `reset_states()` can be
        # called automatically at the start of each epoch
        # or at the start of `evaluate()`.
        # If you don't implement this property, you have to call
        # `reset_states()` yourself at the time of your choosing.
        return [self.loss_tracker]


    def test_step(self, data):
        # Unpack the data
        x, y, sample_weights = data

        loss, total_mnll_loss, total_counts_loss = \
                self._get_loss(x, y, sample_weights, training=False)

        # Compute our own metrics
        self.loss_tracker.update_state(loss)
        return {"loss": self.loss_tracker.result(),
                "batch_loss": loss,
                "profile_predictions_loss": total_mnll_loss, 
                "logcounts_predictions_loss": total_counts_loss}

def get_model(model_path):
    with CustomObjectScope({'MultichannelMultinomialNLL': lambda n='0':n,
        "kb": kb,
        "CustomMeanSquaredError":lambda n='0':n,
        "tf":tf,
        "CustomModel":CustomModel}):
        model = load_model(model_path)
    return model


def random_seq(seqlen):
    return ''.join(random.choices("ACGT", k=seqlen))


def fix_sequence_length(sequence, length):
    """
        Function to check if length of sequence matches specified
        length and then return a sequence that's either padded or
        truncated to match the given length
        Args:
            sequence (str): the input sequence
            length (int): expected length
        Returns:
            str: string of length 'length'
    """

    # check if the sequence is smaller than expected length
    if len(sequence) < length:
        # pad the sequence with 'N's
        sequence += 'N' * (length - len(sequence))
    # check if the sequence is larger than expected length
    elif len(sequence) > length:
    # truncate to expected length
        sequence = sequence[:length]

    return sequence
def one_hot_encode(sequences, seq_length):
    """

       One hot encoding of a list of DNA sequences 

       Args:
           sequences (list): python list of strings of equal length
           seq_length (int): expected length of each sequence in the 
               list

       Returns:
           numpy.ndarray: 
               3-dimension numpy array with shape 
               (len(sequences), len(list_item), 4)
    """

    if len(sequences) == 0:
        logging.error("'sequences' is empty")
        return None

    # First, let's make sure all sequences are of equal length
    sequences = list(map(
        fix_sequence_length, sequences, [seq_length] * len(sequences)))

    # Step 1. convert sequence list into a single string
    _sequences = ''.join(sequences)

    # Step 2. translate the alphabet to a string of digits
    transtab = str.maketrans('ACGTNYRMSWK', '01234444444')    
    sequences_trans = _sequences.translate(transtab)

    # Step 3. convert to list of ints
    int_sequences = list(map(int, sequences_trans))

    # Step 4. one hot encode using int_sequences to index 
    # into an 'encoder' array
    encoder = np.vstack([np.eye(4), np.zeros(4)])
    X = encoder[int_sequences]

    # Step 5. reshape 
    return X.reshape(len(sequences), len(sequences[0]), 4)

def argsparser_input():
    """ Command line arguments for the predict script
    for our wrapper
    """

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--batch-size', '-b', type=int, 
                        help="length of output profile", default=500)
    
    parser.add_argument('--vectorized', '-v', action='store_true', default=False,
                        help="whether to generate per-base signal values")
    
    parser.add_argument('--signal-profiles', '-s', action='store_true', default=False,
                        help="whether to generate per-base signal values")
    
    # path to the pre-trained basepair model
    parser.add_argument('--model', '-m', type=str,
            help="path to the saved model: include one file saved_model.pb, and two directories: assests and variables")

    # input sequence file name
    parser.add_argument('--target-sequence-file', '-f', type=str, required=True,
                        help='the input file containing target sequence. Default to example.tsv')

    parser.add_argument('--output-dir', '-o', type=str, required=True,
                        help="destination directory to store predictions as a "
                        "bigWig file")

    return parser

def random_seq(seqlen):
    return ''.join(random.choices("ACGT", k=seqlen))

#encoded_inserted_sequences = one_hot_encode([random_seq(2114) for i in range(900)],2114)

class NoTracebackException(Exception):
    """
        An exception that when raised results in the error message
        being printed without the traceback
    """
    pass

#option 1
def scalar_prediction_to_profile(batch_size, predictions):
    output_len = 1000
    num_output_tracks = 2
    pred_profiles = np.zeros((batch_size, output_len, num_output_tracks))
    pred_logcounts = np.zeros((batch_size, num_output_tracks))
    logcounts_prediction_arr = np.zeros((batch_size, 1))
    
    for idx in range(batch_size):
        for j in range(num_output_tracks): 
            # combined counts prediction from the count head
            logcounts_prediction = predictions[1][idx] # this is the logsum of both the strands that is 
                                                       # predicted       
            logcounts_prediction_arr[idx]=logcounts_prediction
    
            # predicted profile — for now assuming that there are two strands and only one task
            pred_profile_logits = np.reshape(predictions[0][idx, :, :],[1,output_len*2])
    
            profile_predictions = (np.exp(pred_profile_logits - logsumexp(pred_profile_logits)) * (np.exp(logcounts_prediction)))
    
            pred_profiles[idx, :, j] = np.reshape(profile_predictions,[output_len,2])[:,j]
    
            # counts prediction
            pred_logcounts[idx, j] = np.log(np.sum(np.reshape(profile_predictions, [output_len, 2])[:,j]))

    plus = pred_profiles[:,:,0]
    minus = pred_profiles[:,:,1]
    counts_arr = logcounts_prediction_arr
    return (plus,minus,counts_arr)

#option 2
def vectorized_prediction_to_profile(predictions):
    output_len = 1000
    logits_arr = predictions[0]
    counts_arr = predictions[1]
    pred_profile_logits = np.reshape(logits_arr,[-1,1,output_len*2])
    probVals_array = np.exp(pred_profile_logits-logsumexp(pred_profile_logits,axis=2).reshape([len(logits_arr),1,1]))
    profile_predictions = np.multiply((np.exp(counts_arr)).reshape([len(counts_arr),1,1]),probVals_array)
    plus = np.reshape(profile_predictions,[len(counts_arr),output_len,2])[:,:,0]
    minus = np.reshape(profile_predictions,[len(counts_arr),output_len,2])[:,:,1]

    return (plus,minus,counts_arr)

def extract_total_counts_and_profiles(predictions, vectorized):
    if vectorized:
        plus, minus, logcounts_arr = vectorized_prediction_to_profile(predictions)
    else:
        plus, minus, logcounts_arr = scalar_prediction_to_profile(predictions[0].shape[0], predictions)

    counts = np.exp(logcounts_arr)
    per_base_counts = plus + minus

    return counts, per_base_counts

def write_results_to_tsv(TSV_FILE_OUTPUT, start, end, variant_counts, ref_counts, variant_profiles, ref_profiles):

    for row in range(start, end+1):
        current_idx = row-start
        Sequence_num_tsv = row
        Ref_total_count_N_tsv = ref_counts[current_idx][0]
        Variant_total_count_N_tsv = variant_counts[current_idx][0]
        Ref_signal_values_tsv = ",".join(["%.6f" % e for e in ref_profiles[current_idx].tolist()])
        Variant_signal_values_tsv = ",".join(["%.6f" % e for e in variant_profiles[current_idx].tolist()])
        TSV_FILE_OUTPUT.write(f'{Sequence_num_tsv}\t{Ref_total_count_N_tsv}\t{Variant_total_count_N_tsv}'
                          f'\t{Ref_signal_values_tsv}\t{Variant_signal_values_tsv}\n')

def predict(args, input_data, pred_dir):
    # load the model
    model = get_model(args.model)

    # tsv
    TSV_FILE_OUTPUT = open(pred_dir+'/prediction.tsv', 'w')
    TSV_FILE_HEADER = "Sequence_num\tRef_total_count_N\tVariant_total_count_N\tRef_signal_values\tVariant_signal_values\n"
    TSV_FILE_OUTPUT.write(TSV_FILE_HEADER)

    input_df = pd.read_csv(input_data, sep='\t')
    last_row_num = input_df.shape[0]
    progress = 0
    
    batch_size = args.batch_size
    output_len = 1000
    output_seq_len = 1000
    number_of_strands = 2
    num_output_tracks = 2 # num_output_tracks=2 represents both the strands
    last_processed_row = -1 # last processed row number  
    left_rows = last_row_num # number of left rows
    
    upperbound = ceil(last_row_num/batch_size)
    with tqdm(total=upperbound) as pbar:
        for iteration in range(upperbound):
            # prepare the input sequence in batch
            if left_rows >= batch_size:
                # process next 500 rows
                current_batch_size = batch_size
            else:
                # process all left_rows as the final batch
                current_batch_size = left_rows

            # one hot encoding for variant sequences and reference sequences
            variant_sequences = one_hot_encode([e.upper() for e in input_df.loc[last_processed_row+1:last_processed_row + current_batch_size]['Window_sequence']], 2114)
            reference_sequences = one_hot_encode([e.upper() for e in input_df.loc[last_processed_row+1:last_processed_row + current_batch_size]['Window_reference']], 2114)

            # prediction
            predictions = model.predict([variant_sequences,
                np.zeros(output_seq_len*number_of_strands*variant_sequences.shape[0]).reshape((variant_sequences.shape[0], output_seq_len, number_of_strands)),    
                np.zeros(variant_sequences.shape[0]*number_of_strands).reshape((variant_sequences.shape[0], number_of_strands))])
            
            predictions_ref = model.predict([reference_sequences,
                np.zeros(output_seq_len*number_of_strands*reference_sequences.shape[0]).reshape((reference_sequences.shape[0], output_seq_len, number_of_strands)),    
                np.zeros(reference_sequences.shape[0]*number_of_strands).reshape((reference_sequences.shape[0], number_of_strands))])

            # extract total counts and per-base signal profiles
            variant_counts, variant_profiles = extract_total_counts_and_profiles(predictions, args.vectorized)
            ref_counts, ref_profiles = extract_total_counts_and_profiles(predictions_ref, args.vectorized)

            #if args.signal_profiles
            # write results to tsv
            write_results_to_tsv(TSV_FILE_OUTPUT, last_processed_row+1, last_processed_row+current_batch_size, variant_counts, ref_counts, variant_profiles, ref_profiles)

            # progress report
            last_processed_row += current_batch_size
            left_rows -= current_batch_size
            pbar.update(1)


def predict_main():
    """ The main entry to our wrapper predictor.
    """

    # parse the command-line arguments
    parser = argsparser_input()
    args = parser.parse_args()

    # check if the output directory exists
    if not os.path.exists(args.output_dir):
        logging.error("Directory {} does not exist".format(args.output_dir))

    # make sure that the target sequence file exists
    if not os.path.isfile(args.target_sequence_file):
        raise NoTracebackException(
                "File not found: {} OR you may have accidentally "
                "specified a directory path.".format(args.target_sequence_file))

    if os.path.isdir(args.output_dir):
        pred_dir = args.output_dir
    else:
        logging.error("Directory does not exist {}.".format(args.output_dir))
        return

    input_data = args.target_sequence_file
    
    # predict
    predict(args, input_data, pred_dir)

if __name__ == '__main__':
    predict_main()
