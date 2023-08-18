# HepG2-K562-Toolkits
Toolkits for TF binding prediction, including preprocessing and wrapper for TF basepair models.
Standalone code:
1. merger.py: standardize adjacent variants into a uniform variant
2. SeqEnumerator_fast.py to enumerate sequences with variants in a vcf file
3. statistics.R: post processing code in R to calculate p-values and q-values.

# Step 1. Prepare inputs
driver_snp_ins.sh # for snp and insertion
driver_deletion.sh # for deletion

# Step 2. Run the models on HPC & Post-processing
driver_042523.sh # for snp and insertion
driver_042523_deletion.sh # for deletion
