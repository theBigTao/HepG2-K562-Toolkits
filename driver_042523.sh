module load R/4.2.0

cell_line=files

tfs=(`awk '{print $1}' ${cell_line}`)
input_files=(`awk '{print $NF}' ${cell_line}`)
upbound=`wc -l files | awk '{print $1/5-1}'`
runDate='042523'

submit_jobs() {

    job_submitted=0
    
    for i in `seq 0 ${upbound}`
    do
        id=$((i*5))
        echo "${tfs[$id]} submitted ..."
        input_file=${input_files[$id]}
        output_file=/scratch/users/taowang9/${runDate}/files/${tfs[$id]}/fold4/S4/prediction.tsv
        
        # test whether the input file has been successfully processed.
        input_file_size=`wc -l ${input_file} | awk '{print $1}'`
        output_file_size=0
        if [ -f ${output_file} ]; then
            output_file_size=`wc -l ${output_file} | awk '{print $1}'`
        fi
    
        if [ ${input_file_size} -eq ${output_file_size} ];
        then
            continue
        fi
    
        sbatch --wait predict_single_round4.sh $id ${cell_line} &
        pids[${id}]=$!
    
        job_submitted=`squeue -u taowang9 | wc -l`
        while [ ${job_submitted} -gt 48 ]; # be to precise, the biggest job_submitted could get is 51
        do
            sleep 10
            job_submitted=`squeue -u taowang9 | wc -l`
        done
    
        sleep 10
    done
    
    # wait for all jobs
    for pid in ${pids[@]}
    do
        #echo "wait $pid" 
        wait $pid
    done
}

statistical_test() {
    # pvalues: both pair t-test and paired wilcox rank test
    # effect size: log(median(variants)/median(refs))

    echo "Statistical test w/ both pair t-test and paired wilcox rank test"

    # median
    Rscript statistics.R ${SCRATCH}/042523/files "median_ref.txt" "median_var.txt" "tf1_with_tf2_motif_scrambled_median.tsv"

    # median of median
    Rscript statistics.R ${SCRATCH}/042523/files "median_of_median_ref.txt" "median_of_median_var.txt" "tf1_with_tf2_motif_scrambled_median_of_median.tsv"
}

#submit_jobs
statistical_test
