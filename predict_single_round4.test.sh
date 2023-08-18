#!/bin/bash
#SBATCH --job-name=ChIP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --account=mpsnyder
#SBATCH --time=12:00:00
#SBATCH --gpus 1
#SBATCH --partition=gpu

#srun --nodes=1 --ntasks=1 --cpus-per-task=10 --partition=interactive --account=default --time=24:00:00 --pty /bin/bash
#srun --nodes=1 --ntasks=1 --cpus-per-task=10 -p gpu --gpus 1 --time=24:00:00 --pty /bin/bash
#srun --nodes=1 --ntasks=1 --cpus-per-task=10 -p gpu --account=mpsnyder --gpus 1 --time=24:00:00 --pty /bin/bash

# set up environment for models

module load cuda/11.2.0
module load cudnn/8.1.1.33
source conda/etc/profile.d/conda.sh 

#conda activate py39 # round1

conda activate py36 # round3
module load py-tensorflow/2.4.1_py36 # round3


extract_median_counts_by_median() {
    echo "Extract median of all 5*5 runs"

    work_dir=$1
    for fold in {0..4}
    do
        for batch in {0..4}
        do
            awk -F "\t" 'NR>1 {print $2}' ${work_dir}/fold${fold}/S${batch}/prediction.tsv > b_${batch}_ref
            awk -F "\t" 'NR>1 {print $3}' ${work_dir}/fold${fold}/S${batch}/prediction.tsv > b_${batch}_var
        done

        paste -d "\t" b_0_ref b_1_ref b_2_ref b_3_ref b_4_ref > m_${fold}_ref
        paste -d "\t" b_0_var b_1_var b_2_var b_3_var b_4_var > m_${fold}_var
    done

    paste -d "\t" m_0_ref m_1_ref m_2_ref m_3_ref m_4_ref | awk -F "\t" '{n=split($0, arr, "\t"); asort(arr); if(n%2==0) median=(arr[n/2]+arr[n/2+1])/2; else median=arr[(n+1)/2]; print median}' > ${work_dir}/median_ref.txt
    paste -d "\t" m_0_var m_1_var m_2_var m_3_var m_4_var | awk -F "\t" '{n=split($0, arr, "\t"); asort(arr); if(n%2==0) median=(arr[n/2]+arr[n/2+1])/2; else median=arr[(n+1)/2]; print median}' > ${work_dir}/median_var.txt    
}

extract_median_counts_by_median_of_median() {
    work_dir=$1

    echo "Extract median of medians of total var,ref counts"

    for fold in {0..4}
    do
        for batch in {0..4}
        do
            awk -F "\t" 'NR>1 {print $2}' ${work_dir}/fold${fold}/S${batch}/prediction.tsv > b_${batch}_ref
            awk -F "\t" 'NR>1 {print $3}' ${work_dir}/fold${fold}/S${batch}/prediction.tsv > b_${batch}_var
        done

        paste -d "\t" b_0_ref b_1_ref b_2_ref b_3_ref b_4_ref | awk -F "\t" '{n=split($0, arr, "\t"); asort(arr); if(n%2==0) median=(arr[n/2]+arr[n/2+1])/2; else median=arr[(n+1)/2]; print median}' > m_${fold}_ref
        paste -d "\t" b_0_var b_1_var b_2_var b_3_var b_4_var | awk -F "\t" '{n=split($0, arr, "\t"); asort(arr); if(n%2==0) median=(arr[n/2]+arr[n/2+1])/2; else median=arr[(n+1)/2]; print median}' > m_${fold}_var
    done

    paste -d "\t" m_0_ref m_1_ref m_2_ref m_3_ref m_4_ref | awk -F "\t" '{n=split($0, arr, "\t"); asort(arr); if(n%2==0) median=(arr[n/2]+arr[n/2+1])/2; else median=arr[(n+1)/2]; print median}' > ${work_dir}/median_of_median_ref.txt
    paste -d "\t" m_0_var m_1_var m_2_var m_3_var m_4_var | awk -F "\t" '{n=split($0, arr, "\t"); asort(arr); if(n%2==0) median=(arr[n/2]+arr[n/2+1])/2; else median=arr[(n+1)/2]; print median}' > ${work_dir}/median_of_median_var.txt    
}

predict() {
    cd /home/users/taowang9/ChIP
    i=$1
    cell_line=$2


    tfs=(`awk '{print $1}' ${cell_line} `)
    models=(`awk '{print $2}' ${cell_line} `)
    variant_files=(`awk '{print $3}' ${cell_line} `)
   
    prefix="/home/users/taowang9/ChIP"
    #predict_tool="${prefix}/predict_new.py" # round 1
    predict_tool="${prefix}/predict_round4_042523.py" # round 3

    echo "${tfs[$i]} prediction in progress ..."
    tf=${tfs[$i]}
    
    dir=/scratch/users/taowang9/042523/${cell_line}

    if [ ! -d ${dir}/${tf} ];
    then
        mkdir -p ${dir}/${tf}
    fi
    
    for fold in {0..4}
    do
        if [ ! -d ${dir}/${tf}/fold${fold} ];
        then
            mkdir ${dir}/${tf}/fold${fold}
        fi

        # all the five input batch share the same model
        model=${prefix}/Models_from_Anshul/fold${fold}/${models[$i]}/${models[$i]}_split000

        for batch in {0..4}
        do
            echo "tf ${tf}, fold ${fold}, batch $((batch+1))"
            if [ ! -d ${dir}/${tf}/fold${fold}/S${batch} ];
            then
                mkdir ${dir}/${tf}/fold${fold}/S${batch}
            fi

            output_dir=${dir}/${tf}/fold${fold}/S${batch}
            
            variant_file=${variant_files[i+batch]}
            #echo "python ${predict_tool} --model $model -f ${variant_file} -o ${output_dir} > /dev/null"
            python ${predict_tool} --model $model -f ${variant_file} -o ${output_dir} > /dev/null
        done
    done

    # extract variant/ref counts
    echo "Extract medians"
    extract_median_counts_by_median_of_median ${dir}/${tf}
    extract_median_counts_by_median ${dir}/${tf}
}

predict $1 $2
