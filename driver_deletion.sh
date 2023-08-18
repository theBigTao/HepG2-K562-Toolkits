root_dir=/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/prepareInput
runDate='042523'
output_dir=${root_dir}/${runDate}/TFs_intersect
motif_dir=${root_dir}/${runDate}/TFs_motifs

module_regions_dir=${root_dir}/Modules

cd ${root_dir}

module load bedtools
source ~/.pyenvrc # python3.9

tf_peaks_dir=/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/prepareInput/Shannon
targets_dir="/home/taowang9/ChIP/basepairmodel/basepairmodels/tutorial"

modules=(`awk -F "\t" 'NR > 1 {split($1, a, "_"); print a[1]}' Tao_Binding_Predictions.tsv`)
cells=(`awk -F "\t" 'NR > 1 { if ($1 ~ "_k") {print "K562"} else {print "HepG2"} }' Tao_Binding_Predictions.tsv`)
tfs1=(`awk -F "\t" 'NR > 1 {print $2}' Tao_Binding_Predictions.tsv`)
tfs2=(`awk -F "\t" 'NR > 1 {print $3}' Tao_Binding_Predictions.tsv`)

export motif_dir
export output_dir
export root_dir
export tf_peaks_dir

# Given a tf, return its tf chip ENC number.
find_tf_chip() {
    tf=$1
    cell=$2

    tf_chips=${tf_peaks_dir}/${cell}.Finalcommon.accessionlist.withGene.nofiltering.txt
    tf_chip=`awk -v tf=${tf} -v cell=${cell} ' $0 ~ "^"tf"\t" {print $2}' ${tf_chips}` # exact match
    tf_chip_file=${tf_peaks_dir}/${cell}/${tf_chip}.bed.gz

    #echo "$tf $cell ${tf_chip}" # for debugging
    echo "${tf_chip_file}" # return value
}

intersect_tf1_tf2_module_region_worker() {
    tf1=$1
    tf2=$2
    cell=$3
    
    tf1_chip_file=`find_tf_chip ${tf1} ${cell}`
    tf2_chip_file=`find_tf_chip ${tf2} ${cell}`

    # step 1. Make bed file of TF1 regions that overlap TF2 regions (TF1 AND TF2)
    ## Bedtools intersect -a TF1.bed -b TF2.bed
    bedtools intersect -a ${tf1_chip_file} -b ${tf2_chip_file} -u > ${output_dir}/${tf1}_${tf2}_${cell}_intersect.bed

    ## step 2. Subset bed file from #1 by what intersects with Module regions
    ### Bedtools intersect -a result.from.number1.bed -b module.regions.bed > overlap.bed
    #bedtools intersect -a ${output_dir}/${tf1}_${tf2}_intersect.bed -b ${module_regions_dir}/${cell}/Module_${module}.full.bed -u > ${output_dir}/${tf1}_${tf2}_${cell}_${module}_intersect.bed

    # sort in place
    sort -k1,1V -k2,2n -k3,3n ${output_dir}/${tf1}_${tf2}_${cell}_intersect.bed -o  ${output_dir}/${tf1}_${tf2}_${cell}_intersect.bed
}

intersect_tf1_tf2_module_region() {
    tf1=$1
    tf2=$2
    cell=$3
    echo "$idx $tf1 ${cell} in processing..."
    
    if [ -f ${root_dir}/033023/TFs_intersect/${tf1}_${tf2}_${cell}_intersect.bed ];
    then
        cp ${root_dir}/033023/TFs_intersect/${tf1}_${tf2}_${cell}_intersect.bed ${output_dir}
    else
        echo "Non-existing: ${root_dir}/033023/TFs_intersect/${tf1}_${tf2}_${cell}_intersect.bed"
        intersect_tf1_tf2_module_region_worker ${tf1} ${tf2} ${cell}
    fi
}

export -f find_tf_chip
export -f intersect_tf1_tf2_module_region_worker
export -f intersect_tf1_tf2_module_region

intersect_tf1_tf2_modules() {
    idx=0
    for tf1 in ${tfs1[@]}
    do

        tf2=${tfs2[$idx]}
        cell=${cells[$idx]}
        module=${modules[$idx]}


        # step 1. intersect_tf1_tf2_module_region
        echo "intersect_tf1_tf2_module_region ${tf1} ${tf2} ${cell}"
        idx=$((idx+1))
    done | parallel -j 30
}

extract_motif_for_tf_by_motif_name_worker() {

    # Subset file (https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0 ) for rows w/ motif_name
    tf=$1
    motif_name=$2
    echo "extract_motif_for_tf_by_motif_name ${tf} ${motif_name}"

    tf_motif_file=${motif_dir}/${tf}_motifs.bed
    zcat hg38.all_motifs.v1.0.bed.gz | awk -F "\t" -v motif_name=${motif_name} '$4 == motif_name {print $0}' | sort -k1,1V -k2,2n -k3,3n > ${tf_motif_file}
}

extract_motif_for_tf_by_motif_name() {

    ## Subset file (https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0 ) for rows w/ motif_name
    tf=$1
    motif_name=$2
    
    # try to reuse existing ones from previous run
    if [ -f ${root_dir}/033023/TFs_motifs/${tf}_motifs.bed ];
    then
        cp ${root_dir}/033023/TFs_motifs/${tf}_motifs.bed ${motif_dir}
    else
        echo "Non-existing: ${root_dir}/033023/TFs_motifs/${tf}_motifs.bed"

        # calculate from scratch
        extract_motif_for_tf_by_motif_name_worker ${tf} ${motif_name} 
    fi
}

export -f extract_motif_for_tf_by_motif_name_worker
export -f extract_motif_for_tf_by_motif_name

extract_motif_for_tf() {
    # Step 3. Get motifs from the following files https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/motif_annotations.xlsx 
    tf_cluster_id=$1

    # step 3.2: Find the “Name” in archetypes clusters tab (from above link) for that Cluster_ID
    tf_cluster_name=`awk -F "\t" -v cluster=${tf_cluster_id} '$1 == cluster {print $4}' Archetype_clusters.tsv`

    # step 3.3: Subset file (https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0 )for rows that contain the “Name” from above in column #4
    tf_motif_file=${motif_dir}/${tf_cluster_id}_motifs.bed
    echo "zcat hg38.all_motifs.v1.0.bed.gz | awk -F "\t" -v cluster_name=${tf_cluster_name} '$4 == cluster_name {print $0}' > ${tf_motif_file}"
    zcat hg38.all_motifs.v1.0.bed.gz | awk -F "\t" -v cluster_name=${tf_cluster_name} '$4 == cluster_name {print $0}' > ${tf_motif_file}

    # sort
    sort -k1,1V -k2,2n -k3,3n ${tf_motif_file} > tmp_${tf_motif_file}
    mv tmp_${tf_motif_file} ${tf_motif_file}
}


motif_name="" # global variable, return value of choose_motif_name_for_tf
choose_motif_name_for_tf() {
    # Given a TF, decide the motif name
    tf=$1

    # first try Jaspar motif
    num_rows=`grep -P "\t${tf}_" Motifs.tsv | grep "Jaspar" | wc -l`
    if [ ${num_rows} -eq 1 ];
    then
        motif_name="`grep -P "\t${tf}_" Motifs.tsv | grep "Jaspar" | awk '{print $2}'`"
        return
    elif [ ${num_rows} -gt 1 ];
    then
        echo "WARNING: multiple Jasper databases for ${tf}"
        motif_name="`grep -P "\t${tf}_" Motifs.tsv | grep "Jaspar" | head -1 | awk '{print $2}'`"
        return
    fi

    # If Jaspar does not exist, try Hocomoco “HUMAN” motif
    num_rows=`grep -P "\t${tf}_" Motifs.tsv | grep "HOCOMOCO" | grep "HUMAN" | wc -l`
    if [ ${num_rows} -eq 1 ];
    then
        motif_name="`grep -P "\t${tf}_" Motifs.tsv | grep "HOCOMOCO" | grep "HUMAN" | awk '{print $2}'`"
        return
    elif [ ${num_rows} -gt 1 ];
    then
        echo "WARNING: multiple Jasper databases for ${tf}"
        motif_name="`grep -P "\t${tf}_" Motifs.tsv | grep "HOCOMOCO" | grep "HUMAN" | head -1 | awk '{print $2}'`"
        return
    fi

    # Otherwise, any one
    motif_name="`grep -P "\t${tf}_" Motifs.tsv | head -1 | awk '{print $2}'`"

}

parallel_commands_extract_motif_for_tfs_by_motif_name() {
    declare -A tf_motif_names

    for tf2 in `awk -F "\t" 'NR > 1 {print $3}' Tao_Binding_Predictions.tsv | sort -u`
    do
        choose_motif_name_for_tf ${tf2}
        tf_motif_names[${tf2}]=${motif_name}
    done

    # print out the mapping between tf and motif_name
    #for key in ${!tf_motif_names[@]}
    #do
    #    echo "${key} ${tf_motif_names[$key]}"
    #done

    #exit 0

    ## sequential run
    #for key in ${!tf_motif_names[@]}
    #do
    #    #echo "${key} ${tf_motif_names[$key]}"
    #    extract_motif_for_tf_by_motif_name ${key} ${tf_motif_names[$key]}
    #done

    # parallel run
    for key in ${!tf_motif_names[@]}
    do
        echo "extract_motif_for_tf_by_motif_name ${key} ${tf_motif_names[$key]}"
    done | parallel -j 30
}

extract_motif_for_tfs_missed_by_vierstra() {
    # download tf motif files from Factorbook.
    # tf_motifs_files.tsv is manually curated by Tao 
    # See note: # 12/13/22 in https://docs.google.com/document/d/1ubLjP6WGLonEaiOj8WcyfCOJgmBYI5EQgJJmvJvPco4/edit
    #awk 'NR>1 {print $5}' tf_motifs_files.tsv | xargs -P 4 -n 1 -I {} wget -P TFs_motifs_Tao/ {}

    # decompress tf motifs
    #cd TFs_motifs_Tao
    #ls *.gz | xargs gzip -d
    #cd -
    
    ## extract tf1 tf2 combination file for the 15 TFs that Tao found motifs on Factorbook.
    #awk -F "\t" 'BEGIN{ 
    #line=1; 
    #while( getline < "tf_motifs_files.tsv") { if(line > 1) {tfs[$3]=1} ; line += 1; } } 
    #NR==1 {print $0}
    #NR>1 { if($3 in tfs) {print $0} } ' Tao_Binding_Predictions.121022.txt > Tao_Binding_Predictions.1214.byTao.txt 
    ## extract tf2 motif files
    #local tfs=(`awk -F "\t" 'NR > 1 {print $3}' tf_motifs_files.tsv`)
    #local motif_files=(`awk -F "\t" 'NR > 1 {n=split($5, a, "/"); split(a[n], b, ".gz"); print b[1]}' tf_motifs_files.tsv`)
    #local idx=0
    
    # Note: Tao_Binding_Predictions.1214.byTao.txt contains 317 tf combinations that "Viestra" does not have motifs.
    # On 12/14/2022, Tao prepared the motifs for tf2 in Tao_Binding_Predictions.1214.byTao.txt by manually selecting
    # motif file from FactorBook.
    # On 01/09/2023, Shannon emailed Tao her list of selection for these tf2s: only 6 have motifs from Jasper 
    # and 1 (NFE2L1, same as NRF1 in Viestra). For all of them, no need to rerun intersect_tf1_tf2_modules (step 1&2) since all
    # intersected files have been produced.
    # For NFE2L1, run the following 
    #   awk 'NR == 1 || $3 == "NFE2L1" { gsub("NFE2L1", "NRF1"); print $0; }' Tao_Binding_Predictions.1214.byTao.txt > Tao_Binding_Predictions_NFE2L1.byTao.txt
    #   cp Tao_Binding_Predictions_NFE2L1.byTao.txt Tao_Binding_Predictions.tsv
    #   parallel_commands_extract_motif_for_tfs_by_motif_name
    # For the other 6 TF2s, run the following
    #   head -6 tf_motifs_files_shannon.tsv | awk '{print $2}' | xargs -P 4 -n 1 -I {} wget -P TFs_motifs_Tao/ {} ## download motif bed files
    #   awk 'BEGIN{ while(getline < "tf_motifs_files_shannon.tsv") {tfs[$1]=1} } $3 in tfs || NR == 1 { print $0}' Tao_Binding_Predictions.1214.byTao.txt > Tao_Binding_Predictions_011123.tsv
    #   local tfs=(`head -6 tf_motifs_files_shannon.tsv | awk '{print $1}'`)
    #   local motif_files=(`head -6 tf_motifs_files_shannon.tsv | awk '{n=split($2, a, "/"); print a[n]}'`)

    awk 'BEGIN{ while(getline < "tf_motifs_files_shannon.tsv") {tfs[$1]=1} } $3 in tfs || NR == 1 { print $0}' Tao_Binding_Predictions.1214.byTao.txt > Tao_Binding_Predictions_011123.tsv
    local tfs=(`head -6 tf_motifs_files_shannon.tsv | awk '{print $1}'`)
    local motif_files=(`head -6 tf_motifs_files_shannon.tsv | awk '{n=split($2, a, "/"); print a[n]}'`)
    for tf in ${tfs[@]}
    do
        motif_file=${motif_files[$idx]}
        motif_output=${tf}_motifs.bed

        # sort
        echo "sorting ${tf}, ${idx}"
        sort -k1,1V -k2,2n -k3,3n TFs_motifs_Tao/${motif_file} | uniq > TFs_motifs_Tao/t.${tf}
        mv TFs_motifs_Tao/t.${tf} TFs_motifs_Tao/${motif_file}

        # check their chromesomes:
        #fs=(`ls TFs_motifs_Tao`); for f in ${fs[@]}; do awk '{print $1}' TFs_motifs_Tao/$f | sort -k1,1V | uniq; done  | sort | uniq
        
        # produce motif file
        echo "produce tf2 motif ${tf}, ${idx}"
        awk -v f="TFs_motifs_Tao/${motif_file}" -v tf=${tf} -v contigs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY" 'function load_chr(chr) {
            #print chr
            n=0; file="/home/taowang9/Long-read-RNA/hg38/"chr".fa"; 
            while(getline < file) {
                n++; 
                if(substr($1,1,1)==">"){if(n>1){print chr"\t"seq;} chr=substr($0,2,length($0)-1); seq=""; } 
                else {seq=seq""$0;} 
                }; 
            return (seq) 
        }  
        BEGIN {
            OFS="\t"
            split(contigs, keys, " "); 
            for(i=1; i<=length(keys); i++) { key=keys[i]; chrs[key]=load_chr(key) } 
        }
        {
            seq=chrs[$1]

            # motif files are 0-based, but chr is 1-based in awk
            motif=substr(seq, $2+1, $3-$2)

            if (f ~ /.bed/) {print $1, $2, $3, tf, $5, $6, toupper(motif) }
            else {print $1, $2, $3, tf, $5, $4, toupper(motif)}
        }' TFs_motifs_Tao/${motif_file} > TFs_motifs_Tao/t.${tf}

        mv TFs_motifs_Tao/t.${tf} TFs_motifs_Tao/${motif_output}

        idx=$((idx+1))
    done

    cp TFs_motifs_Tao/* TFs_motifs/
}

parallel_commands_extract_motif_for_tfs() {
    declare -A tf_cluster_ids

    idx=0
    for tf2 in ${tfs2[@]}
    do
        idx=$1
        # step 3.1: Find TF name in “Motifs” tab column named “Motif” and identify corresponding Cluster_ID (in above link) column “Cluster_ID” (column A)
        tf_cluster_id=`awk -v tf=${tf2} 'NR > 1 && $2 ~ tf {print $1}' Motifs.tsv | uniq`
        tf_cluster_ids[${tf_cluster_id}]=1
    done

    rm commands
    for tf_cluster_id in ${!tf_cluster_ids[@]};
    do
        # step 2. extract_motif_for_tf
        #echo "extract_motif_for_tf ${tf_cluster_id}"
        echo "bash cmd6.sh ${tf_cluster_id}" >> commands
    done
}

intersect_motif_with_tf1_tf2_module_region_by_motif_name_parallel() {
    tf1=$1
    tf2=$2
    cell=$3
    module=$4

    # figure out motif file
    motif_file=${motif_dir}/${tf2}_motifs.bed
    
    # tf1_tf2-intersected bed file
    tf1_tf2_module_region_file=${output_dir}/${tf1}_${tf2}_${cell}_intersect.bed

    # keep rows in motif file
    echo "bedtools intersect -a ${motif_file} -b ${tf1_tf2_module_region_file} -u > ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed"
    bedtools intersect -a ${motif_file} -b ${tf1_tf2_module_region_file} -u > ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed

    # remove redundant
    awk '{key=$1"_"$2"_"$3; a[key] += 1; if(a[key] == 1) {print $0}  }' ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed > t.${tf1}_${tf2}_${cell}_${module}
    mv t.${tf1}_${tf2}_${cell}_${module} ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed
    
    ## handle overlapping cases w/ merger
    ## pre-merging statistics
    #awk 'NR > 1 { if ($2 < prev) {printf "%s %s\n", prevLine, $0} } { prev=$3; prevLine=$0;} ' ${motif_dir}/${tf1}_${tf2}_${cell}_intersect_${tf2}_motif_deleted.bed | wc -l
    #wc -l ${motif_dir}/${tf1}_${tf2}_${cell}_intersect_${tf2}_motif_deleted.bed

    # prepare the bed file for merging
    awk 'BEGIN{OFS="\t"} { print $1, $2, $3, $6, $7, "" } ' ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed > t.${tf1}_${tf2}_${cell}_${module}
    mv t.${tf1}_${tf2}_${cell}_${module} ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed
    
    # merge overlapping motifs
    merger="/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/Data_from_Tao/merger.py"
    echo "python ${merger} -g reference/hg38.fa -i ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed -o r.${tf1}_${tf2}_${cell}_${module}.bed -d"
    #python ${merger} -g reference/hg38.fa -i ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed -o r.${tf1}_${tf2}_${cell}_${module}.bed -r -s -d # scramble
    #python ${merger} -g reference/hg38.fa -i ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed -o r.${tf1}_${tf2}_${cell}_${module}.bed -r -d # random
    python ${merger} -g reference/hg38.fa -i ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed -o r.${tf1}_${tf2}_${cell}_${module}.bed -d # deletion
    
    # since the range of $5, i.e., [$1, $2), is formed by the left-most boundary and right-most boundary of overlapping motifs,
    # setting $6 to empty is effectively to delete all of them simultaneously. Combining w/ SeqEnumerator's default setting (no
    # reference sequence is used as one of combinatorial factors in the enumeration), all motifs within the same window size are
    # removed simultaneously.
    #echo "${tf1}_${tf2}_${cell}_${module}"
    #source_fileName="${root_dir}/042523/Input_sequences.scramble/S1/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed"
    #dest_fileName="${root_dir}/042523/Input_sequences.deletion/S1/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed"
    #awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, toupper($5), ""}' ${source_fileName} > ${dest_fileName}
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, toupper($5), ""}' ${fileName} > tmp.${fileName}
    mv tmp.${fileName} ${fileName}

    ## post-merging statistics
    #awk 'NR > 1 { if ($2 < prev) {printf "%s %s\n", prevLine, $0} } { prev=$3; prevLine=$0;} ' ${motif_dir}/${tf1}_${tf2}_${cell}_intersect_${tf2}_motif_deleted.bed | wc -l
    #wc -l ${motif_dir}/${tf1}_${tf2}_${cell}_intersect_${tf2}_motif_deleted.bed
}

export -f intersect_motif_with_tf1_tf2_module_region_by_motif_name_parallel

intersect_motif_with_tf1_tf2_module_region_by_motif_name() {
    tf1=$1
    tf2=$2
    cell=$3
    module=$4

    # figure out motif file
    motif_file=${motif_dir}/${tf2}_motifs.bed
    
    # tf1_tf2-intersected bed file
    tf1_tf2_module_region_file=${output_dir}/${tf1}_${tf2}_${cell}_intersect.bed

    # keep rows in motif file
    echo "bedtools intersect -a ${motif_file} -b ${tf1_tf2_module_region_file} -u > ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed"
    bedtools intersect -a ${motif_file} -b ${tf1_tf2_module_region_file} -u > ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed

    # remove redundant
    awk '{key=$1"_"$2"_"$3; a[key] += 1; if(a[key] == 1) {print $0}  }' ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed > t
    mv t ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed
    
    ## handle overlapping cases w/ merger
    ## pre-merging statistics
    #awk 'NR > 1 { if ($2 < prev) {printf "%s %s\n", prevLine, $0} } { prev=$3; prevLine=$0;} ' ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed | wc -l
    #wc -l ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed

    # prepare the bed file for merging
    awk 'BEGIN{OFS="\t"} { print $1, $2, $3, $6, $7, "" } ' ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed > t
    mv t ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed
    
    # merge overlapping motifs
    merger="/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/Data_from_Tao/merger.py"
    python ${merger} -g reference/hg38.fa -i ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed -o r.bed -r
    
    # since the range of $5, i.e., [$1, $2), is formed by the left-most boundary and right-most boundary of overlapping motifs,
    # setting $6 to empty is effectively to delete all of them simultaneously. Combining w/ SeqEnumerator's default setting (no
    # reference sequence is used as one of combinatorial factors in the enumeration), all motifs within the same window size are
    # removed simultaneously.
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, toupper($5), toupper($6)}' r.bed > ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed

    ## post-merging statistics
    #awk 'NR > 1 { if ($2 < prev) {printf "%s %s\n", prevLine, $0} } { prev=$3; prevLine=$0;} ' ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed | wc -l
    #wc -l ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed
}

intersect_motif_with_tf1_tf2_module_region() {
    tf1=$1
    tf2=$2
    cell=$3
    module=$4

    # figure out cluster id
    tf_cluster_id=`awk -v tf=${tf2} 'NR > 1 && $2 ~ tf {print $1}' Motifs.tsv | uniq`

    # tf1_tf2_module-intersected bed file
    tf1_tf2_module_region_file=${output_dir}/${tf1}_${tf2}_${cell}_${module}_intersect.bed

    # motif file
    motif_file=${motif_dir}/${tf_cluster_id}_motifs.bed

    # keep rows in motif file
    echo "bedtools intersect -a ${motif_file} -b ${tf1_tf2_module_region_file} -u > ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed"
    bedtools intersect -a ${motif_file} -b ${tf1_tf2_module_region_file} -u > ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed

    # remove redundant
    awk '{key=$1"_"$2"_"$3; a[key] += 1; if(a[key] == 1) {print $0}  }' ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed > t
    mv t ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed

    # handle overlapping cases
    merger="/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/Data_from_Tao/merger.py"
    awk 'NR > 1 { if ($2 < prev) {printf "%s %s\n", prevLine, $0} } { prev=$3; prevLine=$0;} ' ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed | wc -l
    wc -l ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed
    
    awk 'BEGIN{OFS="\t"} { print $1, $2, $3, $6, $7, "" } ' ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed > t
    mv t ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed
    
    python ${merger} -g reference/hg38.fa -i ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed -o r.bed
    # since the range of $5, i.e., [$1, $2), is formed by the left-most boundary and right-most boundary of overlapping motifs,
    # setting $6 to empty is effectively to delete all of them simultaneously. Combining w/ SeqEnumerator's default setting (no
    # reference sequence is used as one of combinatorial factors in the enumeration, all motifs within the same window size are
    # removed simultaneously.
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, toupper($5), ""}' r.bed > ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed

    awk 'NR > 1 { if ($2 < prev) {printf "%s %s\n", prevLine, $0} } { prev=$3; prevLine=$0;} ' ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed | wc -l
    wc -l ${motif_dir}/${tf_cluster_id}_${tf1}_${tf2}_${cell}_${module}_intersect.bed

}

# Step 4 (parallel): Delete the motif regions from #3 from the bed file obtained in #2
motif_vcf_for_tfs_parallel() {
    idx=0
    for tf1 in ${tfs1[@]}
    do
        tf2=${tfs2[$idx]}
        cell=${cells[$idx]}
        module=${modules[$idx]}

        # intersect motif w/ tf1_tf2_module-intersected bed file
        ## keep rows in motif bed files
        
        # the first implementation identifies motif file by tf cluster id. 
        # It is no longer valid. Now we are using tf name to identify motif file directly
        #intersect_motif_with_tf1_tf2_module_region ${tf1} ${tf2} ${cell} ${module}

        # 2nd implementation
        echo "intersect_motif_with_tf1_tf2_module_region_by_motif_name_parallel ${tf1} ${tf2} ${cell} ${module}"

        idx=$((idx+1))
    done | parallel -j 30 
}

# Step 4: Delete the motif regions from #3 from the bed file obtained in #2
motif_vcf_for_tfs() {
    idx=0
    for tf1 in ${tfs1[@]}
    do
        tf2=${tfs2[$idx]}
        cell=${cells[$idx]}
        module=${modules[$idx]}

        # intersect motif w/ tf1_tf2_module-intersected bed file
        ## keep rows in motif bed files
        
        # the first implementation identifies motif file by tf cluster id. 
        # It is no longer valid. Now we are using tf name to identify motif file directly
        #intersect_motif_with_tf1_tf2_module_region ${tf1} ${tf2} ${cell} ${module}

        # 2nd implementation
        intersect_motif_with_tf1_tf2_module_region_by_motif_name ${tf1} ${tf2} ${cell} ${module}

        idx=$((idx+1))
    done
}

prepare_sequences_for_ref_variant() {
    tf1=$1
    tf2=$2
    cell=$3
    module=$4

    #tsv_file=${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.vcf.tsv
    #echo "${tsv_file}"
    #return

    ## keep rows in motif bed files
    predictor="/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/prepareInput/K562/SeqEnumerator_fast.py"
    reference_genome="/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/prepareInput/reference/hg38.fa"

    i=1
    #vcf_file=${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.vcf
    bed_file=${root_dir}/042523/Input_sequences.deletion/S${i}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.bed
    vcf_file=${root_dir}/042523/Input_sequences.deletion/S${i}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.vcf

    # prepare the vcf file. since each file has < 5k lines, there is no need to split it for parallel prediction.
    awk -F "\t" 'BEGIN{OFS="\t"} { print $1,$2+1, $4,$5,$6}' ${bed_file} > ${vcf_file} 

    echo "$vcf"
    echo "python ${predictor} -g ${reference_genome} -c 1 -s FORWARD -wd 1056 -wu 1057 -q BOTH --independent no -i ${vcf_file} -o ${motif_dir}/ins_containing.fasta -r report.txt"
    python ${predictor} -g ${reference_genome} -c 1 -s FORWARD -wd 1056 -wu 1057 -q BOTH --independent no -i ${vcf_file} -o ${motif_dir}/ins_containing.fasta -r report.txt
}

export -f prepare_sequences_for_ref_variant

tf_model_prediction_parallel_work_multi_nodes() {
    idx=0
    job_submitted=0
    rm 042523/slurm_outputs/*
    
    for tf1 in ${tfs1[@]}
    do
        tf2=${tfs2[$idx]}
        cell=${cells[$idx]}
        module=${modules[$idx]}

        #tf_cluster_id=`awk -v tf=${tf2} 'NR > 1 && $2 ~ tf {print $1}' Motifs.tsv | uniq`
        echo "Submitted: ${tf1} ${tf2} ${cell} ${module}"
        sbatch --wait batch_prepare_sequences_for_ref_variant_deletion.sh ${tf1} ${tf2} ${cell} ${module} &

        pids[${idx}]=$!
        job_submitted=`squeue -u taowang9 | wc -l`
        while [ ${job_submitted} -gt 50 ];
        do
            job_submitted=`squeue -u taowang9 | wc -l`
        done

        idx=$((idx+1))
    done

    wait # wait for all child processes
}

tf_model_prediction_parallel_work_one_node() {
    idx=0
    for tf1 in ${tfs1[@]}
    do
        tf2=${tfs2[$idx]}
        cell=${cells[$idx]}
        module=${modules[$idx]}
        #tf_cluster_id=`awk -v tf=${tf2} 'NR > 1 && $2 ~ tf {print $1}' Motifs.tsv | uniq`
        
        echo "prepare_sequences_for_ref_variant ${tf1} ${tf2} ${cell} ${module}"

        idx=$((idx+1))
    done | parallel -k -j 30
}

tf_model_prediction_parallel() {
    # single-node parallelization
    #tf_model_prediction_parallel_work_one_node

    ## multi-node parallelization
    tf_model_prediction_parallel_work_multi_nodes

    idx=0
    prefix_sherlock="/home/groups/mpsnyder/taowang9/ChIP/tsv_files"
    rm files
    for tf1 in ${tfs1[@]}
    do
        tf2=${tfs2[$idx]}
        cell=${cells[$idx]}
        module=${modules[$idx]}

        ## produce prediction commands
        # prediction is run on sherlock
        # K562.all and HepG2.all on sherlock was produced on SCG with /home/taowang9/pulsarpy-to-encodedcc/pulsarpy_to_encodedcc/scripts/get_target_from_ENCS.py
        # the input was from Shannon's email
        # select the model according to cell line and TF
        count=`grep "\b${tf1}\b" ${targets_dir}/targets.txt.${cell} | wc -l`
        if [ ${count} -eq 0 ];
        then
            echo "${tf1} does not have TF binding model"
            idx=$((idx+1))
            continue
        else
            i=1
            grep "\b${tf1}\b" ${targets_dir}/targets.txt.${cell} | awk -v module=${module} -v cell=${cell} -v tf2=${tf2} -v file="${prefix_sherlock}/S${i}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.vcf.tsv" '{print module"_"cell"_"$2"_"tf2"\t"$1"\t"file}' >> files
            
            idx=$((idx+1))
        fi
    done

    # filter out tsv files with no contents
    #cat files | awk -F "/" '{print "wc -l TFs_motifs/"$NF}' | bash | awk '{print $1}' > rs
    cat files | awk -F "/" '{print "042523/Input_sequences.deletion/"$(NF-1)"/"$NF}' | parallel -k -j 20 wc -l | awk '{print $1}' > rs
    paste -d "\t" rs files | awk '$1 > 1 {print $2"\t"$3"\t"$4}' > files_with_contents
    mv files_with_contents files
}

# Step 5: Plot TF1 signal at motif-deleted regions and motif-containing regions (violin plot similar to the one you just created for ref and var)
tf_model_prediction() {
    idx=0
    for tf1 in ${tfs1[@]}
    do
        tf2=${tfs2[$idx]}
        cell=${cells[$idx]}
        module=${modules[$idx]}
        #tf_cluster_id=`awk -v tf=${tf2} 'NR > 1 && $2 ~ tf {print $1}' Motifs.tsv | uniq`

        # prepare the vcf file. since each file has < 5k lines, there is no need to split it for parallel prediction.
        awk -F "\t" 'BEGIN{OFS="\t"} { print $1,$2+1, $4,$5,$6}' ${motif_dir}/${tf1}_${tf2}_${cell}_intersect_${tf2}_motif_deleted.bed > ${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.vcf
        
        ## keep rows in motif bed files
        predictor="/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/prepareInput/K562/SeqEnumerator_fast.py"
        reference_genome="/oak/stanford/scg/prj_ENCODE/ChIP/basepairmodel/basepairmodels/tutorial/prepareInput/reference/hg38.fa"
        vcf_file=${motif_dir}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.vcf
        echo "python ${predictor} -g ${reference_genome} -c 1 -s FORWARD -wd 1056 -wu 1057 -q BOTH --independent no -i ${vcf_file} -o ins_containing.fasta -r report.txt"
        break

        ## produce prediction commands
        # prediction is run on sherlock
        # K562.all and HepG2.all on sherlock was produced on SCG with /home/taowang9/pulsarpy-to-encodedcc/pulsarpy_to_encodedcc/scripts/get_target_from_ENCS.py
        # the input was from Shannon's email
        
        # select the model according to cell line and TF
        prefix_sherlock="/home/groups/mpsnyder/taowang9/ChIP/tsv_files"
        grep "\b${tf1}\b" ${targets_dir}/targets.txt.${cell} | awk -v module=${module} -v cell=${cell} -v tf2=${tf2} -v file="${prefix_sherlock}/${tf1}_${tf2}_${cell}_${module}_intersect_${tf2}_motif_deleted.vcf.tsv" '{print module"_"cell"_"$2"_"tf2"\t"$1"\t"file}' >> files

        idx=$((idx+1))
    done
}

extract_ref_variant_read_counts_simplified() {
    #cd ${root_dir}/010523
    cd ${root_dir}/042523

    mkdir ref_variant_reads
    for d in `ls prediction/`;
    do
        if [ -f prediction/${d}/prediction.tsv ];
        then
            awk -F "\t" '{print $2"\t"$3}' prediction/${d}/prediction.tsv > ref_variant_reads/${d}.tsv
        else
            echo "${d}/prediction.tsv does not exist or have contents"
        fi
    done
}

extract_ref_variant_read_counts() {
    cd ${root_dir}/120822

    idx=0
    for tf1 in ${tfs1[@]}
    do
        tf2=${tfs2[$idx]}
        cell=${cells[$idx]}
        module=${modules[$idx]}

        result_tsv_file=${module}_${cell}_${tf1}_${tf2}/prediction.tsv

        awk -F "\t" '{print $2"\t"$3}' ${result_tsv_file} > ref_variant_reads/${module}_${cell}_${tf1}_${tf2}.tsv
        
        idx=$((idx+1))
    done

    cd -
    
}

## step 1 & 2. intersect tf1, tf2
intersect_tf1_tf2_modules ## 042523, step 1 and 2.

## step 3. produce cluster-specific motif files
parallel_commands_extract_motif_for_tfs_by_motif_name ## 042523, step 3, extract motif sequences by motif names


## Step 4. remove tf2's motif from tf1 peak region and prepare sequence files for model prediction.
motif_vcf_for_tfs_parallel ## 042523, step 4. (parallel)

# Step 5: Plot TF1 signal at motif-deleted regions and motif-containing regions (violin plot similar to the one you just created for ref and var)
#tf_model_prediction
tf_model_prediction_parallel ## 042523, step 5.

# Step 6. parallel run TF models on Sherlock
# to run the tf models, copy files, and output files of SeqEnumerator_fast.py to sherlock /home/groups/mpsnyder/taowang9/ChIP/tsv_files
# scp files TFs_motifs/*_motif_deleted.vcf.tsv taowang9@login.sherlock.stanford.edu:/home/users/taowang9/ChIP/tsv_files/
#scp files 033023/TFs_motifs/*_motif_deleted.vcf.tsv taowang9@login.sherlock.stanford.edu:/home/users/taowang9/ChIP/tsv_files/ ## 033023, step 6.

## scp is very slow and cannot skip transferring those transferred
## use rsync instead ## 042523
for i in {1..5}; do  rsync -avP --ignore-existing 042523/Input_sequences/S${i}/*.tsv taowang9@login.sherlock.stanford.edu:/scratch/users/taowang9/042523/tsv_files/S${i}/ ; done ## 042523

## then run bash driver_120722.sh
## the outputs are under /scratch/users/taowang9/120722/files

# Step 7. download outputs from sherlock to laptop and use R for plotting
#mkdir 120822 & cd 120822
#scp -r taowang9@login.sherlock.stanford.edu:/scratch/users/taowang9/120722/files/* .
