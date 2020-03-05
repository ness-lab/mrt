inpath=$1
outpath=$2

mkdir -p ${outpath}

find ${inpath} -name "*.fasta" -print0 | while read -d $'\0' fasta_file; do
#    echo $fasta_file
    base=$(basename $fasta_file .fasta)
    echo "Currently on ${base}"

    time ldhelmet find_confs \
    --num_threads 50 \
    --window_size 50 \
    --output_file ${outpath}/${base}.conf ${fasta_file}

    sleep 1

    echo "table_gen for ${base}"

    time ldhelmet table_gen \
    --num_threads 50 \
    --conf_file ${outpath}/${base}.conf \
    --theta 0.001 \
    --rhos 0.0 0.1 10.0 1.0 100.0 \
    --output_file ${outpath}/${base}.lk > table_gen 2> table_gen2
    
    sleep 1

    rm table_gen*

    sleep 1

    time ldhelmet pade \
    --num_threads 50 \
    --conf_file ${outpath}/${base}.conf \
    --theta 0.001 \
    --output_file ${outpath}/${base}.pade

    sleep 1

    time ldhelmet rjmcmc \
    --num_threads 30 \
    --window_size 50 \
    --seq_file ${fasta_file} \
    --lk_file ${outpath}/${base}.lk \
    --pade_file ${outpath}/${base}.pade \
    --num_iter 1000000 \
    --burn_in 100000 \
    --block_penalty 50 \
    --mut_mat_file ../../resources/mut_mat.txt \
    --output_file ${outpath}/${base}.post

    sleep 1

    time ldhelmet post_to_text \
    --mean \
    --perc 0.025 \
    --perc 0.50 \
    --perc 0.975 \
    --output_file ${outpath}/${base}.txt \
    data/ldhelmet/${outpath}/${base}.post

    sleep 3

    echo "Removing temp files..."
    rm ${outpath}/${base}.conf
    rm ${outpath}/${base}.lk
    rm ${outpath}/${base}.pade
    rm ${outpath}/${base}.post

done 

