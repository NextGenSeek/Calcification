for f1 in *.fastq.gz
do
       echo -n "quality control $f1"
       fastqc $f1
       echo " Done"
done

multiqc .

for f1 in *_1.fastq
do
        echo -n "Quality trimming and trimming adapters... $f1 "
        f2=${f1%%_1.fastq}"_2.fastq"
            ~/TrimGalore-0.6.6/trim_galore --paired $f1 $f2 --cores 8 --dont_gzip
        echo " Done"
done

for f1 in *_val_1.fq
do
        echo -n "mapping to human transcriptome $f1"
        f2=${f1%%_val_1.fq}"_val_2.fq"
        salmon quant --threads 8 -i /PATH/TO/Hsapiens_index -l A -1 $f1 -2 $f2  --validateMappings --minScoreFraction 0.5 -o quants/$f1
        echo " Done"
done
