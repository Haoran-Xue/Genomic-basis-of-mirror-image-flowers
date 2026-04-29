#Wachendorfia paniculata RNA-seq mapping

for f1 in /mnt/barn/hxue/2024-05-06-Wpanic_stranded/01.RawDataLaneMerged/Wp_*_1_val_1.fq.gz

do

f2=$(echo $f1 | sed 's/1_val_1.fq.gz/2_val_2.fq.gz/')

f1u=$(echo $f1 | sed 's/\.gz//')

f2u=$(echo $f2 | sed 's/\.gz//')

sample=$(echo $f1 | sed 's/.*\///' | sed 's/_1_val_1\.fq\.gz//')

echo $sample $f1 $f2 $f1u $f2u

gunzip $f1

gunzip $f2

nice -n 10 /mnt/barn/hxue/programs/hisat2-2.2.1/hisat2 -x /mnt/barn/hxue/2022-11-29-Wpaniculata-RNA-seq/3.HISAT2_alignments/Wp002.hifi -p 10 -1 $f1u -2 $f2u --rna-strandness RF | samtools view -S -b | samtools sort > $sample.HISAT2.Wp002.sorted.bam

samtools index $sample.HISAT2.Wp002.sorted.bam

nohup gzip $f1u &

nohup gzip $f2u &

done



#Wachendorfia multiflora RNA-seq mapping

for f1 in /mnt/barn/hxue/2023-10-10-Wmultiflora-RNA-seq/X204SC23084391-Z01-F001/01.RawData/*val_1.fq.gz

do

f2=$(echo $f1 | sed 's/1_val_1.fq.gz/2_val_2.fq.gz/')

f1u=$(echo $f1 | sed 's/\.gz//')

f2u=$(echo $f2 | sed 's/\.gz//')

sample=$(echo $f1 | sed 's/.*\///' | sed 's/_1_val_1\.fq\.gz//')

echo $sample $f1 $f2 $f1u $f2u

gunzip $f1

gunzip $f2

nice -n 10 /mnt/barn/hxue/programs/hisat2-2.2.1/hisat2 -x /mnt/barn/hxue/2022-11-29-Wpaniculata-RNA-seq/3.HISAT2_alignments/Wp002.hifi -p 10 -1 $f1u -2 $f2u --rna-strandness RF | samtools view -S -b | samtools sort > $sample.HISAT2.Wp002.sorted.bam

samtools index $sample.HISAT2.Wp002.sorted.bam

nohup gzip $f1u &

nohup gzip $f2u &

done



#Barberetta aurea RNA-seq mapping

for f1 in /mnt/barn/hxue/Barberetta-RNA-seq/01.RawData/*_val_1.fq.gz

do

f2=$(echo $f1 | sed 's/1.fq.gz/2.fq.gz/')

sample=$(echo $f1 | sed 's/.*\///' | sed 's/_1\.fq\.gz//')

echo $sample $f1 $f2

nice -n 10 /mnt/barn/hxue/programs/hisat2-2.2.1/hisat2 -x /mnt/barn/hxue/Ba_Ngeli_ONT/1.genome_assembly/Ba_Ngeli_ONT -p 20 -1 $f1 -2 $f2 | samtools view -Sb - | samtools sort -o Ba.$sample.HISAT2.BaN.bam

samtools index Ba.$sample.HISAT2.BaN.bam


done



#StringTie transcript assembly

for f in *.bam
	
do

sample=$(echo $f | sed 's/.*\///' | sed 's/\.bam//')

nice -n 10 nohup /mnt/barn/hxue/programs/stringtie-2.2.1.Linux_x86_64/stringtie $f -p 1 -o $sample.stringtie.gtf &

done



#Merging transcript assemblies of RNA-seq samples

nice -n 10 nohup /mnt/barn/hxue/programs/stringtie-2.2.1.Linux_x86_64/stringtie --merge -o Wp002.WmWp.ST.gtf -F 0 -T 0 -l Wp002WmWpST W*.gtf > stringtie_merge_Wp002WmWp.out &

nice -n 10 nohup /mnt/barn/hxue/programs/stringtie-2.2.1.Linux_x86_64/stringtie --merge -o BaN.ST.gtf -F 0 -T 0 -l BaNST *BaN.stringtie.gtf > stringtie_merge_BaN.out &



#Wachendorfia Salmon RNA-seq read quantification

for f1 in /mnt/barn/hxue/2024-05-06-Wpanic_stranded/01.RawDataLaneMerged/*val_1.fq.gz
do
        f2=$(echo $f1 | sed 's/1_val_1.fq.gz/2_val_2.fq.gz/')

        sample=$(echo $f1 | sed 's/.*\///' | sed 's/_1_val_1\.fq\.gz//')

        echo $sample $f1 $f2

nohup nice -n 10 /mnt/barn/hxue/programs/salmon-1.9.0_linux_x86_64/bin/salmon quant --gcBias -i /mnt/barn/hxue/Wachendorfia_RNA-Seq_analyses/Wp002WmWpStrandedStringtie.transcriptome -l IU \
	-1 $f1 \
	-2 $f2 \
	-p 1 --validateMappings -o quants/$sample\_quant > $sample.salmon.log &

done



#Barberetta Salmon RNA-seq read quantification

for f1 in  /mnt/barn/hxue/sequencing_data/Barberetta-RNA-seq/X208SC24104850-Z01-F014/01.RawData/*_R_val_1.fq.gz
do
        f2=$(echo $f1 | sed 's/_1.fq.gz/_2.fq.gz/')

        sample=$(echo $f1 | sed 's/.*\///' | sed 's/_1\.fq\.gz//')

        echo $sample $f1 $f2

nice -n 10 /mnt/barn/hxue/programs/salmon-1.9.0_linux_x86_64/bin/salmon quant --gcBias -i ./BaN.ST -l IU \
	-1 $f1 \
	-2 $f2 \
	-p 40 --validateMappings -o quants/$sample\_quant > $sample.salmon.log

done
