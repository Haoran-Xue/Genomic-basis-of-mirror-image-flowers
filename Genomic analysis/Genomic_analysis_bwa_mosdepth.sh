## Align Illumina DNA-seq reads to the genome assembly of W. paniculata

bwa index /mnt/barn/hxue/2022-11-24-Wpaniculata-pacbio/1.genome_assembly/Wp002.hifi.bp.p_ctg.fa

for f1 in /mnt/barn/hxue/sequencing_data/W*_1.fq.gz; do
        f2=$(echo $f1 | sed 's/_1\.f/_2.f/')
        f=$(echo $f1 | sed 's/.*\///' )
        echo $f $f1 $f2
        nice -n 10 /mnt/barn/hxue/programs/bwa/bwa mem -Y -t 30 -R  "@RG\tID:$f\tSM:$f\tPL:ILLUMINA" /mnt/barn/hxue/2022-11-24-Wpaniculata-pacbio/1.genome_assembly/Wp002.hifi.bp.p_ctg.fa $f1 $f2 | samtools view -b -h - | samtools sort -o bwa/$f.bwa.Wp002.bam

samtools index bwa/$f.bwa.Wp002.bam

fRM="bwa/${f}.bwa.Wp002.RM.bam"
fR="bwa/${f}.bwa.Wp002.Repeats.bam"
nohup nice -n 10 samtools view bwa/$f.bwa.Wp002.bam -b -h -L /mnt/barn/hxue/2022-11-24-Wpaniculata-pacbio/1.genome_assembly/Wp002.hifi.bp.p_ctg.RM.bed -U $fRM -o $fR &
samtools index $fRM
rm $fR

done



## Align Illumina DNA-seq reads to the genome assembly of W. thyrsiflora

bwa index /mnt/barn/hxue/2023-07-21-Wachendorfia-thyrsiflora-genomes-PacBio/1.genome_assembly/Wt302R.hifi.bp.p_ctg.fa

for f1 in /mnt/barn/hxue/Wthyrsiflora-pool-seq/W*_1.fq.gz; do
        f2=$(echo $f1 | sed 's/_1/_2/')
        fr=$(echo $f1 | sed 's/.*\///' | sed 's/_1.*//' | sed 's/_2.*//')
        echo $f $f1 $f2

        bwa mem -t 44 /mnt/barn/hxue/2023-07-21-Wachendorfia-thyrsiflora-genomes-PacBio/1.genome_assembly/Wt302R.hifi.bp.p_ctg.fa $f1 $f2 | samtools view -b -h - | samtools sort -o bwa/$f.bwa.Wt302R.bam
        samtools index bwa/$f.bwa.Wt302R.bam

fRM="bwa/${f}.bwa.Wt302R.RM.bam"
fR="bwa/${f}.bwa.Wt302R.Repeats.bam"
nohup nice -n 10 samtools view bwa/$f.bwa.Wp002.bam -b -h -L /mnt/barn/hxue/2023-07-21-Wachendorfia-thyrsiflora-genomes-PacBio/1.genome_assembly/Wt302R.hifi.bp.p_ctg.RM.bed -U $fRM -o $fR &
samtools index $fRM
rm $fR

done



## Align Illumina DNA-seq reads to the genome assembly of Barberetta aurea

for f1 in /mnt/barn/hxue/sequencing_data/2023-03-27-Barberetta-pool-seq/X204SC23030985-Z01-F001/01.RawData/*/*_1.fq.gz; do


        f2=$(echo $f1 | sed 's/_1/_2/')
        fr=$(echo $f1 | sed 's/.*\///' | sed 's/_1.*//' )
        echo $fr $f1 $f2
        /mnt/barn/hxue/programs/bwa/bwa mem -Y -t 40 -R  "@RG\tID:$fr\tSM:$fr\tPL:ILLUMINA" /mnt/barn/hxue/Ba_Ngeli_ONT/1.genome_assembly/Ba_Ngeli_ONT.bp.p_ctg.fa $f1 $f2 | samtools view -b -h - | samtools sort -o $fr.bwa.BaN_ONT.bam

samtools index $fr.bwa.BaN_ONT.bam

fRM="bwa/${f}.bwa.BaN_ONT.RM.bam"
fR="bwa/${f}.bwa.BaN_ONT.Repeats.bam"
nohup nice -n 10 samtools view bwa/$f.bwa.BaN_ONT.bam -b -h -L /mnt/barn/hxue/Ba_Ngeli_ONT/1.genome_assembly/Ba_Ngeli_ONT.bp.p_ctg.RM.bed -U $fRM -o $fR &
samtools index $fRM
rm $fR

done



## Read coverage quantification for 50-kb genomic windows (all species)

for f in bwa/*.RM.bam; do
fr=$(echo $f | sed 's/bam/50kb/' | sed 's/\.\././')
nohup /mnt/barn/hxue/programs/mosdepth -t 4 -b 50000 -n -x $fr $f & 
done



## Read coverage quantification for protein coding genes in Wachendorfia species

for bamf in bwa/*.bwa.Wp002.RM.bam; do
        depthf=$(echo $bamf | sed 's/bam/Wp002tx09.gene/')
        depthf2=$(echo $bamf | sed 's/bam/E-locus.gene/')


nohup /mnt/barn/hxue/programs/mosdepth -b /mnt/barn/hxue/2022-11-24-Wpaniculata-pacbio/1.genome_assembly/braker3/Wp002.750m.RM.braker3WmWpStrRNASeq.Wp002tx09/Wp002.750m.RM.braker3WmWpStrRNASeq.Wp002tx09.gene.bed $depthf $bamf &

nohup /mnt/barn/hxue/programs/mosdepth -b E-locus.bed $depthf2 $bamf &

done



## Read coverage quantification for protein coding genes in Barberetta aurea

for bamf in bwa/*.bwa.BaN_ONT.RM.bam; do
        depthf=$(echo $bamf | sed 's/bam/BaNtx1.gene/')
        depthf2=$(echo $bamf | sed 's/bam/R-locus.gene/')

/mnt/barn/hxue/programs/mosdepth -t 10 -b /mnt/barn/hxue/Ba_Ngeli_ONT/1.genome_assembly/braker3/BaN.RM.braker3.BaRNARM.run1/BaN.RM.braker3.BaRNA.RM.gene.bed $depthf $bamf -n -x

/mnt/barn/hxue/programs/mosdepth -t 10 -b BaN.R-locus.bed $depthf2 $bamf -n -x

done
