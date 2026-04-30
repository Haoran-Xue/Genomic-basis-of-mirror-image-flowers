----------------------------------
| Wachendorfia paniculata genome |
----------------------------------

## genome assembly
nohup /path/programs/hifiasm/hifiasm -t 40 -o Wp --hg-size 750m /path/Wp.fastq.gz > Wp.out &
awk '/^S/{print ">"$2"\n"$3}' Wp.bp.p_ctg.gfa | fold > Wp.bp.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' Wp.bp.hap1.p_ctg.gfa | fold > Wp.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' Wp.bp.hap2.p_ctg.gfa | fold > Wp.bp.hap2.p_ctg.fa

## BUSCO assessment
nice -n 10 nohup busco -i Wp.bp.p_ctg.fa -l embryophyta_odb10 -o Wp.bp.p_ctg.busco.embryophyta_odb10 -c 40 -m genome > Wp.bp.p_ctg.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i Wp.bp.hap1.p_ctg.fa -l embryophyta_odb10 -o Wp.bp.hap1.p_ctg.busco.embryophyta_odb10 -c 40 -m genome > Wp.bp.hap1.p_ctg.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i Wp.bp.hap2.p_ctg.fa -l embryophyta_odb10 -o Wp.bp.hap2.p_ctg.busco.embryophyta_odb10 -c 40 -m genome > Wp.bp.hap2.p_ctg.busco.embryophyta_odb10.out &

## QUAST assessment
nohup /path/programs/quast/quast.py Wp.bp.p_ctg.fa -o quast.Wp.bp.p_ctg.out &
nohup /path/programs/quast/quast.py Wp.bp.hap1.p_ctg.fa -o quast.Wp.bp.hap1.p_ctg.out &
nohup /path/programs/quast/quast.py Wp.bp.hap2.p_ctg.fa -o quast.Wp.bp.hap2.p_ctg.out &

## RepeatMasker Step 1: BuildDatabase
nice -n 10 /path/programs/RepeatModeler/BuildDatabase -name Wp ../Wp.bp.p_ctg.fa > Ba.BuildDatabase.log

## RepeatMasker Step 2: RepeatModeler
nice -n 10 /path/programs/RepeatModeler/RepeatModeler \
 -engine ncbi -threads 10 -database Wp > Wp.RepeatModeler.out

## RepeatMasker Step 3: RepeatMasker
mkdir Wp.hardMasked
nice -n 10 /path/programs/RepeatMasker/RepeatMasker -pa 4 -no_is -dir ./Wp.hardMasked/ -lib consensi.fa.classified ../Wp_Ngeli_ONT.bp.p_ctg.fa > Wp.RMhard.out

## Run RepeatMasker for Haplotypes 1 and 2 as well

## BRAKER3 gene annotation
nohup braker.pl --genome=/path/Wp.bp.p_ctg.fa --prot_seq=Viridiplantae.odb11.for_braker.fa --bam=/path/HISAT2.Wp.merged.bam --workingdir=Wp.RM.braker3WmWp --threads 40 --verbosity=4 > Wp.RM.braker3WmWp.run2.log &
nohup braker.pl --genome=/path/Wp.bp.hap1.p_ctg.fa --prot_seq=Viridiplantae.odb11.for_braker.fa --bam=/path/HISAT2.Wp.hap1.merged.bam --workingdir=Wp.hap1.RM.braker3WmWp --threads 40 --verbosity=4 > Wp.hap1.RM.braker3WmWp.run2.log &
nohup braker.pl --genome=/path/Wp.bp.hap2.p_ctg.fa --prot_seq=Viridiplantae.odb11.for_braker.fa --bam=/path/HISAT2.Wp.hap2.merged.bam --workingdir=Wp.hap2.RM.braker3WmWp --threads 40 --verbosity=4 > Wp.hap2.RM.braker3WmWp.run2.log &


----------------------------------
| Wachendorfia thyrsiflora genome |
----------------------------------

## genome assembly
nohup /path/programs/hifiasm/hifiasm -t 40 -o Wt --hg-size 700m /path/Wt.fastq.gz > Wt.out &
awk '/^S/{print ">"$2"\n"$3}' Wt.bp.p_ctg.gfa | fold > Wt.bp.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' Wt.bp.hap1.p_ctg.gfa | fold > Wt.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' Wt.bp.hap2.p_ctg.gfa | fold > Wt.bp.hap2.p_ctg.fa

## BUSCO assessment
nice -n 10 nohup busco -i Wt.bp.p_ctg.fa -l /path/busco_downloads/lineages/embryophyta_odb10 -o Wt.bp.p_ctg.busco.embryophyta_odb10 -c 40 -m genome -f > Wt.bp.p_ctg.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i Wt.bp.hap1.p_ctg.fa -l /path/busco_downloads/lineages/embryophyta_odb10 -o Wt.bp.hap1.p_ctg.busco.embryophyta_odb10 -c 40 -m genome -f > Wt.bp.hap1.p_ctg.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i Wt.bp.hap2.p_ctg.fa -l /path/busco_downloads/lineages/embryophyta_odb10 -o Wt.bp.hap2.p_ctg.busco.embryophyta_odb10 -c 40 -m genome -f > Wt.bp.hap2.p_ctg.busco.embryophyta_odb10.out &

## QUAST assessment
nohup /path/programs/quast/quast.py Wt.bp.p_ctg.fa -o quast.Wt.bp.p_ctg.out &
nohup /path/programs/quast/quast.py Wt.bp.hap1.p_ctg.fa -o quast.Wt.bp.hap1.p_ctg.out &
nohup /path/programs/quast/quast.py Wt.bp.hap2.p_ctg.fa -o quast.Wt.bp.hap2.p_ctg.out &

## RepeatMasker Step 1: BuildDatabase
nice -n 10 /path/programs/RepeatModeler/BuildDatabase -name Wt ../Wp.bp.p_ctg.fa > Wt.BuildDatabase.log

## RepeatMasker Step 2: RepeatModeler
nice -n 10 /path/programs/RepeatModeler/RepeatModeler \
 -engine ncbi -threads 10 -database Wt > Wt.RepeatModeler.out

## RepeatMasker Step 3: RepeatMasker
mkdir Wt.hardMasked
nice -n 10 /path/programs/RepeatMasker/RepeatMasker -pa 4 -no_is -dir ./Wt.hardMasked/ -lib consensi.fa.classified ../Wt_Ngeli_ONT.bp.p_ctg.fa > Wt.RMhard.out

## Run RepeatMasker for Haplotypes 1 and 2 as well

## BRAKER3 gene annotation
nohup nice -n 10 braker.pl --genome=/path/Wt.bp.p_ctg.fa.masked --prot_seq=/mnt/barn/hxue/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=Wt --threads 40 --verbosity=4 > Wt.RM.braker3.log &
nohup nice -n 10 braker.pl --genome=/path/Wt.bp.hap1.p_ctg.fa.masked --prot_seq=/mnt/barn/hxue/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=Wt.hap1.RM --threads 20 --verbosity=4 > Wt.hap1.RM.braker3.log &
nohup nice -n 10 braker.pl --genome=/path/Wt.bp.hap2.p_ctg.fa.masked --prot_seq=/mnt/barn/hxue/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=Wt.hap2.RM --threads 20 --verbosity=4 > Wt.hap2.RM.braker3.log &



---------------------------------
| Barberetta aurea Ngeli genome |
---------------------------------

## genome assembly
nohup /path/programs/hifiasm-0.25.0/hifiasm --ont --chem-c 1 --chem-f 256 --rl-cut 1000 --sc-cut 10 -o Ba_Ngeli_ONT -t 40 /path/sequencing_data/2025-12-11-Baurea_Ngeli_ONT/barberetta_aurea_fastq_pass_con.fastq.gz > Ba_Ngeli_ONT.out &
awk '/^S/{print ">"$2"\n"$3}' Ba_Ngeli_ONT.bp.p_ctg.gfa | fold > Ba_Ngeli_ONT.bp.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' Ba_Ngeli_ONT.bp.hap1.p_ctg.gfa | fold > Ba_Ngeli_ONT.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' Ba_Ngeli_ONT.bp.hap2.p_ctg.gfa | fold > Ba_Ngeli_ONT.bp.hap2.p_ctg.fa

## BUSCO assessment
nice -n 10 nohup busco -i Ba_Ngeli_ONT.bp.p_ctg.fa -l embryophyta_odb10 -o Ba_Ngeli_ONT.busco.embryophyta_odb10 -c 40 -m genome -f > Ba_Ngeli_ONT.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i Ba_Ngeli_ONT.bp.hap1.p_ctg.fa -l embryophyta_odb10 -o Ba_Ngeli_ONT.hap1.busco.embryophyta_odb10 -c 40 -m genome -f > Ba_Ngeli_ONT.hap1.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i Ba_Ngeli_ONT.bp.hap2.p_ctg.fa -l embryophyta_odb10 -o Ba_Ngeli_ONT.hap2.busco.embryophyta_odb10 -c 40 -m genome -f > Ba_Ngeli_ONT.hap2.busco.embryophyta_odb10.out &

## QUAST assessment
nohup /path/programs/quast/quast.py Ba_Ngeli_ONT.bp.p_ctg.fa -o quast.Ba_Ngeli_ONT.bp.p_ctg.out &
nohup /path/programs/quast/quast.py Ba_Ngeli_ONT.bp.hap1.p_ctg.fa -o quast.Ba_Ngeli_ONT.bp.hap1.p_ctg.out &
nohup /path/programs/quast/quast.py Ba_Ngeli_ONT.bp.hap2.p_ctg.fa -o quast.Ba_Ngeli_ONT.bp.hap2.p_ctg.out &

## RepeatMasker Step 1: BuildDatabase
nice -n 10 /path/programs/RepeatModeler/BuildDatabase -name Ba ../Ba_Ngeli_ONT.bp.p_ctg.fa > Ba.BuildDatabase.log

## RepeatMasker Step 2: RepeatModeler
nice -n 10 /path/programs/RepeatModeler/RepeatModeler \
 -engine ncbi -threads 10 -database Ba > Ba.RepeatModeler.out

## RepeatMasker Step 3: RepeatMasker
mkdir Ba.hardMasked
nice -n 10 /path/programs/RepeatMasker/RepeatMasker -pa 4 -no_is -dir ./Ba.hardMasked/ -lib consensi.fa.classified ../Ba_Ngeli_ONT.bp.p_ctg.fa > Ba.RMhard.out

## Run RepeatMasker for Haplotypes 1 and 2 as well

## BRAKER3 gene annotation
nice -n 10 nohup braker.pl --genome=/path/Ba_Ngeli_ONT/Ba_Ngeli_ONT.bp.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=BaN.RM.braker3 --threads 40 --verbosity=4 > BaN.RM.braker3.run1.log &
nice -n 10 nohup braker.pl --genome=/path/Ba_Ngeli_ONT/Ba_Ngeli_ONT.bp.hap1.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=BaN.hap1.RM.braker3 --threads 20 --verbosity=4 > BaN.hap1.RM.braker3.run1.log &
nice -n 10 nohup braker.pl --genome=/path/Ba_Ngeli_ONT/Ba_Ngeli_ONT.bp.hap2.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=BaN.hap2.RM.braker3 --threads 20 --verbosity=4 > BaN.hap2.RM.braker3.run1.log &



------------------------------------
| Barberetta aurea Karkloof genome |
------------------------------------

## genome assembly
nohup /path/programs/hifiasm-0.25.0/hifiasm -o BaK -t 40 BaK.fastq.gz > BaK.out &
awk '/^S/{print ">"$2"\n"$3}' BaK.bp.p_ctg.gfa | fold > BaK.bp.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' BaK.bp.hap1.p_ctg.gfa | fold > BaK.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' BaK.bp.hap2.p_ctg.gfa | fold > BaK.bp.hap2.p_ctg.fa

## BUSCO assessment
nice -n 10 nohup busco -i BaK.bp.p_ctg.fa -l /path/busco_downloads/lineages/embryophyta_odb10 -o BaK.bp.p_ctg.busco.embryophyta_odb10 -c 40 -m genome -f > BaK.bp.p_ctg.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i BaK.bp.hap1.p_ctg.fa -l /path/busco_downloads/lineages/embryophyta_odb10 -o BaK.bp.hap1.p_ctg.busco.embryophyta_odb10 -c 40 -m genome -f > BaK.bp.hap1.p_ctg.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i BaK.bp.hap2.p_ctg.fa -l /path/busco_downloads/lineages/embryophyta_odb10 -o BaK.bp.hap2.p_ctg.busco.embryophyta_odb10 -c 40 -m genome -f > BaK.bp.hap2.p_ctg.busco.embryophyta_odb10.out &

## RepeatMasker Step 1: BuildDatabase
nice -n 10 /path/programs/RepeatModeler/BuildDatabase -name Ba ../BaK.hifi.bp.p_ctg.fa > BaK.BuildDatabase.log

## RepeatMasker Step 2: RepeatModeler
nice -n 10 /path/programs/RepeatModeler/RepeatModeler \
 -engine ncbi -threads 10 -database Ba > BaK.RepeatModeler.out

## RepeatMasker Step 4: RepeatMasker
mkdir Ba.hardMasked
nice -n 10 /path/programs/RepeatMasker/RepeatMasker -pa 4 -no_is -dir ./BaK.hardMasked/ -lib consensi.fa.classified ../BaK.bp.p_ctg.fa > BaK.RMhard.out

## Run RepeatMasker for Haplotypes 1 and 2 as well

## BRAKER3 gene annotation
nice -n 10 nohup braker.pl --genome=/path/BaK.hifi.bp.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=BaK.RM.braker3 --threads 40 --verbosity=4 > BaK.RM.braker3.run1.log &
nice -n 10 nohup braker.pl --genome=/path/BaK.hifi.bp.hap1.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=BaK.hap1.RM.braker3 --threads 20 --verbosity=4 > BaK.hap1.RM.braker3.run1.log &
nice -n 10 nohup braker.pl --genome=/path/BaK.hifi.bp.hap2.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=BaK.hap2.RM.braker3 --threads 20 --verbosity=4 > BaK.hap2.RM.braker3.run1.log &



----------------------------
| Dilatris ixioides genome |
----------------------------

## genome assembly
nice -n 10 nohup /path/programs/hifiasm/hifiasm -o Di -t 20 --hg-size 350m Dix_1_0.hifi.fastq.gz > Di.out &
awk '/^S/{print ">"$2"\n"$3}' Di.bp.p_ctg.gfa | fold > Di.bp.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' Di.bp.hap1.p_ctg.gfa | fold > Di.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' Di.bp.hap2.p_ctg.gfa | fold > Di.bp.hap2.p_ctg.fa

## BUSCO assessment
nice -n 10 nohup busco -i Di.bp.p_ctg.fa -l embryophyta_odb10 -o Di.busco.embryophyta_odb10 -c 40 -m genome -f > Di.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i Di.bp.hap1.p_ctg.fa -l embryophyta_odb10 -o Di.hap1.busco.embryophyta_odb10 -c 40 -m genome -f > Di.hap1.busco.embryophyta_odb10.out &
nice -n 10 nohup busco -i Di.bp.hap2.p_ctg.fa -l embryophyta_odb10 -o Di.hap2.busco.embryophyta_odb10 -c 40 -m genome -f > Di.hap2.busco.embryophyta_odb10.out &

## QUAST assessment
nohup /path/programs/quast/quast.py Di.fa -o Di.quast &
nohup /path/programs/quast/quast.py Di.hap1.fa -o Di.hap1.quast &
nohup /path/programs/quast/quast.py Di.hap2.fa -o Di.hap2.quast &

## RepeatMasker Step 1: BuildDatabase
nice -n 10 /path/programs/RepeatModeler/BuildDatabase -name Di ../Di.bp.p_ctg.fa > Di.BuildDatabase.log

## RepeatMasker Step 2: RepeatModeler
nice -n 10 /path/programs/RepeatModeler/RepeatModeler \
 -engine ncbi -threads 10 -database Di > Di.RepeatModeler.out

## RepeatMasker Step 4: RepeatMasker
mkdir Di.hardMasked
nice -n 10 /path/programs/RepeatMasker/RepeatMasker -pa 4 -no_is -dir ./Di.hardMasked/ -lib consensi.fa.classified ../Di.bp.p_ctg.fa > Di.RMhard.out

## Run RepeatMasker for Haplotypes 1 and 2 as well

## BRAKER3 gene annotation
nice -n 10 nohup braker.pl --genome=/path/Di.bp.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=Di.RM.braker3 --threads 40 --verbosity=4 > Di.RM.braker3.run1.log &
nice -n 10 nohup braker.pl --genome=/path/Di.bp.hap1.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=Di.hap1.RM.braker3 --threads 20 --verbosity=4 > Di.hap1.RM.braker3.run1.log &
nice -n 10 nohup braker.pl --genome=/path/Di.bp.hap2.p_ctg.fa.masked --prot_seq=/path/braker3/Viridiplantae.odb11.for_braker.fa --workingdir=Di.hap2.RM.braker3 --threads 20 --verbosity=4 > Di.hap2.RM.braker3.run1.log &