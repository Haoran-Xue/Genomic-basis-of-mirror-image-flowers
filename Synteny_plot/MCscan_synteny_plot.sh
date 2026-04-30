conda activate jcvi

awk '$3 == "transcript" { 
	print $1"\t"$4"\t"$5"\t"$9"\t0\t"$7 
}' Wp.RM.braker3.gtf > Wp.bed

awk '$3 == "transcript" {
	print $1"\t"$4"\t"$5"\t"$9"\t0\t"$7
}' Wt.RM.braker3.gtf > Wt.bed

awk '$3 == "transcript" {
print $1"\t"$4"\t"$5"\t"$9"\t0\t"$7
}' BaN.RM.braker3.gtf > BaN.bed

awk '$3 == "transcript" { 
    print $1"\t"$4"\t"$5"\t"$9"\t0\t"$7
}' BaK.hap2.RM.braker3.gtf > Ba.bed

awk '$3 == "transcript" {
    print $1"\t"$4"\t"$5"\t"$9"\t0\t"$7
}' Di.RM.braker3.gtf > Di.bed


awk ' $1 == "ptg000001l" ' Wp.bed > WpR.bed
awk ' $1 == "ptg000033l" ' Wp.bed > WpP.bed
awk ' $1 == "ptg000020l" ' Wt.bed > WtR.bed
awk ' $1 == "ptg000024l" ' Wt.bed > WtP.bed
awk ' $1 == "ptg000064l"' BaN.bed > BaNR.bed
awk  '$1 == "ptg000031l"' BaN.bed > BaNP.bed
awk ' $1 == "h2tg000059l" ' BaK.bed > BaKR.bed
awk ' $1 == "h2tg000013l" ' BaK.bed > BaKP.bed
awk ' $1 == "ptg000005l" ' Di.bed > DiP.bed

cut -f4 WpR.bed | sort | uniq > WpR.tx.txt
cut -f4 WpP.bed | sort | uniq > WpP.tx.txt
cut -f4 WtR.bed | sort | uniq > WtR.tx.txt
cut -f4 WtP.bed | sort | uniq > WtP.tx.txt
cut -f4 BaNR.bed | sort | uniq > BaNR.tx.txt
cut -f4 BaNP.bed | sort | uniq > BaNP.tx.txt
cut -f4 BaKR.bed | sort | uniq > BaKR.tx.txt
cut -f4 BaKP.bed | sort | uniq > BaKP.tx.txt
cut -f4 DiP.bed | sort | uniq > DiP.tx.txt

/mnt/barn/hxue/programs/seqtk/seqtk subseq Wp.RM.braker3.cds.fa WpR.tx.txt > WpR.cds
/mnt/barn/hxue/programs/seqtk/seqtk subseq Wp.RM.braker3.cds.fa WpP.tx.txt > WpP.cds
/mnt/barn/hxue/programs/seqtk/seqtk subseq Wt.RM.braker3.cds.fa WtR.tx.txt > WtR.cds
/mnt/barn/hxue/programs/seqtk/seqtk subseq Wt.RM.braker3.cds.fa WtP.tx.txt > WtP.cds
/mnt/barn/hxue/programs/seqtk/seqtk subseq BaN.RM.braker3.cds.fa BaNR.tx.txt > BaNR.cds
/mnt/barn/hxue/programs/seqtk/seqtk subseq BaN.RM.braker3.cds.fa BaNP.tx.txt > BaNP.cds
/mnt/barn/hxue/programs/seqtk/seqtk subseq BaK.hap2.RM.braker3.cds.fa BaKR.tx.txt > BaKR.cds
/mnt/barn/hxue/programs/seqtk/seqtk subseq BaK.hap2.RM.braker3.cds.fa BaKP.tx.txt > BaKP.cds
/mnt/barn/hxue/programs/seqtk/seqtk subseq DiP.RM.braker3.cds.fa DiP.tx.txt > BaNP.cds

## Add MIR156 manually

python -m jcvi.compara.catalog ortholog -n 2 DiP  WpR --cscore=0.5
python -m jcvi.compara.catalog ortholog -n 2 DiP  WtR --cscore=0.5
python -m jcvi.compara.catalog ortholog -n 2 DiP   BaNR --cscore=0.5
python -m jcvi.compara.catalog ortholog -n 2 DiP   BaKR --cscore=0.5

python -m jcvi.compara.catalog ortholog -n 2 DiP  WpP --cscore=0.5
python -m jcvi.compara.catalog ortholog -n 2 DiP  WtP --cscore=0.5
python -m jcvi.compara.catalog ortholog -n 2 DiP   BaNP --cscore=0.5
python -m jcvi.compara.catalog ortholog -n 2 DiP   BaKP --cscore=0.5

python -m jcvi.compara.synteny mcscan DiP.bed DiP.WpR.lifted.anchors --iter=1 -o DiP_WpR.i1.blocks
python -m jcvi.compara.synteny mcscan DiP.bed DiP.WtR.lifted.anchors --iter=1 -o DiP_WtR.i1.blocks
python -m jcvi.compara.synteny mcscan DiP.bed DiP.BaNR.lifted.anchors --iter=1 -o DiP_BaNR.i1.blocks
python -m jcvi.compara.synteny mcscan DiP.bed DiP.BaKR.lifted.anchors --iter=1 -o DiP_BaKR.i1.blocks
python -m jcvi.compara.synteny mcscan DiP.bed DiP.WpP.lifted.anchors --iter=1 -o DiP_WpP.i1.blocks
python -m jcvi.compara.synteny mcscan DiP.bed DiP.WtP.lifted.anchors --iter=1 -o DiP_WtP.i1.blocks
python -m jcvi.compara.synteny mcscan DiP.bed DiP.BaNP.lifted.anchors --iter=1 -o DiP_BaNP.i1.blocks
python -m jcvi.compara.synteny mcscan DiP.bed DiP.BaKP.lifted.anchors --iter=1 -o DiP_BaKP.i1.blocks

python -m jcvi.formats.base join DiP_WpR.i1.blocks DiP_WtR.i1.blocks DiP_BaNR.i1.blocks DiP_BaKR.i1.blocks  DiP_WpP.i1.blocks DiP_WtP.i1.blocks DiP_BaNP.i1.blocks DiP_BaKP.i1.blocks   --noheader | cut -f1,2,4,6,8,10,12,14,16 > Dip_ref.blocks

awk '{printf "%s\t.\t.\t.\t.\t.\t.\t.\t.\n", $4}' DiP.bed > DiP.all.blocks
awk '{printf ".\t%s\t.\t.\t.\t.\t.\t.\t.\n", $4}' WpR.bed > WpR.all.blocks
awk '{printf ".\t.\t%s\t.\t.\t.\t.\t.\t.\n", $4}' WtR.bed > WtR.all.blocks
awk '{printf ".\t.\t.\t%s\t.\t.\t.\t.\t.\n", $4}' BaNR.bed > BaNR.all.blocks
awk '{printf ".\t.\t.\t.\t%s\t.\t.\t.\t.\n", $4}' BaKR.bed > BaKR.all.blocks
awk '{printf ".\t.\t.\t.\t.\t%s\t.\t.\t.\n", $4}' WpP.bed > WpP.all.blocks
awk '{printf ".\t.\t.\t.\t.\t.\t%s\t.\t.\n", $4}' WtP.bed > WtP.all.blocks
awk '{printf ".\t.\t.\t.\t.\t.\t.\t%s\t.\n", $4}' BaNP.bed > BaNP.all.blocks
awk '{printf ".\t.\t.\t.\t.\t.\t.\t.\t%s\n", $4}' BaKP.bed > BaKP.all.blocks

cat  Dip_ref.blocks   DiP.all.blocks WpR.all.blocks   WtR.all.blocks BaNR.all.blocks  BaKR.all.blocks  WpP.all.blocks  WtP.all.blocks BaNP.all.blocks BaKP.all.blocks >  Dip_ref.withnomatch.blocks

python -m jcvi.graphics.synteny Dip_ref.withnomatch.blocks newmega.bed  Dip_ref.all.vertical.micro.layout --outputprefix Microsyn.all.Dip_ref --font Arial 
