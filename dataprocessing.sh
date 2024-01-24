#Hybrid capture lcWGS

#load conda environment
#my conda env: eabbott@ls6.tacc.utexas.edu:/work/05410/eabbott/ls6/conda/envs/mcavsyms
conda activate mcavsyms

cd

#cat lanes ----------------------------------------------
>do_concat
for file in *L001_R1_001.fastq
do echo "cat $file ${file/_L001_R1_001.fastq/}_L002_R1_001.fastq > ${file/_L001_R1_001.fastq/}_R1.fastq" >>do_concat
done

ls6_launcher_creator.py -n do_concat -j do_concat -q normal -N 10 -w 20 -a $allo -t 02:00:00

>do_concat
for file in *L001_R2_001.fastq
do echo "cat $file ${file/_L001_R2_001.fastq/}_L002_R2_001.fastq > ${file/_L001_R2_001.fastq/}_R2.fastq" >>do_concat
done
ls6_launcher_creator.py -n do_concat -j do_concat -q normal -N 2 -w 5 -a $allo -t 00:30:00

#trim adapters ----------------------------------------------
#trim illumina
>tg
for file in *R1.fastq
do echo "trim_galore --illumina --paired \
${file} ${file/R1.fastq}R2.fastq" >>tg
done
ls6_launcher_creator.py -n tg -j tg -a $allo -e $email -q normal -t 04:00:00 -N 8 -w 4

#trim 5' end
>trimpe
for file in *_R1_val_1.fq
do echo "cutadapt \
-g AGATGTGTATAAGAGACAG \
-G AGATGTGTATAAGAGACAG \
-o ${file/_R1_val_1.fq/}_R1.trim \
-p ${file/_R1_val_1.fq/}_R2.trim \
${file} \
$file > ${file}_trimlog.txt" >> trimpe
done

ls6_launcher_creator.py -n trimpe -j trimpe -a $allo -e $email -q normal -t 3:00:00 -N 8 -w 46

#trim 3'end 
>trimpe3
for file in *R1.trim
do echo "cutadapt \
-a CTGTCTCTTATACACATCT \
-A CTGTCTCTTATACACATCT \
-o ${file/R1.trim/}3p_R1.trim \
-p ${file/R1.trim/}3p_R2.trim \
${file} \
$file > ${file}_trimlog.txt" >> trimpe3
done

ls6_launcher_creator.py -n trimpe3 -j trimpe3 -a $allo -e $email -q normal -t 01:00:00 -N 8 -w 48


#check quality with fastqc ----------------------------------------------
mkdir Fastqc_Restults_postTrim/
> runFQC2
for file in *3p*
do echo "fastqc -o /scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/Fastqc_Restults_postTrim -f fastq $file &" >> runFQC2
done

ls6_launcher_creator.py -n runFQC2 -j runFQC2 -a $allo -e $email -q development -t 02:00:00 -N 2 -w 10

#MAPPING READS --------------------------
#must be a concatenated reference including the host and symbiont
export GENOME_FASTA=/work/05410/eabbott/ls6/genomes/mcavzooxcarly/mcav_abcd_ref.fasta

#runIDs to match paired end reads
ls *_3p_R1.trim | awk '{split($0, a, "_3p_R");print a[1]}' > paired_runIDs.txt

>mappe
while read runID
do echo "\
bowtie2 -x $GENOME_FASTA -1 ${runID}_3p_R1.trim -2 ${runID}_3p_R2.trim --local -p 8 -S ${runID/_3p_R1.trim/}.sam" >> mappe
done < paired_runIDs.txt
ls6_launcher_creator.py -n mappe -j mappe -q normal -N 10 -w 10 -a $allo -e $email -t 12:00:00
sbatch mappe.slurm

#SPLIT OUT SYM AND CORAL READS ------------------------------
#split sam files by organism. use sams, not bams, to preserve headers necessary for samtools sort
#create directory for each organism and copy sam files to each
>cpsams
for file in *sam
do echo "cp $file /scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd" >>cpsams
done
ls6_launcher_creator.py -n cpsams -j cpsams -t 09:00:00 -q vm-small -a $allo -e $email -N 1 -w 1

#chr22 is the scaffold containing the entire symbiont reference
>cladsep
for file in *sam
do echo "grep  "chr22" $file > ${file/.sam}.clad.sam" >>cladsep
done
ls6_launcher_creator.py -n cladsep -j cladsep -t 24:00:00 -q vm-small -a $allo -e $email -N 1 -w 1

#host scaffolds have all other prefixes
>mcavsep
for file in *sam
do echo "grep  -v "chr" $file > ${file/.sam}.mcav.sam" >>mcavsep
done
ls6_launcher_creator.py -n mcavsep -j mcavsep -t 24:00:00 -q vm-small -a $allo -e $email -N 1 -w 1
sbatch --dependency=afterok:1405518 mcavsep.slurm 

#sort and index sam files. do this in the host and symbiont directories
>sort
for file in *.sam
do echo "samtools sort -o ${file/.sam/}.bam $file && samtools index -c ${file/.sam}.bam" >>sort
done
ls6_launcher_creator.py -n sort -j sort -t 04:00:00 -q gpu-a100 -a $allo -e $email -N 4 -w 8

#move mcavs to their own directory for angsd
mv *mcav.bam.csi /scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/mcavangsd


#GET RATIO OF ZOOX TO CORAL READS ------------------------
#use original, unsplit sam files
>bestread
for file in *.bam; do
echo "samtools view -bq 1 $file > ${file/.bam}.unique.bam && samtools index -c ${file/.bam}.unique.bam" >> bestread; done

ls6_launcher_creator.py -n bestread -j bestread -t 10:00:00 -q vm-small -a $allo -e $email -N 1 -w 1

for file in *.bam; 
do samtools idxstats $file > ${file/.bam}.stats; done

#scp stats to local
scp -r eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/hybcapJA232326/JA23377/catlanes/st .

############################
### POPULATION STRUCTURE ###
############################
#do this twice, once in symbiont directory and once in coral directory

#run qualRanks to identify low quality samples
FILTERSQ="-uniqueOnly 1 -remove_bads 1 -minMapQ 20"
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo "angsd -b bams -r chr22 -GL 1 $FILTERSQ $TODOQ -P 12 -out dd && Rscript ~/bin/plotQC.R prefix=dd >qualRanks">a0
ls6_launcher_creator.py -j a0 -n a0 -a $allo -e $email -t 00:10:00 -w 1 -q vm-small
sbatch a0.slurm

#scp dd.pdf to look estimate genotyping rate and minInd for ANGSD
scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd/dd.pdf .

#remove low quality samples
#qualRanks should output bams.qc, a list of bams that passed the quality check
#sometimes this doesn't work, so you must delete them manually
#remove samples that are close to 0. Sometimes samples are very close to 1, remove those too
#cat qualRanks | grep -v '05$' | grep -v '04$' | sed 's/\s.*$//' > bams #make sure to remove header

#INITIAL IBS RUN TO IDENTIFY CLONES -------------------------------------
export GRate=0.3 #shoot for 5 coverage as per the qualRanks results
export MI=70 #aim for ~75% of samples. lower MI if the output includes less than 10,000 sites
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat bams.qc | wc -l`; export MI=`echo "($NIND*$GRate+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams -GL 1 $FILTERS0 $TODO0 -P 12 -out myresult && Rscript ~/bin/detect_clones.R bams.qc myresult.ibsMat 0.06">a1
ls6_launcher_creator.py -j a1 -n a1 -a $allo -e $email -t 00:10:00 -w 1 -q vm-small
sbatch a1.slurm

#samples which fall on y-axis near the techincal replicates are the clones.
#use this to set the IBS cutoff. If your cutoff above is too high/low rerun. This will generate bams.nr, a list of bams of unique genotypes

#sometimes IBS may generate NAs. This is usually caused by low quality samples.
#fix NAs by scp'ing the IBS Matrix and bams.qc and follow removeNAs.R
scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd/bams.qc .
scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd/myresult.ibsMat .

#upload fixed matrix from removeNAs.R
scp -v /Users/evelynabbott/Dropbox/project/Projects/hybcap_mcav_clad/hybcap_december/minind_admix_comp/bams_na_rem eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd
scp -v /Users/evelynabbott/Dropbox/project/Projects/hybcap_mcav_clad/hybcap_december/minind_admix_comp/clad_ngs.ibsMat_na_rem eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd


#rerun IBS excluding clones -------------------------
export GRate=0.3
export MI=70
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat bams.nr | wc -l`; export MI=`echo "($NIND*$GRate+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams.nr -GL 1 $FILTERS0 $TODO0 -P 12 -out myresult && Rscript ~/bin/detect_clones.R bams.nr myresult.ibsMat 0.06">a1
ls6_launcher_creator.py -j a1 -n a1 -a $allo -e $email -t 00:10:00 -w 1 -q vm-small
sbatch a1.slurm

#sites remaining: check a1.e* file to make sure at least 4,000 sites remain
scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd/hctree.pdf .

#Run ngs_admix ---------------------------------
#this will generate input for PCA and admixture in PCA_admix.R
#`seq 2 6` is the range of potential population numbers
>ngs
for K in `seq 2 6`
do echo "NGSadmix -likes myresult.beagle.gz -K $K -P 10 -o clad_k${K}" >> ngs
done

#this is fast so you can run this on the login node

#scp .qopt file and bams_na_rem and ibsMat
scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd/bams.nr .
scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd/myresult.ibsMat .
scp -v eabbott@ls6.tacc.utexas.edu:'/scratch/05410/eabbott/JA23419/JA23419/fastqs/cattedf/cladangsd/*qopt' .
