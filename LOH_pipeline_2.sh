export R_LIBS="/home/local/AMC/nithisha/R/R.3.3-library"
set -x -e
#declare arguments
TUMOR_PFX=$1 #tumor prefix
NORMAL_PFX=$2
SAVEROOT=$3 #path to where you want your results saved
DATAPATH=$4 #path to directory containing sample-specific fastq files
CORES=$5

#paths to required programs
BWA="/mnt/disk4/labs/salipante/programs/bwa-0.7.12/bwa"
SAMTOOLS="/mnt/disk4/labs/salipante/programs/samtools-0.1.19/samtools"
SEQUENZAUTILS="/home/local/AMC/nithisha/R/R.3.3-library/sequenza/exec/sequenza-utils.py"
ALIGNMENT_REF_GENOME="/mnt/disk2/com/Genomes/hg19_chr/human_g1k_v37.fasta"  
GC_WINDOW="/mnt/disk4/labs/salipante/stevesal/LOH_calling/scripts/gcwindow/hg19.gc5Base_v4.txt.gz"
SUPER_DEDUPER="/mnt/disk4/labs/salipante/programs/Super-Deduper-master/super_deduper"
#SCORING_SCRIPT="/mnt/disk4/labs/salipante/stevesal/LOH_calling/scripts/sequenza_loss_score_v2.pl"  
SCORING_SCRIPT_PYTHON="/mnt/disk4/labs/salipante/nithisha/LOH_Scoring/nit_code/LOH_score.py"
#SUMMARY_FILE="/mnt/disk4/labs/salipante/nithisha/LOH_Scoring/nit_summary_file.txt"

#create a subdirectory for this tumor specimen
SAVEPATH=$SAVEROOT/$TUMOR_PFX
mkdir -p $SAVEPATH

#find tumor sequences
TUMOR_R1=$(ls $DATAPATH/$TUMOR_PFX/*_1_sequence.txt.gz)
TUMOR_R2=$(ls $DATAPATH/$TUMOR_PFX/*_2_sequence.txt.gz)

#find normal sequences
NORMAL_R1=$(ls $DATAPATH/$NORMAL_PFX/*_1_sequence.txt.gz)
NORMAL_R2=$(ls $DATAPATH/$NORMAL_PFX/*_2_sequence.txt.gz)

#-----------------------------------------------------------------------------------------------
#CHECK FOR MPILEUPS, create them if they dont exist
#-----------------------------------------------------------------------------------------------
#ls $SAVEPATH/$NORMAL_PFX.mpileup.gz
if [ ! -f $SAVEPATH/$NORMAL_PFX.mpileup.gz ] #remember spaces in []
then 
#-----------------------------------------------------------------------------------------------
#preprocess tumor reads 
$SUPER_DEDUPER -1 $TUMOR_R1 -2 $TUMOR_R2 -p $SAVEPATH/$TUMOR_PFX
echo "line 41...Super deduper for tumor done, check for 2 _nodup_PE1/2.fastq files"
gzip $SAVEPATH/${TUMOR_PFX}_nodup_PE1.fastq &
gzip $SAVEPATH/${TUMOR_PFX}_nodup_PE2.fastq &

wait
echo "line 46...gzipping tumor done, check for 2 _nodup_PE1/2.fastq.gz files"

#align tumor
#Align fastq files to gatk reference genome - this finds SA coordinates of input reads (.sai- suffix array indices)
$BWA aln -t $CORES $ALIGNMENT_REF_GENOME $SAVEPATH/${TUMOR_PFX}_nodup_PE1.fastq.gz > $SAVEPATH/$TUMOR_PFX.1.sai 2>> $SAVEPATH/$TUMOR_PFX.log
$BWA aln -t $CORES $ALIGNMENT_REF_GENOME $SAVEPATH/${TUMOR_PFX}_nodup_PE2.fastq.gz > $SAVEPATH/$TUMOR_PFX.2.sai 2>> $SAVEPATH/$TUMOR_PFX.log
echo "line 52...BWA aln for tumors done, check for 2 .sai files"
#question in line 36, why is {} needed?
# .sai to .sam file and convert SA coordinates to chromosomal loci
$BWA sampe -r "@RG\tID:${PFX}\tPL:ILLUMINA\tPU:NA\tLB:null\tSM:${TUMOR_PFX:?}" $ALIGNMENT_REF_GENOME $SAVEPATH/$TUMOR_PFX.1.sai $SAVEPATH/$TUMOR_PFX.2.sai $SAVEPATH/${TUMOR_PFX}_nodup_PE1.fastq.gz $SAVEPATH/${TUMOR_PFX}_nodup_PE2.fastq.gz 2>> $SAVEPATH/$TUMOR_PFX.log >$SAVEPATH/$TUMOR_PFX.sam
echo "line 56... BWA sampe done for tumor, check for 1 .sam file"
#remove large files
rm -f $SAVEPATH/${TUMOR_PFX}_nodup_PE1.fastq.gz
rm -f $SAVEPATH/${TUMOR_PFX}_nodup_PE2.fastq.gz

#use SAMTOOLS to view, sort and index, .sam (sequence alignment/mapping) to .bam
#view - converts .sam to .bam without any action taken
$SAMTOOLS view -u -b -S -F 4 -Q 20 $SAVEPATH/$TUMOR_PFX.sam | $SAMTOOLS sort - $SAVEPATH/$TUMOR_PFX 2>> $SAVEPATH/$TUMOR_PFX.log
# u: uncompressed if pipe is used, b: bam file output, S: checking for compatibility with older SAMTOOLS version
# F: do not output alignments with FLAG, q: skip alignments with MAPQ lesser than
# T: prefix to be added to output

#index BAM file (Index a coordinate-sorted BAM file for fast random access)
$SAMTOOLS index $SAVEPATH/$TUMOR_PFX.bam
echo "line 70... SAMTOOLS view, sort and index done for tumor, check for 1 .bam file"
#remove more intermediate files
rm -f $SAVEPATH/$TUMOR_PFX.1.sai
rm -f $SAVEPATH/$TUMOR_PFX.2.sai
rm -f $SAVEPATH/$TUMOR_PFX.sam

#-----------------------------------------------------------------------------------------------
#preprocess normal reads 
$SUPER_DEDUPER -1 $NORMAL_R1 -2 $NORMAL_R2 -p $SAVEPATH/$NORMAL_PFX

gzip $SAVEPATH/${NORMAL_PFX}_nodup_PE1.fastq &
gzip $SAVEPATH/${NORMAL_PFX}_nodup_PE2.fastq &

wait

#align normal
#Align fastq files to gatk reference genome - this finds SA coordinates of input reads (.sai- suffix array indices)
$BWA aln -t $CORES $ALIGNMENT_REF_GENOME $SAVEPATH/${NORMAL_PFX}_nodup_PE1.fastq.gz > $SAVEPATH/$NORMAL_PFX.1.sai 2>> $SAVEPATH/$NORMAL_PFX.log
$BWA aln -t $CORES $ALIGNMENT_REF_GENOME $SAVEPATH/${NORMAL_PFX}_nodup_PE2.fastq.gz > $SAVEPATH/$NORMAL_PFX.2.sai 2>> $SAVEPATH/$NORMAL_PFX.log

# .sai to .sam file and convert SA coordinates to chromosomal loci
$BWA sampe -r "@RG\tID:${PFX}\tPL:ILLUMINA\tPU:NA\tLB:null\tSM:${NORMAL_PFX:?}" $ALIGNMENT_REF_GENOME $SAVEPATH/$NORMAL_PFX.1.sai $SAVEPATH/$NORMAL_PFX.2.sai $SAVEPATH/${NORMAL_PFX}_nodup_PE1.fastq.gz $SAVEPATH/${NORMAL_PFX}_nodup_PE2.fastq.gz 2>>$SAVEPATH/$NORMAL_PFX.log > $SAVEPATH/$NORMAL_PFX.sam
echo "line 92... BWA done for normal, check for .sam file"
#remove large files
rm -f $SAVEPATH/${NORMAL_PFX}_nodup_PE1.fastq.gz
rm -f $SAVEPATH/${NORMAL_PFX}_nodup_PE2.fastq.gz

#use SAMTOOLS to view, sort and index, .sam (sequence alignment/mapping) to .bam
#view - converts .sam to .bam without any action taken
#QUESTION q or Q?
$SAMTOOLS view -u -b -S -F 4 -Q 20 $SAVEPATH/$NORMAL_PFX.sam | $SAMTOOLS sort - $SAVEPATH/$NORMAL_PFX 2>>$SAVEPATH/$NORMAL_PFX.log
# u: uncompressed if pipe is used, b: bam file output, S: checking for compatibility with older SAMTOOLS version
# F: do not output alignments with FLAG, q: skip alignments with MAPQ lesser than
# T: prefix to be added to output
echo "line 104... SAMTOOLS view and sort done for normal, check for .bam file"
#index BAM file (Index a coordinate-sorted BAM file for fast random access)
$SAMTOOLS index $SAVEPATH/$NORMAL_PFX.bam
echo "line 107... SAMTOOLS index done for normal, check for .bam file"
#remove more intermediate files
rm -f $SAVEPATH/$NORMAL_PFX.1.sai
rm -f $SAVEPATH/$NORMAL_PFX.2.sai
rm -f $SAVEPATH/$NORMAL_PFX.sam

#-----------------------------------------------------------------------------------------------
#CHECK for 2 bams at this stage, if there was a samtool permission denied, segmentation faults (core dumped) error, you will not have 2 bams - one for tumor and one for germline
#at this point just exit and re-run this sample
if [ ! -f $SAVEPATH/$TUMOR_PFX.bam ] && [ ! -f $SAVEPATH/$NORMAL_PFX.bam ]
then
echo 'ERROR ENCOUNTERED- RE_RUN SAMPLE'
exit
fi
#-----------------------------------------------------------------------------------------------
#Make mpileups and gzip them using SAMTOOLS
$SAMTOOLS mpileup -f $ALIGNMENT_REF_GENOME -d 10000 -A -B $SAVEPATH/$TUMOR_PFX.bam | gzip > $SAVEPATH/$TUMOR_PFX.mpileup.gz
$SAMTOOLS mpileup -f $ALIGNMENT_REF_GENOME -d 10000 -A -B $SAVEPATH/$NORMAL_PFX.bam | gzip > $SAVEPATH/$NORMAL_PFX.mpileup.gz
#f: reference genome, d: max depth, A: Do not skip anomalous read pairs in variant calling
#B: Disable probabilistic realignment for the computation of base alignment quality (BAQ). 
#BAQ is the Phred-scaled probability of a read base being misaligned. 
#Applying this option greatly helps to reduce false SNPs caused by misalignments. 

fi
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
#RUN SEQUENZA, PYTHON
#-----------------------------------------------------------------------------------------------
if [ ! -f $SAVEPATH/$TUMOR_PFX.binned.seqz.gz ]
then
python $SEQUENZAUTILS pileup2seqz -gc $GC_WINDOW -n $SAVEPATH/$NORMAL_PFX.mpileup.gz -t $SAVEPATH/$TUMOR_PFX.mpileup.gz | gzip > $SAVEPATH/$TUMOR_PFX.seqz.gz
python $SEQUENZAUTILS seqz-binning -w 50 -s $SAVEPATH/$TUMOR_PFX.seqz.gz | gzip > $SAVEPATH/$TUMOR_PFX.binned.seqz.gz
#gc: gc window, n: normal mpileup, t:tumor mpileup
#w:-w indicate a window size (in bases), to be used for the binning, s:save to file??
fi
echo "line 133... Sequenza in python done!"
#-----------------------------------------------------------------------------------------------
#RUN SEQUENZA, R
#-----------------------------------------------------------------------------------------------
#Create an executable R script, run it and quit it!
echo "library(\"sequenza\")">$SAVEPATH/$TUMOR_PFX.sequenza.r
echo "data.file <- \"$SAVEPATH/$TUMOR_PFX.binned.seqz.gz\"" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
#echo "data.file <- \"$SAVEPATH/$TUMOR_PFX.seqz.gz\"" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "seqz.data <- read.seqz(data.file)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "gc.stats <- gc.sample.stats(data.file)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "test <- sequenza.extract(data.file)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "CP.example <- sequenza.fit(test)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = \"$TUMOR_PFX\", out.dir=\"$SAVEPATH\")" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "cint <- get.ci(CP.example)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r

#Plot cellularity
echo "jpeg(\"$SAVEPATH/$TUMOR_PFX.nitz.cellularity.jpg\")" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "cp.plot(CP.example)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "cp.plot.contours(CP.example, add = TRUE, likThresh=c(0.95))" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "dev.off()" >> $SAVEPATH/$TUMOR_PFX.sequenza.r

#Call CNVs
echo "cellularity <- cint\$max.cellularity" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "ploidy <- cint\$max.ploidy" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "avg.depth.ratio <- mean(test\$gc\$adj[,2])" >> $SAVEPATH/$TUMOR_PFX.sequenza.r

#Save parameters to file to file
echo "cellularity" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "write(cellularity, file = \"$SAVEPATH/$TUMOR_PFX.nitz.cellularity.txt\")" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "write(ploidy, file = \"$SAVEPATH/$TUMOR_PFX.nitz.ploidy.txt\")" >>$SAVEPATH/$TUMOR_PFX.sequenza.r
echo "write(avg.depth.ratio, file = \"$SAVEPATH/$TUMOR_PFX.nitz.ave_depth.txt\")" >> $SAVEPATH/$TUMOR_PFX.sequenza.r

#Detect variant alleles
echo "mut.tab <- na.exclude(do.call(rbind, test\$mutations))" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "mut.alleles <- mufreq.bayes(mufreq = mut.tab\$F, depth.ratio = mut.tab\$adjusted.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r

#Detect CN variation
echo "seg.tab <- na.exclude(do.call(rbind, test\$segments))" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "cn.alleles <- baf.bayes(Bf = seg.tab\$Bf, depth.ratio = seg.tab\$depth.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "seg.tab <- cbind(seg.tab, cn.alleles)" >>$SAVEPATH/$TUMOR_PFX.sequenza.r
echo "seg.tab" >> $SAVEPATH/$TUMOR_PFX.sequenza.r

#write sequenza matrix to file, this will serve as input to loss score script's 2nd arg
echo "write.table(seg.tab, file = \"$SAVEPATH/$TUMOR_PFX.nitz.copynumber_calls.txt\", append = FALSE)" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
#QUESTION: how does it know what $SAVEPATH is in R?
#exit
echo "q()" >> $SAVEPATH/$TUMOR_PFX.sequenza.r
echo "n" >> $SAVEPATH/$TUMOR_PFX.sequenza.r

#execute the R script
R --vanilla < $SAVEPATH/$TUMOR_PFX.sequenza.r

#Calculate loss score
if [ -f $SAVEPATH/$TUMOR_PFX.nitz.copynumber_calls.txt ]
then 
python $SCORING_SCRIPT_PYTHON $SAVEPATH/$TUMOR_PFX.nitz.copynumber_calls.txt $SAVEPATH/$TUMOR_PFX.nitz.score.txt 0.75

#Append some usefull QC factoids to the lossScore output text
echo "" >> $SAVEPATH/$TUMOR_PFX.nitz.score.txt
echo -n "Estimated tumor cellularity: " >> $SAVEPATH/$TUMOR_PFX.nitz.score.txt
cat $SAVEPATH/$TUMOR_PFX.nitz.cellularity.txt >> $SAVEPATH/$TUMOR_PFX.nitz.score.txt
echo -n "Estimated ploidy: " >> $SAVEPATH/$TUMOR_PFX.nitz.score.txt
cat $SAVEPATH/$TUMOR_PFX.nitz.ploidy.txt >> $SAVEPATH/$TUMOR_PFX.nitz.score.txt

else
exit

fi
