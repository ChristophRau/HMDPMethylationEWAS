# The goal of this script is to call methylation on BAM files for the HMDP.

# First, set up a conda environment
# In one terminal:
alias conda='/home/user/anaconda2/bin/conda' # Do this once

# Create conda enironment
conda create --name methHMDP python=2.7
exit

# In a NEW terminal
PATH=/home/user/anaconda2/bin/:$PATH
source activate methHMDP # do this without the conda alias from above

# Install pysam
conda install -c bioconda pysam

# Now, cat all chrs to mm10.fa
cd /media/user/Methylation/chromFa/

for i in *.fa
do
  cat $i >> mm10.fa
done

cd /media/user/Methylation/

# Now, build the RRBS digested genome (MspI)
bs_path=/media/user/Methylation/BSseeker2-master
bs_build=$bs_path/bs_seeker2-build.py
bs_align=$bs_path/bs_seeker2-align.py

mm10_path=/media/user/Methylation/chromFa/mm10.fa

python $bs_build -f $mm10_path -r -l 20 -u 500 --aligner=bowtie2

# Define function to check number of jobs
# Note that alignment needs multiple cores per job
NUM_JOBS=12
function max_job_check {
  while [ `jobs | wc -l` -ge $NUM_JOBS ]
    do
      sleep 10
    done
}


#Here we go, actually doing the aligning
met=/media/user/Methylation/
bs_path=/media/user/Methylation/BSseeker2-master
bs_build=$bs_path/bs_seeker2-build.py
bs_align=$bs_path/bs_seeker2-align.py




for directory in $met/Raw/Unzipped/*
do
  for sample in $directory/*
  do
    sample_pref=$(basename $sample)
    echo 'Working on:' $sample_pref

    cd $sample

    echo 'Sample directory:' >> report2.txt
    echo $sample >> report2.txt
    echo 'Sample name for BAM:' >> report2.txt
    echo $sample_pref >> report2.txt
    echo 'Sample FASTQ (gzipped)' >> report2.txt
    echo *fastq.gz >> report2.txt

    cat *fastq.gz > $sample_pref.fastq.gz

    outdir=$met/Raw_aligned2/$sample_pref
    mkdir $outdir

    max_job_check

    $bs_align -i $sample_pref.fastq.gz -g $met/chromFa/mm10.fa -r -L 20 -U 500 --aligner=bowtie2 -o $outdir/$sample_pref.bam -a $met/adapter.txt &

  done
done

#Sort BAM files
for directory in $met/Raw_aligned/*
do
    sample_pref=$(basename $directory)
    echo 'Working on:' $sample_pref
    
    cd $directory
    
    max_job_check
    samtools sort -o $sample_pref.sorted.bam $sample_pref.bam &
done 


#Quantify Methylation with cgmaptools
for directory in $met/Raw_aligned/*
do
  sample_pref=$(basename $directory)
  echo 'Working on:' $sample_pref

  cd $directory

  max_job_check
  /media/user/Methylation/cgmaptools-master/cgmaptools convert bam2cgmap -b $sample_pref.sorted.bam -g /media/user/Methylation/chromFa/mm10.fa -o /media/user/Methylation/Raw_called/$sample_pref &
done

#Identify SNPs with cgmaptools
for directory in $met/Raw_aligned/*
do
  sample_pref=$(basename $directory)
  echo 'Working on:' $sample_pref

  #outdir=$met/Raw_called/$sample_pref
  #mkdir $outdir

  max_job_check
  /media/user/Methylation/cgmaptools-master/cgmaptools snv -i /media/user/Methylation/Raw_called/$sample_pref.ATCGmap.gz -m bayes -o /media/user/Methylation/Raw_called/SNPs/$sample_pref.SNPs &
done

#cleanup

for directory in $met/Raw_aligned/*
do
  sample_pref=$(basename $directory)
  echo 'Working on:' $sample_pref

  cd $directory

  rm *.bam_sorted.bam

done






