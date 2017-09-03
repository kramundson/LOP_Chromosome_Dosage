Clean LOP analysis:

Experiment description: Potato dihaploids produced at the International Potato Center
(CIP) were produced by pollinating tetraploid landrace Alca Tarma (LOP868) with pollen
from IVP-101 or PL-4 (an improved inducer and hybrid of IVP-101 and IVP-35).
The pollen donor of each haploid is known for some, but not all dihaploids.

From the total progeny, only those lacking the dominant haploid inducer specific
"embryo spot" phenotypic marker were planted. Progeny were then screened by guard cell
chloroplast count, and those presenting more than 8 per guard cell were discarded. Minimum
ten measured guard cells per individual.

This crossing and selection strategy resulted in a population of 168 dihaploids.

The primary aim of this project was to survey the population for overlapping aneuploidy and
paternal genotypes that suggest haploid induction via selective elimination of most, but
not all haploid inducer chromosomes.

Additionally, sequencing data was used to address the following questions:

A. How does the load of deleterious and dysfunctional alleles scale with increasing ploidy
(monoploid -> diploid -> tetraploid)?

B. How does ploidy affect mutational load in core vs. conserved genes?

C. Do genomic regions exhibit significant segregation distortion?

D. What is the landscape of copy number variation in the tetraploid genome?

A brief description of bench work is provided, and is followed by aneuploidy and CNV
analysis. Mutational load analysis is in progress.

###1. Genomic DNA extraction, library prep, and sequencing

Potato genomic DNA isolated from leaf tissue according to Ghislain et al (1999)

Approximately 750ng of DNA sheared to 300bp using Covaris E-220
Settings:
50ul volume
175W peak incident power
10% duty factor
200 cycles per burst
50s treatment time
4C minimum temp
9C max temp

Half of gDNA carried forward for KAPA Hyper (KR0961) library prep using half-scale rxn,
all other steps according to protocol except as noted below:

Custom 8bp single index adaptors used for each library (see LOPHAP barcodes)

Libraries of haploids amplified using half-scale rxn; Alca Tarma, PL-4, and IVP-101 at full-scale.
4 cycles of library amplification for each library unless noted PCR-free below.

Haploids normalized to 1 ng/ul based on SYBR Green I quantification of cleaned post-amplification library
81 haploids pooled per lane and sequenced on Illumina HiSeq 4000 at UC Davis DNA Tech Core.
One lane sequenced per pool of 81.

PL-4 and Alca Tarma subject to double-cut size selection with SeraMag beads.
Parents PL-4 and Alca Tarma pooled 1:1 based on qPCR quant (qPCR and pooling done by DNA tech core)
and sequenced on HiSeq4000 at UCD DNA Tech core. 

One lane for PL-4 and Alca Tarma.
Haploid inducer parent IVP-101 was sequenced 150PE as part of a different sequencing run.
Library of IVP-101 was prepared by David Castillo as previously described.

Six haploids were prepared by Han Tan PCR-free using KAPA Hyper half-scale reactions.
These haploids were sequenced to somewhat higher depth than the others (1-4x coverage).

Five haploids (LOP868-279, LOP868-280, LOP868-281, LOP868-286, LOP868-289) generated < 200,000 reads.
Libraries of all these except LOP868-279 were remade as previously described and sequenced 90nt SR.
LOP868-279 was withheld from analysis.

Sequencing of all libraries except 90nt SR was done with Illumina HiSeq 4000 at UC Davis DNA tech core.
90nt SR sequenced at Vincent Coates Genome Sequencing Laboratory (UC Berkeley).

###2. Demultiplex without quality trim (allprep12.py). Save demultiplexed reads for SRA submission.

Reads downloaded via FTP. Checksums were not provided by either core, but file sizes were
consistent with expectation (commands not shown).

```
# LOP haploids, lane 1
gunzip LOP-1_S45_L004_R1_001.fastq.gz LOP-1_S45_L004_R2_001.fastq.gz
python ./scripts/allprep-12.py -b LOPHAP1_barcodes.txt \
-f LOP-1_S45_L004_R1_001.fastq -i LOP-1_S45_L004_R2_001.fastq -m -D -n \
> LOP1_demult.out 2> LOP1_demult.err


# LOP haploids, lane 2 to demultiplex libraries with rice index collisions.
# The LOP libraries with collisions can be extracted based on their reverse index seq.
gunzip H718P-Amundson_S46_L005_R*
python ./scripts/allprep-12.py -b modified-LOPHAP2_barcodes.txt \
-f H718P-Amundson_S46_L005_R1_001.fastq -i  H718P-Amundson_S46_L005_R2_001.fastq \
-I H718P-Amundson_S46_L005_R3_001.fastq -m -D -n \
> LOP2_demult.out 2> LOP2_demult.err

LOP haploids, lane 2 to get non-collision libraries
python ./scripts/allprep-12.py -b LOPHAP2_barcodes.txt \
-f H718P-Amundson_S46_L005_R1_001.fastq -i H718P-Amundson_S46_L005_R2_001.fastq \
-m -D -n > LOP2_norice_demult.out 2> LOP2_norice_demult.err


# LOP parents, Alca Tarma/PL-4 lane
gunzip KA-1-2_S26_L002_R1_001.fastq.gz KA-1-2_S26_L002_R2_001.fastq.gz KA-1-2_S26_L002_R4_001.fastq.gz
python ./scripts/allprep-12.py -b KA1-KA2-barcodes_fwdonly.txt
-f KA-1-2_S26_L002_R1_001.fastq -r KA-1-2_S26_L002_R4_001.fastq \
-i KA-1-2_S26_L002_R2_001.fastq -m -D -n > KA1_KA2_demult.out 2> KA1_KA2_demult.err


# LOP parents, IVP-101/IVP-35 lane
gunzip H821P-Amudson_S34_L005_R1_001.fastq.gz H821P-Amudson_S34_L005_R2_001.fastq.gz H821P-Amudson_S34_L005_R4_001.fastq.gz
python ./scripts/allprep-12.py -b IVP35-IVP101-barcodes.txt \
-f H821P-Amudson_S34_L005_R1_001.fastq -r H821P-Amudson_S34_L005_R4_001.fastq \
-i H821P-Amudson_S34_L005_R2_001.fastq -m -D -n > IVP_demult.out 2> IVP_demult.err

# Reads from Han (previously demultiplexed from FRAG run)
for i in *_S{35..36}*
do
    cp /cato2pool/backup-share/kramundson/Potato_dosage_backup/data/EXP021/Raw_Data/$i ./
done

for i in *_S{41..48}*
do
    cp /cato2pool/backup-share/kramundson/Potato_dosage_backup/data/EXP021/Raw_Data/$i ./
done
gunzip FRAG*.fastq.gz

# Reads from pre-demultiplexed run by QB3 (Berkeley GSL)
# These reads have the .fq extension but have not been quality-trimmed and are 90nt single end
# For consistency with remaining data, cut these reads to 50nt later.
cp /cato2pool/backup-share/kramundson/Potato_dosage_backup/data/LOP-MM/LOP*.gz ./
gunzip LOP-2*
```

Save md5sums, read counts, and base counts of all untrimmed demultiplexed files before proceeding:

```
rm Tai_*.fq # rice reads pooled with LOP2 lane
rm LOP868_361.fq LOP868_398.fq LOP868_420.fq LOP868_439.fq LOP868_461.fq LOP868_474.fq LOP868_493.fq LOP868_510.fq
md5sum *.fq  FRAG*.fq >> md5sums_demult_LOP.txt
for i in *.fq ; do echo $i >> readcounts_raw_LOP.txt ; wc -l $i | awk '{print $0/4}' >> readcounts_raw_LOP.txt ; done
for i in FRAG* ; do echo $i >> readcounts_raw_LOP.txt ; wc -l $i | awk '{print $0/4}' >> readcounts_raw_LOP.txt ; done
md5sum *.* >> md5_all.txt # once this is done, start compressing .fastq files

# make link with .fq extension to get picked up by basecount.sh
for i in FRAG* ; do BASE=$(basename $i '.fastq') ; ln -s $i $BASE.fq ; done
scripts/basecount.sh
rm FRAG*.fq
gzip *.fastq
```

###3. Quality trim (uses Filter_N_Adapter_Trim_Batchmode_nounzip.py)

This would have taken forever in serial, so I made subdirectories and symlinked
raw .fq files with a .fastq extension so they would be recognized by a custom Python script I
inherited from Han. Trimming was then run behind screens in parallel using this script.

```
mkdir LOP1_trim LOP2_trim IVP_trim KA1_KA2_trim Han_trim LOP_berk_trim

cd LOP1_trim
for i in LOP868_{004..280}.fq ; do BASE=$(basename $i '.fq') ; ln -s ../$i $BASE.fastq ; done
file *.fastq | grep "broken" | cut -d ":" -f 1 | xargs rm
python ../scripts/Filter_N_Adapter_Trim_Batchmode_nounzip.py 40 > LOP1_trim.out 2> LOP1_trim.err

cd ../LOP2_trim
for i in LOP868_{281..538}.fq ; do BASE=$(basename $i '.fq') ; ln -s ../$i $BASE.fastq ; done
for i in S_LOP868_{361..510}.fq ; do BASE=$(basename $i '.fq') ; ln -s ../$i $BASE.fastq ; done
file *.fastq | grep "broken" | cut -d ":" -f 1 | xargs rm
for i in LOP868_{361,398,420,439,461,474,493,510}.fastq ; do rm $i ; done
python ../scripts/Filter_N_Adapter_Trim_Batchmode_nounzip.py 40 > LOP2_trim.out 2> LOP2_trim.err

cd ../IVP_trim
for i in ../IVP*.fq ; do BASE=$(basename $i '.fq') ; ln -s $i $BASE.fastq ; done
python ../scripts/Filter_N_Adapter_Trim_Batchmode_nounzip.py 40 -p > IVP_trim.out 2> IVP_trim.err

cd ../KA1_KA2_trim
for i in ../KA*.fq ; do BASE=$(basename $i '.fq') ; ln -s $i $BASE.fastq ; done
python ../scripts/Filter_N_Adapter_Trim_Batchmode_nounzip.py 40 -p > KA1_KA2_trim.out 2> KA1_KA2_trim.err

cd ../Han_trim
ln -s ../FRAG*.fastq ./
python ../scripts/Filter_N_Adapter_Trim_Batchmode_nounzip.py 40 > Han_trim.out 2> Han_trim.err

cd ../LOP_berk_trim
for i in ../LOP-2* ; do BASE=$(basename $i '.fq'); ln -s $i $BASE.fastq ; done
python ../scripts/Filter_N_Adapter_Trim_Batchmode_nounzip.py 40 > LOP_berk_trim.out 2> LOP_berk_trim.err
```

###4. Prep control dataset for bin by sam (Alca Tarma R1 reads, cut to 50nt)

This will be needed later for bin-by-sam. Here, I uninterleaved paired-end reads from
tetraploid Alca Tarma, hard-cut reads to 50nt with a custom script, and quality trimmed
the hard cut reads using the same script as for the haploids.

```
cd ~/LOP_clean/
mkdir control_reads
cd control_reads
ln -s KA1.fq ./
../scripts/fastUninterleave.sh KA1.fq
rm KA1-2.fq KA1.fq
mv KA1-1.fq KA1-1.fastq
python ../scripts/readcut.py KA1-1.fastq 50
python ../scripts/Filter_N_Adapter_Trim_Batchmode_nounzip.py 40
```

###5. Tidy up before read mapping

```
cd ~/LOP_clean
gzip *.fastq # running
gzip *.fq # running

mkdir demultiplex_only multiplexed_reads all_trim_single all_trim_pair

cd all_trim_single
ln -s ../LOP1_trim/*.fq ./
ln -s ../LOP2_trim/*.fq ./
ln -s ../LOP_berk_trim/*.fq ./
ln -s ../control_reads/*.fq ./
ln -s ../Han_trim/*.fq ./

# play with pooling .fq from multiple libraries before and after mapping
cat LOP-280.fq LOP868_280.fq > LOP868_280_combined.fq
cat LOP-281.fq LOP868_281.fq > LOP868_281_combined.fq
cat LOP-286.fq LOP868_286.fq > LOP868_286_combined.fq
cat LOP-289.fq LOP868_289.fq > LOP868_289_combined.fq

cd ../all_trim_pair
ln -s ../KA1_KA2_trim/*.fq ./
ln -s ../IVP_trim/IVP-101.fq ./
ls -l *.fq > LOP_raw_filesizes.txt


cd ../

mv LOP-1* multiplexed_reads
mv H718P-Amundson_S46_L005_R* multiplexed_reads/
mv H821P* multiplexed_reads/
mv KA-1-2* multiplexed_reads/
mv FRAG* demultiplex_only
mv IVP-101.fq.gz demultiplex_only # moved early
mv IVP-35.fq.gz demultiplex_only # moved early
mv *.fq.gz demultiplex_only # move remaining
```

###6. Read mapping: Paired-end reads for parents. Genome is DM1-3 v404 assembly from Michigan State University, described in Hardigan et al (2016)

```
mkdir genome
cd genome
wget http://solanaceae.plantbiology.msu.edu/data/potato_dm_v404_all_pm_un.fasta.zip
unzip potato_dm_v404_all_pm_un.fasta.zip
bwa index potato_dm_v404_all_pm_un.fasta

cd all_trim_pair
../scripts/bwa-doall-vModules-current.py -d ../genome/potato_dm_v404_all_pm_un.fasta \
-m ps -M -t 12 > paired_reads_map.out 2> paired_reads_map.err

# Overamp was not run for bam files, so I sorted bams for paired end reads
for i in *.bam
do
    $BASE=$(basename $i '.bam') ; echo $BASE ; samtools sort $i $BASE.sorted
done

####SINGLE END READS WERE MAPPED ON CABERNET####

# Uses the same version of bwa 0.7.5a r405
# Slurm sbatch scripts are provided in scripts folder
```

###6.1 Commands used for moving to cabernet, running bwa-doall, and moving back to isner:

```
###DONE IN CABERNET###
cd /share/comailab/
mkdir LOP_haploid_aln
scp -r kramundson@isner.genomecenter.ucdavis.edu:~/LOP_clean/all_trim_single/ ./LOP_haploid_aln/
scp -r kramundson@isner.genomecenter.ucdavis.edu:~/LOP_clean/genome/potato_dm_v404_all_pm_un.fasta
cd LOP_haploid_aln
mkdir map{1..10}

# divide in to 10 groups of 18 .fq files for mapping
for i in {1..10}
do
    ls *.fq | head -n 18 | xargs mv map$i
done

cd ~/scripts
sbatch bwa-doall-LOP-1.slurm
for i in {2..10}
do
    sbatch bwa-doall-LOP-$i.slurm
done

###DONE IN ISNER###
cd ~/LOP_clean/all_trim_single
mkdir fq bai bam sai sam usam
cd bam
scp kramundson@cabernet.genomecenter.ucdavis.edu:/share/comailab/kramundson/LOP_haploid_aln/all_trim_single/map*/*.bam ./


# convert bam to sam
for i in *.bam
do
    BASE=$(basename $i '.bam')
    samtools view -h $i > $BASE.sam
done

# run overamp-5.py on all sam
for i in *.sam
do
    ../../scripts/overamp-5.py -m s -f $i -o u$i
done

# convert usam to ubam and mv to new folder
for i in u*.sam
do
    BASE=$(basename $i '.sam')
    samtools view -bS $i > $BASE.bam
done
mkdir ../ubam
mv u*.bam ../ubam
cd ../ubam

# Merge bam files corresponding to different libraries of same biological sample
samtools merge -r uLOP868_1192_aln.sorted.bam uFRAG01514_S41_L005_R1_001_aln.sorted.bam uFRAG01518_S45_L005_R1_001_aln.sorted.bam &
samtools merge -r uLOP868_1180_aln.sorted.bam uFRAG01515_S42_L005_R1_001_aln.sorted.bam uFRAG01519_S46_L005_R1_001_aln.sorted.bam &
samtools merge -r uLOP868_1190_aln.sorted.bam uFRAG01516_S43_L005_R1_001_aln.sorted.bam uFRAG01520_S47_L005_R1_001_aln.sorted.bam &
samtools merge -r uLOP868_1157_aln.sorted.bam uFRAG01517_S44_L005_R1_001_aln.sorted.bam uFRAG01521_S48_L005_R1_001_aln.sorted.bam &
samtools merge -r uLOP868_280_merged_aln.sorted.bam uLOP868_280_aln.sorted.bam uLOP-280_aln.sorted.bam &
samtools merge -r uLOP868_281_merged_aln.sorted.bam uLOP868_281_aln.sorted.bam uLOP-281_aln.sorted.bam &
samtools merge -r uLOP868_286_merged_aln.sorted.bam uLOP868_286_aln.sorted.bam uLOP-286_aln.sorted.bam &
samtools merge -r uLOP868_289_merged_aln.sorted.bam uLOP868_289_aln.sorted.bam uLOP-289_aln.sorted.bam &
```

###6.5 Read mapping using updated bwa-doall script from Meric (not done, consider doing).

```
cp /isner/share/scripts/bwa-doall-CURRENT.py ~/LOP_clean/scripts/
mkdir genome2
cd genome2
cp ../genome/potato_dm_v404_all_pm_un.fasta ./
```

###7. Call variants (samtools run-mpileup.py pipeline, from Meric)

Dihaploids

```
cd ~/LOP_clean/all_trim_single/ubam/
python ../../scripts/run-mpileup.py -s /usr/bin/samtools -o LOP_haploid_mpup.txt \
-r ../../genome/potato_dm_v404_all_pm_un.fasta \
> LOP_haploid_mpup.out 2> LOP_haploid_mpup.err
```

Parents:

```
cd ~/LOP_clean/all_trim_pair/ubam
python ../../scripts/run-mpileup.py -s /usr/bin/samtools -o LOP_parent_mpup.txt \
-r ../../genome/potato_dm_v404_all_pm_un.fasta \
> LOP_haploid_mpup.out 2> LOP_haploid_mpup.err
```

###8. More tidying: .out and .err to log folder

```
mkdir log
mv *.out *.err log
mv *trim*/*.out log
mv *trim*/*.err log
```