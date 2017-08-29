Clean LOP analysis

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

Haploids normalized to 1 ng/ul based on SYBR Green I quantification of cleaned post-amplification library
81 haploids pooled per lane and sequenced on Illumina HiSeq 4000 at UC Davis DNA Tech Core.
One lane sequenced per pool of 81.

PL-4 and Alca Tarma subject to double-cut size selection with SeraMag beads.
Parents PL-4 and Alca Tarma pooled 1:1 based on qPCR quant (qPCR and pooling done by DNA tech core)
and sequenced on HiSeq4000 at UCD DNA Tech core. 

One lane for PL-4 and Alca Tarma.
Haploid inducer parent IVP-101 was sequenced 150PE as part of a different sequencing run.
Library of IVP-101 was prepared by David Castillo as previously described.

See file md5_all.txt for checksums of all data and analysis files

###2. Demultiplex without quality trim (allprep12.py). Save demultiplexed reads for SRA submission

```
# LOP haploids, lane 1
gunzip LOP-1_S45_L004_R1_001.fastq.gz LOP-1_S45_L004_R2_001.fastq.gz
python ./scripts/allprep-12.py -b LOPHAP1_barcodes.txt \
-f LOP-1_S45_L004_R1_001.fastq -i LOP-1_S45_L004_R2_001.fastq -m -D -n \
> LOP1_demult.out 2> LOP1_demult.err
gzip LOP-1_S45_L004_R1_001.fastq LOP-1_S45_L004_R2_001.fastq

# LOP haploids, lane 2
gunzip H718P-Amundson_S46_L005_R*
python ./scripts/allprep-12.py -b modified-LOPHAP2_barcodes.txt \
-f H718P-Amundson_S46_L005_R1_001.fastq -i LOP-1_S45_L004_R2_001.fastq \
-I H718P-Amundson_S46_L005_R3_001.fastq -m -D -n \
> LOP2_demult.out 2> LOP2_demult.err
gzip H718P-Amundson_S46_L005_R*

# LOP parents, Alca Tarma/PL-4 lane
gunzip KA-1-2_S26_L002_R1_001.fastq.gz KA-1-2_S26_L002_R2_001.fastq.gz KA-1-2_S26_L002_R4_001.fastq.gz
python ./scripts/allprep-12.py -b KA1-KA2-barcodes_fwdonly.txt
-f KA-1-2_S26_L002_R1_001.fastq -r KA-1-2_S26_L002_R4_001.fastq \
-i KA-1-2_S26_L002_R2_001.fastq -m -D -n
gzip KA-1-2_S26_L002_R*.fastq

# LOP parents, IVP-101/IVP-35 lane
gunzip H821P-Amudson_S34_L005_R1_001.fastq.gz H821P-Amudson_S34_L005_R2_001.fastq.gz H821P-Amudson_S34_L005_R4_001.fastq.gz
python ./scripts/allprep-12.py -b IVP35-IVP101-barcodes.txt \
-f H821P-Amudson_S34_L005_R1_001.fastq -r H821P-Amudson_S34_L005_R4_001.fastq \
-i H821P-Amudson_S34_L005_R2_001.fastq -m -D -n
gzip H821P-Amudson_S34_L005_R*.fastq
```

###3. Quality trim (Filter_N_Adapter_Trim_Batchmode_nounzip.py)

