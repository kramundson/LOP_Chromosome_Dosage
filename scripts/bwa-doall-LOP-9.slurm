#!/bin/bash
#
#SBATCH --job-name=LOP-bwa-9
#SBATCH -c 24 # number of cores
#SBATCH -N 1 # ensure that all cores are on one machine
#SBATCH --mem=2G # Memory pool for all cores in MB (see also --mem-per-cpu)
#SBATCH -p gc # partition to submit to
#SBATCH -o LOP-bwa-9.out # file to which STDOUT will be written
#SBATCH -e LOP-bwa-9.err # file to which STDERR will be written
#SBATCH --mail-type=ALL # type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kramundson@ucdavis.edu # Email to which notifications will be sent
#SBATCH --time=3-00:00:00

cd /share/comailab/kramundson/LOP_haploid_aln/all_trim_single/map9
python ~/scripts/bwa-doall-vModules-current.py -d /share/comailab/kramundson/LOP_haploid_aln/genome/potato_dm_v404_all_pm_un.fasta -O -t 24

exit
