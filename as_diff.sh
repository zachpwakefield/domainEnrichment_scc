#!/bin/bash -l

#########################################################################
# Name: Zach Wakefield
# Date: Nov 14, 2023
# Description: protein impact backbone
#########################################################################
# Set SCC project to charge
#$ -P hybrids

#$ -pe omp 28

# Specify hard time limit for the job.
#   The job will be aborted if it runs longer than this time.
#   The default time, also selected here, is 12 hours.
#   You can increase this up to 720:00:00 for single processor jobs but your job will take longer to start.
#$ -l h_rt=5:00:00


# Name job
#$ -N pib

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

# Specify email
#$ -M zachpwakefield@gmail.com

# Join error and output streams in one file
#$ -j y

## Protein Change Bash Script pipeline

## Protein Change Bash Script pipeline
out_directory='/projectnb2/evolution/zwakefield/proteinImpacts/mpc_cpc/'
mkdir $out_directory
## SCC Settings
cd /projectnb2/evolution/zwakefield/proteinImpacts/

module load R/4.1.2
Rscript ./pipe.R '/projectnb2/evolution/zwakefield/proteinImpacts/' 'pc' 'AFE' '3' '3' $out_directory "/projectnb2/evolution/zwakefield/Finished_Projects/CardiacDifferentiation/HIT_Stat_Xingpei_Analysis/HIT_Stat_Analysis/mpc_cpc" '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/4-1744_T35/hit/' '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/5-1748_T35/hit/' '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/6-1749_T35/hit/' '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/7-1744_T21/hit/' '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/8-1748_T21/hit/' '/projectnb/evolution/zwakefield/proteinChange/pipeline_6_13/samples/9-1749_T21/hit/'
#

module load java/16.0.2 gcc/8.3.0 python3/3.8.10

cd
for file in $out_directory/bgoutFast.fa $out_directory/fgoutFast.fa
do
  # ./my_interproscan/interproscan-5.65-97.0/interproscan.sh -i ${file} -f tsv -goterms -T $TMPDIR
  ./my_interproscan/interproscan-5.65-97.0/interproscan.sh -i ${file} -b $out_directory/${file} -goterms -T $TMPDIR
done


module unload java/16.0.2 gcc/8.3.0 python3/3.8.10
module load R/4.1.2

Rscript ./domainEnrichment.R $out_directory
