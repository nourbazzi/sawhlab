#!/bin/bash -l
#SBATCH --job-name=matlab_test
#SBATCH --account=def-hilsawh 
#SBATCH --time=5:00:00        
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem=8000             # adjust this according to the memory requirement per node you need
#SBATCH --mail-user=nat.kwong@mail.utoronto.ca # adjust this to match your email address
#SBATCH --mail-type=ALL

module load matlab/2023b.2
matlab -singleCompThread -batch "dnaMERFISH_decode_new(1,'codesonly.csv')"
