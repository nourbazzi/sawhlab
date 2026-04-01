#!/bin/bash

#SBATCH --time=0:20:00
#SBATCH --account=def-hilsawh
#SBATCH --gpus-per-node=1
#SBATCH  --cpus-per-task=4
#SBATCH --mem=80G



module load python/3.11 gcc/12.3 opencv/4.10.0
source $HOME/cellpose/bin/activate

python -m cellpose --dir /home/bazzi/scratch/sampleTiffs --anisotropy 1.82 --do_3D --diameter 38 --min_size 4000 --verbose --pretrained_model nuclei --save_tif --use_gpu --no_npy --savedir /home/bazzi/scratch/masks
