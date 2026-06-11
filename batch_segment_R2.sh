#!/bin/bash

#SBATCH --time=0:35:00
#SBATCH --account=def-hilsawh
#SBATCH --gpus-per-node=1
#SBATCH  --cpus-per-task=4
#SBATCH --mem=80G



module load python/3.11 gcc/12.3 opencv/4.10.0
source $HOME/cellpose/bin/activate

python -m cellpose --dir /scratch/bazzi/sampleTiffs_R2 --anisotropy 1 --do_3D  --diameter 15 --verbose --pretrained_model cpsam --save_tif --use_gpu --no_npy --flow3D_smooth 0 --cellprob_threshold 0 --savedir /scratch/bazzi/masks_R2_diam15
