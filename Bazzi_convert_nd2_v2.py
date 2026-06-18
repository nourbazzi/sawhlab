import nd2
import numpy as np
import os
os.makedirs('./sequential', exist_ok=True)
f = nd2.ND2File('/scratch/bazzi/20260612_ANS9_mNG_GAM488_R2_denoised/20260612_ANS9_mNG_GAM488_R2_crop_1.nd2')
arr = f.asarray()
num_fov = arr.shape[0]
HYB = 0
CHANNEL_NAMES = ['SD561', 'SD488', 'SD405']
for fov in range(num_fov):
    for ch in range(3):
        stack = arr[fov, :, ch, :, :]
        np.save(f'./sequential/{CHANNEL_NAMES[ch]}_{HYB}_{fov+1}.npy', stack)
        print(f'Saved {CHANNEL_NAMES[ch]} FOV {fov+1}')
f.close()
print('Done')
