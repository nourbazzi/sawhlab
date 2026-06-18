import nd2
import numpy as np
import os
os.makedirs('sequential', exist_ok=True)
f = nd2.ND2File('260320_ANS9_no_antibody_Denoised.nd2')
print('Channel names:', f.metadata.channels)
arr = f.asarray()
print('Shape:', arr.shape)
f.close()
