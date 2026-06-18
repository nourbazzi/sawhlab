import nd2
f = nd2.ND2File('/scratch/bazzi/20260612_ANS9_mNG_GAM488_R1_denoised/20260612_ANS9_mNG_GAM488_R1_crop.nd2')
arr = f.asarray()
print('shape:', arr.shape)
print('axes:', f.sizes)
f.close()
