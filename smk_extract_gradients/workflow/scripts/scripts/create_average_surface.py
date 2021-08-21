import nibabel as nib
import pandas as pd
import numpy as np

gii = nib.load(snakemake.input[0])
arr = gii.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0]
vertices = arr.data
midsurf_concat = np.zeros((
    vertices.shape[0],
    vertices.shape[1],
    len(snakemake.input))
)

for s in range(0,len(snakemake.input)):
    g = nib.load(snakemake.input[s])
    a = g.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0]
    v = a.data
    midsurf_concat[:,:,s] = v

midsurf_avg = np.mean(midsurf_concat,2)

gii.vertices = midsurf_avg
nib.save(gii, snakemake.output.surf)