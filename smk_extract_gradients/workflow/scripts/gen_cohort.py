import os
import pandas as pd

scan = ['sub-' + snakemake.wildcards.subject]
    
if snakemake.params.session != False:
    scan.append('ses-' + snakemake.params.session[0])

# if snakemake.params.run != False:
#     scan.append('run-' + snakemake.params.run[0])

scan.append(str(snakemake.input.rfmri).replace(snakemake.params.fmriprep_path,'').strip('/'))

gen_cohort_header = lambda cohort_row: [ 'id' + str(idx) if idx != len(scan)-1 else 'img' for idx, list_item in enumerate(cohort_row) ]

df_scan = pd.DataFrame([scan], columns = gen_cohort_header(scan))

if os.path.exists(snakemake.output.cohort):
    df_cohort = pd.read_csv(snakemake.input.cohort, index_col=0)
    df_cohort.append(df_scan).to_csv(snakemake.output.cohort,index=0)
else:
    df_scan.to_csv(snakemake.output.cohort,index=0)