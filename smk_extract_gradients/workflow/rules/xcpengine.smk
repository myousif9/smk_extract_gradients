from scripts.utilities import get_friprep_dir, gen_cohort
import os
import pandas as pd

rule gen_cohort:
    input:
        rfmri = config['input_path']['bold_volume'],
    params:
        fmriprep_dir = get_friprep_dir(config['input_path']['bold_volume']), 
        subject = wildcards.subject,
        run =  config['input_lists']['bold_volume']['run']  if 'run' in config['input_lists']['bold_volume'] else False,
        session =  config['input_lists']['bold_volume']['session'] if 'session' in config['input_lists']['bold_volume']['session'] else False
    output:
        cohort = bids(
            root = "work",
            run = '{run}',
            suffix = 'cohort.csv',
            **config['subj_wildcards']
        )
    group: 'subj'
    run:
        scan = [wildcards.subject]
       
        if params.session != False:
            scan.append(params.session)

        if params.run != False:
            scan.append(params.run)
        
        scan.append(str(input.rfmri).replace(params.fmriprep_dir,'').strip('/'))
        
        gen_cohort_header = lambda cohort_row: [ 'id'+str(count) for count in range(0:len(cohort)-1)].append('img')
        
        df_scan = pd.DataFrame([scan], columns = gen_cohort_header(scan))
        
        if os.path.exists(output.cohort):
            df_cohort = pd.read_csv(input.cohort,)
            df_cohort.append(df_scan).to_csv(output.cohort,index=0)
        else:
            df_scan.to_csv(output.cohort,index=0)
        
rule run_xcpengine:
    input:
        rfmri = config['input_path']['bold_volume'],
        cohort = bids(
            root = "work",
            run = '{run}',
            suffix = 'cohort.csv',
            **config['subj_wildcards']
        )
    params:
        fmriprep_dir = get_friprep_dir(config['input_path']['bold_volume']),
        pipeline_design = '../resources/fc-36p.dsn'
    output:
        work_dir = directory(
            bids(
                root = "work",
                dataype = 'xcpengine'
                **config['subj_wildcards']
            )
        ),
        output_dir = directory(
            bids(
                root = "result",
                dataype = 'xcpengine'
                **config['subj_wildcards']
            )
        )
    container: config['singularity']['xcpengine']
    group: 'subj'
    shell:
        """
        xcpengine \
        -d {param.pipeline_design} \
        -c {input.cohort} \
        -r {params.fmriprep_dir} \
        -i {output.work_dir} \
        -o {output.output_dir}
        """