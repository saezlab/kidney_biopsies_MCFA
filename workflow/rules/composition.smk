rule clr:
    input:
        adata = config['input']
    output:
        clr_csv = 'results/composition/clr~{annotation}~{subset}.csv',
        prop_csv = 'results/composition/prop~{annotation}~{subset}.csv',
    resources:
        mem_mb=50000
    singularity:
        config["singularities"]["scanpy"]#'workflow/envs/scanpy.0.0.5.sif'
    script:
        '../scripts/composition/clr.py'