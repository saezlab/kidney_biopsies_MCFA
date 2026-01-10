rule create_mdata: # Done
    input:
        h5ad_input = config['input']
        coldata_path = config['metadata']
    output:
        h5ad_output = 'results/MCFA/adata.h5ad',
        h5mu_output = 'results/MCFA/adata.h5mu'
    resources:
        mem_mb=100000
    wildcard_constraints:
        subset = 'includingKPMP'
    params:
        sample_key = config["sample_key"],
        condition_key = config["condition_key"],
        min_cells = config["min_cells"],
        groupby = config["groupby"],
        sample_subset = config["sample_subset"],
        drop_views = config["drop_views"]
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/create_mdata.py'

rule run_liana:
    input:
        h5ad_input = 'results/MCFA/adata.h5ad'
    output:
        h5ad_output = 'results/MCFA/adata_LIANA.h5ad'
    resources:
        mem_mb=100000
    params:
        sample_key = config["sample_key"],
        groupby = config["groupby"]
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/run_liana.py'

rule create_liana_mdata:
    input:
        h5ad_input = 'results/MCFA/adata_LIANA.h5ad'
    output:
        h5mu_output = 'results/MCFA/adata_LIANA.h5mu'
    resources:
        mem_mb=100000
    params:
        sample_key = config["sample_key"],
        groupby = config["groupby"]
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/create_liana_mdata.py'

rule combine_liana_mdata:
    input:
        h5mu_liana = 'results/MCFA/adata_LIANA.h5mu',
        h5mu_input = 'results/MCFA/adata.h5mu'
    output:
        h5mu_output = 'results/MCFA/adata_combined.h5mu'
    params:
        hvg = 'all'
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/combine_liana_mdata.py'
        
rule run_mofa_subset: # Done
    input:
        h5mu_input = 'results/MCFA/adata_combined.h5mu'
    output:
        model_output = 'results/MCFA/model~{hvg}~{nfactors}.h5ad',
        scores_output = 'results/MCFA/scores~{hvg}~{nfactors}.csv',
        scores_wide_output = 'results/MCFA/widescores~{hvg}~{nfactors}.csv',
        loadings_output = 'results/MCFA/loadings~{hvg}~{nfactors}.csv'
    wildcard_constraints:
        hvg = 'all|hvg'
    params:
        nfactors = lambda wildcards: wildcards.nfactors,
        hvg = lambda wildcards: wildcards.hvg,
        sample_key = config["sample_key"]
    singularity:
        config["singularities"]["mofaxmetacells"]
    script:
        '../scripts/MCFA/run_mofa_subset.py'

rule calculate_r2_subset: # Done, updated 21.02.24
    input:
        model_input = 'results/MCFA/model~{hvg}~{nfactors}.h5ad'
    output:
        variance_explained_by_sample_output = 'results/MCFA/r2_by_sample~{hvg}~{nfactors}.csv',
        variance_explained_total_output = 'results/MCFA/r2_total~{hvg}~{nfactors}.csv'
    params:
        condition_key = config["condition_key"]
    singularity:
        config["singularities"]["scanpy_update"]
    script:
        '../scripts/MCFA/calculate_r2.py'

rule combine_scores:
    input:
        expand('results/MCFA/scores~{{hvg}}~{nfactors}.csv', nfactors=config["nfactors_list"])
    output:
        'results/MCFA/robustness_scores~{hvg}.csv'
    wildcard_constraints:
        hvg = 'all|hvg'
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/combine.py'

rule combine_loadings:
    input:
        expand('results/MCFA/loadings~{{hvg}}~{nfactors}.csv', nfactors=config["nfactors_list"])
    output:
        'results/MCFA/robustness_loadings~{hvg}.csv'
    wildcard_constraints:
        hvg = 'all|hvg'
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/combine.py'

rule combine_r2_by_sample:
    input:
        expand('results/MCFA/r2_by_sample~{{hvg}}~{nfactors}.csv', nfactors=config["nfactors_list"])
    output:
        'results/MCFA/robustness_r2_by_sample~{hvg}.csv'
    wildcard_constraints:
        hvg = 'all|hvg'
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/combine.py'

rule combine_r2_total:
    input:
        expand('results/MCFA/r2_total~{{hvg}}~{nfactors}.csv', nfactors=config["nfactors_list"])
    output:
        'results/MCFA/robustness_r2_total~{hvg}.csv'
    wildcard_constraints:
        hvg = 'all|hvg'
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/combine.py'


rule factor_robustness:
    """Compute factor similarity (product of score and loading correlations)."""
    input:
        scores='results/MCFA/robustness_scores~all.csv',
        loadings='results/MCFA/robustness_loadings~all.csv'
    output:
        factor_robustness='results/MCFA/factor_robustness.csv'
    singularity:
        config["singularities"]["renv"]
    script:
        '../scripts/MCFA/factor_robustness.R'

rule factors_obs:
    input:
        scores_wide = 'results/MCFA/widescores~all~5.csv',
        metadata_csv = 'data/metadata/patient_metadata_EDA.csv',
        clr = 'results/composition/clr~main_class.csv'
    output:
        factors_obs_csv = 'results/MCFA/factors_obs~all~5.csv'
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/factors_obs.py'

rule collectri:
    input:
        loadings = 'results/MCFA/loadings~all~5.csv'
    output:
        acts = 'results/MCFA/acts_collectri.csv'
    singularity:
        config["singularities"]["scanpy_update"]
    script:
        '../scripts/MCFA/collectri.py'

rule progeny:
    input:
        loadings = 'results/MCFA/loadings~all~5.csv'
    output:
        acts = 'results/MCFA/acts_progeny.csv'
    singularity:
        config["singularities"]["scanpy_update"]
    script:
        '../scripts/MCFA/progeny.py'

rule factor_enrichment:
    input:
        loadings = 'results/MCFA/loadings~all~5.csv'
    output:
        enrichment_output = 'results/MCFA/hallmark_Factor2.csv'
    params:
        gene_sets = config["gene_sets"]
    singularity:
        config["singularities"]["scanpy"]
    script:
        '../scripts/MCFA/factor_enrichment.py'

rule metadata_associations:
    input:
        factors_obs = 'results/MCFA/factors_obs~all~5.csv'
    output:
        metadata_pearson_correlations = 'results/MCFA/metadata_pearson_correlations.csv',
        metadata_anova = 'results/MCFA/metadata_anova.csv',
        clr_pearson_correlations = 'results/MCFA/clr_pearson_correlations.csv'
    singularity:
        config["singularities"]["renv"]
    script:
        '../scripts/MCFA/metadata_associations.R'