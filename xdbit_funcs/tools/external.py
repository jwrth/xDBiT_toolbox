
def run_mofa(adata, output_file="ST_model.hdf5", 
             n_factors=10, features_subset="highly_variable", 
             obsm_key="spatial",
             seed=2021, frac_inducing=0.5):
    '''
    Run mefisto/mofa on adata.
    '''

    from mofapy2.run.entry_point import entry_point
    
    # Set up mofa options
    ent = entry_point()
    ent.set_data_options(use_float32=True)
    ent.set_data_from_anndata(adata, features_subset=features_subset)
    ent.set_model_options(factors=n_factors)
    ent.set_train_options()
    ent.set_train_options(seed=seed)
    
    ent.set_covariates([adata.obsm[obsm_key]], covariates_names=["imagerow", "imagecol"])
    ent.set_smooth_options(sparseGP=True, frac_inducing=frac_inducing,
                           start_opt=10, opt_freq=10)
    
    # build, run and save
    ent.build()
    ent.run()
    ent.save(output_file) #hdf5 file
