#!/usr/bin/env nextflow

params.rna_data = '/data/projects/punim1597/Data/cite-seq_datasets/hao_data/hao_RNA.prepd.h5ad'
params.adt_data = '/data/projects/punim1597/Data/cite-seq_datasets/hao_data/hao_ADT.prepd.h5ad'
// more params?

params.num_alphas = 20
params.num_etas = 10
params.range_alphas = [0,0.3]
params.range_etas = [0.00001,0.001]

params.data_subset = 0.2
params.outdir = 'grid_search_out'
params.adt_list = ['CD8','CD45RA','CD16','CD11c','CD14','CD19','CD34']


process PREP_DATA {
    conda = '/home/danrawlinson/.conda/envs/scvi'
    memory '50 GB'
    time '10min'
    publishDir params.outdir

    input:
        path rna_data
        path adt_data
        val data_subset

    output:
        tuple path("train_rna.h5ad"),path("train_adt.h5ad"),path("test_rna.h5ad"),path("test_adt.h5ad")

    """
    #!/usr/bin/env python

    import scanpy as sp
    import random
    
    import importlib.util
    spec = importlib.util.spec_from_file_location("CitePred_functions", "/data/gpfs/projects/punim1597/Projects/CITE-seq/CitePred_functions.py")
    CitePred_functions = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CitePred_functions)


    hao_rna = sp.read_h5ad("$rna_data")
    hao_adt = sp.read_h5ad("$adt_data")

    hao_rna.raw = hao_rna
    hao_adt.raw = hao_adt

    random.seed(42)
    sampled_names = random.sample(hao_rna.obs_names.tolist(), round($data_subset * hao_rna.n_obs))

    sampled_rna = hao_rna[sampled_names,:]
    sampled_adt = hao_adt[sampled_names,:]

    data_splits = CitePred_functions.split_AnnData(sampled_rna, sampled_adt, renormalise = True, test_portion = 0.2)

    #train partition
    split_RNA_train = data_splits['train'][0].copy()
    split_ADT_train = data_splits['train'][1].copy()

    sp.pp.log1p(split_RNA_train)
    sp.pp.highly_variable_genes(split_RNA_train, n_top_genes = 10000)

    sp.pp.scale(split_RNA_train)
    sp.pp.scale(split_ADT_train)

    #test partition
    split_RNA_test = data_splits['test'][0].copy()
    split_ADT_test = data_splits['test'][1].copy()

    sp.pp.log1p(split_RNA_test)

    sp.pp.scale(split_RNA_test)
    sp.pp.scale(split_ADT_test)

    split_RNA_train.write_h5ad("train_rna.h5ad")
    split_ADT_train.write_h5ad("train_adt.h5ad")
    split_RNA_test.write_h5ad("test_rna.h5ad")
    split_ADT_test.write_h5ad("test_adt.h5ad")
    """
}

process ALPHA_GRID { //Just spit out all the different options and use the .combine() channel operator
    executor 'local'
    input:
        val num_alphas
        val range_alphas
    
    output:
        stdout
    script:
        first = range_alphas[0]
        last = range_alphas[1]
    """
    incr=`echo "scale=2;($last-$first)/($num_alphas-1)" | bc `
    seq $first \$incr $last
    """
}

 process ETA_GRID {
    executor 'local'
    input:
        val num_etas
        val range_etas

    output:
        stdout

    script:
        first = range_etas[0]
        last = range_etas[1]
    """
    incr=`echo "scale=6;($last-$first)/($num_etas-1)" | bc`
    seq $first \$incr $last 
    """
 }

process SGD_REGRESS {
    time '2hr'
    memory '50 GB'
    conda '/home/danrawlinson/.conda/envs/scvi'

    input:
        path adata_ready_to_model
        tuple val(alpha), val(eta)
    output:
        path trained_model

    """
    #!/usr/bin/env python

    import scanpy as sp
    from sklearn.linear_model import SGDRegressor
    
    SGDmodels = {}
    lasCV_scores = {}
    chunksize = 1000

    #couple of ideas:
        #run parameter search on smaller subset of this data
        #use snakemake to do grid search
        #add cross-fold here. so that each data point is being tested on a model that hasn't used it for training. The loss values of each prediction are just concatenated to give the overall accuracy.
            #cross-fold means I might have to change my DATA_PREP process
            #some options here https://scikit-learn.org/stable/modules/cross_validation.html

    for a in [adts_in_cbmc[0]]:
        SGDmodels[a] = SGDRegressor(loss='squared_error', penalty = 'l1', alpha = $alpha, learning_rate='constant', eta0=$eta)
        for i in range(0, hao_RNA_train.shape[0], chunksize):
            SGDmodels[a].partial_fit(X = hao_RNA_train[i:i+chunksize, hao_RNA_train.var['highly_variable']].X, y = hao_ADT_train[i:i+1000,a].X.ravel())
        #for i, x in enumerate(hao_RNA_train[:,hao_RNA_train.var['highly_variable']]): #one sample at a time for .partial_fit()
        #    SGDmodels[a].partial_fit(X = x.X.reshape(1,-1), y = hao_ADT_train[i,a].X.ravel())
        lasCV_scores[a] = SGDmodels[a].score( X =  hao_RNA_train[:,hao_RNA_train.var['highly_variable']].to_df(), y =  hao_ADT_train.to_df()[a])
        print(a, lasCV_scores[a])


    """
}
/*

*/



workflow { //define channels in here
alphas_ch = ALPHA_GRID(params.num_alphas, params.range_alphas).splitText().toFloat()
etas_ch = ETA_GRID(params.num_etas, params.range_etas).splitText().toFloat()
alphas_ch.combine(etas_ch)

prep_ch = PREP_DATA(params.rna_data, params.adt_data, params.data_subset)

}
