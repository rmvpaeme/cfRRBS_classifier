#%%
## Import libraries
import os
import numpy as np
import pandas as pd
import argparse
import deconvolve_resids
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

#%%
parser = argparse.ArgumentParser(
    """
    runMeth_atlas.py

    Runs metatlas on a testmatrix and reference matrix.
    Files will be written to ./classifySamples/output/classification

    With the default parameters (iterations = 1), the replicates are collapsed with mean/median (depending on the -c parameter).
    If iterations are > 1, a beta distribution is calculated from every CpG/region, and in every iteration the methylation (beta) values are random samples from this beta distribution, which gives an uncertainty around the tumor classification call (plot_error.R). 

    Example:
    runMeth_atlas.py -a /folder/test -b /folder/ref -p outprefix -n cfDNA,WBC
    """,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-a', '--test', help = "prefix + beta.txt.gz from makeTest.py", default = None, required=True)
parser.add_argument('-b', '--reference',  help = "output from makeTrain.py", default = None, required=True)
parser.add_argument('-n', '--normal', help = "comma-separated list of labels of the normal tissues in the reference dataset to exclude from tumor assignment (after the deconvolution step) e.g. cfdna,wbc", default = None, action=SplitArgs)
parser.add_argument('-x', '--exclude', help = "comma-separated list with tumor entities to exclude from the reference dataset (before the deconvolution step)", default = None, action=SplitArgs)
parser.add_argument('-c', '--collapse', help = "collapse replicate columns with mean or median", default = "median", required = False, choices=['mean', "median"])
parser.add_argument('-p', '--outprefix', help = "prefix for output files", default = None, required = True)
parser.add_argument('-f', '--fs', help = "enable feature selection by selecting the top n hyper and hypomethylated regions per entity", action='store_true', required = False)
parser.add_argument('-m', '--mrc', help = "top n hyper and hypomethylated regions will be selected with feature selection", default=100, required = False)
parser.add_argument('-i', '--iter', help = "number of iterations (every iteration samples random from the beta distribution for every CpG)", default=1, required = False)

args = parser.parse_args()

inputFile_samples = args.test
inputFile_tumor_atlas = args.reference
normal = args.normal
collapse = args.collapse
prefix = args.outprefix
feature_selection = args.fs
exclude = args.exclude
mrc = int(args.mrc)
iter = int(args.iter)

## Set the working directory
#os.chdir('/Users/rmvpaeme/Repos/cfRRBS_classifier_v0.5')

# Input files (= outputfiles from MakeTrainTest.py)
#inputFile_tumor_atlas = "./classifySamples/resources/train_plasma_WGBS_as_normal_MANUSCRIPT.gz"
#inputFile_samples = "./classifySamples/resources/20190323_test_beta_plasma_MANUSCRIPT.gz"
#prefix = "test"
#normal = "normal"
# python runMeth_atlas.py -a /Users/rmvpaeme/Repos/RVP2101_deconvolutionBenchmark/methAtlas/output/test_beta_benchmark_hg38.gz -b /Users/rmvpaeme/Repos/RVP2101_deconvolutionBenchmark/methAtlas/output/train_methatlas_benchmark_hg38.gz -n normal,wbc -p examplerun -i 50
# python runMeth_atlas.py -a ./classifySamples/resources/20190323_test_beta_plasma_MANUSCRIPT.gz -b ./classifySamples/resources/train_plasma_WGBS_as_normal_MANUSCRIPT.gz -n normal,wbc -p examplerun
#%%
outputFile_atlas = prefix + "_train_methatlas.csv"
outputFile_samples = prefix + "_test_methatlas.csv"
outDir = "./classifySamples/output/classification"

if not os.path.exists(outDir):
    os.makedirs(outDir)

classificationResults = prefix + "_classificationResults_methAtlas.csv"

print("""
        Running deconvolution with NNLS (meth_atlas)...
            - feature selection: %s 
            - reference matrix: %s
                %s will be excluded
                saving collapsed reference matrix as %s/%s
                replicates will be collapse with the %s method
            - test matrix: %s
                saving reformated test matrix as %s/%s
            - full deconvolution results will be saved to %s/%s
            - tumor prediction will be saved to %s/%s (highest fraction after removal of %s)
        """ % (feature_selection,
                inputFile_tumor_atlas, exclude, outDir, outputFile_atlas, 
                collapse,
                inputFile_samples, outDir,outputFile_samples,
                outDir,outputFile_samples.split('.')[0] + "_deconv_output.csv",
                outDir,classificationResults,
                normal))
#%%

for i in range(1, iter+1, 1):
    num = str(i)
    print("Iteration: %s" % str(i))
    outputFile_atlas = prefix + num +"_train_methatlas.csv"
    outputFile_samples = prefix + num + "_test_methatlas.csv"
    classificationResults = prefix + num + "_classificationResults_methAtlas.csv"

    df_tumor = pd.read_csv(inputFile_tumor_atlas, sep="\t", header = None, index_col = None, compression = "gzip")
    df_tumor = df_tumor.transpose()
    df_tumor.columns = df_tumor.iloc[0]
    df_tumor = df_tumor.reindex(df_tumor.index.drop(0))
    df_tumor = df_tumor.astype('float64')

    if iter > 1:
        appended_data = []
        for column in df_tumor.columns.unique():
            bootstrap = df_tumor[column]

            # fit with scipy, currently not used
            #def calcalpha(x):
            #    if x.isnull().all():
            #        fit = [np.nan, np.nan]
            #    else:
            #        fit = scipy.stats.beta.fit(x,floc=-0.5,fscale=1.5)
            #    return [fit[0],fit[1]]
            #params = bootstrap.progress_apply(lambda row: calcalpha(row.dropna()), axis = 1)

            m=bootstrap.mean(axis = 1)
            v=bootstrap.var(axis = 1)

            alpha=-m*(v+m*m-m)/v
            alpha[alpha <= 0] = np.nan

            beta=(m-1)*(v+m*m-m)/v
            beta[beta <= 0] = np.nan

            # uncomment below to get a sanity check on the fit of the beta distribution
            # for cpg in np.random.randint(1,14000, 10):
            #     import statsmodels.api as sm
            #     import pylab
            #     from scipy.special import gamma
            #     import scipy as sp
            #     import matplotlib.pyplot as plt

            #     fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(5, 3))
            #     a = alpha[cpg]
            #     b = beta[cpg]
            #     print("alpha = %f; beta = %f, cpg = %i, tissue = %s" % (a, b, cpg, column))
            #     def pdf(a, b, x):
            #         Z = gamma(a + b) / gamma(a) / gamma(b)
            #         return x**(a - 1) * (1 - x)**(b - 1) / Z
            #     x = np.arange(0.0, 1.01, 0.01)
            #     y = pdf(a, b, x)
            #     axes[0].plot(x, y)
            #     axes[1].hist(bootstrap.loc[cpg], x)
            #     fig.tight_layout()
            #     fig.savefig('beta_check_%s_cpg_%s.png' % (column, str(cpg)), transparent = False)

            sampled = pd.DataFrame(np.random.beta(alpha,beta))
            sampled = sampled.rename(columns = {0: column})
            sampled['IlmnID'] = sampled.index
            sampled['IlmnID'] = 'cg' + sampled['IlmnID'].astype(str)
            sampled = sampled.set_index("IlmnID")
            sampled = sampled.clip(lower=0)
            sampled = sampled.clip(upper=1)
            appended_data.append(sampled)

        df_tumor = pd.concat(appended_data, axis = 1)
        df_tumor.to_csv("%s/%s" % (outDir,outputFile_atlas), header=True, index = True, sep=',', mode = 'w')
 
    else:
        if collapse == "median":
            df_tumor = df_tumor.groupby(by=df_tumor.columns, axis=1).median()
        elif collapse == "mean":
            df_tumor = df_tumor.groupby(by=df_tumor.columns, axis=1).mean()
        df_tumor['IlmnID'] = df_tumor.index
        df_tumor['IlmnID'] = 'cg' + df_tumor['IlmnID'].astype(str)
        df_tumor = df_tumor.set_index("IlmnID")
        if exclude is not None:
            df_tumor = df_tumor.drop(exclude, axis = 1) # confounding entities

        if feature_selection == True:
            topup = (df_tumor.div(df_tumor.sum(axis=1), axis=0))
            markers_list = []
            for column in topup:
                markers_df = topup[[column]].sort_values(by=column, ascending = False).head(mrc)
                markers_list.append(list(markers_df.index.map(str)))
            topdown = ((1-df_tumor).div((1-df_tumor).sum(axis=1), axis=0))
            for column in topdown:
                markers_df = topdown[[column]].sort_values(by=column, ascending = False).head(mrc)
                markers_list.append(list(markers_df.index.map(str)))
            markers_list = list(set([item for sublist in markers_list for item in sublist]))
            markers_df = pd.DataFrame(markers_list)
            markers_df.columns = ["IlmnID"]
            markers_df = markers_df.set_index("IlmnID")
            df_tumor = pd.merge(df_tumor, markers_df, left_index=True, right_index=True, how='right')

        df_tumor.to_csv("%s/%s" % (outDir,outputFile_atlas), header=True, index = True, sep=',', mode = 'w')
        # df_tumor = pd.read_csv("./train_methatlas_plasma_manuscript.csv", sep=",", header = 0, index_col = 0) # uncomment this to reproduce the manuscript results

    #%%
    df_samples = pd.read_csv(inputFile_samples, sep="\t", header = None, index_col = None, compression = "gzip")
    df_samples = df_samples.transpose()
    df_samples.columns = df_samples.iloc[0]
    df_samples = df_samples.reindex(df_samples.index.drop(0))
    df_samples = df_samples.astype('float64')
    df_samples['IlmnID'] = df_samples.index
    df_samples['IlmnID'] = 'cg' + df_samples['IlmnID'].astype(str)
    df_samples = df_samples.set_index("IlmnID")
    df_samples = df_samples.reindex(sorted(df_samples.columns), axis=1)
    df_samples.to_csv("%s/%s" % (outDir,outputFile_samples), header=True, index = True, sep=',', mode = 'w')

    deconvolve_resids.Deconvolve(atlas_path="%s/%s" % (outDir,outputFile_atlas),samp_path="%s/%s" % (outDir,outputFile_samples),out_dir=outDir,plot=False,resid=True).run()

    results = pd.read_csv("%s/%s" % (outDir,outputFile_samples.split('.')[0] + "_deconv_output.csv"), sep=",", header = 0, index_col = 0)
    if normal is not None:
        results = results.drop(normal)
    results = results.replace(0, np.nan)
    results.sort_values(by=list(results.columns.values), inplace=True)

    dict = {}
    for column in results:
        dict[column] = [results[column].idxmax(), results[column].max()]
        results_tumor = pd.DataFrame.from_dict(dict, orient = "index")
    results_tumor = results_tumor.rename(columns={0: "Classification", 1: "Tumor burden"})
    results_tumor = results_tumor.replace(np.nan, 0)
    results_tumor = results_tumor.sort_index()
    results_tumor["iteration"] = num
    results_tumor.to_csv("%s/%s" % (outDir,classificationResults), header=True, index = True, sep=',', mode = 'w')

print("Done.")
