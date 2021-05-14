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

    Example:
    MakeTest.py -a /folder/test -b /folder/ref -p outprefix -n cfDNA,WBC
    """)

parser.add_argument('-a', '--test', help = "prefix + beta.txt.gz from makeTest.py", default = None, required=True)
parser.add_argument('-b', '--reference',  help = "output from makeTrain.py", default = None, required=True)
parser.add_argument('-n', '--normal', help = "comma-separated list of labels of the normal tissues in the reference dataset to exclude from tumor assignment e.g. cfdna,wbc", default = None, action=SplitArgs)
parser.add_argument('-p', '--outprefix', help = "prefix for output files", default = None, required = True)
args = parser.parse_args()

inputFile_samples = args.test
inputFile_tumor_atlas = args.reference
normal = args.normal
prefix = args.outprefix

## Set the working directory
#os.chdir('/Users/rmvpaeme/Repos/cfRRBS_paper/code/')

# Input files (= outputfiles from MakeTrainTest.py)
#inputFile_tumor_atlas = "./classifySamples/resources/train_plasma_WGBS_as_normal_MANUSCRIPT.gz"
#inputFile_samples = "./classifySamples/resources/20190323_test_beta_plasma_MANUSCRIPT.gz"
#prefix = "test"
#normal = "normal"
#python runMeth_atlas.py -a ./classifySamples/resources/20190323_test_beta_plasma_MANUSCRIPT.gz -b ./classifySamples/resources/train_plasma_WGBS_as_normal_MANUSCRIPT.gz -n normal -p examplerun
#%%
outputFile_atlas = prefix + "_train_methatlas.csv"
outputFile_samples = prefix + "_test_methatlas.csv"
outDir = "./classifySamples/output/classification"

if not os.path.exists(outDir):
    os.makedirs(outDir)

classificationResults = prefix + "_classificationResults_methAtlas.csv"

print("""
        Running deconvolution with NNLS (meth_atlas)...
            - reference matrix: %s
                saving collapsed reference matrix as %s/%s
            - test matrix: %s
                saving reformated test matrix as %s/%s
            - full deconvolution results available at %s/%s
            - tumor prediction available at %s/%s (removal of %s)
        """ % (inputFile_tumor_atlas, outDir, outputFile_atlas, 
                inputFile_samples, outDir,outputFile_samples,
                outDir,outputFile_samples.split('.')[0] + "_deconv_output.csv",
                outDir,classificationResults,
                normal))
#%%
df_tumor = pd.read_csv(inputFile_tumor_atlas, sep="\t", header = None, index_col = None, compression = "gzip")
df_tumor = df_tumor.transpose()
df_tumor.columns = df_tumor.iloc[0]
df_tumor = df_tumor.reindex(df_tumor.index.drop(0))
df_tumor = df_tumor.astype('float64')
df_tumor = df_tumor.groupby(by=df_tumor.columns, axis=1).median()
df_tumor['IlmnID'] = df_tumor.index
df_tumor['IlmnID'] = 'cg' + df_tumor['IlmnID'].astype(str)
df_tumor = df_tumor.set_index("IlmnID")

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
results_tumor.to_csv("%s/%s" % (outDir,classificationResults), header=True, index = True, sep=',', mode = 'w')
