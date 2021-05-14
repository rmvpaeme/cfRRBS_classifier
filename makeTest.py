#%%
## Import libraries
import glob
import os
import numpy as np
import pandas as pd
import subprocess
import argparse
import multiprocess
from multiprocess import Manager, Pool
import imports_cfRRBS_classifier as cfRRBS
cpuCount = (multiprocess.cpu_count() - 2)

class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

# testargs
#folder = "./classifySamples/testfiles"
#regions = "./classifySamples/resources/RRBS_450k_intersectClusters.tsv"
#prefix = "testscript"
#cutoff = 30

#%%

parser = argparse.ArgumentParser(
    """
    MakeTest.py

    Creates a testmatrix to use for NNLS with regions as defined by the user. Takes Bismark coverage files as input.

    Example:
    MakeTest.py -n /folder/ -r regions.tsv -o testmatrix.txt
    """)

parser.add_argument('-f', '--folder', default = None, required=True)

parser.add_argument('-r', '--regions', default = None, required=True)
parser.add_argument('-c', '--cutoff', default = 30)

parser.add_argument('-p', '--outprefix', default = None)
args = parser.parse_args()

folder = args.folder
regions = args.regions
cutoff = args.cutoff
prefix = args.outprefix

print("""Running makeTest.py
            - Bismark coverage folder: %s
                Files found: %s
            - Regions: %s
            - Cutoff: %i
            - output: %s
    """ % (folder, glob.glob(os.path.join(folder, "*.cov.gz")), regions, cutoff, prefix))

#%%
# Folder to store intermediate files
tmp_folder = "/tmp/"

clusters = cfRRBS.import_clusters(regions)[0]
clusterFile = cfRRBS.import_clusters(regions)[1]

# Name of output
testMethyName = prefix + "_methy"
testDepthName = prefix + "_depth"
testBetaName = prefix + "_beta"

# Folder to store intermediate files
tmp_folder = "/tmp/"

## The location of the cfRRBS files, after running the preprocessing pipeline
test_files = glob.glob(os.path.join(folder, "*.cov.gz"))

#%%
def import_test_files(x):
        # Goal: to obtain one file, containing all the test files, where the first column is all the samples and the rest of the column either the beta values, total depth or # methylated reads for that cluster.
        # 1. Read in the bismark coverage file and convert them to a sort-of bed file, so that it can be manipulated with bedtools.
        file = x
        file_name = os.path.splitext(os.path.basename(file))[0]
        df = pd.read_csv(file, sep="\t",usecols=[0,1,2,3,4,5], header=None, compression = "gzip", low_memory=False)
        df[3] = df[3]/100    # From methylation percentage to methylation ratio
        df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

        # 2. Intersect between the cluster file and the bismark coverage file
        outfile = open(tmp_folder + '%s_intersect.txt' % file_name, 'w')
        print("     Running bedtools intersect on %s.txt..." % file_name)
        arg = "bedtools intersect -wb -b %s/%s.txt -a %s%s.txt" % (tmp_folder, clusterFile, tmp_folder, file_name)
        arg = arg.split()
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()
        df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols=[6,7,8,3,4,5,9], header=None, low_memory=False) # The previous step shuffles the column order, so this step rearranges the column order
        df = df[[6,7,8,3,4,5,9]] # chr, start, stop, beta value, count methylated, count unmethylated, cluster number
        df.to_csv(tmp_folder + "%s_reordered.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

        # 3. Group all the rows that are within a cluster, and get the sum of the methylated and unmethylated values. In addition, get the row with the CpG cluster number for future indexing.
        arg = "bedtools groupby -i %s%s_reordered.txt -g 1-3,7 -c 5,6 -o sum" % (tmp_folder, file_name)
        arg = arg.split()
        outfile = open(tmp_folder + '%s_clustered.txt' % file_name, 'w')
        print("     Running bedtools groupby on %s.txt..." %file_name)
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        # 4. Remove all clusters that have less than 30 reads
        df = pd.read_csv(tmp_folder + '%s_clustered.txt' % file_name, sep="\t", header=None, index_col = 3, low_memory=False)
        df.index.name = None    # Remove index.name for consistency
        df[6] = df[4]/(df[4] + df[5])   # Get beta value per cluster
        df = df[[0,1,2,6,4,5]] # Reorder the columns in chr, start, stop, beta value, no methylated, no unmethylated
        df.sort_values(by=[0,1,2], inplace=True) # Sort by chromosome
        df[7] = df[4] + df[5]   # Get total depth (=methylated + unmethylated count)
        df[[7,6,4,5]] = df[[7,6,4,5]].mask(df[7] < cutoff)  # Mark all clusters lower than 30 reads with NA
        print("The amount of clusters in %s remaining after removing NA: %s" % (file_name, len(df[7].replace(np.nan, 'NA').apply(pd.to_numeric, errors='coerce').dropna())))
        df = df.replace(np.nan, 'NA', regex=True) # Replace numpy NaN with NA string
        print("     Extracting %s and %s from file %s.txt..." % (testMethyName, testDepthName, file_name))

        # 5. Add the first file to the testMethy_list ( = number of methylated CpGs per cluster)
        testMethy_df = df
        testMethy_df.columns = [0,1,2,6, "%s" % file_name, 5,7]
        testMethy_df = testMethy_df.drop([0,1,2,6,5,7], axis = 1).astype(str)
        testMethy_df[file_name] = testMethy_df[file_name].apply(lambda x: x.split('.')[0]) # Make integer of float numbers, but as the dataframe is astype(str), we need to do this with a lambda function
        testMethy_list.append(testMethy_df)

        # Identical to testMethy
        testDepth_df = df
        testDepth_df.columns = [0,1,2,5,3,4,"%s" % file_name]
        testDepth_df = testDepth_df.drop([0,1,2,3,4,5], axis = 1).astype(str)
        testDepth_df[file_name] = testDepth_df[file_name].apply(lambda x: x.split('.')[0])
        testDepth_list.append(testDepth_df)

        # Make a new variable for visualisation that contains the beta values of the clusters
        testBeta_df = df
        testBeta_df.columns = [0,1,2, "%s" % file_name, 4, 5,7]
        testBeta_df = testBeta_df.drop([0,1,2,4,5,7], axis = 1).astype(str)
        testBeta_list.append(testBeta_df)

        os.remove(tmp_folder + '%s_intersect.txt' % file_name)
        os.remove(tmp_folder + '%s_clustered.txt' % file_name)
        os.remove(tmp_folder + "%s.txt" % file_name)

#%%
with Manager() as manager:
    # Define empty lists
    testMethy_list = manager.list()
    testDepth_list = manager.list()
    testBeta_list  = manager.list()

    pool = Pool(cpuCount)  # Parallelisation function

    pool.map(import_test_files, test_files)    # Import the files in parallel
    pool.terminate()
    
    print("Merging all %s in one file..." % testMethyName)  # Merge the testMethy_list in a pandas dataframe, merging the same indices
    testMethy = pd.concat(testMethy_list, axis = 1)
    testMethy = pd.merge(clusters, testMethy, how = "left", left_index=True, right_index=True)    # Merge the pandas df with the clusters, leaving NA values for clusters that were not covered.
    testMethy = testMethy.transpose().fillna('NA')   # Transpose and Fill NaN with NA string
    testMethy.to_csv("./classifySamples/output/%s" % testMethyName, header=None,sep='\t', mode = 'w')

    print("Merging all %s in one file..." % testDepthName)
    testDepth = pd.concat(testDepth_list, axis = 1)
    testDepth = pd.merge(clusters, testDepth, how = "left", left_index=True, right_index=True)
    testDepth = testDepth.transpose().fillna('NA')
    testDepth.to_csv("./classifySamples/output/%s" % testDepthName, header=None,sep='\t', mode = 'w')

    testBeta = pd.concat(testBeta_list, axis = 1)
    testBeta = pd.merge(clusters, testBeta, how = "left", left_index=True, right_index=True)
    testBeta = testBeta.transpose().fillna('NA')
    print("Writing to disk...")
    testBeta.to_csv("./classifySamples/output/%s" % testBetaName, header=None,sep='\t', mode = 'w')

    testMethy_rmNA = testMethy.apply(pd.to_numeric, errors='coerce').dropna(axis=1)
    testDepth_rmNA = testDepth.apply(pd.to_numeric, errors='coerce').dropna(axis=1)
    print("The number of rows in the %s is: %i" %  (regions,clusters.shape[0])) 
    print("The number of columns in the %s file is: %i" % (testMethyName,testMethy.shape[1]))
    print("The number of columns in the %s file is: %i" % (testDepthName,testDepth.shape[1]))

# %%
