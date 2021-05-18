#%%
## Import libraries
import glob
import os
import tempfile
import numpy as np
import pandas as pd
import subprocess
import argparse
import multiprocess
from multiprocess import Manager, Pool
import imports_cfRRBS_classifier as cfRRBS
import csv  
cpuCount = (multiprocess.cpu_count() - 2)

class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))

# Folder to store intermediate files
tempdir = (tempfile.TemporaryDirectory())
tmp_folder = os.path.join(tempdir.name)
# testargs
#ngsfolder = ["/Users/rmvpaeme/Repos/2003_CelFiE/NBL_reference_set/NBL/","/Users/rmvpaeme/Repos/2003_CelFiE/NBL_reference_set/cfDNA"]
#ngslabels = ["NBL","cfDNA"]

#infiniumfolder= ["/Users/rmvpaeme/Repos/cfRRBS_classifier_v0.2/classifySamples/train/infinium/EWS/","/Users/rmvpaeme/Repos/cfRRBS_classifier_v0.2/classifySamples/train/infinium/EWS2/"]
#infiniumlabels=["EWS","EWS2"]

#regions = "./classifySamples/resources/RRBS_450k_intersectClusters.tsv"
#infiannot = "./classifySamples/resources/HumanMethylation450_15017482_v1-2.csv.gz"
#epicannot = "./classifySamples/resources/MethylationEPIC_v-1-0_B4.csv.gz"
#output = "./testscript.tsv"
#python MakeTrain.py -n /Users/rmvpaeme/Repos/2003_CelFiE/NBL_reference_set/NBL/","/Users/rmvpaeme/Repos/2003_CelFiE/NBL_reference_set/cfDNA -a NBL,cfDNA -r ./classifySamples/resources/RRBS_450k_intersectClusters.tsv -o test -t methatlas
#%%

parser = argparse.ArgumentParser(
    """
    MakeTrain.py

    Creates a reference matrix to use for NNLS or CelFIE with regions as defined by the user. Takes bismark coverage files, Infinium HM450K and methylationEPIC as input.
    Saves the reference matrix in ./classifySamples/output/ folder.

    Example:
    MakeTrain.py -n /folder/NBL,/folder/cfDNA -a NBL,cfDNA -r regions.tsv -o reference.txt
    """, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-n', '--ngsfolder', help = "comma-separated list of location of the folder that contains bismark coverage files in cov.gz format (e.g. /path/to/folder/tumor1,/path/to/folder/tumor2). All cov.gz files in this folder will be added to the reference dataset.", 
                    default = None, action=SplitArgs)
parser.add_argument('-a', '--ngslabels', help = "comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)", default = None, action=SplitArgs)

parser.add_argument('-i', '--infiniumfolder', help = "comma-separated list of location of the folder that contains HM450K files in .txt format (e.g. /path/to/folder/tumor1,/path/to/folder/tumor2). All .txt files in this folder will be added to the reference dataset. The files should be tab-separated and contain the cg-identifiers in the first column and the beta-values in the second column.", default = None, action=SplitArgs)
parser.add_argument('-b', '--infiniumlabels', help = "comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)", default = None, action=SplitArgs)

parser.add_argument('-e', '--epicfolder', help = "comma-separated list of location of the folder that contains MethylationEPIC files in .txt format (e.g. /path/to/folder/tumor1,/path/to/folder/tumor2). All .txt files in this folder will be added to the reference dataset. The files should be tab-separated and contain the cg-identifiers in the first column and the beta-values in the second column.", default = None, action=SplitArgs)
parser.add_argument('-d', '--epiclabels', help = "comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)", default = None, action=SplitArgs)

parser.add_argument('-r', '--regions', required = True, default = None, help = "tab-separated file contain the regions of interest, containing 4 columns with chrom\tstart\tstop\tclusterID")
parser.add_argument('-c', '--cutoff', default = 30, help = "all clusters with reads below this threshold will be marked as NA.")

parser.add_argument('-y', '--infiannot', help = "annotation file of HM450K in csv.gz format, see Illuina website.", default = "./classifySamples/resources/HumanMethylation450_15017482_v1-2.csv.gz")
parser.add_argument('-z', '--epicannot', help = "annotation file of MethylationEPIC in csv.gz format, see Illumina website.", default = "./classifySamples/resources/MethylationEPIC_v-1-0_B4.csv.gz")

parser.add_argument('-t', '--type', choices=['celfie', "celfie_individ_cpg", 'methatlas'], required = True, help = "Make reference for celfie or meth_atlas (=NNLS). Celfie only supports bismark coverage files. Default = 'methatlas'.", default = "methatlas")

parser.add_argument('-o', '--output', required = True, help = "the reference matrix will be saved as this file in ./classifySamples/output/", default = None)
args = parser.parse_args()

ngsfolder = args.ngsfolder
ngslabels = args.ngslabels

infiniumfolder = args.infiniumfolder
infiniumlabels = args.infiniumlabels

epicfolder = args.epicfolder
epiclabels = args.epiclabels

regions = args.regions
infiannot = args.infiannot
epicannot = args.epicannot

cutoff = int(args.cutoff)

trainFileName = args.output

if not os.path.exists("./classifySamples/output/"):
    os.makedirs("./classifySamples/output/")

trainFileName = "./classifySamples/output/" + trainFileName

type = args.type

print("""Running makeTrain.py for %s
            - Bismark coverage folders: %s
                will be labeled in the order %s
            - Infinium 450K folders: %s
                will be labeled in the order %s
            - HumanmethylationEPIC folders: %s
                will be labeled in the order: %s
            - Regions: %s
            - Cutoff: %i
            - HM450K annotation file: %s
            - MethylationEPIC annotation file: %s
            - output: %s,
            - temp directory: %s
    """ % (type, ngsfolder,ngslabels,
            infiniumfolder, infiniumlabels,
            epicfolder, epiclabels,
            regions, cutoff,
            infiannot, epicannot,
            trainFileName, tmp_folder))
#%%
clusters = cfRRBS.import_clusters(regions, tmp_folder)[0]
clusterFile = cfRRBS.import_clusters(regions, tmp_folder)[1]
if infiniumfolder is not None:  
    array450k = cfRRBS.import_450k(infiannot)
if epicfolder is not None:
    array850k = cfRRBS.import_450k(epicannot)


#%%
def getAvg(x):
    ## This function gets the average beta value in a cluster, or writes NA if more than half are not available.
    if isinstance(x, float):
        return x
    elif len(x) == 0:
        x = "NA"
    else:
        line_values = []
        countNA = 0
        line_values = x.split(',')
        line_values = [i.strip(' ') for i in line_values]
        line_values = list(filter(None, line_values))
        countTot = len(line_values)
        for value in line_values:
            if value == 'NA':
                countNA = countNA +1
        if countTot == 0:
            x = "NA"
        elif countNA/countTot >= 0.5:
            x = "NA"
        else:
            line_values_rmNA = filter(lambda a: a != 'NA', line_values)
            line_values_rmNA_list = list(line_values_rmNA)
            calcmean = np.array(line_values_rmNA_list).astype(float)
            x = np.mean(calcmean)
        return x

def generateTrain_Infinium(label, file_name):
    # The input for this function is an ordered infinium 450k file with the order chr - start - stop - beta value. If it doesn't have this structure, some preprocessing needs to be done.
    outfile = open(tmp_folder + "%s_intersect.txt" % file_name, 'w')
    print("     Running bedtools intersect on %s.txt..." % file_name)
    proc = subprocess.Popen(args=["bedtools", "intersect", "-b", tmp_folder + "%s.txt" % clusterFile, "-a", tmp_folder + "%s.txt" % file_name, "-wb"], stdout=outfile, universal_newlines=True).wait()

    df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols=[6,3,4,5,7], header=None, low_memory=False)
    df = df[[4,5,6,3,7]]
    df[3] = df[3].replace(np.nan, 'NA', regex=True)
    df.to_csv(tmp_folder + "%s_reordered.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

    # Group the total and methylated reads per cluster
    arg = "bedtools groupby -i %s%s_reordered.txt -g 1-3,5 -c 4 -o collapse" % (tmp_folder, file_name)
    arg = arg.split()
    outfile = open(tmp_folder + '%s_clustered.txt' % file_name, 'w')
    print("     Running bedtools groupby on %s.txt..." %file_name)
    proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

    df = pd.read_csv(tmp_folder + '%s_clustered.txt' % file_name, sep="\t", header=None, index_col=3 )
    df[4] = df[4].astype(str)
    df = df.groupby([df.index,0,1,2])[4].apply(','.join) # BEDtools groupby doesnt always make the index unique for some reason, this fixes this.
    df = df.reset_index().set_index(3)
    df = df[~df.index.duplicated(keep='first')] # This shouldn't be necessary anymore, but keep it here as a double check
    df.index.name = None
    df.sort_values(by=[0,1,2], inplace=True)
    df[4] = df[4].apply(getAvg)
    df.columns = [0,1,2,label]
    df = df.drop([0,1,2], axis = 1)

    os.remove(tmp_folder + '%s_intersect.txt' % file_name)
    os.remove(tmp_folder + '%s_clustered.txt' % file_name)
    os.remove(tmp_folder + "%s.txt" % file_name)

    return df

def generateTrain_NGS(inputfile, label, file_name):
    # The structure of this function is very similar to import_test_files()
    df = pd.read_csv(inputfile, sep="\t",usecols=[0,1,2,3,4,5], header=None, compression="gzip", low_memory=False)
    df[3] = df[3]/100
    df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

    outfile = open(tmp_folder + '%s_intersect.txt' % file_name, 'w')
    print("     Running bedtools intersect on %s.txt..." % file_name)
    arg = "bedtools intersect -wb -b %s%s.txt -a %s%s.txt" % (tmp_folder, clusterFile, tmp_folder, file_name)
    arg = arg.split()
    proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

    df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols=[6,7,8,3,4,5,9], header=None, low_memory=False )
    df = df[[6,7,8,3,4,5,9]]
    df.to_csv(tmp_folder + "%s_reordered.txt" % file_name, header=None, index=None, sep='\t', mode = 'w')

    arg = "bedtools groupby -i %s%s_reordered.txt -g 1-3,7 -c 5,6 -o sum" % (tmp_folder, file_name)
    arg = arg.split()
    outfile = open(tmp_folder + '%s_clustered.txt' % file_name, 'w')
    print("     Running bedtools groupby on %s.txt..." %file_name)
    proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

    df = pd.read_csv(tmp_folder + '%s_clustered.txt' % file_name, sep="\t", header=None, index_col = 3, low_memory=False )
    df.index.name = None

    df[6] = df[4]/(df[4] + df[5])  # Get beta value per cluster

    df = df[[0,1,2,6,4,5]]         # Reorder
    df.sort_values(by=[0,1,2], inplace=True)

    df[7] = df[4] + df[5]         # Get total depth

    df[[7,6,4,5]] = df[[7,6,4,5]].mask(df[7] < cutoff)         # Mark all clusters lower than cutoff reads with NA
    df = df.replace(np.nan, 'NA', regex=True)
    df.columns = [0,1,2,label,4,5,7]
    df = df.drop([0,1,2,4,5,7], axis = 1)

    os.remove(tmp_folder + '%s_intersect.txt' % file_name)
    os.remove(tmp_folder + '%s_clustered.txt' % file_name)
    os.remove(tmp_folder + "%s.txt" % file_name)

    return df

def generateTrain_celfie(inputfile, label, file_name):
    #inputfile = "/Users/rmvpaeme/Repos/2003_CelFiE/NBL_reference_set/NBL/DNA044134_S32.cov.gz"
    #inputfile = './NBL/DNA044134_S32.cov.gz'
    #file_name = 'a'
    df = pd.read_csv(inputfile, sep="\t",usecols=[0,1,2,3,4,5], header=None, compression = "gzip", low_memory=False)
    df[3] = df[3]/100
    df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

    outfile = open(tmp_folder + '%s_intersect.txt' % file_name, 'w')
    print("     Running bedtools intersect on %s.txt..." % file_name)
    arg = "bedtools intersect -wb -b %s%s.txt -a %s%s.txt" % (tmp_folder, clusterFile, tmp_folder, file_name)
    arg = arg.split()
    proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

    if type == "celfie_individ_cpg":
        df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols = [0,1,2,3,4,5,9], header=None )
        df = df[[0,1,2,3,4,5,9]]
        df.to_csv(tmp_folder + "%s_reordered.txt" % file_name, header=None, index=None, sep='\t', mode = 'w')

        df = pd.read_csv(tmp_folder + "%s_reordered.txt" % file_name, sep="\t", header=None)

        df.sort_values(by=[0,1,2], inplace=True)

        df[7] = df[4] + df[5]         # Get total depth

        df["index"] = df[0].astype(str)+'_'+df[1].astype(str)+'_'+df[2].astype(str)
        df.index = df["index"]
        df.index.name = None
        df = df.drop([0,1,2,3,5,6, "index"], axis = 1)

        os.remove(tmp_folder + '%s_intersect.txt' % file_name)
        os.remove(tmp_folder + "%s.txt" % file_name)

        return df
    else:
        df = pd.read_csv(tmp_folder + '%s_intersect.txt' % file_name, sep="\t", usecols=[6,7,8,3,4,5,9], header=None )
        df = df[[6,7,8,3,4,5,9]]
        df.to_csv(tmp_folder + "%s_reordered.txt" % file_name, header=None, index=None, sep='\t', mode = 'w')

        arg = "bedtools groupby -i %s%s_reordered.txt -g 1-3,7 -c 5,6 -o sum" % (tmp_folder, file_name)
        arg = arg.split()
        outfile = open(tmp_folder + '%s_clustered.txt' % file_name, 'w')
        print("     Running bedtools groupby on %s.txt..." %file_name)
        proc = subprocess.Popen(args=arg, stdout=outfile, universal_newlines=True).wait()

        df = pd.read_csv(tmp_folder + '%s_clustered.txt' % file_name, sep="\t", header=None, index_col = 3 )
        df.index.name = None

        df[6] = df[4]/(df[4] + df[5])  # Get beta value per cluster

        df = df[[0,1,2,6,4,5]]         # Reorder
        df.sort_values(by=[0,1,2], inplace=True)

        df[7] = df[4] + df[5]         # Get total depth

        df[[7,6,4,5]] = df[[7,6,4,5]].mask(df[7] < 30)         # Mark all clusters lower than 30 reads with NA

        df = df.drop([0,1,2,5,6], axis = 1)

        os.remove(tmp_folder + '%s_intersect.txt' % file_name)
        os.remove(tmp_folder + "%s_clustered.txt" % file_name)
        os.remove(tmp_folder + "%s.txt" % file_name)

        return df
#%%
if type == "methatlas":
    with Manager() as manager:
        #Define empty list
        trainFile_list = manager.list()

        if ngsfolder is not None:
            ngsindex = 0
            for folder in ngsfolder:
                files = glob.glob(os.path.join(str(folder), "*.gz"))  
                labels = ngslabels[ngsindex]

                def import_NGS_train(x):
                    file = x
                    #file = '/Users/rmvpaeme/Repos/2003_CelFiE/NBL_reference_set/cfDNA/DNA050873_S4_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz'
                    file_name = os.path.splitext(os.path.basename(file))[0]
                    df = generateTrain_NGS(inputfile = file, label = labels, file_name = file_name)
                    trainFile_list.append(df)

                print("Running on %s which should contain bismark coverage files... " % folder) 
                print("Found files! %s " % files) 
                print("Labeling these files as %s" % labels)
                pool = Pool(cpuCount)
                pool.map(import_NGS_train, files)
                pool.terminate()
                ngsindex = ngsindex + 1
        
        if infiniumfolder is not None:
            infiniumindex = 0 
            for folder in infiniumfolder:
                files = glob.glob(os.path.join(str(folder), "*.txt"))  
                labels = infiniumlabels[infiniumindex]

                def import_HM450K_train(x):
                    file = x
                    file_name = os.path.splitext(os.path.basename(file))[0]
                    ## Open the file
                    df = pd.read_csv(file, sep="\t", header=None, index_col=0, names = ["Beta_Value"], low_memory=False)
                    ## Add the chromosomal position to the sample
                    df = pd.merge(array450k, df, how = "inner", left_index=True, right_index=True)
                    ## Add a stop and reorder the columns
                    df["MAPINFO_Stop"] = df["MAPINFO"]
                    df = df[["CHR", "MAPINFO", "MAPINFO_Stop", "Beta_Value"]]
                    df.sort_values(by = ["CHR", "MAPINFO"], inplace=True)
                    df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

                    df = generateTrain_Infinium(label = labels, file_name = file_name)
                    trainFile_list.append(df)

                print("Running on %s which should contain Infinium HM450K files (tab separated with cg and beta value)... " % folder) 
                print("Found files! %s " % files) 
                print("Labeling these files as %s" % labels)

                pool = Pool(cpuCount)
                pool.map(import_HM450K_train, files)
                pool.terminate()
                infiniumindex = infiniumindex + 1

        if epicfolder is not None:
            epicindex = 0
            for folder in epicfolder:
                files = glob.glob(os.path.join(str(folder), "*.txt"))  
                labels = epiclabels[epicindex]

                def import_epic_train(x):
                    file = x
                    file_name = os.path.splitext(os.path.basename(file))[0]
                    ## Open the file
                    df = pd.read_csv(file, sep="\t", header=None, index_col=0, names = ["Beta_Value"], low_memory=False)
                    ## Add the chromosomal position to the sample
                    df = pd.merge(array850k, df, how = "inner", left_index=True, right_index=True)
                    ## Add a stop and reorder the columns
                    df["MAPINFO_Stop"] = df["MAPINFO"]
                    df = df[["CHR", "MAPINFO", "MAPINFO_Stop", "Beta_Value"]]
                    df.sort_values(by = ["CHR", "MAPINFO"], inplace=True)
                    df.to_csv(tmp_folder + "%s.txt" % file_name , header=None, index=None, sep='\t', mode = 'w')

                    df = generateTrain_Infinium(label = labels, file_name = file_name)
                    trainFile_list.append(df)

                print("Running on %s which should contain Infinium HumanMethylationEPIC files (tab separated with cg and beta value)... " % folder) 
                print("Found files! %s " % files) 
                print("Labeling these files as %s" % labels)

                pool = Pool(cpuCount)
                pool.map(import_epic_train, files)
                pool.terminate()
                epicindex = epicindex + 1 

        # Generate full matrix from list
        trainFile = pd.concat(trainFile_list, axis = 1)
        # Make sure that the file contains all the clusters
        trainFile = pd.merge(clusters, trainFile, how = "left", left_index=True, right_index=True)
        trainFile = trainFile.transpose().fillna('NA')
        trainFile.to_csv(trainFileName + ".tsv.gz", header=None, sep='\t', mode = 'w', compression = "gzip")
        trainFile_rmNA = trainFile.select_dtypes(include=['float64'])
        print("The number of rows in the %s is: %i" %  (regions,clusters.shape[0])) 
        print("The number of columns in the %s file is: %i" %  (trainFileName,trainFile.shape[1]))
        print("The number of columns in the %s file after removing all NA values is: %i" %  (trainFileName,trainFile_rmNA.shape[1]))

elif (type == "celfie" or type == "celfie_individ_cpg"):
    with Manager() as manager:
        #Define empty list
        trainFile_list = manager.list()

        if ngsfolder is not None:
            ngsindex = 0
            for folder in ngsfolder:
                tumorGroup_list =  manager.list()
               # tumorGroup_list = 
                files = glob.glob(os.path.join(str(folder), "*.gz"))  
                labels = ngslabels[ngsindex]

                def import_celfie_NGS_train(x):
                    file = x
                    file_name = os.path.splitext(os.path.basename(file))[0]
                    df = generateTrain_celfie(inputfile = file, label = labels, file_name = file_name)
                    tumorGroup_list.append(df)
                
                print("Running on %s which should contain bismark coverage files... " % folder) 
                print("Found files! %s " % files) 
                print("Labeling these files as %s" % labels)

                pool = Pool(cpuCount)
                pool.map(import_celfie_NGS_train, files)
                pool.terminate()
                ngsindex = ngsindex + 1 

                tumorGroup = pd.concat(tumorGroup_list, axis = 1)
                tumorGroup[8] = tumorGroup[4].mean(axis= 1)
                tumorGroup[9] = tumorGroup[7].mean(axis= 1)
                tumorGroup = tumorGroup.drop([4,7], axis = 1)
                trainFile_list.append(tumorGroup)

        # Generate full matrix from list
        trainFile = pd.concat(trainFile_list, axis = 1)
        if type == "celfie_individ_cpg":
            a = trainFile
            a["index"] = a.index
            a[['chr','start', 'stop']] = a["index"].str.split('_',expand=True)
            a = a.drop(["index"], axis = 1 )
            def move_column_inplace(df, col, pos):
                col = df.pop(col)
                df.insert(pos, col.name, col)
            move_column_inplace(a, "chr", 0)
            move_column_inplace(a, "start", 1)
            move_column_inplace(a, "stop", 2)
            a.sort_values(by=["chr","start","stop"], inplace=True)
            a = a.fillna(np.NaN)
            a.to_csv(trainFileName + "_celfie.tsv.gz", header=None, sep='\t', mode = 'w', index = False, na_rep = np.NaN)
        else:
            # Make sure that the file contains all the clusters
            trainFile = pd.merge(clusters, trainFile, how = "left", left_index=True, right_index=True)
            trainFile_rmNA = trainFile.dropna(axis = 0)
            trainFile_rmNA.to_csv(trainFileName + "_celfie.tsv.gz", header=None, sep='\t', mode = 'w', index = False, na_rep = np.NaN)
            ngslabels.append("unknown")
            with open(trainFileName + '_celfie_samplekey.tsv', 'w', newline='') as f_output:
                tsv_output = csv.writer(f_output, delimiter='\t')
                tsv_output.writerow(ngslabels)
print("Cleaning up temporary directory...")
tempdir.cleanup()
print("Done.")# %%
