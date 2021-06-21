#imports
import pandas as pd
import os
# Folder to store intermediate files

def import_clusters(regions, tempdir, type = "methatlas"):
    # The file containing the features (= the intersect between HM450K data and RRBS data, see GitHub README)
    if (type == "methatlas"):
        clusters = pd.read_csv(regions, sep="\t",usecols=[0,1,2], skiprows=[0], header=None, index_col=None)
    elif (type == "cancerdetector"):
        clusters = pd.read_csv(regions, sep="\t",usecols=[0,1,2,3], skiprows=[0], header=None, index_col=3)        
    clusters = clusters[clusters[0] != "chrY"]
    clusters = clusters[clusters[0] != "chrX"]
    clusters = clusters[clusters[0] != "Y"]
    clusters = clusters[clusters[0] != "X"]
    clusters[3] = clusters.index
    clusterFile = str(os.path.splitext(os.path.basename(regions))[0]) + "_tmp"
    clusters.to_csv(tempdir + "%s.txt" % clusterFile, header=None, index=None, sep='\t', mode = 'w')
    if (type == "celfie" or type == "celfie_individ_cpg"):
        clusters = clusters.drop([3], axis = 1)
        clusters.columns = ["CHROM", "START", "END"] # Use empty index to later extract all the clusters from, so that every sample has the same number of clusters
    else:
        clusters = clusters.drop([0,1,2,3], axis = 1) # Use empty index to later extract all the clusters from, so that every sample has the same number of clusters
    return clusters,str(clusterFile)

def import_450k(infiannot, annotbuild):
    if annotbuild == "hg19":
        array450k = pd.read_csv(infiannot, dtype={"CHR": str}, header = 7, usecols = (0,10,11,12), index_col="IlmnID", compression = "gzip")
        array450k = array450k.dropna()
        array450k[['MAPINFO', 'Genome_Build']] = array450k[['MAPINFO', 'Genome_Build']].astype(int)
        array450k = array450k[array450k['Genome_Build'] == 37] # Extract locations with genome build GRCh37
        array450k["MAPINFO"] = array450k["MAPINFO"] - 1
        array450k = array450k.drop(['Genome_Build'], axis = 1)
        array450k[['CHR', 'MAPINFO']] = array450k[['CHR', 'MAPINFO']].astype(str)
        array450k.index.name = None
    elif annotbuild == "hg38":
        array450k = pd.read_csv(infiannot, compression = "gzip", sep = "\t", index_col="probeID", usecols = (0,1,4)) # is 0-based
        array450k = array450k.dropna()
        array450k.columns = ["CHR", "MAPINFO"]
        array450k["MAPINFO"] = array450k["MAPINFO"].astype(int) 
        array450k[['CHR', 'MAPINFO']] = array450k[['CHR', 'MAPINFO']].astype(str)
        array450k.index.name = None
    return array450k

def import_850k(epicannot, annotbuild):
    # Load MethylationEPIC reference file
    if annotbuild == "hg19":
        array850k = pd.read_csv(epicannot, dtype={"CHR": str}, header = 7, usecols = (0,10,11,12), index_col="IlmnID", compression = "gzip")
        array850k = array850k.dropna()
        array850k[['MAPINFO', 'Genome_Build']] = array850k[['MAPINFO', 'Genome_Build']].astype(int)
        array850k = array850k[array850k['Genome_Build'] == 37] # Extract locations with genome build GRCh37
        array850k = array850k.drop(['Genome_Build'], axis = 1)
        array850k["MAPINFO"] = array850k["MAPINFO"] - 1 # make 0-based
        array850k[['CHR', 'MAPINFO']] = array850k[['CHR', 'MAPINFO']].astype(str)
        array850k.index.name = None
    elif annotbuild == "hg38":
        array850k = pd.read_csv(epicannot, compression = "gzip", sep = "\t", index_col="probeID", usecols = (0,1,4)) # is 0-based
        array850k = array850k.dropna()
        array850k.columns = ["CHR", "MAPINFO"]
        array850k["MAPINFO"] = array850k["MAPINFO"].astype(int) 
        array850k[['CHR', 'MAPINFO']] = array850k[['CHR', 'MAPINFO']].astype(str)
        array850k.index.name = None
    return array850k