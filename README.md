# Summary
This repository contains the scripts to perform deconvolution on cf-RRBS data. 
1. Create a reference dataset compatible with the algorithm = `makeTrain.py -h`, supports both creation of the reference dataset for [meth_atlas](https://github.com/nloyfer/meth_atlas) (`makeTrain.py --t meth_atlas`) and [CelFIE](https://github.com/rmvpaeme/celfie/blob/master/celfiePipeline.snakefile) (`makeTrain.py --t celfie`).
1. Create a matrix containing the samples on which deconvolution will be performed.
    - meth_atlas: `makeTest.py -h`
    - CelFIE: see snakemake pipeline at https://github.com/rmvpaeme/celfie/blob/master/celfiePipeline.snakefile
1. Perform the deconvolution and obtain results for further downstream processing.
    - meth_atlas: `runMeth_atlas.py -h`
    - CelFIE: see snakemake pipeline at https://github.com/rmvpaeme/celfie/blob/master/celfiePipeline.snakefile

# Installation
## Dependencies

```
python==3.9.2
pandas==1.2.4
numpy==1.20.2
scipy==1.6.3
multiprocess==0.70.11.1
bedtools==2.30.0
```
## Installation with Conda
```
git clone https://github.com/rmvpaeme/cfRRBS_classifier
cd cfRRBS_classifier
conda create --name cfRRBS_classifier python=3.9.2
conda activate cfRRBS_classifier
conda install --file requirements.txt
conda install -c bioconda bedtools 
mkdir -p ./classifySamples/resources && cp RRBS_450k_intersectClusters.tsv ./classifySamples/resources 
```

Furthermore, you need the Infinium HM450K and methylationEPIC annotation files if you are going to use these datatypes in the reference dataset.

```
# these are in hg19
mkdir -p ./classifySamples/resources

wget -qO - ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv | gzip -c > HumanMethylation450_15017482_v1-2.csv.gz && mv HumanMethylation450_15017482_v1-2.csv.gz ./classifySamples/resources

wget -q -O tmp.zip https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip && unzip tmp.zip && rm tmp.zip && gzip MethylationEPIC_v-1-0_B4.csv && mv MethylationEPIC_v-1-0_B4.csv.gz ./classifySamples/resources

# these are in GRCh38
# if you use these annotation files, you have to add --hg38 to makeTrain.py and specify their locations with --infiannot and --epicannot

mkdir -p ./classifySamples/resources

wget http://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.tsv.gz && mv HM450.hg38.manifest.tsv.gz ./classifySamples/resources

wget http://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.tsv.gz && mv EPIC.hg38.manifest.tsv.gz ./classifySamples/resources
```
## Quick start
```
# make the reference dataset (only needs to happen once or after including other samples)
$ python makeTrain.py -n ./classifySamples/train/examples/NGS/NBL/,./classifySamples/train/examples/NGS/cfDNA \
                    -a NBL,cfDNA \
                    -r ./classifySamples/resources/RRBS_450k_intersectClusters.tsv \
                    -o test.txt -t methatlas

# make the test matrix
$ python makeTest.py -f ./classifySamples/testfiles/examples/ \
                    -r ./classifySamples/resources/RRBS_450k_intersectClusters.tsv \
                    -c 30 -p outtest

# run the deconvolution with meth-atlas
python runMeth_atlas.py -a ./classifySamples/resources/20190323_test_beta_plasma_MANUSCRIPT.gz \
                        -b ./classifySamples/resources/train_plasma_WGBS_as_normal_MANUSCRIPT.gz \
                        -n normal,wbc \
                        -p examplerun
```

# Changelog
- 2020-05-31: added array annotation GRCh38, made the Illumina annotation files 0-based.
