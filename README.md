# Summary
This repository contains the scripts to perform deconvolution on cf-RRBS data with NNLS (wrapper provided by [@nloyfer](https://github.com/nloyfer)) or CeLFIE (algorithm by [@christacaggiano](https://github.com/christacaggiano)). 
1. Create a reference dataset compatible with the algorithm = `makeTrain.py -h`, supports both creation of the reference dataset for [meth_atlas](https://github.com/nloyfer/meth_atlas) (`makeTrain.py --t meth_atlas`) and [CelFIE](https://github.com/rmvpaeme/celfie/blob/master/celfiePipeline.snakefile) (`makeTrain.py --t celfie`).
1. Create a matrix containing the samples on which deconvolution will be performed.
    - meth_atlas: `makeTest.py -h`
    - CelFIE: see snakemake pipeline at https://github.com/rmvpaeme/celfie/blob/master/celfiePipeline.snakefile
1. Perform the deconvolution and obtain results for further downstream processing.
    - [meth_atlas by nloyfer](https://github.com/nloyfer/meth_atlas): `runMeth_atlas.py -h`
    - CelFIE: see snakemake pipeline at https://github.com/rmvpaeme/celfie/blob/master/celfiePipeline.snakefile
1. (optional) By specifying the number of iterations in `runMeth_atlas.py` a confidence around the final tumor call can be obtained. 
    - If the number of iterations = 1 (default), replicates in the reference dataset (e.g. 20x normal tissue, 30x tumor tissue) are collapsed into 1 column for each tissue type and the mean/median is used for further deconvolution.
    - If the number of iterations > 1, a beta distribution is fitted for every CpG/region. Replicates in the reference dataset are collapsed by random sampling from this beta distribution, so the reference dataset is slightly different in every itation. 
    
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

For R, tested with
```
Rmisc_1.5
tidyverse_1.3.1
reshape2_1.4.4
data.table_1.14.0
optparse_1.6.6 
```
## Installation with Conda

```bash
git clone https://github.com/rmvpaeme/cfRRBS_classifier
cd cfRRBS_classifier
conda create --name cfRRBS_classifier python=3.9.2
conda activate cfRRBS_classifier
conda install --file requirements.txt
conda install -c bioconda bedtools 
mkdir -p ./classifySamples/resources && cp RRBS_450k_intersectClusters.tsv ./classifySamples/resources 
```

Furthermore, you need the Infinium HM450K and methylationEPIC annotation files if you are going to use these datatypes in the reference dataset.

This does not install the R dependencies.

```bash
# these are in hg19
mkdir -p ./classifySamples/resources

wget -qO - ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv | gzip -c > HumanMethylation450_15017482_v1-2.csv.gz && mv HumanMethylation450_15017482_v1-2.csv.gz ./classifySamples/resources

wget -q -O tmp.zip https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip && unzip tmp.zip && rm tmp.zip && gzip MethylationEPIC_v-1-0_B4.csv && mv MethylationEPIC_v-1-0_B4.csv.gz ./classifySamples/resources

# annotation files in GRCh38 are available at ./classifySamples/resources
# if you use these annotation files, you have to add --hg38 to makeTrain.py and specify their locations with --infiannot and --epicannot
# you can also download them directly from http://zhouserver.research.chop.edu/InfiniumAnnotation/

wget http://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/HM450/HM450.hg38.manifest.tsv.gz

wget http://zhouserver.research.chop.edu/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.tsv.gz 
```

Generating the reference dataset (makeTrain.py) takes Bismark coverage files and Infinium data as input. Bismark coverage files can be used without any modifications. If you use Infinium data, you first need to preprocess your data so that you have 1 column with cg-IDs and a second column with beta values (according to this [example](https://github.com/rmvpaeme/cfRRBS_classifier/blob/main/classifySamples/train/examples/infinium/EWS/GSM2357802.txt)). You can use a simple for loop in bash for this

```bash
for file in *txt; do
# remove the tail -n+2 step if the files are already headerless
tail -n+2 ${file} | cut -f x,y  > ${file}_cut.txt
done
```

Where `x,y` are the columns to extract.


## Quick start
```bash
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

# if deconvolution is run with more than 1 iteration (e.g. --iter 50)
Rscript plot_error.R -n normal,wbc 
```

## Features
### makeTrain.py
```
python makeTrain.py -h
usage:
    MakeTrain.py

    Creates a reference matrix to use for NNLS or CelFIE with regions as defined by the user. Takes bismark coverage files, Infinium HM450K and methylationEPIC as input.
    Saves the reference matrix in ./classifySamples/output/ folder.

    Example:
    MakeTrain.py -n /folder/NBL,/folder/cfDNA -a NBL,cfDNA -r regions.tsv -o reference.txt

       [-h] [-n NGSFOLDER] [-a NGSLABELS] [-i INFINIUMFOLDER] [-b INFINIUMLABELS] [-e EPICFOLDER] [-d EPICLABELS] -r REGIONS [-c CUTOFF] [-y INFIANNOT] [-z EPICANNOT] [--annotbuild {hg19,hg38}] -t
       {celfie,celfie_individ_cpg,methatlas} -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -n NGSFOLDER, --ngsfolder NGSFOLDER
                        comma-separated list of location of the folder that contains bismark coverage files in cov.gz format (e.g. /path/to/folder/tumor1,/path/to/folder/tumor2). All cov.gz files in this folder will
                        be added to the reference dataset. (default: None)
  -a NGSLABELS, --ngslabels NGSLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2) (default: None)
  -i INFINIUMFOLDER, --infiniumfolder INFINIUMFOLDER
                        comma-separated list of location of the folder that contains HM450K files in .txt format (e.g. /path/to/folder/tumor1,/path/to/folder/tumor2). All .txt files in this folder will be added to
                        the reference dataset. The files should be headerless, tab-separated and contain the cg-identifiers in the first column and the beta-values in the second column. (default: None)
  -b INFINIUMLABELS, --infiniumlabels INFINIUMLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2) (default: None)
  -e EPICFOLDER, --epicfolder EPICFOLDER
                        comma-separated list of location of the folder that contains MethylationEPIC files in .txt format (e.g. /path/to/folder/tumor1,/path/to/folder/tumor2). All .txt files in this folder will be
                        added to the reference dataset. The files should be headerless, tab-separated and contain the cg-identifiers in the first column and the beta-values in the second column. (default: None)
  -d EPICLABELS, --epiclabels EPICLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2) (default: None)
  -r REGIONS, --regions REGIONS
                        tab-separated file contain the regions of interest, containing 4 columns with chrom start stop clusterID (default: None)
  -c CUTOFF, --cutoff CUTOFF
                        all clusters with reads below this threshold will be marked as NA. (default: 30)
  -y INFIANNOT, --infiannot INFIANNOT
                        annotation file of HM450K in csv.gz format, see README and Illuina website. (default:
                        /Users/rmvpaeme/Repos/cfRRBS_classifier_v0.5/classifySamples/resources/HumanMethylation450_15017482_v1-2.csv.gz)
  -z EPICANNOT, --epicannot EPICANNOT
                        annotation file of MethylationEPIC in csv.gz format, see README and Illumina website. (default:
                        /Users/rmvpaeme/Repos/cfRRBS_classifier_v0.5/classifySamples/resources/MethylationEPIC_v-1-0_B4.csv.gz)
  --annotbuild {hg19,hg38}
                        Reference genome for the HM450K/EPIC annotation files. (default: hg19)
  -t {celfie,celfie_individ_cpg,methatlas}, --type {celfie,celfie_individ_cpg,methatlas}
                        Make reference for celfie or meth_atlas (=NNLS). Celfie only supports bismark coverage files. Default = 'methatlas'. (default: methatlas)
  -o OUTPUT, --output OUTPUT
                        the reference matrix will be saved as this file in ./classifySamples/output/ (default: None)

```

### makeTest.py

```
python makeTest.py -h
usage:
    MakeTest.py

    Creates a testmatrix to use for NNLS with regions as defined by the user. Takes bismark coverage files as input.
    Outputs three files in the folder ./classifySamples/output/:
        1. prefix + "_depth.tsv.gz" = contains the total depth per region, per sample
        2. prefix + "_methy.tsv.gz" = contains the number of methylated reads per region, per sample
        3. prefix + "_beta.tsv.gz" = contains beta value per region, per sample.

    Example:
    python makeTest.py -f ./folder/ -r RRBS_450k_intersectClusters.tsv -c 30 -p test

       [-h] -f FOLDER -r REGIONS [-c CUTOFF] [-p OUTPREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -f FOLDER, --folder FOLDER
                        folder containg bismark cov.gz files that will be added to the test matrix. Automatically detects all .cov.gz files. (default: None)
  -r REGIONS, --regions REGIONS
                        tab-separated file contain the regions of interest, containing 4 columns with chrom start stop clusterID. Should be identical to the reference dataset from makeTrain.py. (default: None)
  -c CUTOFF, --cutoff CUTOFF
                        all clusters with reads below this threshold will be marked as NA. (default: 30)
  -p OUTPREFIX, --outprefix OUTPREFIX
                        prefix for the depth, methy and beta file. (default: None)
```

### runMeth_atlas.py

```
python runMeth_atlas.py -h
usage: 
    runMeth_atlas.py

    Runs metatlas on a testmatrix and reference matrix.
    Files will be written to ./classifySamples/output/classification

    With the default parameters (iterations = 1), the replicates are collapsed with mean/median (depending on the -c parameter).
    If iterations are > 1, a beta distribution is calculated from every CpG/region, and in every iteration the methylation (beta) values are random samples from this beta distribution, which gives an uncertainty around the tumor classification call (plot_error.R). 

    Example:
    runMeth_atlas.py -a /folder/test -b /folder/ref -p outprefix -n cfDNA,WBC
    
       [-h] -a TEST -b REFERENCE [-n NORMAL] [-x EXCLUDE] [-c {mean,median}] -p OUTPREFIX [-f] [-m MRC] [-i ITER]

optional arguments:
  -h, --help            show this help message and exit
  -a TEST, --test TEST  prefix + beta.txt.gz from makeTest.py (default: None)
  -b REFERENCE, --reference REFERENCE
                        output from makeTrain.py (default: None)
  -n NORMAL, --normal NORMAL
                        comma-separated list of labels of the normal tissues in the reference dataset to exclude from tumor assignment (after the deconvolution step) e.g. cfdna,wbc (default: None)
  -x EXCLUDE, --exclude EXCLUDE
                        comma-separated list with tumor entities to exclude from the reference dataset (before the deconvolution step) (default: None)
  -c {mean,median}, --collapse {mean,median}
                        collapse replicate columns with mean or median (default: median)
  -p OUTPREFIX, --outprefix OUTPREFIX
                        prefix for output files (default: None)
  -f, --fs              enable feature selection by selecting the top n hyper and hypomethylated regions per entity (default: False)
  -m MRC, --mrc MRC     top n hyper and hypomethylated regions will be selected with feature selection (default: 100)
  -i ITER, --iter ITER  number of iterations (every iteration samples random from the beta distribution for every CpG) (default: 1)
```

### plot_error.R

If `runMeth_atlas.py` is run with more than 1 iteration (e.g. `--iter 50`), a new reference is made in every iteration by sampling from the beta distribution for every CpG/region.

The results can be visualised with `plot_error.R`. Example figures are available in `./classifySamples/output/plots/examples/`

- `countplot.png` shows the majority vote in every iteration (so the tumor class that is the highest in every iteration).
- `errorplot.png` shows the mean +/- SD.

```
Usage: plot_error.R [options]


Options:
        -d CHARACTER, --directory=CHARACTER
                directory for *deconv_output.csv files

        -o CHARACTER, --outdir=CHARACTER
                output directory for figures [default= ./classifySamples/output/plots/]

        -n CHARACTER, --normal=CHARACTER
                comma-separated list of labels of the normal tissues in the reference dataset to exclude from tumor assignment e.g. cfdna,wbc

        -h, --help
                Show this help message and exit
```
# Changelog
- 2021-06-16: added iterations
- 2021-06-03: added "filename" argument for --labels.
- 2021-05-31: added array annotation GRCh38, made the Illumina annotation files 0-based.
