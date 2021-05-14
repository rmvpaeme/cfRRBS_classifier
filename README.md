# Installation
## Dependencies

```
python==3.9.2
pandas==1.2.4
numpy==1.20.2
scipy==1.6.3
multiprocess==0.70.11.1
```
## Installation with Conda
```
git clone https://github.com/rmvpaeme/cfRRBS_classifier
conda create --name cfRRBS_classifier python=3.9.2
conda activate cfRRBS_classifier
conda install --file requirements.txt
```

## Quick start
```
# make the reference dataset (only needs to happen once or after including other samples)
$ python makeTrain.py -n ./classifySamples/train/examples/NGS/NBL/,./classifySamples/train/examples/NGS/cfDNA \
                    -a NBL,cfDNA \
                    -r ./classifySamples/resources/RRBS_450k_intersectClusters.tsv \
                    -o test.txt
                    
# make the test matrix
$ python makeTest.py -f ./classifySamples/testfiles/examples/ \
                    -r ./classifySamples/resources/RRBS_450k_intersectClusters.tsv \
                    -c 30 -p outtest -t nnls

# run the deconvolution with meth-atlas
python runMeth_atlas.py -a ./classifySamples/resources/20190323_test_beta_plasma_MANUSCRIPT.gz \
                        -b ./classifySamples/resources/train_plasma_WGBS_as_normal_MANUSCRIPT.gz \
                        -n normal,wbc \
                        -p examplerun
```

# Manual
## makeTest.py
Builds the sample matrix. 

Available arguments:
```
 python makeTrain.py -h
usage: 
    MakeTrain.py

    Creates a reference matrix to use for NNLS or CelFIE with regions as defined by the user. Takes bismark coverage files, Infinium HM450K and methylationEPIC as input.

    Example:
    MakeTrain.py -n /folder/NBL,/folder/cfDNA -a NBL,cfDNA -r regions.tsv -o reference.txt
    
       [-h] [-n NGSFOLDER] [-a NGSLABELS] [-i INFINIUMFOLDER] [-b INFINIUMLABELS] [-e EPICFOLDER] [-d EPICLABELS] -r REGIONS [-c CUTOFF] -y INFIANNOT -z
       EPICANNOT -t {celfie,celfie_individ_cpg,methatlas} -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -n NGSFOLDER, --ngsfolder NGSFOLDER
                        comma-separated list of location of the folder that contains bismark coverage files in cov.gz format (e.g.
                        /path/to/folder/tumor1,/path/to/folder/tumor2). All cov.gz files in this folder will be added to the reference dataset.
  -a NGSLABELS, --ngslabels NGSLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)
  -i INFINIUMFOLDER, --infiniumfolder INFINIUMFOLDER
                        comma-separated list of location of the folder that contains HM450K files in .txt format (e.g.
                        /path/to/folder/tumor1,/path/to/folder/tumor2). All .txt files in this folder will be added to the reference dataset. The files should be
                        tab-separated and contain the cg-identifiers in the first column and the beta-values in the second column.
  -b INFINIUMLABELS, --infiniumlabels INFINIUMLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)
  -e EPICFOLDER, --epicfolder EPICFOLDER
                        comma-separated list of location of the folder that contains MethylationEPIC files in .txt format (e.g.
                        /path/to/folder/tumor1,/path/to/folder/tumor2). All .txt files in this folder will be added to the reference dataset. The files should be
                        tab-separated and contain the cg-identifiers in the first column and the beta-values in the second column.
  -d EPICLABELS, --epiclabels EPICLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)
  -r REGIONS, --regions REGIONS
                        tab-separated file contain the regions of interest, containing 4 columns with chrom start stop clusterID
  -c CUTOFF, --cutoff CUTOFF
                        all clusters with reads below this threshold will be marked as NA.
  -y INFIANNOT, --infiannot INFIANNOT
                        annotation file of HM450K in csv.gz format, see Illuina website.
  -z EPICANNOT, --epicannot EPICANNOT
                        annotation file of MethylationEPIC in csv.gz format, see Illumina website.
  -t {celfie,celfie_individ_cpg,methatlas}, --type {celfie,celfie_individ_cpg,methatlas}
                        make reference for celfie or meth_atlas (=NNLS). Celfie only supports bismark coverage files.
  -o OUTPUT, --output OUTPUT
                        reference matrix will be saved as this file
```

## makeTrain.py
Builds the reference matrix. 

```
python makeTrain.py -h
usage: 
    MakeTrain.py

    Creates a reference matrix to use for NNLS with regions as defined by the user. Takes NGS, Infinium HM450K and methylationEPIC as input.

    Example:
    MakeTrain.py -n /folder/NBL,/folder/cfDNA -a NBL,cfDNA -r regions.tsv -o reference.txt
    
       [-h] [-n NGSFOLDER] [-a NGSLABELS] [-i INFINIUMFOLDER] [-b INFINIUMLABELS] [-e EPICFOLDER] [-d EPICLABELS] -r REGIONS [-c CUTOFF] -y INFIANNOT -z
       EPICANNOT -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -n NGSFOLDER, --ngsfolder NGSFOLDER
                        comma-separated list of location of the folder that contains bismark coverage files in cov.gz format (e.g.
                        /path/to/folder/tumor1,/path/to/folder/tumor2). All cov.gz files in this folder will be added to the reference dataset.
  -a NGSLABELS, --ngslabels NGSLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)
  -i INFINIUMFOLDER, --infiniumfolder INFINIUMFOLDER
                        comma-separated list of location of the folder that contains HM450K files in .txt format (e.g.
                        /path/to/folder/tumor1,/path/to/folder/tumor2). All .txt files in this folder will be added to the reference dataset. The files should be
                        tab-separated and contain the cg-identifiers in the first column and the beta-values in the second column.
  -b INFINIUMLABELS, --infiniumlabels INFINIUMLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)
  -e EPICFOLDER, --epicfolder EPICFOLDER
                        comma-separated list of location of the folder that contains MethylationEPIC files in .txt format (e.g.
                        /path/to/folder/tumor1,/path/to/folder/tumor2). All .txt files in this folder will be added to the reference dataset. The files should be
                        tab-separated and contain the cg-identifiers in the first column and the beta-values in the second column.
  -d EPICLABELS, --epiclabels EPICLABELS
                        comma-separated list of labels corresponding to the folders (e.g. tumor1,tumor2)
  -r REGIONS, --regions REGIONS
                        tab-separated file contain the regions of interest, containing 4 columns with chrom start stop clusterID
  -c CUTOFF, --cutoff CUTOFF
                        all clusters with reads below this threshold will be marked as NA.
  -y INFIANNOT, --infiannot INFIANNOT
                        annotation file of HM450K in csv.gz format, see Illuina website.
  -z EPICANNOT, --epicannot EPICANNOT
                        annotation file of MethylationEPIC in csv.gz format, see Illumina website.
  -o OUTPUT, --output OUTPUT
                        reference matrix will be saved as this file
```

## runMeth_atlas.py
Runs deconvolution on the output of the previous two steps. 

```
python runMeth_atlas.py -h
usage: 
    runMeth_atlas.py

    Runs metatlas on a testmatrix and reference matrix.
    Files will be written to ./classifySamples/output/classification

    Example:
    MakeTest.py -a /folder/test -b /folder/ref -p outprefix -n cfDNA,WBC
    
       [-h] -a TEST -b REFERENCE [-n NORMAL] -p OUTPREFIX

optional arguments:
  -h, --help            show this help message and exit
  -a TEST, --test TEST  prefix + beta.txt.gz from makeTest.py
  -b REFERENCE, --reference REFERENCE
                        output from makeTrain.py
  -n NORMAL, --normal NORMAL
                        comma-separated list of labels of the normal tissues in the reference dataset to exclude from tumor assignment e.g. cfdna,wbc
  -p OUTPREFIX, --outprefix OUTPREFIX
                        prefix for output files
```