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

## Running on example files

Making the reference matrix:
```
python makeTrain.py -n ./classifySamples/train/examples/NGS/NBL/,./classifySamples/train/examples/NGS/cfDNA \
                    -a NBL,cfDNA \
                    -r ./classifySamples/resources/RRBS_450k_intersectClusters.tsv \
                    -o test.txt
```

Making the sample matrix:
```
python makeTest.py -f ./classifySamples/testfiles/examples/ \
        -r ./classifySamples/resources/RRBS_450k_intersectClusters.tsv \
         -c 30 -p outtest
```
