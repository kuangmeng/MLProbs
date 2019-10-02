# MLProbs

## Dependency Libraries
PnpProbs: https://github.com/ytye/PnpProbs

QuickProbs: https://github.com/refresh-bio/QuickProbs

Scikit-Learn: https://scikit-learn.org/stable/


## Installation
Please compile our provided PnpProbs in `./base_MSA/` path to `altered_pnpprobs` and save it in `./base_MSA/PnpProbs/` path.

Please compile our provided QuickProbs in `./realign/` to `quickprobs` and save it in `./realign/QuickProbs/bin/` path.

Use `miniconda` or `pip` to install `sk-learn`, `subprocess`, `joblib`, `xlrd`, `xlwt`.


## Run
Run on all the benchmarks (The benchmarks are in `./bench_all/` folder):

```
cd ./path/to/MLProbs/
python script.py
```

Run on a single file:

```
cd ./path/to/MLProbs/
python MLProbs.py input.file.fasta -b [output.fasta.msa]
```
or

```
cd ./path/to/MLProbs/
python MLProbs.py input.file.fasta -s [output.fasta.msa]
```

## Citation
