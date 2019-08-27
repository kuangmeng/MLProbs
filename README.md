# MLProbs

## Dependency Libraries
PnpProbs: https://github.com/ytye/PnpProbs

QuickProbs: https://github.com/refresh-bio/QuickProbs

Scikit-Learn: https://scikit-learn.org/stable/


## Installation
Please compile PnpProbs to `altered_pnpprobs` and move it to `./base_MSA/pnpprobs/` path.

Please compile QuickProbs to `quickprobs` and move it to `./realign/` path.

Use "miniconda" or "pip" to install `sk-learn`, `subprocess`, `joblib`, `xlrd`, `xlwt` and so on.


## Run
Run on all the benchmarks:

```
cd ./path/to/MLProbs/
python script.py
```

Run on a single file:

```
cd ./path/to/MLProbs/
python MLProbs.py input.file.fasta -b [output.fasta.msa]
```

## Citation
