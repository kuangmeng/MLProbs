# Title

MLProbs: A Data-Centric Pipeline for better Multiple Sequence Alignment

## Dependency Libraries

PnpProbs: https://github.com/ytye/PnpProbs

QuickProbs: https://github.com/refresh-bio/QuickProbs

Scikit-Learn: https://scikit-learn.org/stable/


## Installation
Please compile our provided $\mathcal C_{P, NP}^{Aln}$ in `./base_MSA/` path to `c_p_np_aln` and save it in `./base_MSA/C_P_NP_Aln/` path.

Please compile our provided QuickProbs in `./realign/` to `quickprobs` and save it in `./realign/QuickProbs/bin/` path.

Use `miniconda` or `pip` to install `sk-learn`, `subprocess`, `joblib`, `xlrd`, `xlwt` or `pip install -r requirements.txt` to install them.


## Run
Run on all the benchmarks (The benchmarks are in `./TEST/bench/` folder):

```
cd ./path/to/MLProbs/TEST/
python script.py
```

Run on a single file:

```
cd ./path/to/MLProbs/
python MLProbs.py input.file.fasta [output.fasta.msa]
```

**If it doesn't work, please move the content in `./TEST/` directory to the `./` directory. **

## Citation