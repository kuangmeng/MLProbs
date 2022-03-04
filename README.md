# MLProbs: A Data-Centric Pipeline for better Multiple Sequence Alignment

## Evaluation
Run on all the benchmarks (in `./TEST/` folder):

```
cd ./path/to/MLProbs/
python script.py
```

Run on a single file:

```
cd ./path/to/MLProbs/
python MLProbs.py input.file.fasta [output.fasta.msa]
```

## Installation
Please follow the `README` file in `./baseMSA/C_P_NP_Aln/` to compile our provided <img src="./images/classifier.png" width = "30" height = "25" alt="$\mathcal C_{P, NP}^{Aln}$" align=center /> in `./baseMSA/` folder to `c_p_np_aln` and save it to `./baseMSA/C_P_NP_Aln/`.

Please follow the `quickprobs-manual.pdf` to compile our provided QuickProbs in `./realign/QuickProbs/` to `quickprobs` and save it to `./realign/QuickProbs/bin/`.

Use `miniconda` or `pip` to install `sk-learn`, `subprocess`, `joblib`, `xlrd`, `xlwt` or `pip install -r requirements.txt` to install them.


## Main Dependency Libraries

PnpProbs: https://github.com/ytye/PnpProbs

QuickProbs: https://github.com/refresh-bio/QuickProbs

Scikit-Learn: https://scikit-learn.org/stable/ (recommend version: 0.21.3)

Others: `subprocess`, `joblib`, `xlrd`, `xlwt` (or just `pip install -r requirements.txt`)


## Citation
```
@article{9705141,  
    author={Kuang, Mengmeng and Zhang, Yong and Lam, Tak-Wah and Ting, Hing-Fung},  
    journal={IEEE/ACM Transactions on Computational Biology and Bioinformatics},   
    title={MLProbs: A Data-Centric Pipeline for better Multiple Sequence Alignment},   
    year={2022},  
    doi={10.1109/TCBB.2022.3148382}
}
```
