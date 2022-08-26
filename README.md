## QATEN: High-accuracy protein model quality assessment using attention graph neural networks
 

### Dependencies

- Dependencies can be installed using the `requirements.txt` file.

### Datasets:

- We use DeepAccNet and GNNRefine datasets to train the QATEN, and use [CASP14](http://predictioncenter.org/download_area/CASP14/server_predictions/) and [RCSB](https://www.rcsb.org/) datasets for local(residue) and global Quality Assessment predictions. You can download them through [links](http://www.csbio.sjtu.edu.cn/bioinf/QATEN/).

### Demo pictures:
- We compare QATEN with confidence predictor in AlphaFold2 and evaluation module in RosettaFold on test datasets. The demo pictures are in our network `http://www.csbio.sjtu.edu.cn/bioinf/QATEN/`.

### To start a quick testing:

1) Install all the requirements by executing `pip install -r requirements.txt.`

2) Install required processing library by executing `sh preprocess.sh` which clones and installs [this](https://github.com/gjoni/mylddt) github repository.

3) move all the pdb files in `data` folder.

4) Quick run:
  ```shell
  python test.py
  ```
  Once successfully run, this creates feature files under the path `./data/`. Global score and local scores are predicted under the path `./result/`.

Please cite the following paper if you use this code in your work.
```bibtex
@article {Sanyal2020.04.06.028266,
	author = {Pei-Dong Zhang, Chun-Qiu Xia and Hong-Bin Shen},
	title = {QATEN: High-accuracy protein model quality assessment using attention graph neural networks},
	year = {2022},
}
```
For any clarification, comments, or suggestions please create an issue or contact [Pei-Dong Zhang](cq-zhang-2016@sjtu.edu.cn).
