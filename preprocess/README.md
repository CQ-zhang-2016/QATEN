# mylddt

### Requirements

`gcc` version `6.1.0` and above (tested on Linux platform)

### Download & compile
```
git clone https://github.com/gjoni/mylddt
cd mylddt
make
```

### Usage

```
Usage:   ./get_features [-option] [argument]

Options:  -i input.pdb  \       # (required) input PDB file
          -r reference.pdb \    # (optional) reference PDB file
          -j output.json \      # (optional) atomic features
          -p output.pdb \       # (optional) cleaned model
          -d DMAX \             # (999.99) contact distance
          -t TOPN \             # (100) top contacts to save
          -v VERBOSITY \        # (1) verbosity level

```

#### Extract features from a PDB structure using 8A distance cutoff
```
./get_features -i example/tag0001.al.pdb -j example/tag0001.al.json -d 8.0
```

#### Extract features and compare to the native structure
```
./get_features -i example/tag0001.al.pdb -r example/native.pdb -j example/tag0001.al.json
```

#### Extract features for all files in a folder
Unpack example PDBs
```
tar xf init.tar.gz
```

Do feature extraction
```
./preprocessor.sh ./init ./json
```

### Troubleshooting

If ```make``` does not work:
```
cd src
g++ -Wall -Wno-unused-result -pedantic -O3 -mtune=native -std=c++11 *.cpp -o ../get_features
```

