# Sidescan

## Installing

```
$ git clone git@github.com:omixlab/sidescan.git
$ cd sidescan
$ conda env create
$ conda activate sidescan
```

## Saving
```
git add .
git commit -m 'descrição'
git push -u origin 'nome da branch'
```

## Download datasets

Downloads data from the SIDERS (side effects database) to the `data/` directory.

```
$ sidescan-download data/
```

## Preprocess datasets

```
$ sidescan-preprocess data/
```

## Train model

```
$ sidescan-train data/
```


## Predict side effects for molecule

```
$ sidescan-search --input molecule.sdf --model data/models.pickle --output molecule.json
```

## Run server

```
$ sidescan-server
```