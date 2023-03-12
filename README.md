# Sidescan

## Installing

```
$ git clone git@github.com:omixlab/sidescan.git
$ cd sidescan
$ conda env create
$ conda activate sidescan
```

## Git
```
git add .
git commit -m 'descrição'
git push -u origin 'nome da branch'
```

## ???

```
$sudo apt-get install libxrender1
```

## Download datasets

Downloads, preprocess and train data from the SIDERS (side effects database) to the `data/` directory.

```
$ sidescan-download data/
$ sidescan-preprocess data/
$ sidescan-train data/
```


## Predict side effects for molecule

```
$ sidescan-search --input molecule.sdf --model data/models.pickle --output molecule.json
```

## Redis (celery?)
```
$ sudo apt-get install redis-server
$ sudo service redis-server start
```

## Run server

```
$ sidescan-server run
```