# Sidescan

## Installing

```
$ git clone git@github.com:omixlab/sidescan.git
$ cd sidescan
$ conda env create
$ conda activate sidescan
```

## Requirements

```
$ sudo apt-get install libxrender1
$ sudo apt-get install redis-server
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

## Setting venv variables on .env file
```
MAIL_USERNAME=<email>
MAIL_PASSWORD=<password>
```

## Run celery
```
$ celery -A sidescan.worker.broker worker
```

## Run server
```
$ sidescan-server
```



## Git
```
git add .
git commit -m 'descrição'
git push -u origin 'nome da branch'
```