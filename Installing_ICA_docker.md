# Software installation - option 2 (Docker)

## Step 1: Install the docker image

```sh
docker load < IlluminaConnectedAnnotations-3.22.0-0-gc13dcb61-net6.0-amd64-docker.tar.gz
```

## Step 2: Install a wrapper script for the annotator

**run_ica_docker.sh:**

```sh
#!/bin/bash
N_DATA_DIR=$HOME/ICA/Data
GENOME=GRCh37
VCF=$(realpath $1)
BASE=${VCF%%.vcf.gz}
echo "*********************************************************"
echo "*  Annotating $(basename $VCF)"
echo "*********************************************************"
docker run -it --rm -v $HOME:$HOME illumina-connected-annotations:3.22.0 Nirvana \
-c   $N_DATA_DIR/Cache/ \
--sd $N_DATA_DIR/SupplementaryAnnotation/$GENOME \
-r   $N_DATA_DIR/References/Homo_sapiens.$GENOME.Nirvana.dat \
-i   $VCF \
-o   $BASE \
2>&1 \
| tee $(dirname $BASE)/nirvana.$(basename $BASE).log
```

## Step 3: Install annotation data sources

Download the annotation data from Illumina using Illumina's Downloader tool as shown below. Illumina's Downloader won't include COSMIC because COSMIC requires a license for commercial use. For non-commercial use, or when having a license, you can download COSMIC's original data files and transform them to the format required by _Illumina Connected Annotations_ with Illumina's SAUtils program, which is included in their software bundle. You can use the script below for this.

This will download about 77 GB, make sure you have enough disk space in your local environment. If you need only one genome version, use "--ga GRCh37" or "--ga GRCh38" instead of "--ga all".

```bash
#!/bin/bash
N_DATA_DIR=$HOME/ICA-3.22.0/Data
mkdir -p $N_DATA_DIR
docker run -it --rm -v $HOME:$HOME illumina-connected-annotations:3.22.0 Downloader \
     --ga all \
     -o $N_DATA_DIR
```

## Step 4 (optional): Download test VCF file

```bash
curl -O https://illumina.github.io/ICADocumentation/files/HiSeq.10000.vcf.gz
```

## Step 5 (optional): Run ICA on the test VCF file

Running the pipeline on the test VCF file takes about 9 seconds.

```bash
run_ica_docker.sh HiSeq.10000.vcf.gz
```

## Step 6 (optional): Add COSMIC as additional annotation

Here is how to create ICA annotation files for COSMIC based on the original COSMIC files. This needs to be done for each genome version separately.

```bash
#!/bin/bash

##### MAKE SURE TO UPDATE THIS AS REQUIRED #####
COSMIC_VERSION=v94
COSMIC_DATE=2021-04-29
GRCH=GRCh37
################################################

cd $HOME/ICA
TMPDIR=$HOME/tmp/create_cosmicdb.$$
OUTDIR=Data/SupplementaryAnnotation/$GRCH
COSMIC_DIR=$HOME/Cosmic/$COSMIC_VERSION/cosmic/$GRCH/

mkdir -p $TMPDIR

ln -s $COSMIC_DIR/CosmicCodingMuts.vcf.gz $TMPDIR
ln -s $COSMIC_DIR/CosmicMutantExportCensus.tsv.gz $TMPDIR

cat >$TMPDIR/CosmicCodingMuts.vcf.gz.version  <<EOT
NAME=cosmic
VERSION=$COSMIC_VERSION
DATE=$COSMIC_DATE
DESCRIPTION=Cosmic Mutant Export Census
EOT

mkdir -p $OUTDIR
dotnet bin/Release/net6.0/SAUtils.dll Cosmic \
        -r Data/References/Homo_sapiens.$GRCH.Nirvana.dat \
        -i $TMPDIR/CosmicCodingMuts.vcf.gz \
        -t $TMPDIR/CosmicMutantExportCensus.tsv.gz \
        -o $OUTDIR

rm -rf $TMPDIR

```

