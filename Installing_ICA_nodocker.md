# Software installation - option 1 (no Docker)

## Step 1: Install the dotnet 6.0 SDK

```sh
wget https://packages.microsoft.com/config/ubuntu/18.04/packages-microsoft-prod.deb \
    -O packages-microsoft-prod.deb
sudo dpkg -i packages-microsoft-prod.deb
rm packages-microsoft-prod.deb
sudo apt-get update; \
  sudo apt-get install -y apt-transport-https && \
  sudo apt-get update && \
  sudo apt-get install -y dotnet-sdk-6.0
```

## Step 2: Install the _Illumina Connected Annotations_ binaries

```sh
N_DIR=$HOME/ICA
unzip -d $N_DIR IlluminaConnectedAnnotations-3.22.0-0-gc13dcb61-net6.0.zip
```

## Step 3: Install a wrapper script for the annotator

**run_ica.sh:**

```sh
#!/bin/bash
NDIR=$HOME/ICA
GENOME=GRCh37
VCF=$1
BASE=$(basename $VCF .vcf)
echo "*********************************************************"
echo "*  Annotating $VCF"
echo "*********************************************************"
dotnet $NDIR/Nirvana.dll \
-c   $NDIR/Data/Cache/ \
--sd $NDIR/Data/SupplementaryAnnotation/$GENOME \
-r   $NDIR/Data/References/Homo_sapiens.$GENOME.Nirvana.dat \
-i   $VCF \
-o   $BASE \
2>&1 \
| tee nirvana.$BASE.log
```

## Step 4: Install annotation data sources

Download the annotation data from Illumina using Illumina's Downloader tool as shown below. Illumina's Downloader won't include COSMIC because COSMIC requires a license for commercial use. For non-commercial use, or when having a license, you can download COSMIC's original data files and transform them to the format required by _Illumina Connected Annotations_ with Illumina's SAUtils program, which is included in their software bundle. You can use the script below for this.

This will download about 77 GB, make sure you have enough disk space in your local environment. If you need only one genome version, use "--ga GRCh37" or "--ga GRCh38" instead of "--ga all".

```bash
#!/bin/bash
N_DATA_DIR=$HOME/ICA/Data
mkdir -p $N_DATA_DIR
docker run -it --rm -v $HOME:$HOME illumina-connected-annotations:3.22.0 Downloader \
     --ga all \
     -o $N_DATA_DIR
```

## Step 5 (optional): Download test VCF file

```bash
curl -O https://illumina.github.io/ICADocumentation/files/HiSeq.10000.vcf.gz
```

## Step 6 (optional): Run ICA on the test VCF file

Running the pipeline on the test VCF file takes about 9 seconds.

```bash
run_ica.sh HiSeq.10000.vcf.gz
```
## Step 7 (optional): Add COSMIC as additional annotation

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