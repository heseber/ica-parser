# Running ICA for many samples

## Creating multi-sample VCF

First, make sure that VCF files are compressed with bgzip and indexed with tabix:

```sh
for i in *.vcf.gz; do gunzip $i; bgzip ${i%%.gz}; done
for i in *.vcf.gz; do tabix -p vcf $i; done
```

Then merge the single-sample VCF files to multi-sample VCF files with `bcftools`. This only works for a maximum of about 1000 samples because Linux has a maximum of 1024 open files per shell by default, so if the total number of samples is larger, then multiple multi-sample-VCF files need to be generated, or the limit needs to be increased, e.g., to 4096 with `ulimin -n 4096`. Splitting into parts may be a good idea anyway because the parts can be run in parallel and the results can be merged afterwards. However, don't use too many batches per machine because then the I/O is the limit. For example, running 8 parts in parallel on a 8-core machine is not faster than running 4 parts.

```sh
bcftools merge -Oz -o ../part1.vcf.gz $(/bin/ls *.vcf.gz|head -n 1000)
bcftools merge -Oz -o ../part2.vcf.gz $(/bin/ls *.vcf.gz|tail -n +1001|head -n 1000)
bcftools merge -Oz -o ../part3.vcf.gz $(/bin/ls *.vcf.gz|tail -n +2001|head -n 1000)
# etc.
```

This takes about 10-12 minutes for each batch of about 1000 samples.

## Run ICA on each of the multi-sample VCF files

```sh
#!/bin/bash
NDIR=$HOME/ICA
GENOME=GRCh37
for VCF in $@
do
    BASE=${VCF%%.vcf.gz}
    if [ -f $BASE.json.gz ]; then continue; fi
    echo "*********************************************************"
    echo "*  Annotating $VCF"
    echo "*********************************************************"
    dotnet $NDIR/bin/Release/net6.0/Nirvana.dll \
        -c   $NDIR/Data/Cache/$GENOME/Both \
        --sd $NDIR/Data/SupplementaryAnnotation/$GENOME \
        -r   $NDIR/Data/References/Homo_sapiens.$GENOME.Nirvana.dat \
        -i   $VCF \
        -o   $BASE \
        2>&1 \
        | tee $BASE.log
done
```

This takes about 40-50 minutes for each multi-sample VCF file.

## Split the multi-sample JSON files

```Python
import icaparser as icap
icap.split_multi_sample_json_file("part1.json.gz", "SINGLE_SAMPLE_JSON_DIR")
icap.split_multi_sample_json_file("part2.json.gz", "SINGLE_SAMPLE_JSON_DIR")
```

This takes about 25 minutes for each multi-sample JSON file.

## Summary of runtimes

For 1000 single-sample VCF files, we get approximately these runtimes:

- 00:12 for merging VCF files to a multi-sample VCF file
- 00:50 for annotating the multi-sample VCF file with ICA
- 00:25 for splitting the multi-sample JSON file into single-sample JSON files

This is about 1.5 hours total runtime. If we ran ICA on each VCF file for the 1000 samples separately, the runtime would be about 15 hours, i.e., 10 times as long.

