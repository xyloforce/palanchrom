# Tips

- all scripts are in the `/Scripts` subfolder

# Convert dat files to bed containg human regions only

Use `datToBed.py` for *each species*. Inputs : folder where dat files are, output name.

Outputs one bed containing human intervals.

# Intersect all files to get human common regions

Use `intersectAll.sh` *one time*. Inputs : folder with outputs from previous scripts are, output folder.

Outputs multiple sorted bed and one intersected file.

# Update species intervals according to human ones

Use `updateInterval.py` for *each species*. Inputs : folder where dat files are, intersected file from previous step, output name.

Outputs one bed by species with all interval common accross species.

# Filter barriers to keep ones separated by more than 1kb

Use `filterInterNIEBS.R` *one time*. Inputs : interbarrier file, output name.

Outputs bed with one line by valid barrier.

# Intersect chimp intervals with barriers

Just use bedtools *one time* :

```
bedtools intersect -sorted -a {previous file} -b <(sort -k1,1 -k2,2n {bed for chimp generated before}) > {output name}
```

**WARNING : will report incomplete barriers if they arent included in integrality in the alignement**

# Get fasta seq for common intervals

Get files on UCSC :
```
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'/ -P ../Data
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.fa.gz' -P ../Data
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/gorGor4/bigZips/gorGor4.fa.gz' -P ../Data
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/ponAbe2/bigZips/ponAbe2.fa.gz' -P ../Data
gunzip ../Data/*.fa.gz
```

Add file for hg38 from the disk bc it isnt available on ucsc (why ??)

Then intersect for each species:
```
bedtools getfasta -fi {input fasta} -bed {common bed generated before} > {output name}
```

# Check sequence length


