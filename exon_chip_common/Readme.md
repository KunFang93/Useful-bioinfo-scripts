# Finding Exon and Chip overlap information
## Programming language
python
## Required Packages
pandas  
pybedtools
## Usage
```
python find_chip_exon.py meta_info.xlsx exon_expand_length(FLOAT) cutoff(FLOAT,[0,1])

Parameters:
  exon_expand_length: control how long the exon region extent from its ends
  cutoff: minimum fraction of the bedtools intersection (bedtools intersect -f)
```
