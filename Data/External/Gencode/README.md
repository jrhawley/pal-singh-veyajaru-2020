# File descriptions

This GENCODE annotation file (v19) was downloaded on 2018-12-10 to annotate the hg19-aligned differentially accessible sites, directly from GENCODE's website.
[Download link](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz).

This was filtered for genes only, via

```shell
awk -v OFS="\t" '{if ($3 == "gene") print $1, $4 - 1, $5, ".", $6, $7, $9}' gencode.v19.annotation.gff3 > gencode.v19.annotation.genes.bed
```