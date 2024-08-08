# Plots

TRACS includes several plots that can be used to further investigate the relationship between samples. 
These can be generated using the `tracs plot` command, specifying the plot type using the `--type` parameter.


### Heatmap

A heatmap provides an efficient way to compare pairwise transmission distances between multiple samples simultaneously. The following script accepts a CSV file produced by the [tracs distance](distance.md) command and generates a heatmap of the distances specified by the `--column-name` parameter.

For example, to plot the filtered SNP distances between pairs of samples you can run

```
tracs plot --type heatmap -i example_dist.csv -p output_prefix  --column-name "filtered SNP distance"
```

Distances can also be filtered prior to plotting using the `--threshold` parameter. This command will result in the following plot

![heatmap showing the pairwise SNP distance relationships of 3 samples](/_figures/heatmap.png)

### Minority variant line plot

This plot aims to determine if a suspected transmission event is driven by a strain that exists at a minority frequency. It requires pileup files generated for at least two samples using the [tracs align](alignment.md) command.

The plot displays the frequency of variants found in both samples, focusing on those that are at a minority frequency in at least one sample. Variants below the frequency specified by the `--min-freq` parameter are filtered out.

Large bands of variants changing in frequency can indicate a strain that has shifted in prevalence between the two samples. Multiple bands may suggest the presence of more than one strain in at least one of the samples.

```
tracs plot --type line -i sampleA_pileup.txt.gz sampleB_pileup.txt.gz -p output_prefix --min-freq 0.05 
```

![Line plot showing the change in frequency of minority variants between two samples](/_figures/line.png)


### Interactive scatter plot of variants

To investigate variants and compare their frequency and location in two samples, the `scatter` plot type can be used. This generates an interactive HTML plot that allows users to zoom into areas of interest. It requires pileup files generated for at least two samples using the [tracs align](alignment.md) command.
Variants below the frequency specified by the `--min-freq` parameter are filtered out.

```
tracs plot --type scatter -i sampleA_pileup.txt.gz sampleB_pileup.txt.gz -p output_prefix --min-freq 0.01
```

### Options

```
tracs plot [-h] -i INPUT_FILES [INPUT_FILES ...] -p OUTPUT_FILE --type {scatter,line,heatmap} [--min-freq MIN_FREQ] [--either-strand] [--contigs CONTIGS [CONTIGS ...]] [--column-name COLUMN_NAME] [--threshold THRESHOLD] [--alpha ALPHA] [--height HEIGHT]
                      [--width WIDTH]

Generates plots from a pileup file.

options:
  -h, --help            show this help message and exit

Input/output:
  -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        path to query signature
  -p OUTPUT_FILE, --prefix OUTPUT_FILE
                        prefix of output file
  --type {scatter,line,heatmap}
                        Type of plot (scatter, line, heatmap)

Pileup options:
  --min-freq MIN_FREQ   minimum frequency to include a variant (default=0.0)
  --either-strand       turns off the requirement that a variant is supported by both strands
  --contigs CONTIGS [CONTIGS ...]
                        contigs for plotting (default=All)

Transmission distance options:
  --column-name COLUMN_NAME
                        Column name in distance matrix to use (default='SNP distance')
  --threshold THRESHOLD
                        threshold to filter transmission distances (default=None)

Plot options:
  --alpha ALPHA         alpha value for plotting (default=0.1)
  --height HEIGHT       height value for plotting (default=7)
  --width WIDTH         width value for plotting (default=10)
```