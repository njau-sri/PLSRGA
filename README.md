# PLSRGA

An R package for QTL mapping in hybrid population (diallel cross and NCII design) using partial least square regression (PLSR) and genetic algorithm (GA).

## Install dependencies

pls

https://github.com/khliland/pls

https://cran.r-project.org/web/packages/pls/index.html

GALGO

https://github.com/vtrevino/GALGO

## Load PLSRGA

    source("PLSRGA.r")

## Perform analysis

    PLSRGA("Result", "genotype.csv", "phenotype.csv")

Arguments

filepath: a string, the scheduled working path

mrkfile: a string, the dataset of genotype, which include the absolute path and the genotype filename.

traitfile: a string, the dataset of phenotype, which include the absolute path and the genotype filename.

Exp.des: a string, the experiment design, which is ‘NCII’ or ‘diallel’.

singmrk: a logical value, TRUE or FALSE, which is used to set whether single marker analysis filtering is required for genotype. Default to TRUE.

para: a logical value, TRUE or FALSE, which is used to set whether the GA algorithm running in parallel. Default to FALSE.

chromSize: an int, Specify the GA chromosome size (the number of variables/genes to be included in a model). Default to 5.

max_gnr: an int, Maximum number of generations. Default to 1000.

pop_size: an int, Specify the number of chromosomes per niche. Default to 50.

goal_fit: a float, Specify the desired fitness value (fraction of correct classification). Defaults to 0.90.

save_freq: an int, How often the “current” solutions are saved. Defaults to 50.

## Genotype data file format

The data source is a text file with tab delimited or csv file with commas separated according your favoriate. The expected file format is marker genotypes in rows and samples in columns. The first columns must be marker names identifier, accession number, or anything to distinguish uniquely the markers. The first row must contain the sample names (identifier of parents using for diallel or NCII), again unique values. In NCII design, the second row is the class description for each parents group, in other words, the second row contain the identifier that were used for distinguish male and female parents.

## Phenotype data file format

Phenotype were ordered by column, in other words, one trait was included in one column. The most important issue is the rank of each row. In marker dataset, the column ranked as m1, m2, m3,... f1, f2, f3,..., 'mi' and 'fi' denote male and female parents respectively.

