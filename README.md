# PEATH: Probabilistic Evolutionary Algorithm with Toggling for Haplotyping

PEATH is a novel SIH algorithm based on the estimation of distribution algorithm (EDA).
It implementes the algorithm proposed in:
```
J.C. Na et al., PEATH: Single Individual Haplotyping by Probabilistic Evolutionary Algorithm with Toggling for Haplotyping,
submitted to Bioinformatics.
```

## Compiling PEATH

After downloading the CPP code, complie PATH by using

```
g++ PEATH.cpp -o PEATH
```

## Running PEATH

To run PEATH, use the following command:

```
./PEATH <input_file> <output_file>
```

<input_file> is an input matrix for sequence reads and
<outputfile> contains phased haplotype.

ex) ./PEATH chr1.matrix.SORTED chr1.phased
(Assuming that working directory contains PEATH(executable file) and chr1.matrix.SORTED(input file))

## Data sets

We used three data sets for experiments.
1. HuRef dataset (Levy et al. 2007) which has been the most widely used in SIH related articles.
2. Fosmid dataset (Duitama et al. 2012) which has been widely used to assess and compare SIH algorithms.
3. Simulated dataset which was generated based on Fosmid data.
