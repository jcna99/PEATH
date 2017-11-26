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
g++ PEATH.cpp -o PEATH -O2 -std=c++11
```

## Running PEATH

To run PEATH, use the following command:

```
./PEATH <input_file> <output_file> (param)
```

<input_file> is an input matrix for sequence reads and
<output_file> contains phased haplotype.

(param) is an optional parameter for time/accuracy tradeoff which is a positive integer (default: 50).

```
ex) ./PEATH chr1.matrix.SORTED chr1.haplo
```

## Evaluation

To evaluate a haplotype, use the evaluation code.

```
g++ eval.cpp -o eval -O2
./eval <matrix_file> <answer_file> <haplo_file> <output_file>
```

<matrix_file> is a matrix for sequence reads,
<answer_file> is an answer file,
<haplo_file> is a haplotype file,
<output_file> is an output file containing evaluation results.

(mask_file) is a mask haplo file which is optional.

```
ex) ./eval chr1.matrix.SORTED chr1.valid.master chr1.haplo chr1.eval
```

## Data sets

We used three data sets for experiments.
1. Fosmid dataset (Duitama et al. 2012) which has been widely used to assess and compare SIH algorithms.
2. Simulated dataset which was generated based on Fosmid data.
3. HuRef dataset (Levy et al. 2007) which has been the most widely used in SIH related articles.

