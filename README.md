# Umi-collapser
## Description
Umi-collapser identifies families of molecules in scRNA-seq experiments
originating from the same original molecule and merges them into a single output 
read, calling consensus bases for any bases covered by multiple reads
 
## Consensus base calling
Two methods of consensus bases calling are currently supported: calculation of
 the posterior probabilities for each base and simple majority voting.
 
## Inputs and Output
Umi-collapser accepts a BAM files with reads tagged by cellular and molecular 
barcodes and outputs a new BAM file where reads from the same family are collapsed.
Reads that do not form families of more than one are output with no modification.

## Example Usage
```

```