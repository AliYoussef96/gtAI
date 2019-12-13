# gtAI_pkg

## Python Support

Python >=3.7 is required.

## Dependencies

1. Biopython

2. pandas

3. urllib

4. numpy

5. gaft

6. CAI

7. scipy

## Installation Instructions

**Using pip**

```python
pip install gtAI_pkg
```

## Contribution Guidelines

For bugs and suggestions, the most effective way is by raising an issue on the github issue tracker. 
Github allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.

If you wish to contribute some changes to the code then you should submit a [pull request](https://github.com/AliYoussef96/gtAI_pkg/pulls)
How to create a Pull Request? [documentation on pull requests](https://help.github.com/en/articles/about-pull-requests)

## Usage

```python
gtai_analysis(main_fasta, GtRNA, genetic_code_number, size_pop, generation_number=50, bacteria=False)
```

Where:

```python

main_fasta (str): A main fasta file contains genes want to be analyzed (CDS)
GtRNA (dict): the tRNA genes count
ref_fasta (str): A reference genes with the highest gene expression in a genome (CDS)
genetic_code_number (int): default = 1, The Genetic Codes number described by NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
size_pop (int): A parameter for the genetic algorithm to identify the population size containing the possible solutions to optimize Sij-values
generation_number (int): A parameter for the genetic algorithm to identify the generation number
bacteria (bool): True If the tested organism is prokaryotic or archaea, else equal to False ( default = False )


```
Returns:

```python
df_tai (dataframe): Contains each gene id and its gtAI value 
final_dict_wi (dict): Contains each codon and its absolute adaptiveness value
rel_values (dict): Contains each codon and its relative adaptiveness values
```
## Easy use

We provide easy GUI software to run this package.

Just install python3 and download the package using pip.
```python
pip install gtAI_pkg
```
Download this zip file, uncompress it and run the script named gtAI_GUI.py by any way you would like.

Run it on any operating system and enjoy easily using our package.

## Citation

