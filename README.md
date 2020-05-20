# Global tRNA Adaptation index (gtAI)

**To estimate the tRNA adaptation index (tAI).**

- For more information about the gtAI: 

- For more information about the tAI: [Mario dos Reis et. al.,](https://academic.oup.com/nar/article/32/17/5036/1333956).

## Python Support

Python >=3.7 is required.

## Dependencies

1. Biopython

2. pandas

3. urllib

4. numpy

5. gaft

## Installation Instructions

**Using pip**

```python
pip install gtAI
```

## Contribution Guidelines

Contributions to the software are welcome

For bugs and suggestions, the most effective way is by raising an issue on the github issue tracker. 
Github allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.

If you wish to contribute some changes to the code then you should submit a [pull request](https://github.com/AliYoussef96/gtAI_pkg/pulls)
How to create a Pull Request? [documentation on pull requests](https://help.github.com/en/articles/about-pull-requests)

## Usage

```python
from Run_gtAI import gtai_analysis

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

## Output

## Example


## API Documentation

You can access the API documentation from here: [gtAI Documentation]()


## Citation

