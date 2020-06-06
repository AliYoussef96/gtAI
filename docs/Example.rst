Example
========

1- Import gtAI functions.
	
.. code-block:: python

	from gtAI import Run_gtAI
	
	from gtAI import gtAI 
	
2- In this example, we will use `Saccharomyces cerevisiae S288C <https://www.ncbi.nlm.nih.gov/genome/browse/#!/eukaryotes/15/Saccharomyces%20cerevisiae%20S288c>`_ coding sequences.

3- Prepare the tRNA gene copy number of the tested genome.

The user has two options;  a) input the tRNA gene copy number as python dictionary or, b) using GtRNAdb() function, the user can get it automatically from the GtRNA database, using the link to the tested genome (In our case Saccharomyces cerevisiae S288C).
Or by tRNADB_CE() function to get the tRNA gene copy number from tRNADB_CE database using also the link to the tested genome.

In this example, the second option (b) will be used.

.. code-block:: python

	url_GtRNAdb = "http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/"
	
	#### From GtRNAdb
	
	GtRNA = gtAI.GtRNAdb(url_GtRNAdb)

**for more infromation about GtRNAdb() as well as tRNADB_CE(); go to API part.**

4- Parameter settings for gtai_analysis() function.

.. code-block:: python

	genetic_code_number = 1
	ref_fasta = ""
	bacteria = False
	size_pop = 60
	generation_number = 100

**for more infromation about gtai_analysis() and the parameters; go to API part.**


5- Run gtAI.

.. code-block:: python

	df_tai , final_dict_wi, rel_values = Run_gtAI.gtai_analysis(main_fasta,GtRNA,genetic_code_number,bacteria=bacteria, size_pop=size_pop,generation_number=generation_number)


Returns:

.. code-block:: python

	df_tai (dataframe): Contains each gene id and its gtAI value 
	final_dict_wi (dict): Contains each codon and its absolute adaptiveness value
	rel_values (dict): Contains each codon and its relative adaptiveness values
	
	
6- To save the gtAI result as a CSV file.

.. code-block:: python

	import pandas as pd
	df_tai.to_csv("test.csv", header=True)


`Output example <https://github.com/AliYoussef96/gtAI/blob/master/Saccharomyces%20cerevisiae%20S288c.csv>`_
