There are two versions of the softwawre that both implement the independent site IGC expansion.
The first version assumes one single duplication event without loss which the second version does not.

## New version

I'll start with how to use the new version.  For independent site IGC model (IS-IGC), please refer to the [IS-IGC folder](https://github.com/xji3/IGCexpansion/tree/master/tutorials/IS_IGC) for example input files with a python script that imports the package and runs the analysis.

#### Input Files
This section describes the input files and their supposed format for running the analysis.The table below summarizes all input files (as the variable name in the [script](https://github.com/xji3/IGCexpansion/tree/master/tutorials/IS_IGC/Run_IS_IGC.py)) followed with their more detailed descriptions.

| File | Description |
|-------------|:-------|
| [alignment file](#alignment) | Multiple sequence alignment file |
| [newicktree](#newick) | Specis tree file stored in the Newick format |
| [DupLosList](#DupLosList)| A file describes the duplication loss events on the branches of the newick tree |
| [gene\_to\_orlg_file](#gene_to_orlg) | A file describes the orthologous group each gene belongs to |
| [seq\_index\_file](#seq_index) | A file describes corresponding sequence position of each column in the alignment |

##### <a name='alignment'>alignment file</a>
The alignment file is in fasta format.
However, the name of the gene follows a rather weird rule:
it is combined by the species name and gene name with a double underscore **"\__"** connecting them. For example, **Human** **EDN** gene would be denoted as **"Human\__EDN"** this way.  The double score is used to allow species / gene names to have single scores. For example, a **Tree\_Shrew** **EDN** gene would be **Tree\_Shrew__EDN**
##### <a name='newick'>newick tree</a>
##### <a name='DupLosList'>DupLosList</a>
##### <a name='gene_to_orlg'>gene_to_orlg_file</a>
##### <a name='seq_index'>seq_index_file</a>

