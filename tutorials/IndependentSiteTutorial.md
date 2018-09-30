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
it is combined by the species name and gene name with a double underscore **"\__"** connecting them. For example, **Human** **EDN** gene would be denoted as **"Human\__EDN"** this way.  The double score is used to allow species / gene names to have single scores. For example, a **Tree\_Shrew** **EDN** gene would be **Tree\_Shrew__EDN**.  

Please use this [example alignment file](https://github.com/xji3/IGCexpansion/tree/master/tutorials/IS_IGC/EDN_ECP_Cleaned_NewFormat.fasta) as a reference.


##### <a name='newick'>newick tree</a>

The tree file describes the species tree in newick format.

However, unlike usual newick files, all internal nodes have to be named in this tree file for later uses to define the duplication loss history along this species tree (to point events to the right branch). A convention I have been using is to name them by "N" + number, e.g. N0, N1... I usually start the numbering from the root to the tip.

Please use this [example tree file](https://github.com/xji3/IGCexpansion/tree/master/tutorials/IS_IGC/EDN_ECP_tree.newick) as a reference.

##### <a name='DupLosList'>DupLosList</a>
##### <a name='gene_to_orlg'>gene_to_orlg_file</a>
##### <a name='seq_index'>seq_index_file</a>

