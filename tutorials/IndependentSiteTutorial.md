There are two versions of the softwawre that both implement the independent site IGC expansion.
The first version assumes one single duplication event without loss which the second version does not.

## New version

I'll start with how to use the new version.  For independent site IGC model (IS-IGC), please refer to the [IS-IGC folder](/IS_IGC) for example input files with a python script that imports the package and runs the analysis.

#### Input Files
This section describes the input files and their supposed format for running the analysis.The table below summarizes all input files (as the variable name in the [script](/IS_IGC/Run_IS_IGC.py)) followed with their more detailed description.

| File | Description |
|-------------|:-------|
| alignment file | Multiple sequence alignment file |
| newicktree | Specis tree file stored in the Newick format |
| DupLosList| A file describes the duplication loss events on the branches of the newick tree |
| gene_to_orlg_file | A file describes the orthologous group each gene belongs to |
| seq_index_file | A file describes corresponding sequence position of each column in the alignment |
