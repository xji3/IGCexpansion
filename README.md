# IGCexpansion
[![Build Status](https://travis-ci.com/xji3/IGCexpansion.svg?branch=master)](https://travis-ci.com/xji3/IGCexpansion)

IGC expansion development folder

##### Dependent software:

[jsonctmctree package](http://jsonctmctree.readthedocs.org/en/latest/) (powerful likelihood  calculation
engine by Alex Griffing)

[Biopython](http://biopython.org/wiki/Biopython)
[networkx](https://networkx.github.io/)
[numpy](https://numpy.org/)
[scipy](https://www.scipy.org/)

(you could install them by

`
pip install --user Biopython networkx numpy scipy
`)


#### Coding Language

Python 3.6 or higher


#### Preparation

*Mac OS / Linux*

1. To install python packages, you need to use [pip](https://pip.pypa.io/en/stable/installing/) (package management).

2. You might need to install [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

*Windows*

1. I recommand using [Anaconda](https://www.anaconda.com/products/individual#windows) on Windows that ships with pip functionality and many more useful features.

2. You still need to install [git](https://git-scm.com/download/win)


#### Install python packages

1. Install jsonctmctree package by Alex Griffing (slightly modified for version updates):

`
pip install --user git+https://github.com/xji3/jsonctmctree.git
`

2. Install IGCexpansion package:

`
pip install --user git+https://github.com/xji3/IGCexpansion.git
`

3. Similarly install any other python packages (they should have been installed with IGCexpansion already)

`
pip install --user networkx
`

`
pip install --user Biopython
`

To uninstall:

`
pip uninstall IGCexpansion
`

##### Getting a local copy of the package

`
git clone https://github.com/xji3/IGCexpansion.git
`

You can now run the tutorial file or edit it to perform analyses.

`
cd IGCexpansion/tutorials/IS_IGC
`

`
python Run_IS_IGC.py
`


##### Tutorials

For independent site IGC (IS-IGC) analyses, please refer to this [tutorial](tutorials/IndependentSiteTutorial.md).  
There are two versions of the software that implement the same IS-IGC approach.  The difference between them is the flexibility of considering different duplication/loss histories where the first version assumes one single duplication event without loss which the second version does not.