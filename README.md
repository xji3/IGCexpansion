# IGCexpansion
IGC expansion development folder

##### Dependent software: 

[jsonctmctree package](http://jsonctmctree.readthedocs.org/en/latest/) (powerful likelihood  calculation 
engine by Alex Griffing)

[Biopython](http://biopython.org/wiki/Biopython) (you could install it by `pip install --user Biopython`)

[networkx](https://networkx.github.io/) version 1.1 `pip install --user networkx==1.11 --upgrade` 

#### Coding Language
Python 2.7

#### Installatoin

0. To install python packages, you need to use [pip](https://pip.pypa.io/en/stable/installing/) (package management). 

1. Install jsonctmctree package:
	
	`
	pip install --user git+https://github.com/argriffing/jsonctmctree.git
	`

2. Install IGCexpansion package:
	
	`
	git clone https://github.com/xji3/IGCexpansion.git
	`
	
	`
	cd IGCexpansion
	`
	
	`
	pip install --user .
	`  _(preferred)_
	
	or
	
	`
	python setup.py install
	` _(hard to uninstall)_  


3. Similarly install any other python packages
	
	`
	pip install --user networkx==1.11 --upgrade
	`
	
	`
	pip install --user Biopython
	`


4. edit `tutorials/IS_IGC/Run_IS_IGC.py` to perform analyses.


To uninstall:
	
	pip uninstall IGCexpansion


##### Tutorials

For independent site IGC (IS-IGC) analyses, please refer to this [tutorial](tutorials/IndependentSiteTutorial.md).  
There are two versions of the software that implement the same IS-IGC approach.  The difference between them is the flexibility of considering different duplication/loss histories where the first version assumes one single duplication event without loss which the second version does not.