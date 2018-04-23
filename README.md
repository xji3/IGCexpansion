# IGCexpansion
IGC expansion development folder

To install:

`
pip install --user git+https://github.com/xji3/IGCexpansion.git
`

To uninstall:

`
pip uninstall IGCexpansion
`

##### Dependence

The code broke with the newest networkx package update that produce following error:

AttributeError: 'DiDegreeView' object has no attribute 'items'

Unfortunately, the likelihood package also depends on networkx such that reverting back to an older version is necessary before modifying it which probably will never happen...

`
pip install --user networkx==1.11 --upgrade
` 
