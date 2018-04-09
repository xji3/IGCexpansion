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

This should be fixed now. If you are still seeing it, try:

`
pip install --user networkx==1.11 --upgrade
` 
