# Trying to create a GUI for my IGCexpansion package
# Xiang ji
# xji3@ncsu.edu
from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'qt4' #'wx'
from traits.api import HasTraits, Str, Int
from traitsui.api import View, Item
import traitsui

# This should be the Model part in MVC structure
class SimpleTest(HasTraits):
    Alignment = Str

view_test = View(Item(name = 'Alignment'))

if __name__ == '__main__':
    test = SimpleTest()
    test.configure_traits()#view = view_test)
