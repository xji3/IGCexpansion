# Trying to create a GUI for my IGCexpansion package
# Xiang ji
# xji3@ncsu.edu
from traits.etsconfig.api import *
ETSConfig.toolkit = 'qt4' #'wx'
from traitsui.api import *
from traits.api import *
import os
#from mpl_figure_editor import MPLFigureEditor

class TextDisplay(HasTraits):
    string = String()

    view = View( Item('string', show_label = False, springy = True, style = 'custom'))

# This should be the Model part in MVC structure
class Alignment(HasTraits):
    alignment_file = File
    display = Instance(TextDisplay)
    #Align   = Button()
    _updated = Bool(False)

    view = View( Item('alignment_file', show_label = False))#, Item('Align', show_label = False))

    def _alignment_file_changed(self):
        if os.path.isfile(self.alignment_file):
            self.display.string = 'Found alignment file. \n' + self.display.string
        else:
            self.display.string = 'Alignment file does not exist. \n' + self.display.string

# model enumeration list
class SupportedModels(HasTraits):
    model = Enum('MG94', 'HKY',)
    force = Bool
    clock = Bool

    display = Instance(TextDisplay)

    view = View( Group(Item(name = 'model',
                            label = 'Implemented Models'),
                       Item(name = 'force',
                            style = 'simple',
                            label = 'without tau parameter'),
                       Item(name = 'clock',
                            style = 'simple',
                            label = 'strict molecular clock'),
                       style = 'simple')
                 )

    def _model_changes(self):
        self.display.string = 'Model set to ' + self.model + '\n' + self.display.string

    def _force_changes(self):
        self.display.string = 'Constrain tau to be 0. \n' + self.display.string

    def _clock_changes(self):
        self.display.string = 'Assume strict molecular clock. \n' + self.display.string


class TC_Handler(Handler):

    def setattr(self, info, object, name, value):
        Handler.setattr(self, info, object, name, value)
        info.object._updated = True

    def object__updated_changed(self, info):
        if info.initialized:
            info.ui.title += "*"

class Tree(HasTraits):
    tree_file = File
    display = Instance(TextDisplay)
    view = View(Item('tree_file', show_label = True))

    def _tree_file_changed(self):
        if os.path.isfile(self.tree_file):
            self.display.string = 'Found tree file. \n' + self.display.string
        else:
            self.display.string = 'Alignment file does not exist. \n' + self.display.string

class Paralog(HasTraits):
    paralog1 = Str
    paralog2 = Str
    view = View(Item('paralog1', show_label = False),
                Item('paralog2', show_label = False))

class MainWindow(HasTraits):
    display = Instance(TextDisplay)

    alignment = Instance(Alignment)

    model = Instance(SupportedModels)

    tree = Instance(Tree)

    paralog = Instance(Paralog)

    def _alignment_default(self):
        return Alignment(display = self.display)

##    def _paralog_changed(self):
##        if not (self.paralog.paralog1 in self.alignment.alignment_file and
##                self.paralog.paralog2 in self.alignment.alignment_file):
##            self.display.string = 'Please make sure the alignment file is for this paralog. \n' + self.display.string
##                

    view = View(HSplit(Item('display', dock = 'vertical'),
                       Group(Item('alignment', style = 'custom'),
                             Item('_'),
                             Item('model', style = 'custom'),
                             Item('_'),
                             Item('tree', style = 'custom'),
                             Item('_'),
                             Item('paralog', style = 'custom')
                             ),
                       show_labels = True),
                height = 0.15, width = 0.25,
                buttons = [OKButton, CancelButton],
                title = 'IGC GUI',
                style = 'custom', resizable=True)


if __name__ == '__main__':
##    container = Container(camera=Camera(), display=TextDisplay())
##    container.configure_traits()
    test = SupportedModels()
##    test.configure_traits()

    test = MainWindow(model = SupportedModels(), alignment = Alignment(), display = TextDisplay(), tree = Tree(),
                      paralog = Paralog())
    test.configure_traits()
