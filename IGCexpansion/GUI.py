# Trying to create a GUI for my IGCexpansion package
# Xiang ji
# xji3@ncsu.edu
from traits.etsconfig.api import *
ETSConfig.toolkit = 'qt4' #'wx'
from traitsui.api import *
from traits.api import *

class TextDisplay(HasTraits):
    string = String()

    view = View( Item('string', show_label = False, springy = True, style = 'custom'))

# This should be the Model part in MVC structure
class Alignment(HasTraits):
    Alignment_file = File
    display = Instance(TextDisplay)
    #Align   = Button()
    _updated = Bool(False)

   # Display specification (one Item per editor style):
    file_group = Group(
        Item( 'file_name', style = 'simple',   label = 'Simple' ),
        Item( '_' ),
        Item( 'file_name', style = 'custom',   label = 'Custom' ),
        Item( '_' ),
        Item( 'file_name', style = 'text',     label = 'Text' ),
        Item( '_' ),
        Item( 'file_name', style = 'readonly', label = 'ReadOnly' )
    )

    view = View( Item('Alignment_file', show_label = False))#, Item('Align', show_label = False))

    def _Alignment_file_changed(self):
        self.Alignment = '  '
        
    
class TC_Handler(Handler):

    def setattr(self, info, object, name, value):
        Handler.setattr(self, info, object, name, value)
        info.object._updated = True

    def object__updated_changed(self, info):
        if info.initialized:
            info.ui.title += "*"



class MainWindow(HasTraits):
    display = Instance(TextDisplay, ())

    alignment = Instance(Alignment)

    def _alignment_default(self):
        return Alignment(display = self.display)
    view = View('display', 'alignment',
                buttons = [OKButton, CancelButton],
                title = 'IGC GUI',
                style = 'custom', resizable=True)

view_test = View(Item(name = 'Alignment_file'),
                 title = 'Alter Title',
                 handler = TC_Handler(),
                 #buttons = [OKButton, CancelButton, Align])
                 buttons = ['OK', 'Cancel'])

if __name__ == '__main__':
    test = MainWindow()
    test.configure_traits()
