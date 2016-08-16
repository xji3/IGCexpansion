# Trying to create a GUI for my IGCexpansion package
# Xiang ji
# xji3@ncsu.edu
from traits.etsconfig.api import *
from threading import Thread
ETSConfig.toolkit = 'qt4' #'wx'
from traitsui.api import *
from traits.api import *
import os
from IGCexpansion.CodonGeneconv import ReCodonGeneconv
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

##    def _alignment_file_changed(self):
##        if os.path.isfile(self.alignment_file):
##            self.display.string = 'Found alignment file. \n' + self.display.string
##        else:
##            self.display.string = 'Alignment file does not exist. \n' + self.display.string

# model enumeration list
class SupportedModels(HasTraits):
    display = Instance(TextDisplay)
    #display = TextDisplay()
    model = Enum('MG94', 'HKY',)
    force = Bool
    clock = Bool

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
    display = Instance(TextDisplay)
    #display = TextDisplay()
    tree_file = File
    view = View(Item('tree_file', show_label = False))

##    def _tree_file_changed(self):
##        if os.path.isfile(self.tree_file):
##            self.display.string = 'Found tree file. \n' + self.display.string
##        else:
##            self.display.string = 'Alignment file does not exist. \n' + self.display.string

class Paralog(HasTraits):
    paralog1 = Str
    paralog2 = Str
    view = View(Item('paralog1', show_label = False),
                Item('paralog2', show_label = False))

class Save(HasTraits):
    save_file = File
    display = Instance(TextDisplay)
    #display = TextDisplay()
    view = View(Item('save_file', show_label = False),
                      
                style = 'simple')

##    def _save_file_changed(self):
##        if os.path.isfile(self.save_file):
##            self.display.string = 'Found save file. Will overwrite it. \n' + self.display.string
##        else:
##            self.display.string = 'Save file does not exist. Will create it. \n' + self.display.string

class Summary(HasTraits):
    summary_folder = File
    #display = Instance(TextDisplay)
    display = TextDisplay()
    view = View(Item('summary_folder', show_label = False),
                style = 'simple')


##    def _summary_folder_changed(self):
##        if os.path.isdir(self.summary_folder):
##            self.display.string = 'Found summary folder. \n' + self.display.string
##        else:
##            self.display.string = 'Summary folder does not exist.\n' + self.display.string

class ApplySettings(Thread):
    def run(self):
        mle_paralog = [str(self.paralog.paralog1), str(self.paralog.paralog2)]
        mle_model = self.model.model
        mle_force = self.model.force
        mle_clock = self.model.clock
        mle_newicktree = self.tree.tree_file
        if mle_force:
            if mle_model == 'MG94':
                Force = {5:0.0}
            elif mle_model == 'HKY':
                Force = {4:0.0}
        else:
            Force = None
        mle_alignment_file = self.alignment.alignment_file
        mle_save_name = self.save.save_file
        mle_summary_folder = self.summary.summary_folder
        summary_name_appendix = '_summary.txt'
        if mle_clock:
            summary_name_appendix = '_clock' + summary_name_appendix
        else:
            summary_name_appendix = '_nonclock' + summary_name_appendix
        if mle_force:
            summary_name_appendix = '_force' + summary_name_appendix
        else:
            summary_name_appendix = '_noforce' + summary_name_appendix
        mle_summary_file = mle_summary_folder + '/' + '_'.join(mle_paralog) + '_' + mle_model + summary_name_appendix
        mle_sitewise_summary_file = mle_summary_folder + '/' + '_'.join(mle_paralog) + '_' + mle_model + summary_name_appendix.replace('_summary.txt', '_sitewise_summary.txt')

        self.test = ReCodonGeneconv( mle_newicktree, mle_alignment_file, mle_paralog, Model = mle_model, Force = Force, clock = mle_clock, save_name = mle_save_name)
        

class MainWindow(HasTraits):

    display = Instance(TextDisplay, ())
    #display = TextDisplay()

    results_string = String()

    alignment = Instance(Alignment)

    model = Instance(SupportedModels)

    tree = Instance(Tree)

    paralog = Instance(Paralog)

    save = Instance(Save)

    summary = Instance(Summary)

    apply_thread = Instance(ApplySettings)
    apply_button = Button()

    test = None

    

       
    def _alignment_default(self):
        return Alignment(display = self.display)

    def _tree_default(self):
        return Tree(display = self.display)
    
##    def _paralog_changed(self):
##        if not (self.paralog.paralog1 in self.alignment.alignment_file and
##                self.paralog.paralog2 in self.alignment.alignment_file):
##            self.display.string = 'Please make sure the alignment file is for this paralog. \n' + self.display.string
##                
       
    
    view = View(#HSplit(Item('display', dock = 'vertical'),
                       Group(Item('alignment', style = 'custom'),
                             Item('_'),
                             Item('model', style = 'custom'),
                             Item('_'),
                             Item('tree', style = 'custom', label = 'Tree file'),
                             Item('paralog', style = 'custom'),
                             Item('_'),
                             Item('save', style = 'custom', label = 'Save file'),
                             Item('_'),
                             Item('summary', style = 'custom', label = 'Summary Folder'),
                             #Item('_'),
                             #Item('apply_button', label = 'Apply Settings')
                             
                       show_labels = True),
                height = 0.20, width = 0.25,
                buttons = [OKButton, CancelButton],
                title = 'IGC GUI',
                style = 'custom', resizable=True)

    def add_line(self, string):
        """ Adds a line to the textbox display.
        """
        self.results_string = (string + '\n' + self.results_string)[:1000]

    def _apply_button_fired(self):
        if self.apply_thread and self.apply_thread.isAlive():
            self.apply_thread.wats_abort = True
        else:
            self.apply_thread = ApplySettings()
            self.apply_thread.wants_abort = False
            self.apply_thread.display = self.display
            self.apply_thread.paralog = self.paralog
            self.apply_thread.model   = self.model
            self.apply_thread.tree    = self.tree
            self.apply_thread.alignment = self.alignment
            self.apply_thread.save    = self.save
            self.apply_thread.summary = self.summary
            self.apply_thread.test    = self.test
            self.apply_thread.start()

def run_mle(test):
    mle_paralog = [str(test.paralog.paralog1), str(test.paralog.paralog2)]
    mle_model = test.model.model
    mle_force = test.model.force
    mle_clock = test.model.clock
    mle_newicktree = test.tree.tree_file
    if mle_force:
        if mle_model == 'MG94':
            Force = {5:0.0}
        elif mle_model == 'HKY':
            Force = {4:0.0}
    else:
        Force = None
    mle_alignment_file = test.alignment.alignment_file
    mle_save_name = test.save.save_file
    mle_summary_folder = test.summary.summary_folder
    summary_name_appendix = '_summary.txt'
    if mle_clock:
        summary_name_appendix = '_clock' + summary_name_appendix
    else:
        summary_name_appendix = '_nonclock' + summary_name_appendix
    if mle_force:
        summary_name_appendix = '_force' + summary_name_appendix
    else:
        summary_name_appendix = '_noforce' + summary_name_appendix
    mle_summary_file = mle_summary_folder + '/' + mle_model + '_' + '_'.join(mle_paralog) + summary_name_appendix
    mle_sitewise_summary_file = mle_summary_folder + '/' + mle_model + '_' + '_'.join(mle_paralog) + summary_name_appendix.replace('_summary.txt', '_sitewise_summary.txt')

    mle_test = ReCodonGeneconv( mle_newicktree, mle_alignment_file, mle_paralog, Model = mle_model, Force = Force, clock = mle_clock, save_name = mle_save_name)
    mle_test.get_mle(True, True, 0, 'BFGS')
    mle_test.get_individual_summary(summary_path = mle_summary_folder, file_name = mle_summary_file)
    mle_test.get_SitewisePosteriorSummary(summary_path = mle_summary_folder, file_name = mle_sitewise_summary_file)
    return mle_test
        

if __name__ == '__main__':
    test = MainWindow(model = SupportedModels(), alignment = Alignment(), display = TextDisplay(), tree = Tree(),
                      paralog = Paralog(), save = Save(), summary = Summary())

    test.configure_traits()
    mle_test = run_mle(test)
    
