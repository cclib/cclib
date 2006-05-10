import sys
import thread
from qt import QProgressDialog,QString

class QtProgress(QProgressDialog):
    
    def __init__(self,parent):
        
        QProgressDialog.__init__(self,parent,"progress",True)
        
        self.nstep=0
        self.text=None
        self.oldprogress=0
        self.progress=0
        self.calls=0
        
        self.setCaption("Progress...")
        
    def initialize(self,nstep,text=None):
        
        self.nstep=nstep
        self.text=text
        self.setTotalSteps(nstep)
        if text:
          self.setLabelText(text)
        self.setProgress(1)
        #sys.stdout.write("\n")
        
    def update(self,step,text=None):

        self.setLabelText(text)
        self.setProgress(step)    
            
        return

