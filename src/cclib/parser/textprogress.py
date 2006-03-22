import sys
import thread

class TextProgress:
    
    def __init__(self):
        
        self.nstep=0
        self.text=None
        self.oldprogress=0
        self.progress=0
        self.calls=0
        
    def initialize(self,nstep,text=None):
        
        self.nstep=float(nstep)
        self.text=text

        #sys.stdout.write("\n")
        
    def update(self,step,text=None):

        self.progress = int(step*100/self.nstep)
        
        if self.progress/2>=self.oldprogress/2+1 or self.text!=text: #just went through at least an interval of ten, ie. from 39 to 41, so update
            
            str="\r["
            prog=self.progress/10
            str+=prog*"="+(10-prog)*"-"
            str+="] %3i"%(self.progress)+"%"
            
            if text:
                str+="    Parsing "+text

            sys.stdout.write("\r"+70*" ")
            sys.stdout.flush()
            sys.stdout.write(str)
            sys.stdout.flush()
            self.oldprogress=self.progress
            
            if(self.progress>=100 and text=="Done"):
                print " "
            
            
        return

