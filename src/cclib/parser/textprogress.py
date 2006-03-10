
class TextProgress:
    
    def __init__(self):
        
        self.nstep=0
        self.text=None
        self.oldprogress=0
        self.progress=0
        
    def initialize(self,nstep,text=None):
        
        self.nstep=float(nstep)
        self.text=text
        
    def update(self,step):
        
        self.progress = int(step/self.nstep*100)
        #print step,self.nstep
        
        if self.progress/10==self.oldprogress/10+1: #just went through an interval of ten, ie. from 39 to 41, so update
            
            str="["
            for i in range(self.progress/10):
                str+="="
            for i in range(self.progress/10,10):
                str+="-"
                
            str+="]"
            
            print str
            
        self.oldprogress=self.progress
        return

