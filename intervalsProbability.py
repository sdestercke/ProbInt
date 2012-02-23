import numpy as np

class intervalsProbability:

    def __init__(self,lproba):
        if lproba.__class__.__name__!='ndarray':
            raise Exception('Expecting a numpy array as argument')
        if np.size(lproba[:,1])!=2:
            raise Exception('Array should contain two rows: top for upper prob, bottom for lower prob')
        self.lproba=lproba
        self.nbDecision=np.size(lproba[0])
        if np.all(lproba[0]>lproba[1])!=1:
            raise Exception('Some upper bounds lower than lower bounds')

    def isProper(self):
        if np.sum(self.lproba[1,:])<=1 and np.sum(self.lproba[0,:])>=1:
            return 1
        else:
            return 0

    def getLowerProbability(self,subset):
        if subset.__class__.__name__!='ndarray':
            raise Exception('Expecting a numpy array as argument')
        if np.size(subset)!=self.nbDecision:
            raise Exception('Subset incompatible with the frame size')
        if self.isReachable()==0:
            self.setReachableProbability()
        else:
            lowerProbability=max(np.sum(self.lproba[1,subset[:]==1]),1-np.sum(self.lproba[0,subset[:]==0]))
            return lowerProbability

    def getUpperProbability(self,subset):
        if subset.__class__.__name__!='ndarray':
            raise Exception('Expecting a numpy array as argument')
        if np.size(subset)!=self.nbDecision:
            raise Exception('Subset incompatible with the frame size')
        if self.isReachable()==0:
            self.setReachableProbability()
        else:
            upperProbability=min(np.sum(self.lproba[0,subset[:]==1]),1-np.sum(self.lproba[1,subset[:]==0]))
            return upperProbability

    def isReachable(self):
        for i in range(self.nbDecision):
            subset=np.ones(self.nbDecision)
            subset[i]=0
            if self.lproba[0,i] + np.sum(self.lproba[1,subset[:]==1])  > 1.0:
                return 0
            if self.lproba[1,i] + np.sum(self.lproba[0,subset[:]==1]) < 1.0:
                return 0
        return 1

    def setReachableProbability(self):
        if self.isProper()==1:
            lreachableProba=np.zeros((2,self.nbDecision))
            for i in range(self.nbDecision):
                subset=np.ones(self.nbDecision)
                subset[i]=0
                lb=max(self.lproba[1,i],1-np.sum(self.lproba[0,subset[:]==1]))
                ub=min(self.lproba[0,i],1-np.sum(self.lproba[1,subset[:]==1]))
                lreachableProba[1,i]=lb
                lreachableProba[0,i]=ub
            self.lproba[:]=lreachableProba[:]
        else:
            raise Exception('intervals inducing empty set: operation not possible')

    def printProbability(self):
        str1,str2="upper bound |","lower bound |"
        str3=" "*13;
        i=0
        for interval in range(self.nbDecision):
            str3+="   y%d "%i
            str1+=" %.3f" % self.lproba[0,interval]
            str2+=" %.3f" % self.lproba[1,interval]
            i+=1
        print str3
        print " "*11, "-"*20
        print str1
        print str2


if __name__=='__main__':
    lproba =np.array([[0.4,0.4,0.5],[0.2,0,0.2]])
    essai=intervalsProbability(lproba)
    print essai.isProper()
    print essai.getUpperProbability(np.array([1,1,1]))
    essai.printProbability()
    essai.setReachableProbability()
    print essai.getUpperProbability([1,0,1])
    essai.printProbability()
