import numpy as np

class intervalsProbability:

    def __init__(self,lproba):
        self.lproba=lproba
        self.nbDecision=len(lproba)

    def isProper(self):
        sumLowerBound=0
        sumUpperBound=0
        for interval in self.lproba:
            sumLowerBound+=interval[0]
            sumUpperBound+=interval[1]
        if sumLowerBound<=1 and sumUpperBound>=1:
            return 1
        else:
            return 0

    def getLowerProbability(self,subset):
        if len(subset)!=self.nbDecision:
            raise Exception('Subset incompatible with the probabilities')
        else:
            lowerProbability=0
            for i in range(self.nbDecision):
                if subset[i]!=0:
                    lowerProbability+=self.lproba[i][0]
            return lowerProbability

    def getUpperProbability(self,subset):
        if len(subset)!=self.nbDecision:
            raise Exception('Subset incompatible with the probabilities')
        else:
            upperProbability=0
            for i in range(self.nbDecision):
                if subset[i]!=0:
                    upperProbability+=self.lproba[i][1]
            return upperProbability


    def setReachableProbability(self):
        if self.isProper()==1:
            self.lreachableProba=[]
            for i in range(0,self.nbDecision):
                subset=np.ones(self.nbDecision)
                subset[i]=0
                lb=max(lproba[i][0],1-self.getUpperProbability(subset))
                ub=min(lproba[i][1],1-self.getLowerProbability(subset))
                self.lreachableProba.append([lb,ub])
            self.lproba=self.lreachableProba

    def printProbability(self):
        str1,str2="lower bound |","upper bound |"
        str3=" "*13;
        i=0
        for interval in self.lproba:
            str3+="   y%d "%i
            str1+=" %.3f" % interval[0]
            str2+=" %.3f" % interval[1]
            i+=1
        print str3
        print " "*11, "-"*20
        print str1
        print str2


if __name__=='__main__':
    #lproba=[[0.1,0.2],[0.3,0.7],[0.5,0.9]]
    lproba=[[0.3,0.7],[0.5,0.9]]
    essai=intervalsProbability(lproba)
    print essai.isProper()
    print essai.getUpperProbability([1,1])
    essai.printProbability()
    essai.setReachableProbability()
    essai.printProbability()
