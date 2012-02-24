import numpy as np

def getMaxCoherentIntervals(setOfInt):
    """Find and return the maximal subsets of coherent intervals from a list of intervals 

    Argument:
    setOfInt -- a 2xn array containing upper (1st row) and lower (2nd row) bounds of intervals    
    """
    # separate and setting up indices of sorted interval bounds
    upper=setOfInt[0,:]
    lower=setOfInt[1,:]
    indup=upper.argsort()
    indlow=lower.argsort()
    sortup=np.sort(upper)
    sortlow=np.sort(lower)
    # initializing variables for MCS detection
    MCSlist=[]
    currentMCS=[]
    ind_up=0
    ind_low=0 
    for i in range(2*upper.size-1):
        # deals with possible equality of next value
        if sortup[ind_up] == sortlow[ind_low]:
            currentMCS.append(indlow[ind_low])
            ind_low=ind_low+1
            # no more lower bounds encountered (reach end of upper)
            if ind_low == upper.size-1:
                print 'last lower bound reached, adding last MCS'
                MCSlist.append(currentMCS[:])
                print MCSlist
                break
            MCSlist.append(currentMCS[:])
        # first case: the next value is a lower bound
        elif sortlow[ind_low] < sortup[ind_up]:
            currentMCS.append(indlow[ind_low])
            print 'next value is lower bound, adding'
            print indlow[ind_low]
            print 'current MCS'
            print currentMCS
            ind_low=ind_low+1
            # no more lower bounds encountered (reach end of upper)
            if ind_low == upper.size-1:
                print 'last lower bound reached, adding last MCS'
                MCSlist.append(currentMCS[:])
                print MCSlist
                break
            # case where a lower bound is followed by an upper bound
            if min(sortlow[ind_low],sortup[ind_up]) < sortlow[ind_low]:
                print 'next value gonna be upper bound'
                print 'adding to List'
                print currentMCS
                MCSlist.append(currentMCS[:])
                print MCSlist
        # second case: the next value is an upper bound
        else:
            print 'removing element from MCS'
            print indup[ind_up]
            currentMCS.remove(indup[ind_up])
            ind_up=ind_up+1
    return MCSlist
        

class intervalsProbability:
    """Class of probability intervals: upper and lower prob. bounds on singletons

    Argument:
    lproba -- a 2xn array containing upper (1st row) and lower (2nd row) probabilistic bounds
    
    """
    
    def __init__(self,lproba):
        if lproba.__class__.__name__ != 'ndarray':
            raise Exception('Expecting a numpy array as argument')
        if lproba[:,1].size != 2:
            raise Exception('Array should contain two rows: top for upper prob, bottom for lower prob')
        if lproba.ndim != 2:
            raise Exception('Bad dimension of array: should contain 2 dimensions')
        self.lproba=lproba
        self.nbDecision=lproba[0].size
        if np.all(lproba[0] >=lproba[1]) != 1:
            raise Exception('Some upper bounds lower than lower bounds')

    def isProper(self):
        """Check if probability intervals induce a non-empty probability set. 
        
        Return 0 (empty) or 1 (non-empty).
        
        """
        if self.lproba[1,:].sum()<=1 and self.lproba[0,:].sum()>=1:
            return 1
        else:
            return 0

    def getLowerProbability(self,subset):
        """Compute lower probability of an event expressed in binary code. 
        
        Argument:
        subset -- a 1xn vector containing 1 for elements in the event, 0 otherwise.
        
        Return lower probability value.
        
        """
        if subset.__class__.__name__!='ndarray':
            raise Exception('Expecting a numpy array as argument')
        if subset.size != self.nbDecision:
            raise Exception('Subset incompatible with the frame size')
        if self.isReachable()==0:
            self.setReachableProbability()
        else:
            lowerProbability=max(self.lproba[1,subset[:]==1].sum(),1-self.lproba[0,subset[:]==0].sum())
            return lowerProbability

    def getUpperProbability(self,subset):
        """Compute upper probability of an event expressed in binary code. 
        
        Argument:
        subset -- a 1xn vector containing 1 for elements in the event, 0 otherwise.
        
        Return upper probability value.
        
        """    
        if subset.__class__.__name__!='ndarray':
            raise Exception('Expecting a numpy array as argument')
        if subset.size != self.nbDecision:
            raise Exception('Subset incompatible with the frame size')
        if self.isReachable()==0:
            self.setReachableProbability()
        else:
            upperProbability=min(self.lproba[0,subset[:]==1].sum(),1-self.lproba[1,subset[:]==0].sum())
            return upperProbability

    def isReachable(self):
        for i in range(self.nbDecision):
            subset=np.ones(self.nbDecision)
            subset[i]=0
            if self.lproba[0,i] + self.lproba[1,subset[:]==1].sum()  > 1.0:
                return 0
            if self.lproba[1,i] + self.lproba[0,subset[:]==1].sum() < 1.0:
                return 0
        return 1

    def setReachableProbability(self):
        if self.isProper()==1:
            lreachableProba=np.zeros((2,self.nbDecision))
            for i in range(self.nbDecision):
                subset=np.ones(self.nbDecision)
                subset[i]=0
                lb=max(self.lproba[1,i],1-self.lproba[0,subset[:]==1].sum())
                ub=min(self.lproba[0,i],1-self.lproba[1,subset[:]==1].sum())
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
        
class setOfIntProba:
    """Class to handle sets of Int Proba

    Argument:
    intlist -- a mx2xn array containing upper (1st row) and lower (2nd row) bounds
               of the m probability intervals
               
               dim 1: index of probability set
               dim 2: lower or upper prob bounds
               dim 3: values of bounds on each element
    """
    
    def __init__(self,intlist):
        if intlist.__class__.__name__ != 'ndarray':
            raise Exception('Expecting a numpy array as argument')
        if intlist.ndim != 3:
            raise Exception('Expecting a 3-dimensional array')
        self.intlist=intlist
        self.nbProbInt=intlist[:,0,0].size
        self.nbDecision=intlist[0,0].size
        
    def areCompatible(self):
        """Check whether the set of probability intervals are compatible, i.e., if the conjunction is non-empty.
        
        Return 1 if non-empty, 0 if empty
        """
        min=self.intlist[:,1,:].max(axis=0)
        max=self.intlist[:,0,:].min(axis=0)
        if min.sum() <= 1 and max.sum() >= 1:
            return 1
        else:
            return 0
            
    def conjunction(self):
        """Perform a conjunctive merging of the set of probability intervals
        
        Return a possibly non-proper intervalsProbability class object.
        """
        fusedproba=np.zeros((2,self.nbDecision))
        for i in range(self.nbDecision):
            subset=np.ones(self.nbDecision)
            subset[i]=0
            lb=max(self.intlist[:,1,i].max(),1-self.intlist[:,0,subset[:]==1].min(axis=0).sum())
            ub=min(self.intlist[:,0,i].min(),1-self.intlist[:,1,subset[:]==1].max(axis=0).sum())
            fusedproba[1,i]=lb
            fusedproba[0,i]=ub
        result=intervalsProbability(fusedproba)
        return result
        
    def disjunction(self):
        """Perform a disjunctive merging of the set of probability intervals
        
        Return an intervalsProbability class object.
        """
        fusedproba=np.zeros((2,self.nbDecision))
        for i in range(self.nbDecision):
            subset=np.ones(self.nbDecision)
            subset[i]=0
            lb=self.intlist[:,1,i].min()
            ub=self.intlist[:,0,i].max()
            fusedproba[1,i]=lb
            fusedproba[0,i]=ub
        result=intervalsProbability(fusedproba)
        return result        

if __name__=='__main__':
    lproba =np.array([[0.4,0.4,0.5],[0.2,0,0.2]])
    essai=intervalsProbability(lproba)
    print essai.isProper()
    print essai.getUpperProbability(np.array([1,1,1]))
    essai.printProbability()
    essai.setReachableProbability()
    print essai.getUpperProbability(np.array[1,0,1])
    essai.printProbability()
    
    setproba=np.array([[[0.6,0.4,0.5],[0.2,0.3,0.2]],[[0.7,0.3,0.6],[0.2,0.,0.2]],
                    [[0.7,0.3,0.6],[0.2,0.,0.2]],[[0.6,0.4,0.5],[0.2,0.3,0.2]]])
    test=Intprob.setOfIntProba(setproba)
    test.areCompatible()
    test.conjunction()

