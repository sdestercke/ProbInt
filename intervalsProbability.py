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
            # no more lower bounds encountered (reach end of lower)
            if ind_low == upper.size:
                MCSlist.append(currentMCS[:])
                break
            MCSlist.append(currentMCS[:])
        # first case: the next value is a lower bound
        elif sortlow[ind_low] < sortup[ind_up]:
            currentMCS.append(indlow[ind_low])
            ind_low=ind_low+1
            # no more lower bounds encountered (reach end of lower)
            if ind_low == upper.size:
                MCSlist.append(currentMCS[:])
                break
            # case where a lower bound is followed by an upper bound
            if min(sortlow[ind_low],sortup[ind_up]) < sortlow[ind_low]:
                MCSlist.append(currentMCS[:])
        # second case: the next value is an upper bound
        else:
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
        upperProbability=min(self.lproba[0,subset[:]==1].sum(),1-self.lproba[1,subset[:]==0].sum())
        return upperProbability

    def isReachable(self):
        """Check if the probability intervals are reachable (are coherent / correspond to tightest possible 
        constraints) 
        
        Return a 0/1 value (1: are reachable).
        
        """    
        for i in range(self.nbDecision):
            subset=np.ones(self.nbDecision)
            subset[i]=0
            if self.lproba[0,i] + self.lproba[1,subset[:]==1].sum()  > 1.0:
                return 0
            if self.lproba[1,i] + self.lproba[0,subset[:]==1].sum() < 1.0:
                return 0
        return 1

    def setReachableProbability(self):
        """Make the bounds reachable and return them. 
        """    
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
            
    def nc_maximin_decision(self):
        """Return the maximin classification decision (nc: no costs)
        """
        if self.isReachable()==0:
            self.setReachableProbability()
        return self.lproba[1,:].argmax()
        
    def nc_maximax_decision(self):
        """Return the maximax classification decision (nc: no costs)
        """
        if self.isReachable()==0:
            self.setReachableProbability()
        return self.lproba[0,:].argmax()
        
    def nc_maximal_decision(self):
        """Return the classification decisions that are maximal (nc: no costs)
        """
        if self.isReachable()==0:
            self.setReachableProbability()
        maximality_classe=np.ones(nb_classes)
        for i in range(nb_classes):
            for j in range(nb_classes):
                if i != j and maximality_classe[i] == 1 and maximality_classe[j] == 1:
                    if -self.lproba[0,j]+self.lproba[1,i] > 0:
                        maximality_classe[j]=0
        return 0

    def printProbability(self):
        """Print the current bounds 
        """  
        str1,str2="upper bound |","lower bound |"
        str3=" "*13;
        i=0
        for interval in range(self.nbDecision):
            str3+="   y%d " %i
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
        comp=1
        min=self.intlist[:,1,:].max(axis=0)
        max=self.intlist[:,0,:].min(axis=0)
        if min.sum() > 1 or max.sum() < 1:
            comp=0
        for i in range(self.nbDecision):
            if min[i] > max[i]:
                comp=0
        return comp
            
    def conjunction(self):
        """Perform a conjunctive merging of the set of probability intervals
        
        Return a possibly non-proper intervalsProbability class object.
        """
        if self.areCompatible() == 0:
            raise Exception('Probability intervals not compatible, conjunction empty') 
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
        
    def getalmostMCS(self):
        """Internal function to get almost MCS probInt, in order to fusion them.
        
        Return the set of 'almost' MCS at the end"""
        #Initialize MCS as all probability intervals
        listofMCS=[range(self.nbProbInt)]
        for j in range(self.nbDecision):
            temp_list=[]
            for i in range(len(listofMCS)):
                test=setOfIntProba(self.intlist[listofMCS[i],:,:])
                if test.areCompatible() == 1:
                    temp_list.append(listofMCS[i][:])
                else:
                    MCS=getMaxCoherentIntervals(self.intlist[listofMCS[i][:],:,j].transpose())
                    for l in range(len(MCS)):
                        sub_MCS=[]
                        for k in MCS[l]:
                            sub_MCS.append(listofMCS[i][k])
                        temp_list.append(sub_MCS[:])
            listofMCS=temp_list[:]
        listofMCS=[]
        # elmminating redundant MCS
        for elt in temp_list:
            try:
                ind=listofMCS.index(elt)
            except:
                listofMCS.append(elt)
        return listofMCS
            
        
    def almostMCScomb(self):
        """get a list of 'almost' MCS and perform a combination according to it.
        
        Return a proper probability intervals"""
        list=self.getalmostMCS()
        nbMCS=len(list)
        conj=[]
        setofdisj=[]
        for i in range(nbMCS):
            setofprob=setOfIntProba(self.intlist[list[i],:,:])
            if setofprob.areCompatible() == 0:
                setofprob.discountnoncomp
            conj=setofprob.conjunction()
            resconj=np.array([conj.lproba[0,:],conj.lproba[1,:]])
            setofdisj.append(resconj)
        setofprob2=setOfIntProba(np.array(setofdisj))
        return setofprob2.disjunction()
            
        
    def mostMCSconj(self):
        """get a list of 'almost' MCS and perform a conjunctive combination on the MCS
        counting the most elements (in case of ties, first one is chosen)
        
        Return a proper probability intervals
        """
        list=self.getalmostMCS()
        nbMCS=len(list)
        nbsetinMCS=np.zeros(nbMCS)
        for i in range(nbMCS):
            nbsetinMCS[i]=len(list[i])
        setofprob=setOfIntProba(self.intlist[list[nbsetinMCS.argmax()],:,:])
        if setofprob.areCompatible() == 0:
            setofprob.discountnoncomp
        return setofprob.conjunction()
        
        
    def discountnoncomp(self):
        """return set of discounted probability intervals if they are not compatible
        """
        if self.areCompatible() == 0:
            if min.sum() - 1 > 0:
                epsilon_l=1./min.sum()
            if max.sum() - 1 < 0:
                epsilon_u=(1.-self.nbDecision)/(max.sum()-self.nbDecision)
        discount=min(epsilon_l,epsilon_u)
        self.intlist[:,1,:]=self.intlist[:,1,:]*discount
        self.intlist[:,0,:]=self.intlist[:,0,:]*discount+(1-discount)

if __name__=='__main__':
    lproba =np.array([[0.4,0.4,0.5],[0.2,0,0.2]])
    essai=intervalsProbability(lproba)
    print essai.isProper()
    print essai.getUpperProbability(np.array([1,1,1]))
    essai.printProbability()
    essai.setReachableProbability()
    print essai.getUpperProbability(np.array([1,0,1]))
    essai.printProbability()
    
    setproba=np.array([[[0.6,0.5,0.2],[0.4,0.3,0.]],[[0.55,0.55,0.2],[0.35,0.35,0.]],
                    [[0.5,0.2,0.6],[0.3,0.,0.4]],[[0.35,0.6,0.35],[0.15,0.4,0.15]]])
    test=setOfIntProba(setproba)
    test.areCompatible()
    test.getalmostMCS()

