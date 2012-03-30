import numpy as np
import Orange
import orngTree, orngEnsemble
from intervalsProbability import *
import random
import time as t

def test_simpleTree(training, test):
    """Function that takes a training and test data sets, build a tree and return decision
    """
    start = t.time()
    tree_learn = Orange.classification.tree.TreeLearner(minExamples=2, mForPrunning=2, 
                                                        sameMajorityPruning=True, name='tree')
    single_tree = tree_learn(training)

    prec=0
    for i in range(len(test)):
        if single_tree(test[i])==test[i].getclass():
            prec+=1.
    acc=prec/len(test)

    return acc,t.time()-start

def test_forestFusion(training,test,combMethod,nbTree=15):
    """Function that takes a training and test data sets, build forests and return decisions
    """
    start = t.time()

    s=4
    accuracy=0.
    set_accuracy=0.
    disc_accuracy=0.
    nb_classes=len(test.domain.class_var.values)
    tree_learn = Orange.classification.tree.TreeLearner(minExamples=2, mForPrunning=2, 
                            sameMajorityPruning=True, name='tree')
    forest = Orange.ensemble.forest.RandomForestLearner(trees=nbTree, base_learner=tree_learn,rand=random.Random(0))
    result = forest(training)
 
    for j in range(len(test)):
        setofprob=[]
        for i in range(len(result.classifiers)):
            low=np.zeros(nb_classes)
            up=np.zeros(nb_classes)
            answer=result.classifiers[i].descender(result.classifiers[i].tree,test[j])
            divide=sum(answer[0].distribution)+s
            for k in range(nb_classes):
                low[k]=(answer[0].distribution[k])/divide
                up[k]=(answer[0].distribution[k]+s)/divide
            prob=np.array([up,low])
            setofprob.append(prob[:])
        
        resultingset=setOfIntProba(np.array(setofprob))
        #resultingcomb=resultingset.almostMCScomb()
        #resultingcomb=resultingset.mostMCSconj()
        #resultingcomb=resultingset.bestfirstMCS(1)
        #resultingcomb=resultingset.meanfirstMCSweighted(5)
        resultingcomb=resultingset.runCombination(combMethod)
        decision=resultingcomb.nc_hurwicz_decision(0.5)   
        if test[j].getclass()==test.domain.class_var.values[decision]:
            accuracy=accuracy+1

        decision_max=resultingcomb.nc_maximal_decision() 

        true_class=np.zeros(nb_classes)
        for k in range(nb_classes):
            if test[j].getclass()==test.domain.class_var.values[k]:
                true_class[k]=1
   
        if any(np.minimum(decision_max,true_class)==1):
            set_accuracy=set_accuracy+1
            disc_accuracy=disc_accuracy+(1./decision_max.sum())
          

    accuracy=accuracy/len(test)
    set_accuracy=set_accuracy/len(test)
    disc_accuracy=disc_accuracy/len(test)
    
    return accuracy, set_accuracy, disc_accuracy, t.time()-start



def test_forestVote(training,test,nbTree=15):
    """Function that takes a training and test data sets, build forests and return decisions
    """
    start = t.time()

    accuracy=0.
    nb_classes=len(test.domain.class_var.values)
    tree_learn = Orange.classification.tree.TreeLearner(minExamples=2, mForPrunning=2, 
                            sameMajorityPruning=True, name='tree')
    forest = Orange.ensemble.forest.RandomForestLearner(trees=nbTree, base_learner=tree_learn,rand=random.Random(0))    
    result = forest(training)

    for i in range(len(test)):
        y=result(test[i],Orange.classification.Classifier.GetValue)
        if y==test[i].getclass():
            accuracy+=1
    
    accuracy=accuracy/len(test)

    return accuracy, t.time()-start



if __name__=='__main__':
    import orange
    #data = orange.ExampleTable('bupa.tab')
    #data=Orange.data.Table("satimage/satimage")
    #data=Orange.data.Table("segment/segment")
    #data=Orange.data.Table("audiology") # data.domain, len(data)=> nb instance, len(data.domain)=> nb attribut
    #data=Orange.data.Table("phoneme.tab")
    #data=Orange.data.Table("wine")
    data=Orange.data.Table("zoo")
    #print len(data), len(data.domain), len(data.domain.class_var.values)
    indices = Orange.data.sample.SubsetIndices2(p0=0.25) # 
    ind=indices(data)
    iristr = data.select(ind, 0) # ind=>vect de bool de taille du nb d'instance, 0=> prendre tous les indices
    iristst = data.select(ind, 1) # 1=> prendre l'inverse des indices

    nbTree=[1, 20, 50, 100, 150, 200]
    prec=np.zeros(len(nbTree))
    accuracy=np.zeros(len(nbTree))
 
    for iTree in range(0,len(nbTree)):
        prec,tpsst=test_simpleTree(iristr,iristst)
        [accAlmost, set_accAlmost, disc_accAlmost,tpsffAlmost]=test_forestFusion(iristr,iristst,"almostMCScomb",nbTree[iTree])
        [accConj, set_accConj, disc_accConj,tpsffConj]=test_forestFusion(iristr,iristst,"mostMCSconj",nbTree[iTree])
        [accBFirst, set_accBFirst, disc_accBFirst,tpsffBFirst]=test_forestFusion(iristr,iristst,"bestfirstMCS",nbTree[iTree])
        [accMFirst, set_accMFirst, disc_accMFirst,tpsffMFirst]=test_forestFusion(iristr,iristst,"meanfirstMCSweighted",nbTree[iTree])
        accVote,tpsfv=test_forestVote(iristr,iristst,nbTree[iTree])

        print " & %d" %nbTree[iTree] ," & %.2f " %prec, " & %.2f " %accVote, " & %.2f " %disc_accConj, " & %.2f "  %set_accConj, " & %.2f" %accConj, " & %.2f " %disc_accAlmost, " & %.2f "  %set_accAlmost, " & %.2f" %accAlmost, " & %.2f " %disc_accBFirst, " & %.2f "  %set_accBFirst, " & %.2f" %accBFirst, " & %.2f " %disc_accMFirst, " & %.2f "  %set_accMFirst, " & %.2f" %accMFirst, "\\\\"
