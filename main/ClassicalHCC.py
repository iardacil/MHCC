'''
Class for HCC.
Two methods are available:
(mode=4) default method by Ichino based on concept sizes in feature space for interval valued data
(mode=0) modified version of previous for histogram valued data that considers histogram shape in more detail than just span/area in feature space
@author: Kadri Umbleja
'''

from main.Common import SymbolicObject


#mode=4 rectangle/concept size method by Ichino
#mode=0 histogram version by Umbleja
def HCC(aList2,nrf,nre,till,mode=4,Fmin=[],Fmax=[]):
    aList=[]
    #make copy of symbolic objects as they will be destroyed durning HCC
    for i in range(nre):
        aList.append(SymbolicObject(aList2[i].hist,aList2[i].name,aList2[i].types,id=aList2[i].id))
            
    #for normalization if needed
    feature_len=[]
    if(len(Fmin)==0):
        for i in range(nrf):
            Fmin.append(0)
    if(len(Fmax)==0):
        for i in range(nrf):
            Fmax.append(1)
    cnt=nrf        
    for i in range(nrf):
        feature_len.append(Fmax[i]-Fmin[i])

    #dissimilarity matrix
    diss=[ [ 0 for y in range( nre ) ] for x in range( nre ) ]
    
    while(len(aList)>till):
        curMin=float("inf")
        minj=-1
        mini=-1
        minW =None
        curr_size=len(aList)

        for i in range(curr_size):
            for j in range(curr_size):
                if(i>j):
                    if(mode==0): #histogram method
                        b = aList[i].combine(aList[j])
                        sum = b.calcCompactness(cnt)
                        if(curMin >sum):
                            curMin=sum
                            minW=b
                            mini=i
                            minj=j
                    elif(mode==4): #ichino method
                        b = aList[i].combineR(aList[j])
                        sum=b.get_avg_Span()
                        if(curMin >sum):
                            curMin=sum
                            minW=b
                            mini=i
                            minj=j

        print(str(curMin)+","+minW.getName())

        #update dissimilarity matrix
        for i in range(len(aList[mini].list)):
            for j in range(len(aList[minj].list)):
                k = aList[mini].list[i].id;
                m = aList[minj].list[j].id;
                diss[k][m]=curMin
                diss[m][k]=curMin
        #update objects
        aList.append(minW)
        aList.pop(mini)
        aList.pop(minj)


    return diss
