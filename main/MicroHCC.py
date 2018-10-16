'''
Class for Microscopib Hierarchical Conceptual Clustering
@author: Kadri Umbleja
'''


#microscopic similiarities
def MHCC(aList2,nre,nrf,objects,quantiles,till=1):
    nrq=len(quantiles)
    
    pp=[] #placeholder for quantile points
    kk=[] #placeholder for names
    nn=[] #element index holder for dissimilarity matrix
    
    #dissimilarity matrix
    diss=[ [ 0 for y in range( nre ) ] for x in range( nre ) ]
    #initialise placeholders
    for i in range(nre):
        kk.append(objects[i])
        nn.append([i])
    
    #min/max values for quantiles
    minQ=[[ float("inf") for y in range(nrq ) ] for y in range( nrf )]
    maxQ=[[ float("-inf") for y in range( nrq) ]  for y in range( nrf )]
    
    #initialize quantile values and min/max values
    for i in range(nre):
        ll=[]
        for f in range(nrf):
            a=[]
            for j in range(nrq):
                val=aList2[i].getHistogram(f).get_quantile(quantiles[j])
                if(val<minQ[f][j]):
                    minQ[f][j]=val
                if(val>maxQ[f][j]):
                    maxQ[f][j]=val    
                b=[val]
                a.append([val])
            ll.append(a)
        pp.append(ll)
  
    #clustering
    while (len(pp)>till):

        mm=float("inf")
        ind1=-1
        ind2=-1
        
        #over all elements left
        for i in range(len(pp)):
            for j in range(len(pp)):
                if(i>j):
                    d=0
                    nr=0
                    #over all features
                    for f in range(nrf):
                        #over all quantiles
                        for n in range(nrq):
                            if(minQ[f][n]!=maxQ[f][n]):                
                                mn=min(pp[i][f][n]+pp[j][f][n])
                                mx=max(pp[i][f][n]+pp[j][f][n])
                                d+=(mx-mn)/(maxQ[f][n]-minQ[f][n])
                                #for normalization later
                                nr+=1
                    #check if best result          
                    if(d/(nr)<mm):
                        mm=d/(nr)
                        ind1=i
                        ind2=j
        #merge features
        for f in range(nrf):
            for n in range(len(quantiles)):
                pp[ind2][f][n]+=pp[ind1][f][n]
        #update dissimilarity matrix
        for i in range(len(nn[ind1])):
            for j in range(len(nn[ind2])):
                diss[nn[ind1][i]][nn[ind2][j]]=mm
                diss[nn[ind2][j]][nn[ind1][i]]=mm
        
        #new name for concept
        kk[ind2]="("+str(kk[ind2])+"-"+str(kk[ind1])+")"
        nn[ind2]+=nn[ind1]
        #remove not needed 
        pp.pop(ind1)
        kk.pop(ind1)
        nn.pop(ind1)

    #return dissimilarity matrix
    return diss