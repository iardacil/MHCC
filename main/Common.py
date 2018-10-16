'''
Class for handeling symbolic data
SymbolicObject consists of Histograms
*there are methods for both SymbolicObject and Histogram that are not needed for current task
@author: Kadri Umbleja
'''

#Creates symbolic objects and normalizes data
def normalize(X,objects,titles,type,mode=0):
    nre= len(X)
    nrf = len(titles)
    aList=[]
    for i in range(nre):
        aList.append(SymbolicObject(X[i],objects[i],type,id=i))

    
    Fmin=[float("inf") for i in range(nrf)]
    Fmax=[float("-inf") for i in range(nrf)]
    for j in range(nrf):
        for i in range(nre):
            #by default this method is used - empty bins are ignored in normalization
            if(mode==0):
                if(aList[i].getHistogram(j).get_nonZeroS()<Fmin[j]):
                    Fmin[j]=aList[i].getHistogram(j).get_nonZeroS()
                if(aList[i].getHistogram(j).get_nonZeroE()>Fmax[j]):
                    Fmax[j]=aList[i].getHistogram(j).get_nonZeroE()
            #even if bins have probability 0, they are considered
            elif(mode==1):
                if(aList[i].getHistogram(j).start[0]<Fmin[j]):
                    Fmin[j]=aList[i].getHistogram(j).start[0]
                if(aList[i].getHistogram(j).end[aList[i].getHistogram(j).nrofBins-1]>Fmax[j]):
                    Fmax[j]=aList[i].getHistogram(j).end[aList[i].getHistogram(j).nrofBins-1]
    aList2=[]

    #create symbolic objects
    for i in range(nre):
        aList2.append(SymbolicObject(X[i],objects[i],type,Fmin,Fmax,id=i))
        
    return aList2,Fmin,Fmax

class SymbolicObject(object):
    def __init__(self,hist,name,types=[],Fmin=[],Fmax=[],id=0,nr=1,A = None,B = None,C=[],cluster=0,l=[]):
        self.name=name
        self.hist=hist
        self.id=id
        self.nr=nr
        self.list=[]
        self.types=types
        self.histograms=[]
        self.cluster=cluster
        self.l=l
        self.A=A
        self.B=B
        self.Fmin=Fmin
        self.Fmax=Fmax
        if(A is not None and B is not None):
            self.list.extend(A.getList())        
            self.list.extend(B.getList())
            self.nr=A.getNr()+B.getNr() 

        if(len(self.list)==0):
            self.list.append(self)
        abi=hist.split("},")
        self.nrf=len(abi)
        if not types:
            for i in range (self.nrf):
                types.append(0)
        normalize=1
        if not Fmin:
            normalize=0
            for i in range (self.nrf):
                Fmin.append(0)
                Fmax.append(1)
        if(len(C)==0):
            for i in range(self.nrf):
                C.append(1)
        self.C=C
        for i in range (self.nrf):
            self.histograms.append(Histogram(abi[i],id,e=0,type=types[i],nr=self.nr,Fmin=Fmin[i],Fmax=Fmax[i]))
        if(normalize==1):
            self.hist=""
            for n in range(self.nrf):
                self.hist+=self.histograms[n].getHist()+","
            self.hist=self.hist[:-1]
        
    def set_cluster(self,cluster):
        self.cluster=cluster
    def get_cluster(self):
        return self.cluster            
            
    def getList(self):
        return self.list;
    def getNr(self):
        return self.nr;
    def getHist(self):
        return self.hist;
    def getHistRounded(self,nrR=0):
        sttr=""
        for i in range(self.nrf):
            sttr+=self.histograms[i].getHistRounded(nrR)+","
        return sttr[:-1]
    def getName(self):
        return self.name;
    def getHistogram(self,nr):
        if nr<self.nrf:
            return self.histograms[nr]
        else:
            return None

    def get_avg_Span(self,l=None):
        if(l==None):
            l=[]
            for i in range(self.nrf):
                l.append(1)
        a=0
        cnt=0
        for i in range(self.nrf):
            if(l[i]==1):
                a+=self.histograms[i].getSpan()
                cnt+=1
        return a/cnt
    def get_avg_compactnessR(self,l=None):
        if(l==None or l==[]):
            l=[]
            for i in range(self.nrf):
                l.append(1)

        cnt=0
        a=0
        for i in range(self.nrf):
            if(l[i]==1):
                a+=self.histograms[i].compactnessR()
                cnt+=1
        return a/cnt
    def get_avg_compactness(self):
        a=0
        for i in range(self.nrf):
            a+=self.histograms[i].get_compactness()
        return a/self.nrf
    def calcCompactness(self,l=[],cnt=0):

        if(len(l)==0):
            for i in range(self.nrf):
                l.append(1)
        cnt=0        
        for i in range(self.nrf):
            if(l[i]==1):
                cnt+=1
        total=0
        for i in range(len(self.list)):
            a=self.list[i]
            sum=0
            for j in range(self.nrf):
                if(l[j]==1):
                    sum+=self.histograms[j].dissimilarity(a.histograms[j])
#                 print(self.histograms[j].dissimilarity(a.histograms[j]))
            total+=sum
        return  total/(self.nr*cnt)

            
    def meetR(self,H):
        meet=""
        name="("+self.getName()+"-"+H.getName()+")"
        for i in range(self.nrf):
            meet+=self.getHistogram(i).meetR(H.getHistogram(i))+","
        return SymbolicObject(meet[:-1],name,self.types,cluster=self.cluster)
    def meet(self,H):
        meet=""
        name="("+self.getName()+"-"+H.getName()+")"
        for i in range(self.nrf):
            meet+=self.getHistogram(i).meet(H.getHistogram(i))+","
        return SymbolicObject(meet[:-1],name,self.types,cluster=self.cluster)
    
    def joinR(self,H):
        join=""
        name="("+self.getName()+"-"+H.getName()+")"
        for i in range(self.nrf):
            join+=self.getHistogram(i).joinR(H.getHistogram(i))+","
        ll=[]
#         for i in range(self.nrf):
#             if(self.l[i]==0 and H.l[i]==0):
#                 ll.append(0)
#             else:
#                 ll.append(1)
        return SymbolicObject(join[:-1],name,self.types,cluster=self.cluster,l=ll)
    
    def join(self,H):
        join=""
        name="("+self.getName()+"-"+H.getName()+")"
        for i in range(self.nrf):
            join+=self.getHistogram(i).join(H.getHistogram(i))+","
        return SymbolicObject(join[:-1],name,self.types,cluster=self.cluster)
        
    def compactnessR(self):
        comp=0
        for i in range(self.nrf):
            comp+=self.getHistogram(i).compactnessR()
        return comp/self.nrf
    def compactness(self):
        comp=0
        for i in range(self.nrf):
            comp+=self.getHistogram(i).compactness()
        return comp/self.nrf
    
    def purenessR(self,H):
        comp=0
        for i in range(self.nrf):
            comp+=self.getHistogram(i).purenessR(H.getHistogram(i))
        return comp/self.nrf
    

    def combine(self,O):
        name="("+str(self.getName())+"-"+str(O.getName())+")"
        nH=""
        
        for i in range(self.nrf):
            nH+=self.getHistogram(i).join(O.getHistogram(i))+","
        
        return SymbolicObject( nH[:-1],name,self.types,A = self,B = O)
    
    def combineR(self,O):
        name="("+self.getName()+"-"+O.getName()+")"
        nH=""
        for i in range(self.nrf):
            s1=self.getHistogram(i)
            s2=O.getHistogram(i)
            nH+=s1.joinR(s2)+","
#             if(s1.type==0):
#                 nH+="{["+str(min(s1.get_nonZeroS(),s2.get_nonZeroS()))+","+str(max(s1.get_nonZeroE(),s2.get_nonZeroE()))+"]1},";
#             elif(s1.type==1):
        
        return SymbolicObject( nH[:-1],name,self.types,A = self,B = O,cluster=self.cluster)

    
    
class Histogram(object):
   
    def __init__(self, hist,i=0,e=0,type=0,nr=1,Fmin=0,Fmax=1):
        
        self.hist=hist+"}"
        self.splitHist(Fmin,Fmax)
        self.type=type
        self.i=i
        self.e=e
        self.nr=nr
        self.nrofBins=len(self.data)
        
        if(Fmin!=0 or Fmax!=1):
            self.hist="{"
            for n in range(self.nrofBins):
                self.hist+="["+str(self.get_start(n))+","+str(self.get_end(n))+"]"+str(self.get_data(n))+";"
            self.hist=self.hist[:-1]+"}"

    
    def getHist(self):
        return self.hist
    def getHistRounded(self,nrR=0):
        sttr="{"
        for i in range(self.nrofBins):
            sttr+="["+str(round(self.start[i],nrR))+","+str(round(self.end[i],nrR))+"]"+str(round(self.data[i],nrR))+";"
        sttr=sttr[:-1]+"}"
        return sttr
    def getSpan(self):
        
        if(self.type==0):
            if(self.get_nonZeroE()==self.get_nonZeroS()):
                return 0.0000000001
            return self.get_nonZeroE()-self.get_nonZeroS()
        elif(self.type==1):
            return self.getSize()
            
    def getSize(self):
        if(self.nrofBins==0):
            return 0
        elif(self.get_nonZeroS()==self.get_nonZeroE()):
            return 0
        else:
            prob_t=0;
            size=0;
            if(self.type==1):
                for i in range(self.nrofBins):
                    size=size+self.get_data(i)*(self.get_end(i)-self.get_start(i))
            else:
                for i in range(self.nrofBins):
                    prob_t=prob_t+self.get_data(i)
                    size=size+prob_t*(self.get_end(i)-self.get_start(i))
        return size

    def dissimilarity(self,H):
        return Histogram(self.join(H),self.i).getSize()-Histogram(self.meet(H),self.i).getSize()
    def purenessR(self,H):
        meet=Histogram(self.meetR(H),self.i).getSize()
        join=Histogram(self.joinR(H),self.i).getSize()
        if(join==0):
            return 0
        return meet/join
    def pureness(self,H):
        meet=Histogram(self.meet(H),self.i).getSize()
        join=Histogram(self.join(H),self.i).getSize()
        return meet/join
    def compactnessR(self):
        return self.getSpan()
    def compactness(self):
        return self.get_compactness()

    def meetR (self, H):
        if(self.type==0):
            if(max(self.get_nonZeroS(), H.get_nonZeroS())<=min(self.get_nonZeroE(), H.get_nonZeroE())):
                if((self.get_nonZeroS()==self.get_nonZeroE() or H.get_nonZeroS()==H.get_nonZeroE()) and ((self.get_nonZeroE() == H.get_nonZeroE()) or (self.get_nonZeroS()== H.get_nonZeroS()))):

                    return "{[0,0]0}"
                else:
                    return "{["+str(max(self.get_nonZeroS(), H.get_nonZeroS()))+","+str(min(self.get_nonZeroE(), H.get_nonZeroE()))+"]1}"
            else:
                return "{[0,0]0}"
        elif(self.type==1):
            sttr="{"
            for i in range(len(self.data)):
                if(self.data[i]==1 and H.data[i]==1):
                    sttr+="["+str(self.start[i])+","+str(self.end[i])+"]1;"
                else:
                    sttr+="["+str(self.start[i])+","+str(self.end[i])+"]0;"
            sttr=sttr[:-1]+"}"
            if(sttr=="}"):
                return "{[0,0]0}"
            else:
                return sttr
    def meet(self,H):
        nSt="{"
        i=0
        j=0
        start_c=min(self.get_start(i),H.get_start(j))
        if(start_c==self.get_start(i)):
            i=i+1
        if(start_c==H.get_start(j)):
            j=j+1
        while(i<=self.getNrofBins() and j<=H.getNrofBins()):
            end_c=min(self.next(i),H.next(j));
            prob=0
            if(self.type==0):
                prob=self.probability_min_uni(i,j,start_c,end_c,H)
            elif(self.type==1): # categorical value
                if(self.get_data(i-1)==1 and H.get_data(j-1)==1):
                    prob=1;
            nSt=nSt+"["+str(start_c)+","+str(end_c)+"]"+str(prob)+";"
            if(end_c==self.next(i)): 
                i=i+1
            if(end_c==H.next(j)):
                j=j+1

            start_c=end_c;
        if(i>self.getNrofBins()):
            while(j<=H.getNrofBins()):
                end_c=H.next(j)
                prob=0;
                nSt=nSt+"["+str(start_c)+","+str(end_c)+"]"+str(prob)+";"
                start_c=end_c; 
                j=j+1
        else:
            while(i<=self.getNrofBins()):
                end_c=self.next(i)
                prob=0;
                nSt=nSt+"["+str(start_c)+","+str(end_c)+"]"+str(prob)+";"
                start_c=end_c; 
                i=i+1
        return nSt[:-1]+"}";
    
    def joinR (self, H):
        if(self.type==0):
            return "{["+str(min(self.get_nonZeroS(), H.get_nonZeroS()))+","+str(max(self.get_nonZeroE(), H.get_nonZeroE()))+"]1}"
        elif(self.type==1):
            sttr="{"
            for i in range(len(self.data)):
                if(self.data[i]==1 or H.data[i]==1):
                    sttr+="["+str(self.start[i])+","+str(self.end[i])+"]1;"
                else:
                    sttr+="["+str(self.start[i])+","+str(self.end[i])+"]0;"
            sttr=sttr[:-1]+"}"
            if(sttr=="}"):
                return "{[0,0]0}"
            else:
                return sttr

    
    def join(self,H,nr1=0,nr2=0):
        nSt="{"
        i=0
        j=0
        if(nr1==0):
            nr1=self.nr
        if(nr2==0):
            nr2=H.nr
        start_c=min(self.get_start(i),H.get_start(j))
        if(start_c==self.get_start(i)):
            i=i+1
        if(start_c==H.get_start(j)):
            j=j+1
            
        while(i<=self.getNrofBins() and j<=H.getNrofBins()):
            end_c=min(self.next(i),H.next(j));
            prob=0
            if(self.type==0):
                prob=self.probability_norm_uni(i,j,start_c,end_c,H,nr1,nr2)
            elif(self.type==1): # categorical value
                if(self.get_data(i-1)==1 or H.get_data(j-1)==1):
                    prob=1;
     
            nSt=nSt+"["+str(start_c)+","+str(end_c)+"]"+str(prob)+";"
            if(end_c==self.next(i)): 
                i=i+1
            if(end_c==H.next(j)):
                j=j+1

            start_c=end_c;
        
        if(i>self.getNrofBins()):
            while(j<=H.getNrofBins()):
                end_c=H.next(j)
                prob=self.probability_norm_uni(i,j,start_c,end_c,H,nr1,nr2);
                nSt=nSt+"["+str(start_c)+","+str(end_c)+"]"+str(prob)+";"
                start_c=end_c; 
                j=j+1
        else:
            while(i<=self.getNrofBins()):
                end_c=self.next(i)
                prob=self.probability_norm_uni(i,j,start_c,end_c,H,nr1,nr2);
                nSt=nSt+"["+str(start_c)+","+str(end_c)+"]"+str(prob)+";"
                start_c=end_c; 
                i=i+1
        return nSt[:-1]+"}";
    
    def probability_min_uni(self,i,j,start_c,end_c,H): 
        if(i-1<0 or i>self.getNrofBins()):
            return 0;
        elif(j-1<0 or j>H.getNrofBins()):
             return 0  
        else:
            if(i-1<0):
                return 0
            elif(j-1<0):
                return 0  
            elif(end_c==start_c):
                return 0      
            elif(self.get_end(i-1)>=end_c and self.get_start(i-1)<=start_c and H.get_end(j-1)>=end_c and H.get_start(j-1)<=start_c):
                return min(self.get_data(i-1)*(end_c-start_c)/(self.get_end(i-1)-self.get_start(i-1)),
                        H.get_data(j-1)*(end_c-start_c)/(H.get_end(j-1)-H.get_start(j-1)))
            else:
                return 0;        
    def probability_norm_uni(self,i,j,start_c,end_c,H,nr1,nr2):

        if(i-1<0 or i>self.getNrofBins()): #first has no bins
            if(j-1<0 or j>H.getNrofBins()): #second has no bins
                return 0
            else:
                if(H.get_end(j-1)>=end_c and H.get_start(j-1)<=start_c): #second has a bin
                    if(end_c==start_c):
                        return (nr2*H.get_data(j-1))/(nr1+nr2)
                    else:
                        return ((end_c-start_c)*nr2*H.get_data(j-1)/(H.get_end(j-1)-H.get_start(j-1)))/(nr1+nr2);
                else:
                    return 0
        elif(j-1<0 or j>H.getNrofBins()): #second has no bins
            if(i-1<0 or i>self.getNrofBins()): #first has no bins
                return 0
            else:
                if(self.get_end(i-1)>=end_c and self.get_start(i-1)<=start_c): #second has a bin
                    if(end_c==start_c):
                        return (nr1*self.get_data(i-1))/(nr1+nr2)
                    else:
                        return ((end_c-start_c)*nr1*self.get_data(i-1)/(self.get_end(i-1)-self.get_start(i-1)))/(nr1+nr2);
                else:
                    return 0
        else:
            if(end_c==start_c):
                return 1
            elif(self.get_end(i-1)>=end_c and self.get_start(i-1)<=start_c and H.get_end(j-1)>=end_c and H.get_start(j-1)<=start_c):
#                 print((end_c-start_c)*self.nr*self.get_data(i-1)/(self.get_end(i-1)-self.get_start(i-1)))
#                 print((end_c-start_c)*H.nr*H.get_data(j-1)/(H.get_end(j-1)-H.get_start(j-1)))
                return ((end_c-start_c)*nr1*self.get_data(i-1)/(self.get_end(i-1)-self.get_start(i-1))+(end_c-start_c)*nr2*H.get_data(j-1)/(H.get_end(j-1)-H.get_start(j-1)))/(nr1+nr2);
            elif(self.get_end(i-1)>=end_c and self.get_start(i-1)<=start_c):
                return ((end_c-start_c)*nr1*self.get_data(i-1)/(self.get_end(i-1)-self.get_start(i-1)))/(nr1+nr2);
            elif(H.get_end(j-1)>=end_c and  H.get_start(j-1)<=start_c):
                return ((end_c-start_c)*nr2*H.get_data(j-1)/(H.get_end(j-1)-H.get_start(j-1)))/(nr1+nr2);
            else:
                return 0;
    
    def next(self,index):        
        if(index<self.getNrofBins()): #valid index
            return self.get_start(index)
        elif(index==self.getNrofBins()):
            return self.get_end(self.getNrofBins()-1)
        else:
            return -1

    def splitHist(self,Fmin,Fmax):

        self.data=[]
        self.start=[]
        self.end=[]
        abi=self.hist.replace("{", "").replace("}", "").replace("[","").split(";")
        conf=Fmax-Fmin
        for n in range(len(abi)):
            b=abi[n].split("]")
            c=b[0].split(",")
            self.data.append(float(b[1]))
            if(Fmin==0 and Fmax==1):
                self.start.append(float(c[0]))
                self.end.append(float(c[1]))
            else:
                self.start.append(max(0,(float(c[0])-Fmin)/conf))
                self.end.append(min(1,(float(c[1])-Fmin)/conf))
            
    def getNrofBins(self):
        return self.nrofBins

    def get_nonZeroS(self):
        for i in range(self.nrofBins):
            if(self.data[i]>0):
                return self.start[i]
        return 0;
        
    def get_nonZeroE(self):
        for i in reversed(range(self.nrofBins)):
            if(self.data[i]>0):
                    return self.end[i]
        return 0;
    
    def get_data(self,j):
        return self.data[j]
    
    def get_start(self,j):
        return self.start[j]
    
    def get_end(self,j):
        return self.end[j]
    
    def get_type(self):
        return self.type
    
    def get_median(self,type=None):
        if(type==None):
            type=self.type
        if(type==1):
            return ((self.get_nonZeroE()+self.get_nonZeroS())/2)
        else:
            total=0
            for i in range(self.nrofBins):
                total+=self.get_data(i)
            if(total==0):
                return 0
            return self.get_quantile(total/2)
    def get_compactness(self,feature_len=1):
        so_far=0
        for i in range(self.nrofBins):
            so_far+=self.get_data(i)*(self.get_end(i)-self.get_start(i))
        return so_far/feature_len
    
    def get_quantile(self,perc):
        so_far=0
        tt=0
        if(perc==0):
            return self.get_nonZeroS();
        if(perc==1):
            return self.get_nonZeroE();

        for i in range(self.nrofBins):
            tt=tt+self.get_data(i)

        for i in range(self.nrofBins):
            if so_far+(self.get_data(i)/tt)<perc:
                so_far=so_far+(self.get_data(i)/tt)
            else:
                return (perc-so_far)*(self.get_end(i)-self.get_start(i))/(self.get_data(i)/tt)+self.get_start(i)
        return 0
    def get_cummulative(self,x):
        so_far=0
        i=0
        while self.get_start(i)<x:
            if(self.get_end(i)<x):
                so_far+=self.get_data(i)
            else:
                so_far+=self.get_data(i)*(x-self.get_start(i))/(self.get_end(i)-self.get_start(i))
            i=i+1
        return so_far
        
   
    