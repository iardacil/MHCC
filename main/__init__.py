'''
Main Class
Two datasets Oils and Hardwood are included
Two methods - classical HCC for symbolic data and MHCC can be called
Dendograms can be drawn according to dissimilarity matrixes produced by clustering

@author: Kadri Umbleja
'''
import numpy as np
from main.Common import Histogram,SymbolicObject,normalize
from main.ClassicalHCC import HCC
from main.MicroHCC import MHCC
from main.dendogram import drawDendowithDiss


def callMethod(X,objects,titles,type):
    nre= len(X)
    nrf = len(titles)
    
    quantiles=[0,0.1,0.25,0.5,0.75,0.9,1]
    nrq=len(quantiles)
    
    #transform data to symbolic object, normalize
    aList2,Fmin,Fmax=normalize(X,objects,titles,type)
    
    diss=MHCC(aList2,nre,nrf,objects,quantiles)
    drawDendowithDiss(diss,objects)
    
    diss=HCC(aList2,nrf,nre,1,mode=4)
    drawDendowithDiss(diss,objects)

if __name__ == '__main__':
    
#     objects=['Linsead','Perilla','Cotton','Sesame','Camellia','Olive','Beef','Hog']
#     titles = ['Specific Gravity','Freezing Point','Iodine Value','Saponification value','Major Fatty Acids']
#     type=[0,0,0,0,0]
#     X = np.loadtxt("data/oils.txt", dtype='str', delimiter='\n')


    objects=['ACER_EAST','ACER_WEST','ALNUS_EAST','ALNUS_WEST','FRAXINUS_EAST','FRAXINUS_WEST','JUGLANS_EAST','JUGLANS_WEST','QUERCUS_EAST','QUERCUS_WEST']
    titles=['ANNT','JANT','JULT','ANNP','JANP','JULP','GDC5','MITM']
    type=[0,0,0,0,0,0,0,0,0]
    X = np.loadtxt("data/hardwood.txt", dtype='str', delimiter='\n')
    
    
    callMethod(X,objects,titles,type)