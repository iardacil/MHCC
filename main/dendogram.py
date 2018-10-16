'''
Class for drawing dendrograms based on dissimilarity matrix
@author: Kadri Umbleja
'''

from main.ClassicalHCC import HCC
import scipy as sp
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

#class that draws dendograms based on HCC
#mode defines metric 0-histogram, 4-interval
def drawDendo(aList,nre,nrf,objects,titles,type,mode=4,till=1,l=[]):
    

    diss,extra=HCC(aList,nrf,nre,till,mode=mode,Fmin=[],Fmax=[],l=l)

    drawDendowithDiss(diss,objects)

def drawDendowithDiss(diss,objects):
    

    diss=sp.array(diss)
    linkage_matrix = linkage(squareform(diss), "complete")
    dendrogram(linkage_matrix, labels=objects,leaf_rotation=90.,leaf_font_size=8.,color_threshold=0, above_threshold_color='black')
    plt.subplots_adjust(bottom=0.3)
    plt.savefig('graphs/dendo.png', format='png', dpi=300)
    plt.show()