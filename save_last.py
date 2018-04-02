#!/usr/bin/python3

import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys                

def saveLastStep(fname):
    searchObj = re.search( r'(.*)_(.*)x(.*)_T(.*).sav', fname, 0)
    name = searchObj.group(1) 
    M = int(searchObj.group(2))
    N = int(searchObj.group(3))
    T = int(searchObj.group(4))
    # lecture binaire suivant l'ordre de stockage des
    # éléments d'un tableau multi-dimensionnel en C : 
    x = np.fromfile(fname,np.float64,T*N*N,"")
    x.shape = (T,N,N)
    ims = []
    ims.append((plt.pcolormesh(x[T-1,:,:]
                                   , norm=plt.Normalize(-5, 15)
                                   ),))
    plt.subplots_adjust(bottom=0, top=1, left=0, right=1)
    plt.axis('off')
    plt.savefig(name + "_" + str(M) + "x" + str(N) + "_T" + str(T) + ".png")
                                   

if __name__ == "__main__":
    if (len(sys.argv) != 2):
        print("Usage example: ./save_last.py shalw_1024x1024_T20.sav ")
    else:
        saveLastStep(sys.argv[1])



