#!/usr/bin/python3

import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys                

def displayFile(fname):
    searchObj = re.search( r'(.*)_(.*)x(.*)_T(.*).sav', fname, 0)
    name = searchObj.group(1) 
    M = int(searchObj.group(2))
    N = int(searchObj.group(3))
    T = int(searchObj.group(4))
    # lecture binaire suivant l'ordre de stockage des
    # éléments d'un tableau multi-dimensionnel en C : 
    x = np.fromfile(fname,np.float64,T*N*N,"")
    x.shape = (T,N,N)
    #np.set_printoptions(precision=15)
    #print(x[0,:,:])
    fig = plt.figure()
    plt.xlabel('x')
    plt.xlim(0, N)
    plt.ylabel('y')
    plt.ylim(0, M)
    plt.gca().invert_yaxis()
    plt.title(fname)
    #ax = fig.add_subplot(111)
    ims = []
    for t in range(T):
        ###
        #t = ax.annotate("toto",(1,1)) # add text
        #ims.append((plt.pcolormesh(x[t,:,:]),))

        ### pcolormesh pour grands tableaux 
        ims.append((plt.pcolormesh(x[t,:,:]
                                   , norm=plt.Normalize(-5, 15)
#                                   , norm=plt.Normalize(-5, 5)
                                   ),))
       

    im_ani = animation.ArtistAnimation(fig, ims, interval=1,
                                       #repeat_delay=3000,
                                       #blit=True
                                       repeat=False )
    plt.show()
    
##    for t in range(T):
##        plt.title('t={}'.format(t))
##        plt.imshow(x[t,:,:])
##        plt.gca().invert_yaxis()
##    plt.show()



if __name__ == "__main__":
    if (len(sys.argv) != 2):
        print("Usage example: ./visu.py shalw_1024x1024_T20.sav ")
    else:
        displayFile(sys.argv[1])



