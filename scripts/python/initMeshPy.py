# main for initMeshPy()

# import importlib
import matplotlib.tri as tri
import numpy as np
import os
import scipy.io as sio
import sys

HOME  = os.environ['HOME']                                              # just for when rifutils is more evolved - now getTri(),tsearch() embedded here
pypath = os.path.join(HOME,"docs","ingit","riftutils")
sys.path.append(os.path.join(pypath,"utils"))

import riftutils

# importlib.reload(riftutils)

asMilamin = True # false for kinedyn node indexing

# from riftutils import getTri

def main():
    fnamei = sys.argv[1]                                       #fnamei = "initFluidPy_06242_i.mat"
    verify = bool(int(sys.argv[2]))                            #verify=True
    env = sio.loadmat(fnamei, verify_compressed_data_integrity=verify)

    hydroelboo = env['HMESH']['hydroelboo'][0][0].flatten() # uint8 [nel]
    nel        = env['HMESH']['nel'][0][0].flatten()        # uint16     '<u2'
    GCOORD     = env['HMESH']['HCOORD'][0][0]               # [2,nnod7]
    EL2NOD     = env['HMESH']['EL2HNO'][0][0]               # int32 [6,nel] '<u4')
    
    HCObdi3, HCObdi6 = riftutils.getBoundaryNodes(EL2NOD, GCOORD, asMilamin=asMilamin)          # fig, ax = plt.subplots(); ax.triplot(trimesh, lw=0.1); fig.show()
                                                                                                # ax.plot(GCOORD[0,GCObdi3], GCOORD[1, GCObdi3], 'o', color='blue', markersize=10)
                                                                                                # ax.plot(GCOORD[0,GCObdi6], GCOORD[1, GCObdi6], 'o', color='orange', markersize=5)
    fnameo = fnamei.replace("_i","_o")
    sio.savemat(fnameo,{'HCObdi3':HCObdi3+1,'HCObdi6':HCObdi6+1})                               # map back into matlab indices & export

if __name__ == "__main__":
    main()
