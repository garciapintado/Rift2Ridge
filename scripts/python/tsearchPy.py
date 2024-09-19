# main for stand-alone python tsearch()

# env = sio.loadmat(os.path.join(dsn,fname), squeeze_me=True)  

import matplotlib.tri as tri
import os
import scipy.io as sio
import sys

#HOME  = os.environ['HOME']                                              # just for when rifutils is more evolved - now getTri(),tsearch() embedded here
#pypath = os.path.join(HOME,"docs","ingit","riftutils")
#sys.path.append(os.path.join(pypath,"utils"))

def getTri(env):
    """
       get triangular mesh as <class 'matplotlib.tri.triangulation.Triangulation'>
    """
    nvertx = 3
    EL2NOD_vert = env['ELEM2NODE'][0:nvertx,:]
    nel = EL2NOD_vert.shape[1]
    nvnodes = EL2NOD_vert.max()
    GCOORDV = env['GCOORD'][:,0:nvnodes]
    EL2NOD_mesh = EL2NOD_vert.transpose() - 1 # .shape = (nel, 3) :: -1 for phyton indexing
    trimesh = tri.Triangulation(GCOORDV[0,:], GCOORDV[1,:], triangles=EL2NOD_mesh)  # <class 'matplotlib.tri.triangulation.Triangulation'> ::
    meshana = tri.TriAnalyzer(trimesh)                                              #  .x.size==nvnodes, .y.size=nvnodes, .triangles.shape=(nel,3)
    meshmsk = meshana.get_flat_tri_mask()
    trimesh.set_mask(meshmsk)
    return trimesh

def tsearch(env, x, y):
    """
       locates elements in triangular mesh surrounding a set of locations
       x :: [m], x-coordinates, where m is the number of requested samples
       y :: [m], y-coordinates, where m is the number of requested samples
       example:

       import time
       xlim = [env['GCOORD'][0,:].min(),env['GCOORD'][0,:].max()]
       ylim = [env['GCOORD'][1,:].min(),env['GCOORD'][1,:].max()]
       m = 10000
       xo = xlim[0] + np.random.rand(m) * (xlim[1] - xlim[0])
       yo = ylim[0] + np.random.rand(m) * (ylim[1] - ylim[0])
       tic = time.time(); xyo_el = tsearch(env, xo,yo); print(time.time() - tic)  
    """
    trimesh = getTri(env)                                                           # <class 'matplotlib.tri.triangulation.Triangulation'>
    trifinder = trimesh.get_trifinder()                                             #  <class 'matplotlib.tri.trifinder.TrapezoidMapTriFinder'>
    return trifinder(x,y)

def main():
    fnamei = sys.argv[1]                                       #fnamei = "tsearchPy_01711.mat"
    verify = bool(int(sys.argv[2]))
    env = sio.loadmat(fnamei, verify_compressed_data_integrity=verify)
    x = env['xy'][0,:]
    y = env['xy'][1,:]
    
    xyel = tsearch(env, x, y)                                                       # xyel.dtype == dtype('int32') | assumes {"ELEM2NODE","GCOORD"} names for triangulation
    fnameo = fnamei.replace("_i","_o")
    
    #xyel.tofile(fnameo)
    #f = open(fnameo,"wb")
    #f.write(xyel.tobytes())
    #f.close()
    sio.savemat(fnameo,{'xyel':xyel})

if __name__ == "__main__":
    main()
