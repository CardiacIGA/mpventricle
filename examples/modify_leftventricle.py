# -*- coding: utf-8 -*-
"""
Created on Fri May 20 21:51:17 2022

@author: s146407
"""

import numpy, random, os
from vtnurbs.construct import MultipatchSolid
from splipy import BSplineBasis, Surface
import matplotlib.pyplot as plt
from utils.save_load_data import load_pickle

directC       = os.path.realpath(os.path.dirname(__file__)) # Current directory    
folder_pickle = 'output/pickle/'
folder_json   = 'output/json/'
folder_txt    = 'output/txt/'
folder_vtk    = 'output/vtk/'

folder_vtk    = os.path.join(directC, folder_vtk)
folder_pickle = os.path.join(directC, folder_pickle)


## Modify existing multipatch Left-ventricle geometry      
surfaceLV_outer = load_pickle(folder_pickle+'LV_left_outer')
surfaceLV_inner = load_pickle(folder_pickle+'LV_left_inner')

# Surface to be modified
surface = surfaceLV_outer
nrefine = 3
surface.refine(nrefine)
surfaceLV_inner.refine(nrefine)

# Example to be run
RunExample = 4

###########################################################
## Example 1: Modify cps of the first patch (patchID=0)  ##
###########################################################
if RunExample == 1:
    patchID     = 1
    cps_patchID = surface.cps(patchID)       # Select the cps of that patch     
    cps_patchID *= 1.2                       # Multiply some cps by a random factor
    w_patchID   = surface.weights(patchID)   # Select the weights of that patch     
    w_patchID   *= 1.2                       # Multiply some weights by a random factor
    ## Note! Adding knots is not yet supported

    # After modifying, update the entire datastructure with its patch-interfacedependencies
    order = 0,1,2,3,4
    surface.update(*order)


###########################################################
## Example 2: Modify cps of all patches randomly (a)     ##
###########################################################
elif RunExample == 2:
    patchID     = 1
    for ipatch in range(surface.nrpatches()):
        random1 = random.uniform(0.8,1.2) 
        random2 = random.uniform(0.9,1.1)
        
        cps_ipatch = surface.cps(ipatch)
        w_ipatch   = surface.weights(ipatch)
        cps_ipatch *=random1       # Select the cps of that patch     
        w_ipatch   *=random2   # Select the weights of that patch     
        

    # After modifying, update the entire datastructure with its patch-interfacedependencies
    order = 0,1,2,3,4
    surface.update(*order)


###########################################################
## Example 3: Modify cps of all patches randomly (b)     ##
###########################################################
elif RunExample == 3:
    CPS = surface._cps
    W   = surface._w

    cps_mod = CPS.copy()*10
    w_mod   = W.copy()*0.9

    surface._cps = cps_mod
    surface._w   = w_mod 

    cps_ipatch    = surface.cps(0)
    cps_ipatch[2] = numpy.array([ 1, 0, 0.2 ])
    # After modifying, update the entire datastructure with its patch-interfacedependencies
    order = 0,1,2,3,4
    surface.update(*order)


##########################################
## Example 4: Convert a patch to Splipy ##
##########################################
elif RunExample == 4:
    patchID = 1
    cps     = surface.cps(patchID)
    w       = surface.weights(patchID)
    knotvec = surface.knotvec(patchID) # <-- This cannot be modified, only knotval() and knotmult()

    control_net    = numpy.concatenate([ cps*w[:,numpy.newaxis], w[:,numpy.newaxis] ], axis = 1) # Homogenous coordinates
    spline_order   = 2
    basis_u        = BSplineBasis(order=spline_order+1, knots=knotvec[0])
    basis_v        = BSplineBasis(order=spline_order+1, knots=knotvec[1])
    Splipy_surface = Surface(basis_u, basis_v, control_net, rational=True)

    ###################
    # Post-processing #
    ###################

    def plot_surface(*surface,control_net=False):
        u    = numpy.linspace(0,1,20) # 150 visualization points on our parametric domain [0,4]
        v    = numpy.linspace(0,1,20) # 150 visualization points on our parametric domain [0,4]
        
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(projection='3d')
        for isurf in surface:
            surf = isurf.evaluate(u,v)    # compute (x,y)-coordinates of the curve evaluation
            ax.plot_wireframe(surf[:,:,0], surf[:,:,1], surf[:,:,2])
            if control_net:
                ax.plot_wireframe(isurf[:,:,0]/isurf[:,:,3], isurf[:,:,1]/isurf[:,:,3], isurf[:,:,2]/isurf[:,:,3], linestyles='--', color='r') 
                ax.plot(isurf[:,:,0].ravel()/isurf[:,:,3].ravel(), isurf[:,:,1].ravel()/isurf[:,:,3].ravel(), isurf[:,:,2].ravel()/isurf[:,:,3].ravel(), 'bs')
        ax.view_init(20, 45)
        ax.set_title('Surface')
        plt.show()
    
    plot_surface(Splipy_surface, control_net=True)   # <- Insert multiple surfaces to plot in same plot


if RunExample !=4:
    ## Construct the solid LV 
    solid   = MultipatchSolid(surfaceLV_outer, surfaceLV_inner)
    solid.save_vtk(folder_vtk+'Modified_LeftventricleGeom', nrsamples=8)





