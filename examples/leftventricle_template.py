# -*- coding: utf-8 -*-
"""
Created on Tue May 24 14:43:03 2022

@author: s146407
"""

import numpy, os
from vtnurbs.construct import CardiacGeometry, MultipatchSolid
from utils.save_load_data import load_pickle, save_pickle

    
directC       = os.path.realpath(os.path.dirname(__file__)) # Current directory    
folder_pickle = 'output/pickle/'
folder_json   = 'output/json/'
folder_txt    = 'output/txt/'
folder_vtk    = 'output/vtk/'
folder_stl    = 'output/stl/'

folder_vtk    = os.path.join(directC, folder_vtk)
folder_pickle = os.path.join(directC, folder_pickle)
folder_json   = os.path.join(directC, folder_json)
folder_txt    = os.path.join(directC, folder_txt)
folder_stl    = os.path.join(directC, folder_stl)




## Left-ventricle geometry example
mm=1e-3;ml=1e-6  
# These volumetric values is what we are aiming for when we have the actual NURBS geometry
Vpap = 4.0 *ml
Vlv0 = 40.0*ml
Vw   = 140.*ml
Vi   = Vlv0 + Vpap
Ve   = Vw - Vpap

# According to literature ():
C  = 43.0*mm # Focal length
H  = 24.0*mm # Truncation height
ξi = 0.37129680875745236 #0.371296808757 #0.3694447442932904 # Transmural value endocardium
ξo = 0.6783556518284651  #0.678355651828 #0.6754874350057226 # Transmural value epicardium  
Rleft_inner  = C*numpy.array([numpy.sinh(ξi),numpy.sinh(ξi),numpy.cosh(ξi)])
Rleft_outer  = C*numpy.array([numpy.sinh(ξo),numpy.sinh(ξo),numpy.cosh(ξo)])


LV_outer = {'Rx': Rleft_outer[0]   , 'Ry': Rleft_outer[1]  , 'Rz': Rleft_outer[2], 'H': H}
LV_inner = {'Rx': Rleft_inner[0]   , 'Ry': Rleft_inner[1]  , 'Rz': Rleft_inner[2]} # H can be given, but not required if it is in LV_inner

# Constructs the left-ventricle inner and outer surfaces
leftventricle  = CardiacGeometry('left-ventricle', LVi=LV_inner, LVo=LV_outer) 

# Generate the surfaces based on the provided input
surfaceLV_outer = leftventricle.generate_surface('left-outer') 
surfaceLV_inner = leftventricle.generate_surface('left-inner')

## Save/load surface data (does not use the provided input), useful when you want to play with the cps and not run the template algorithm
# save_pickle(surfaceLV_outer, filename=folder_pickle+'LV_left_outer')
# save_pickle(surfaceLV_inner, filename=folder_pickle+'LV_left_inner')
# surfaceLV_outer = load_pickle(folder_pickle+'LV_left_outer')
# surfaceLV_inner = load_pickle(folder_pickle+'LV_left_inner')

## Generate the solid
solid = MultipatchSolid(surfaceLV_outer, surfaceLV_inner)
solid.get_volume()

## Postprocessing
#solid.save_vtk(folder_vtk+'LeftventricleGeom', boundary_only=True, nrsamples=50)#, boundary=True)
#solid.save_stl(folder_stl+'LeftventricleGeom', nrsamples=60)#, boundary=True)
#solid.save_pickle(folder_pickle+'LV_GEOMETRY_DATA',boundary=solid.boundary_names()) # <- This data is used for Nutils simulations
#solid.save_json(folder_json+'LV_GEOMETRY_DATA', boundary=solid.boundary_names())
#solid.save_txt(folder_txt+'LV_GEOMETRY_DATA', boundary=solid.boundary_names())
#surfaceLV_inner.save_controlnet(folder_vtk+'Inner surface controlnet')
#surfaceLV_outer.save_controlnet(folder_vtk+'Outer surface controlnet')

# Reload when desired
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_pickle(folder_pickle+'LV_GEOMETRY_DATA')
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_json(folder_json+'LV_GEOMETRY_DATA')
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_txt(folder_txt+'LV_GEOMETRY_DATA')