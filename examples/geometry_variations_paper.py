# -*- coding: utf-8 -*-
"""
Created on Tue May 24 14:43:03 2022

@author: s146407
"""

import numpy, os
from vtnurbs.construct import CardiacGeometry, MultipatchSolid
from utils.save_load_data import load_pickle, save_pickle    

directC       = os.path.realpath(os.path.dirname(__file__)) # Current directory    
folder_pickle = 'output (geom variations)/pickle/'
folder_json   = 'output (geom variations)/json/'
folder_txt    = 'output (geom variations)/txt/'
folder_vtk    = 'output (geom variations)/vtk/'
folder_vtk    = os.path.join(directC, folder_vtk)
folder_pickle = os.path.join(directC, folder_pickle)
folder_json   = os.path.join(directC, folder_json)
folder_txt    = os.path.join(directC, folder_txt)


## Bi-ventricle geometry  
# Three variations have been created: REF (reference), THK (thick), and LONG (elongated)
varType = "REF" #"REF", "LONG", "THK"

mm=1e-3;ml=1e-6  
name = varType
C    = 43.0*mm
H    = 24.0*mm

# Settings
paramREF  = {'ξi':0.371296808757, 'ξo':0.678355651828,
             'fx':0.9789, 'fy':0.659, 'dZ': 0.005, 'dRright': 0.005 ,
             'fL': 1, 'fR': 1}
paramLONG = {'ξi':0.371296808757, 'ξo':0.678355651828,
             'fx':0.9789, 'fy':0.659, 'dZ': 0.005, 'dRright': 0.005,
             'fL': 1.2, 'fR': 0.981}
paramTHK  = {'ξi':0.371296808757, 'ξo':0.678355651828,
             'fx':0.6, 'fy':0.55, 'dZ': 0.01, 'dRright': 0.008 ,
             'fL': 1, 'fR': 1}
param    = {'REF':paramREF,'LONG':paramLONG,'THK':paramTHK}


## According to literature (Same as Benchmark FEniCS):
ξi = param[name]['ξi'] 
ξo = param[name]['ξo']     
fL = param[name]['fL']
fR = param[name]['fR']
H  = H*fL # scale the truncated height

Rleft_inner  = C*numpy.array([numpy.sinh(ξi)/(fL**0.5),numpy.sinh(ξi)/(fL**0.5),fL*numpy.cosh(ξi)])
Rleft_outer  = C*numpy.array([numpy.sinh(ξo)/(fL**0.5),numpy.sinh(ξo)/(fL**0.5),fL*numpy.cosh(ξo)])

Xcenter       = -0.792*Rleft_outer[0] if name == 'THK' else -Rleft_inner[0]
RV_wall       = param[name]['dRright'] 
dZ_right      = param[name]['dZ'] 
Rright_inner  = numpy.array([param[name]['fx']*Rleft_outer[0]*fR  + RV_wall,
                             param[name]['fy']*Rleft_outer[1]*fR  + RV_wall,
                             abs(Rleft_outer[2]*numpy.cos( 0.85*numpy.pi ))*fR - dZ_right ])
Rright_outer  = Rright_inner + RV_wall
# Rright_outer  = numpy.array([C*numpy.sinh(ξo) + RV_wall,0.7*C*numpy.sinh(ξo)+RV_wall,abs(Rleft_outer[-1]*numpy.cos( 0.85*numpy.pi )) - dZ_right ])
# Rright_inner  = Rright_outer - RV_wall

LV_outer = {'Rx': Rleft_outer[0] , 'Ry': Rleft_outer[1] , 'Rz': Rleft_outer[2] , 'H': H, 'Rxs': 0.8}
LV_inner = {'Rx': Rleft_inner[0] , 'Ry': Rleft_inner[1] , 'Rz': Rleft_inner[2] , 'Rxs': 0.8} # H can be given, but not required if it is in LV_inner 
RV_inner = {'Rx': Rright_inner[0], 'Ry': Rright_inner[1], 'Rz': Rright_inner[2]}
RV_outer = {'Rx': Rright_outer[0], 'Ry': Rright_outer[1], 'Rz': Rright_outer[2], 'Origin': numpy.array([Xcenter,0,0])} # H can be given, but not required if it is in LV_inner

##----------------------------------------------------------------------------

## Constructs the bi-ventricle inner and outer surfaces data structures
biventricle    = CardiacGeometry('bi-ventricle',   LVi=LV_inner, LVo=LV_outer, RVi=RV_inner, RVo=RV_outer) # Constructs the bi-ventricle inner and outer surfaces    

## Generate the surfaces
surfaceLV_outer = biventricle.generate_surface('left-outer')
surfaceLV_inner = biventricle.generate_surface('left-inner')
surfaceRV_outer = biventricle.generate_surface('right-outer')
surfaceRV_inner = biventricle.generate_surface('right-inner')

## Save/load data, useful when you want to play with the cps and not run the template algorithm
save_pickle(surfaceLV_outer, filename=folder_pickle+f'BV_left_{name}_outer')
save_pickle(surfaceLV_inner, filename=folder_pickle+f'BV_left_{name}_inner')
save_pickle(surfaceRV_outer, filename=folder_pickle+f'BV_right_{name}_outer')
save_pickle(surfaceRV_inner, filename=folder_pickle+f'BV_right_{name}_inner')
# surfaceLV_outer = load_pickle(folder_pickle+f'BV_left_{name}_outer')
# surfaceLV_inner = load_pickle(folder_pickle+f'BV_left_{name}_inner')
# surfaceRV_outer = load_pickle(folder_pickle+f'BV_right_{name}_outer')
# surfaceRV_inner = load_pickle(folder_pickle+f'BV_right_{name}_inner')

# Construct the solid
solid = MultipatchSolid(surfaceLV_outer, surfaceLV_inner, surfaceRV_outer, surfaceRV_inner)

# Postprocessing
solid.get_volume(shift = H, wall=True) # Determine the cavity volume
#solid.save_stl(folder_vtk+'BiventricleGeom_3Dprint', boundary_names=solid.boundary_names(), nrsamples=35)
solid.save_vtk(folder_vtk+f'BiventricleGeom_paper_{name}', boundary_only = True, nrsamples = 25)#, patch=9)#, boundary=True)
solid.save_pickle(folder_pickle+f'BV_GEOMETRY_DATA_PAPER_{name}', boundary=solid.boundary_names())
#solid.save_json(folder_json+'BV_GEOMETRY_DATA', boundary=solid.boundary_names())
#solid.save_txt(folder_txt+'BV_GEOMETRY_DATA', boundary=solid.boundary_names())
# nrefine=1
# solid.save_mesh_vtk(folder_vtk+f'BV_mesh_lines_{name}_nrefine{nrefine}', nrefine=nrefine, nsample = 10)


# Reload when desired
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_pickle(folder_pickle+'LV_GEOMETRY_DATA')
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_json(folder_json+'LV_GEOMETRY_DATA')
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_txt(folder_txt+'LV_GEOMETRY_DATA')


## Compute epi- and endocard surface area
topo, geom = solid.get_topo_geom()
from nutils import function
ns   = function.Namespace()
ns.x = geom
Aepi  = topo.boundary['epi'].integrate('d:x'@ns, degree=10)*1e4
Aendo = topo.boundary['endo_l'].integrate('d:x'@ns, degree=10)*1e4
print("Epicardium surface area:         {:.2f} [cm2]".format(Aepi))
print("Endocardium (left) surface area: {:.2f} [cm2]".format(Aendo))
##-----------------------------------
