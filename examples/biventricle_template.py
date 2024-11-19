import numpy, os
from vtnurbs.construct import CardiacGeometry, MultipatchSolid
from utils.save_load_data import load_pickle, save_pickle

    
directC       = os.path.realpath(os.path.dirname(__file__)) # Current directory    
folder_pickle = 'output/pickle/'
folder_json   = 'output/json/'
folder_txt    = 'output/txt/'
folder_vtk    = 'output/vtk/'

folder_vtk    = os.path.join(directC, folder_vtk)
folder_pickle = os.path.join(directC, folder_pickle)
folder_json   = os.path.join(directC, folder_json)
folder_txt    = os.path.join(directC, folder_txt)


## Bi-ventricle geometry  
mm=1e-3;ml=1e-6  
# Geometrical input
C  = 43.0*mm # Focal length
H  = 24.8*mm # Truncation height
ξi = 0.3694447442932904 # Transmural value endocardium
ξo = 0.6754874350057226 # Transmural value epicardium 

dr    = 1 #1.1
dr2   = 1 #1.1 
#Rlong = numpy.array([ 1/dr, 1/dr, dr**2 ])

Rleft_inner  = C*numpy.array([numpy.sinh(ξi),numpy.sinh(ξi),numpy.cosh(ξi)])
Rleft_outer  = C*numpy.array([numpy.sinh(ξo),numpy.sinh(ξo),numpy.cosh(ξo)])

RV_wall       = 0.005 
dZ_right      = 0.005
Rright_inner  = numpy.array([C*numpy.sinh(ξo) + RV_wall,0.7*C*numpy.sinh(ξo)+RV_wall,abs(C*numpy.cosh(ξo)*numpy.cos( 0.85*numpy.pi )) - RV_wall ])
Rright_outer  = Rright_inner + RV_wall
# Rright_outer  = numpy.array([C*numpy.sinh(ξo) + RV_wall,0.7*C*numpy.sinh(ξo)+RV_wall,abs(Rleft_outer[-1]*numpy.cos( 0.85*numpy.pi )) - dZ_right ])
# Rright_inner  = Rright_outer - RV_wall

LV_outer = {'Rx': C*numpy.sinh(ξo)/dr, 'Rxs': 0.8, 'Ry':  C*numpy.sinh(ξo)/dr, 'Rz': C*numpy.cosh(ξo)*dr**2, 'H': H}
LV_inner = {'Rx': C*numpy.sinh(ξi)/dr, 'Rxs': 0.8, 'Ry':  C*numpy.sinh(ξi)/dr, 'Rz': C*numpy.cosh(ξi)*dr**2} # H can be given, but not required if it is in LV_inner 
RV_inner = {'Rx': Rright_inner[0]/dr2, 'Ry': Rright_inner[1]/dr2, 'Rz': Rright_inner[2]*dr**2}
RV_outer = {'Rx': Rright_outer[0]/dr2, 'Ry': Rright_outer[1]/dr2, 'Rz': Rright_outer[2]*dr**2, 'Origin': numpy.array([-Rleft_inner[0],0.,0])} # H can be given, but not required if it is in LV_inner

##----------------------------------------------------------------------------

## Constructs the bi-ventricle inner and outer surfaces data structures
biventricle    = CardiacGeometry('bi-ventricle',   LVi=LV_inner, LVo=LV_outer, RVi=RV_inner, RVo=RV_outer) # Constructs the bi-ventricle inner and outer surfaces    

## Generate the surfaces based on the provided input
surfaceLV_outer = biventricle.generate_surface('left-outer')
surfaceLV_inner = biventricle.generate_surface('left-inner')
surfaceRV_outer = biventricle.generate_surface('right-outer')
surfaceRV_inner = biventricle.generate_surface('right-inner')

## Save/load surface data (does not use the provided input), useful when you want to play with the cps and not run the template algorithm
# save_pickle(surfaceLV_outer, filename=folder_pickle+'BV_left_outer')
# save_pickle(surfaceLV_inner, filename=folder_pickle+'BV_left_inner')
# save_pickle(surfaceRV_outer, filename=folder_pickle+'BV_right_outer')
# save_pickle(surfaceRV_inner, filename=folder_pickle+'BV_right_inner')
# surfaceLV_outer = load_pickle(folder_pickle+'BV_left_outer')
# surfaceLV_inner = load_pickle(folder_pickle+'BV_left_inner')
# surfaceRV_outer = load_pickle(folder_pickle+'BV_right_outer')
# surfaceRV_inner = load_pickle(folder_pickle+'BV_right_inner')

## Construct the solid
solid = MultipatchSolid(surfaceLV_outer, surfaceLV_inner, surfaceRV_outer, surfaceRV_inner)

## Postprocessing
solid.get_volume(shift = H) # Determine the cavity volume
#solid.save_stl(folder_vtk+'BiventricleGeom', boundary_names=solid.boundary_names())
solid.save_vtk(folder_vtk+'BiventricleGeom')#, boundary=True)
solid.save_pickle(folder_pickle+'BV_GEOMETRY_DATA', boundary=solid.boundary_names())
#solid.save_json(folder_json+'BV_GEOMETRY_DATA', boundary=solid.boundary_names())
#solid.save_txt(folder_txt+'BV_GEOMETRY_DATA', boundary=solid.boundary_names())

# Reload when desired
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_pickle(folder_pickle+'LV_GEOMETRY_DATA')
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_json(folder_json+'LV_GEOMETRY_DATA')
# cps, w, patchverts, patchcon, bnelems, knotval, knotmult, boundaries = solid.load_txt(folder_txt+'LV_GEOMETRY_DATA')