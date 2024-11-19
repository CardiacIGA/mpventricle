# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 16:59:06 2022

@author: s146407
"""
from nutils import mesh, function, export, elementseq, transformseq, topology
import matplotlib.pyplot as plt
from .surface_generator import SurfaceEllipsoid
from splipy import BSplineBasis, Surface
import pickle, treelog, numpy, stl, scipy, copy
_ = numpy.newaxis


def save_pickle(*data, filename : str = 'filename'):
    '''save_pickle

    Save a new 'filename.pickle' file.
    '''        
    with open(filename + '.pickle', 'wb') as f: 
         pickle.dump(data,f) 
    return

def load_pickle(filename : str = 'filename'):
    '''load_pickle

    Load an existing 'filename.pickle' file and output its results.
    '''
    with open(filename + '.pickle', 'rb') as f: 
      return pickle.load(f) 
      
#%
class CardiacGeometry:
    '''Class:
        
        Constructs a single left-ventricle geometry, or a complete bi-ventricle geometry 
        
        Input:
            
            Requires Rlv_x, Rlv_y, Rlv_z, R....
        
    '''    
    def __init__(self, geom_type: str, LVi: dict = {}, LVo: dict = {}, RVi: dict = {}, RVo: dict = {}, unit : str = 'm'): # Check and store the input parameters correctly for further use
        
        self.geom_type   = geom_type
        biventricle      = True if geom_type == 'bi-ventricle' else False # useful boolean 
        self.biventricle = biventricle
        
        ## Check if type is correct
        for i, iname in zip((LVi,LVo,RVi,RVo),('LVi','LVo','RVi','RVo')):
            assert type(i) == dict, "'{}' is not a dictionary!".format(i)
            for req_keys in ('Rx','Ry','Rz'):
                if bool(i): # Check if it is not an empty dict
                    assert req_keys in i, "The '{}' dictionary, is missing '{}' as an input".format(iname,req_keys)
        ##---------------------------- 
        
        
        ## Check if LVi and LVo are given, otherwise raise error
        if geom_type == 'bi-ventricle':
            if not bool(RVi) or not bool(RVo) or not bool(LVo) or not bool(LVi):
                raise ValueError("NotSupported: You should always specify the following 4 together: LVo, LVi, RVo, and  RVi")
            # if not bool(RVi) or not bool(RVo) or not bool(LVi):
            #     raise ValueError("NotSupported: You should always specify RVi and RVi when providing LVo")                
            # if bool(RVi) or bool(RVo) and  (not bool(LVo) or not bool(LVi)):
            #     raise ValueError("NotSupported: You may not specify right-ventricle input data, without outer wall left-ventricle input data")           
            
            # Assert origin
            assert 'Origin' in RVo or RVi, "The 'Origin' (numpy array) of the right-ventricle is not specified in 'RVi' or 'RVo'"
            
            # Filter out C (origin of the right-ventricle ellips) input if they are given in RVi AND RVo
            self.Origin = RVo['Origin'] if 'Origin' in RVo else RVi['Origin'] # Store Origin, depending in which dict it is given
            if 'Origin' in RVi and 'Origin' in RVo: 
                print("WarningMessage: Ellips Origin given twice '{}' and '{}', using value '{}' in remainder".format(RVo['Origin'],RVi['Origin'],RVo['Origin']))

        elif geom_type == 'left-ventricle':
            if bool(RVi) or bool(RVo):
                print("Specified right-ventricle input data is neglected because the geometry type is 'left-ventricle'")
            if not bool(LVi) or not bool(LVo):
                raise ValueError("Specify the left-ventricle input data, LVi and LVo, because you specified the geometry type = 'left-ventricle'")
        else: # Unknown geom_type
            raise ValueError("Specified geometry type '{}' unknown, pick either 'left-ventricle' or 'bi-ventricle'".format(geom_type))
        
        assert 'H' in LVo or LVi, "Truncation height 'H' not specified in input 'LVi' or 'LVo'"

        # Filter out H (truncation height) input if they are given in LVi AND LVo
        #self.H = LVo['H'] if 'H' in LVo else LVi['H'] # Store H, depending in which dict it is given
        if 'H' in LVi and 'H' in LVo: 
            print("WarningMessage: Truncation height given twice '{}' and '{}', using value '{}' in remainder".format(LVo['H'],LVi['H'],LVo['H']))
        ##---------------------------- 
        
               
        ## Store the asserted input in the desired dictionaries
        self.ventricle = dict(left=dict(), right=dict()) # Initialize data-type for storing data
        for iname, idata in zip( ('outer','inner'), (LVo, LVi)):
            trunch = True if iname == 'outer' else False
            self.ventricle['left'][iname] = StoredInput(idata, TruncH = trunch).convert(unit)
        if biventricle:
            for iname, idata in zip( ('outer','inner'), (RVo, RVi)):
                origin = True if iname == 'outer' else False
                self.ventricle['right'][iname] = StoredInput(idata, Origin = origin).convert(unit) # Store and convert to SI units




        ## Generate patch vertices
        self.patchverts, self.patchconnect = self.generate_vertices(biventricle) # return patchvertex coordinates and connectivity
        return
    
    def assert_input(self,geom_type): # asserts whether the provided inputdata is physically possible
        
        # left-ventricle asserts
        for i in ('Rx','Ry','Rz'):
            assert getattr(self.LV['inner'],i) < getattr(self.LV['outer'],i), "Radius of LV inner wall '{}' is larger than outer".format(i)
        assert self.LV['inner'].Rz and self.LV['outer'].Rz < self.LV['outer'].H, "Longitudinal LV radius 'Rz' is smaller than truncation height 'H'"
        
        
        # bi-ventricle asserts
        if geom_type == 'bi-ventricle':
            for i in ('Rx','Ry','Rz'):
                assert getattr(self.RV['inner'],i) < getattr(self.RV['outer'],i), "Radius of RV inner wall '{}' is larger than outer".format(i)
                if i != 'Rx': # Rx is allowed to be greater    
                    assert getattr(self.RV['outer'],i) < getattr(self.LV['outer'],i), "Radius of RV outer wall '{}' is larger than LV outer wall".format(i)
            assert ( self.RV['outer'].Rz + self.RV['outer'].Origin[-1] ) and ( self.RV['inner'].Rz + self.RV['outer'].Origin[-1] )  < self.LV['outer'].H, "Longitudinal LV radius 'Rz + Origin_z' is smaller than truncation height 'H' (incl. origin shift)"
            
            # 
            assert self.RV['outer'].Origin[1] > 1e-10, "NotSupported: origin of the right-ventricle can not be shifted in y-direction!"
            assert (self.RV['outer'].Origin[0] - self.RV['outer'].Rx ) < 0.8*self.LV['outer'].Rx, "Radius of the RV (Origin_x - Rx) is too small compared to LV radius Rx" # 0.8 is just a factor to make sure it is sufficiently smaller
            assert self.RV['outer'].Origin[0] < 0, "Origin_x of the RV should be negative (shifted to the left), but is positive"
            
        return
    
    
    def generate_surface(self, surface, pnelems_i : dict = {}, loftd_i : dict = {}, lnelems_i : dict = {}): # Return cps, w and nelems of that surface
        
        pnelems, loftd, lnelems, coons = self.generate_elem_distr(self.biventricle, pnelems_i,loftd_i,lnelems_i) # Number of elements per patch, Lofting direction per patch, lofting nelems per patch
        
        side, surf = surface.split('-') # split the name
        R          = [self.ventricle[side][surf].__getattribute__(r) for r in ("Rx","Ry","Rz")] # Radii 
        patchverts = self.patchverts[side][surf]                                                # Patch vertices
        patchconn  = self.patchconnect[side]                                                    # Patch connectivity
        boundcons  = self.generate_constrains(self.biventricle)[side][surf] if self.biventricle else self.generate_constrains(self.biventricle)                                   # Patch-boundary constrain type

        if self.biventricle:
            if side == 'left':
                popkeys = ['patch{}'.format(i) for i in range(11,15)]
                for delkey in popkeys:
                    pnelems.pop(delkey)
                    loftd.pop(delkey)
                    lnelems.pop(delkey)
                if surf == 'outer':    
                    coons  = {2: True, 6: True, 9: True}   
                origin = numpy.array([0,0,0])
            else:
                popkeys = ['patch{}'.format(i) for i in range(11)]
                newkeys = ['patch{}'.format(i) for i in range(4)]
                oldkeys = ['patch{}'.format(i) for i in range(11,15)]
                for delkey in popkeys:
                    pnelems.pop(delkey)
                    loftd.pop(delkey)
                    lnelems.pop(delkey)
                for newkey, oldkey in zip(newkeys, oldkeys):    
                    pnelems[newkey] = pnelems.pop(oldkey)
                    loftd[newkey]   = loftd.pop(oldkey)
                    lnelems[newkey] = lnelems.pop(oldkey)
                origin = self.ventricle['right']['outer'].Origin    
        else: # It is the left-ventricle
            origin = numpy.array([0,0,0])
            
            
        # Constrain the intersection curve for the right-ventricle geometry (constrained to the already calculated left-ventricle constrain curve)
        if self.biventricle:
            if side == 'right':
                assert hasattr(self, 'bcons_right'), "Left-ventricle outer surface should be determined first before any right-ventricle urfaces can be computed (shared intersection curves)"
                bcons_right = self.bcons_right[surf] 
            else:
                bcons_right = None  
        else:
            bcons_right = None  
        
        
        # Construct the surface
        SurfaceObj = SurfaceEllipsoid(R, patchverts, patchconn, bcons=boundcons, pnelems=pnelems, lnelems=lnelems, loftd=loftd, coons=coons, bcons_cpsw=bcons_right, S=origin)  
        cps, w     = SurfaceObj.get_cps_w()
        nelems     = SurfaceObj.nelems
        
        # save the intersection curve when the left outer surface is determined
        if self.biventricle:
            if side == 'left' and surf == 'outer':
               self.bcons_right   = dict()
               self.bcons_right['outer'] = self.get_intersection_cps_w('outer', SurfaceObj)
               self.bcons_right['inner'] = self.get_intersection_cps_w('inner', SurfaceObj)
        
        return MultipatchSurface(self.geom_type, surface, patchverts, patchconn, cps, w, nelems, pnelems) # Save as a MultipatchSurface class        
    
    
    
    ##-----------------------------## 
    ## General-ventricle functions ##
    ##-----------------------------##    
    def generate_vertices(self,biventricle): # Generates vertex locations, based on the geometry type given
        ## Left-ventricle patch vertices
        connectivity = self.generate_connectivity(biventricle) # dictionary: {left: array, right: array}
        patchverts   = dict(left=dict(), right=dict())
        if biventricle:
            patchvertices = self.Biventricle_vertices()
            i = 0
            for side in ('left','right'):
                for surf in ('outer','inner'):
                    patchverts[side][surf] = patchvertices[i] 
                    i += 1
        else: # it is left-ventricle
            patchverts['left']['outer'] = self.Leftventricle_vertices('outer')
            patchverts['left']['inner'] = self.Leftventricle_vertices('inner')
        return patchverts, connectivity

    
    
    def generate_connectivity(self,biventricle): # Return the connectivity of the entire cardiac geometry = fixed for all input
        connectivity = dict(left=None, right=None)
        if biventricle:
            connectivity['left'] = numpy.array([[ 0, 1, 2, 3], # patch 0
                                                [ 2, 3, 4, 5], # patch 1
                                                [ 4, 5, 6, 7], # patch 2
                                                [ 6, 7, 8, 9], # patch 3
                                                [10,11, 8, 9], # patch 4
                                                [12,13,10,11], # patch 5
                                                [14,15,12,13], # patch 6
                                                [0,  1,14,15], # patch 7
                                                [13, 7,11, 9], # patch 8
                                                [15, 5,13, 7], # patch 9
                                                [ 1, 3,15, 5]])# patch 10
            connectivity['right'] = numpy.array([[0,1,2,3], # patch 0 (11 for bi-ventr)
                                                 [4,5,2,3], # patch 1 (12 for bi-ventr)
                                                 [6,7,4,5], # patch 2 (13 for bi-ventr)
                                                 [7,1,5,3]])# patch 3 (14 for bi-ventr)
        else: # it is the left-ventricle
            connectivity['left'] = numpy.array([[0,1,2,3],
                                                [2,3,4,5],
                                                [6,7,4,5],
                                                [0,1,6,7],
                                                [1,7,3,5]])        
        return connectivity
    
    
    # generate the constrain dictionaries for the specific geometry
    def generate_constrains(self,biventricle):
        if biventricle:
            H      = self.ventricle['left']['outer'].H
            Origin = self.ventricle['right']['outer'].Origin
            Rright_o = [self.ventricle['right']['outer'].__getattribute__(r) for r in ("Rx","Ry","Rz")]
            Rright_i = [self.ventricle['right']['inner'].__getattribute__(r) for r in ("Rx","Ry","Rz")]
            
            boundcons_left   = dict.fromkeys([tuple(i) for i in [(0,2),(2,4),(4,6),(6,8),(10,8),(12,10),(14,12),(0,14)]], ('plane',numpy.array([ 0, 0, 1, H ]))) # Save top boundary connectivity, should be constrained, for left-ventricle
            boundcons_right  = dict.fromkeys([tuple(i) for i in [(0,2),(4,2),(6,4)]], ('plane',numpy.array([ 0, 0, 1, H ]))) # Save top boundary connectivity, should be constrained, for right-ventricle
            
            # intersection boundaries here...
            boundcons_left_intersec_o = dict.fromkeys([tuple(i) for i in [(4,5),(15,5),(14,15)]], ( 'ellips',Origin,Rright_o) ) # Save top boundary connectivity, should be constrained
            boundcons_left_intersec_i = dict.fromkeys([tuple(i) for i in [(6,7),(13,7),(12,13)]], ( 'ellips',Origin,Rright_i) ) # Save top boundary connectivity, should be constrained
            
            # TODO
            # bound_ellips_left = {(1,3): ('default',)}
            # boundcons_left.update(bound_ellips_left)
            # boundcons_right.update(bound_ellips_left)
            
            
            #bound_ellips_left_i = {(4,5): ('default',),(6,7): ('default',),(14,15): ('default',),(12,13): ('default',),
            #                       (15,5): ('default',),(13,7): ('default',)}
            # boundcons_right_intersec_o = dict.fromkeys([tuple(i) for i in [(0,1),(6,7),(7,1)]], ( 'ellips',Origin,Rright_o) ) # Save top boundary connectivity, should be constrained
            # boundcons_right_intersec_i = dict.fromkeys([tuple(i) for i in [(6,7),(13,7),(12,13)]], ( 'ellips',Origin,Rright_i) ) # Save top boundary connectivity, should be constrained
            # boundcons_right.update(boundcons_left_intersec_o)
            # boundcons_right.update(boundcons_left_intersec_i)
            
            boundcons  = dict(left  = dict(outer={**boundcons_left, **boundcons_left_intersec_o, **boundcons_left_intersec_i}, inner = boundcons_left ), #{**boundcons_left, **bound_ellips_left_i}
                              right = dict(outer=boundcons_right                                                             , inner = boundcons_right))
        else: # it is left-ventricle
            boundcons   = dict.fromkeys([tuple(i) for i in [(0,2),(2,4),(6,4),(0,6)]], ('plane',numpy.array([ 0, 0, 1, self.ventricle['left']['outer'].H ]))) # Save top boundary connectivity, should be constrained
            #boundcons = dict()
        
        return boundcons 
    
    
    def generate_elem_distr(self, biventricle, pnelems_i, loftd_i, lnelems_i):
        # Provide pre-determined element distribution, updated by an input and reverified for the existing elements
        # The pre-determined values are chosen such that they give the most desired result in general cases
        if biventricle:
            circ_x  = 1 # number of elements in circumferential direction for patches along the x-axis (5 in total)
            long    = 1 # number of elements in longitudinal direction for most of the patches, starting at patch0
            circ_y  = 1 # number of elements in circumferential direction for patches along the y-axis (3 in total)
            long_s  = 1 # number of elements in longitudinal direction for the patches that form the septum (excluding the middle one)
            right_n = 1 # number of elements of the right-ventricle in its thickness direction
            right_y = 1 # number of elements in circumferential direction for patches along the y-axis of the right-ventricle (3 in total)
            pnelems = {'patch0' : [long,circ_x],  'patch1' : [long,circ_y],  'patch2' : [long,right_n], 'patch3' : [long,long_s],   'patch4': [long,circ_x],
                       'patch5' : [long,long_s],  'patch6' : [long,right_n], 'patch7' : [long,circ_y],  'patch8' : [circ_x,long_s], 'patch9': [circ_x,right_n], 'patch10': [circ_x,circ_y],
                       'patch11': [long,right_y], 'patch12': [long,circ_x],  'patch13': [long,right_y], 'patch14': [circ_x,right_y]} # General element distribution per patch (should match with multipatches)
            loftd   = {'patch0' :  1, 'patch1' :  1, 'patch2' :  1, 'patch3' :  0, 'patch4' :  1,  'patch5':  0,
                       'patch6' :  1, 'patch7' :  1, 'patch8' :  0, 'patch9' :  0, 'patch10':  0,
                       'patch11':  1, 'patch12':  1, 'patch13':  1, 'patch14':  0}          # lofting direction, 0 or 1
            lnelems = {'patch0' :  1, 'patch1' :  1, 'patch2' :  1, 'patch3' :  1, 'patch4' :  1,  'patch5':  1,
                       'patch6' :  1, 'patch7' :  1, 'patch8' :  1, 'patch9' :  1, 'patch10':  1,
                       'patch11':  1, 'patch12':  1, 'patch13':  1, 'patch14':  1}   # nr of elements in lofting direction, value should be below pnelems
            coons   = {None: False}
            
            # # Check for the updated values, change the corresponding patches
            # patch_dep1 = {'patch0': 1,'patch2': 1,'patch4': 1} # Patch dependencies and corresponding direction on which it dependents
            # patch_dep2 = {'patch1': 1,'patch3': 1,'patch4': 0}
            # patch_dep3 = {'patch0': 0,'patch1': 0,'patch2': 0,'patch3': 0} # List patches which are coupled/dependent on eachother given the v-direction (circumferential). i.e., when 1 of them changes in nelems, all the rest should change as well!
            # patch_dep = [patch_dep1,patch_dep2,patch_dep3] # COmbine dependencies
            
            # # Update the nelems distribution, given an input, while using the dependencies
            # for pkey, pitem in pnelems_i.items():
                
            #     for patchd in patch_dep:
            #         if pkey in patchd.keys():
            #             for dkey, ditem in patchd.items():
            #                 pnelems[dkey][ditem] = pitem[ditem]
              
            # loftd.update(loftd_i)
            # lnelems.update(lnelems_i)
                
            for i in lnelems.keys(): 
                assert abs(loftd[i]) <= 1, "Lofting direciton can only have a value of 0 or 1"  
                assert lnelems[i] <= pnelems[i][loftd[i]], "Number of lofting elements is greater than max elements specified for {} in {}-direction".format(i,loftd[i])            
        else: # it is left-ventricle only
            long    = 1
            pnelems = {'patch0': [1,long], 'patch1': [1,long], 'patch2': [1,long], 'patch3': [1,long], 'patch4': [long,long]} # General element distribution per patch (should match with multipatches)
            loftd   = {'patch0':  1      , 'patch1':  1      , 'patch2':  1      , 'patch3':  1      , 'patch4':  1}          # lofting direction, 0 or 1
            lnelems = {'patch0':  1      , 'patch1':  1      , 'patch2':  1      , 'patch3':  1      , 'patch4':  1}          # nr of elements in lofting direction, value should be below pnelems
            coons   = {None: False} # specify which surfaces should be constructed using coons interpolation
            # Check for the updated values, change the corresponding patches
            patch_dep1 = {'patch0': 1,'patch2': 1,'patch4': 1} # Patch dependencies and corresponding direction on which it dependents
            patch_dep2 = {'patch1': 1,'patch3': 1,'patch4': 0}
            patch_dep3 = {'patch0': 0,'patch1': 0,'patch2': 0,'patch3': 0} # List patches which are coupled/dependent on eachother given the v-direction (circumferential). i.e., when 1 of them changes in nelems, all the rest should change as well!
            patch_dep = [patch_dep1,patch_dep2,patch_dep3] # COmbine dependencies
            
            # Update the nelems distribution, given an input, while using the dependencies
            for pkey, pitem in pnelems_i.items():
                
                for patchd in patch_dep:
                    if pkey in patchd.keys():
                        for dkey, ditem in patchd.items():
                            pnelems[dkey][ditem] = pitem[ditem]
              
            loftd.update(loftd_i)
            lnelems.update(lnelems_i)
                
            for i in lnelems.keys(): 
                assert abs(loftd[i]) <= 1, "Lofting direction can only have a value of 0 or 1"  
                assert lnelems[i] <= pnelems[i][loftd[i]], "Number of lofting elements is greater than max elements specified for {} in {}-direction".format(i,loftd[i])            

        return pnelems, loftd, lnelems, coons    
    
    def get_intersection_cps_w(self,surf, LV_surf): # surf = inner, outer 
        bcons_cpsw = {}
        Rright     = numpy.array([self.ventricle['right'][surf].__getattribute__(r) for r in ("Rx","Ry","Rz")])
        Rleft      = numpy.array([self.ventricle['left']['outer'].__getattribute__(r) for r in ("Rx","Ry","Rz")])
        origin     = self.ventricle['right']['outer'].Origin
        idx_left   = dict( outer=[(4,5),(15,5),(14,15)], inner=[(6,7),(13,7),(12,13)])
        for i, j in zip(idx_left[surf],[(0,1),(7,1),(6,7)]):
            W   = LV_surf.LHSb[i][1]
            CPS = LV_surf.LHSb[i][0].reshape(-1,3)/W[:,numpy.newaxis] # Actual coord in unit sphere of left ventr
            Xellips_left  = CPS.reshape(-1,3)*Rleft - origin # convert to left-ventricle ellips from unit sphere and shift it by the origin
            Xellips_right_unit = Xellips_left/Rright             # scale it down to the unitsphere, but this time by the right-ventricle radii
            Xunit = Xellips_right_unit
            bcons_cpsw[j] = (Xunit*W[:,numpy.newaxis], W)
        return bcons_cpsw
    
    
    
    
    
    
    ##--------------------------## 
    ## Left-ventricle functions ##
    ##--------------------------##
    def Leftventricle_vertices(self,wall): # wall = 'inner' or 'outer'
        H       = self.ventricle['left']['outer'].H  
        param   = self.ventricle['left'][wall] 
        φ       = numpy.pi/2
        θtrunc  = self.get_θ_from_Z( H, param.Rz )  # Truncated angle theta
        θbottom = 0.5*(θtrunc + numpy.pi)           # Angle at which bottom patch vertices are located (all 4)
        
        patchverts_sph = numpy.array([[θtrunc ,  0.0],
                                      [θbottom,  0.0],
                                      [θtrunc ,    φ],
                                      [θbottom,    φ],
                                      [θtrunc ,  2*φ],
                                      [θbottom,  2*φ],
                                      [θtrunc ,  3*φ],
                                      [θbottom,  3*φ]])
        patchverts_cart = self.to_cartesian(patchverts_sph, numpy.array([param.Rx, param.Ry, param.Rz]) ) # Convert to cartesian coordinates     
        return patchverts_cart
    
    
    ##------------------------## 
    ## Bi-ventricle functions ##
    ##------------------------##
    def Biventricle_vertices(self,): # wall = 'inner' or 'outer', vertices for the surfaces
    
        # Useful quantities
        H       = self.ventricle['left']['outer'].H  
        origin  = self.ventricle['right']['outer'].Origin
        θtrunc_left  = self.get_θ_from_Z( H, self.ventricle['left']['outer'].Rz )  # Truncated angle theta
        θtrunc_right_outer  = self.get_θ_from_Z( H + origin[2], self.ventricle['right']['outer'].Rz )
        θtrunc_right_inner  = self.get_θ_from_Z( H + origin[2], self.ventricle['right']['inner'].Rz )
        
        θbottom_left  = 0.5*(θtrunc_left + numpy.pi)           # Angle at which bottom patch vertices are located (all 4)
        #θbottom_right = 0.5*(θtrunc_right_outer + numpy.pi)
        
        # Unpack the radii for simplicity
        Rleft_i   = numpy.array([self.ventricle['left']['inner'].__getattribute__(r) for r in ("Rx","Ry","Rz")] )
        Rleft_o   = numpy.array([self.ventricle['left']['outer'].__getattribute__(r) for r in ("Rx","Ry","Rz")] )
        Rright_i  = numpy.array([self.ventricle['right']['inner'].__getattribute__(r) for r in ("Rx","Ry","Rz")]) 
        Rright_o  = numpy.array([self.ventricle['right']['outer'].__getattribute__(r) for r in ("Rx","Ry","Rz")]) 
        
        
        
        
        
        
        ##################################################
        ##### Determine the left-ventricle vertices ######
        ##################################################
        
        # 1) Determine the 4 (but 8 in reality, mirrored) intersection points between the left- and right-ventricles
        Zinter_outer_apex = self.get_Z_from_θ(0.6*(numpy.pi+θtrunc_left),[Rright_o[2],Rright_o[1]],True_θ=False) # z-coordinate at which the vertex is defined at the intersection between the ellipsoids 
        Zinter_inner_apex = self.get_Z_from_θ(0.575*(numpy.pi+θtrunc_left),[Rright_i[2],Rright_i[1]],True_θ=False)
        # Zinter_outer_apex = self.get_Z_from_θ(0.55*(numpy.pi+θtrunc_left),[Rright_o[2],Rright_o[1]],True_θ=False) # z-coordinate at which the vertex is defined at the intersection between the ellipsoids 
        # Zinter_inner_apex = self.get_Z_from_θ(0.5*(numpy.pi+θtrunc_left),[Rright_i[2],Rright_i[1]],True_θ=False)        
        Xinter_outerouter_base =  self.get_intersection(Rleft_o, Rright_o, OriginRight=origin, Z=H) # has size (2x3), because it is mirrored relative to x-axis
        Xinter_outerouter_apex =  self.get_intersection(Rleft_o, Rright_o, OriginRight=origin, Z=Zinter_outer_apex)
        Xinter_innerouter_base =  self.get_intersection(Rleft_o, Rright_i, OriginRight=origin, Z=H)
        Xinter_innerouter_apex =  self.get_intersection(Rleft_o, Rright_i, OriginRight=origin, Z=Zinter_inner_apex)
        
        # 2) Define the 4 Septum coordinates
        w = 0.25 # factor that scales the z-location
        Zseptum_mid   = (1-w)*Zinter_inner_apex + w*H#0.5*(Zinter_inner_apex + H)
        θseptum_outer = self.get_θ_from_Z(Zseptum_mid, Rleft_o[2])
        φinter_inner  = numpy.arctan2( Xinter_innerouter_base[0,1], abs(Xinter_innerouter_base[0,0]) ) # angle of the intersection coordinate with inner right-ventricle and left-ventricle wall
        φseptum       = numpy.pi - (1/2)*φinter_inner
        Xseptum_outer_ypos  = self.to_cartesian(numpy.array([[θtrunc_left,  φseptum], [θseptum_outer,  φseptum+(1/4)*φinter_inner]]), Rleft_o)
        Xseptum_outer_yneg  = self.to_cartesian(numpy.array([[θtrunc_left, -φseptum], [θseptum_outer, -φseptum-(1/4)*φinter_inner]]), Rleft_o)
                
        
        # 3) Define the remaining 4 LV coordinates (opposite to right-ventricle)
        φleft_outer_a  = numpy.arctan2( Xinter_outerouter_base[0,1], abs(Xinter_outerouter_base[0,0]) ) # angle of the intersection coordinate with inner right-ventricle and left-ventricle wall
        φleft_outer    = (1/3)*(numpy.pi -  φleft_outer_a)
        Xleft_outer_ypos = self.to_cartesian(numpy.array([[θtrunc_left,  φleft_outer], [θbottom_left,  φleft_outer]]), Rleft_o)
        Xleft_outer_yneg = self.to_cartesian(numpy.array([[θtrunc_left, -φleft_outer], [θbottom_left, -φleft_outer]]), Rleft_o)


        # 4) Outer patchvertices with correct ordering (corresponding the connectivity array)
        patchverts_left_outer_cartesian = numpy.concatenate([  Xleft_outer_yneg,
                                                               Xleft_outer_ypos,
                                                              [Xinter_outerouter_base[0]],
                                                              [Xinter_outerouter_apex[0]],
                                                              [Xinter_innerouter_base[0]],
                                                              [Xinter_innerouter_apex[0]],
                                                               Xseptum_outer_ypos,
                                                               Xseptum_outer_yneg,
                                                              [Xinter_innerouter_base[1]],
                                                              [Xinter_innerouter_apex[1]],
                                                              [Xinter_outerouter_base[1]],
                                                              [Xinter_outerouter_apex[1]]], axis=0)
        # determine left-ventricle inner vertices based on outer
        idx_base = [0,2,4,6,8,10,12,14] # indices that correspond to the base vertices (should get a different angle theta)
        patchverts_left_inner_ellips             = self.convert_to_spherical(patchverts_left_outer_cartesian, Rleft_o)
        patchverts_left_inner_ellips[idx_base,0] = self.get_θ_from_Z( H, self.ventricle['left']['inner'].Rz )
        
        
        
        
        ## Optional: 
            # a) set the z-coordinate of the inner wall apex intersection points equal to the z-coordinate of the outer wall (gives better shape of element/patch)
        fac1 = Rleft_i[2]/Rright_i[2]# factor that ensures the following coordinates scale accordingly when the lv inner radius is reduced
        fac2 = Rleft_i[2]/Rright_o[2]
        Z_leftinner = Zinter_inner_apex if fac1*1.1 > 1 else -Rleft_i[2] + 0.5*(Rleft_o[2]-Rleft_i[2])
        Z_leftouter = Zinter_outer_apex if fac2*1.1 > 1 else -Rleft_i[2] + 0.3*(Rleft_o[2]-Rleft_i[2])
        #patchverts_left_inner_ellips[[7,13],0] = self.get_θ_from_Z( Z_leftinner, self.ventricle['left']['inner'].Rz ) 
        patchverts_left_inner_ellips[[5,15],0] = self.get_θ_from_Z( Z_leftouter, self.ventricle['left']['inner'].Rz ) 
        
            # 2b) Recalculate Zinter_inner_apex, using vectors to determine ideal location of vertex: Xinter_innerouter_apex
            # create a plane that goes through the origin and Xinter_outerouter_apex, Xseptum_outer_ypos and find the intersection between this plane and the 'left- right (inner) ventricle intersection curve' (latter is unknown, so we adopt an iterative method)
        mapping = ( Rleft_i / Rleft_o ) #  mapper to the inner surface
        nplane  = numpy.cross( self.to_cartesian(patchverts_left_inner_ellips[[5]], Rleft_i), Xseptum_outer_ypos[1]*mapping )# normal of plane
        θi      = self.get_θ_from_Z( Z_leftinner, self.ventricle['left']['inner'].Rz )
        φi      = numpy.pi - numpy.arctan2( Xinter_innerouter_apex[0,1], abs(Xinter_innerouter_apex[0,0]) )
        
        δφ = dφ = 0.01*φi# step-size
        pos     = True
        while True:
            X    = self.to_cartesian(numpy.array([[ θi, φi + δφ ]]), Rleft_i) #self.get_intersection(Rleft_o, Rright_i, OriginRight=origin, Z=Zinter_inner_apex+δZ)[0]*mapping
            dist = numpy.dot(numpy.squeeze(X),numpy.squeeze(nplane)) # distance of the point X to the plane, should be 0 ideally
            if dist > 0:# or sign == -1*sign: # ensure we are not drifting away but approaching the plane + check if the sign switches (which means we overshot the ideal value and should stop)
               if not pos:
                   break
               δφ += dφ 
               pos  = True
            else: 
               δφ -= dφ
               pos  = False  
        patchverts_left_inner_ellips[[7,13]] = numpy.array([[ θi,   φi + δφ ],
                                                            [ θi, -(φi + δφ) ]]) 
        
        
        ## Convert LV inner spherical coordinates to cartesian
        patchverts_left_inner_cartesian     = self.to_cartesian(patchverts_left_inner_ellips, Rleft_i) # convert back to cartesian, but with different inner radius
        

        
        
        #####################################################
        ##### Determine the right-ventricle coordinates #####
        #####################################################
        
        
        # 1) Determine spherical coordinates (angles), based on intersection coordinates 
        # idx_inter = [4,5] # indices of the intersection vertices for y>0
        # θφright_inter_outer = self.convert_to_spherical(patchverts_left_outer_cartesian[idx_inter,:] - origin, Rright_o) # Determine the angles of intersection, when viewed from the right-ventricle
        # φright_outer_ypos   = 0.5*(numpy.pi + θφright_inter_outer[0,-1])  # use this angle to determine remaining cps angles
        idx_inter = [6,7]
        θφright_inter_inner = self.convert_to_spherical(patchverts_left_outer_cartesian[idx_inter,:] - origin, Rright_o) # Determine the angles of intersection, when viewed from the right-ventricle
        φright_outer_ypos   = 0.5*(numpy.pi + θφright_inter_inner[0,-1])  # use this angle to determine remaining cps angles
        

        # 2) Determine the right-ventricle outer vertices and enforce the calculated left-ventricle coordinates
        θbottom_right = self.get_θ_from_Z( Xinter_outerouter_apex[0][2], self.ventricle['right']['outer'].Rz )
        patchverts_right_outer_spherical = numpy.array([ [θtrunc_right_outer,  φright_outer_ypos],
                                                         [θbottom_right     ,  φright_outer_ypos],
                                                         [θtrunc_right_outer, -φright_outer_ypos],
                                                         [θbottom_right     , -φright_outer_ypos] ])
        patchverts_right_outer_cartesian = self.to_cartesian(patchverts_right_outer_spherical, Rright_o) + origin
        patchverts_right_outer_cartesian = numpy.concatenate([  [Xinter_outerouter_base[0]],  # Match specific vertices (intersection verts)
                                                                [Xinter_outerouter_apex[0]],
                                                                patchverts_right_outer_cartesian,
                                                                [Xinter_outerouter_base[1]],  # Match specific vertices (intersection verts)
                                                                [Xinter_outerouter_apex[1]],], axis=0) 
        
        # 3) Copy the outer vertices and modify for the inner vertices
        patchverts_right_inner_spherical   = patchverts_right_outer_spherical.copy()
        patchverts_right_inner_spherical[[0,2],0] = θtrunc_right_inner 

        patchverts_right_inner_cartesian = self.to_cartesian(patchverts_right_inner_spherical, Rright_i) + origin
        patchverts_right_inner_cartesian = numpy.concatenate([  [Xinter_innerouter_base[0]],  # Match specific vertices (intersection verts)
                                                                [Xinter_innerouter_apex[0]], 
                                                                patchverts_right_inner_cartesian,
                                                                [Xinter_innerouter_base[1]],  # Match specific vertices (intersection verts)
                                                                [Xinter_innerouter_apex[1]],], axis=0)        
                                       
        return patchverts_left_outer_cartesian, patchverts_left_inner_cartesian, patchverts_right_outer_cartesian, patchverts_right_inner_cartesian
  
    
    
    ##----------------------##
    ## Conversion functions ##
    ##----------------------##
    def get_Z_from_θ(self,θ,R,True_θ=True):
        if True_θ: # actual ellips angle θ
            return R*numpy.cos(θ)
        else:
            sign = -1 if abs(θ) > 0.5*numpy.pi else 1
            θ = numpy.pi - abs(θ) if sign==-1 else abs(θ)
            t = numpy.arctan( (R[0]/R[1])*numpy.tan( θ ) ) # t is the unit-sphere angle, not scaled with R
            t = t+numpy.pi/2 if θ > numpy.pi/2 else t
            return sign*R[0]*numpy.cos(t)
    def get_θ_from_Z(self,Z,R):
        return numpy.arccos(Z/R)
    def to_cartesian(self, X, R):
        #print(X)
        return numpy.array([R[0]*numpy.sin(X[:,0])*numpy.cos(X[:,1]),
                            R[1]*numpy.sin(X[:,0])*numpy.sin(X[:,1]),
                            R[2]*numpy.cos(X[:,0])]).T 
    def convert_to_spherical(self,Xcart, R): # convert cartesian to unitsphere coordinates
        θ = self.get_θ_from_Z(Xcart[:,-1].ravel(),R[2])
        φ = numpy.zeros(len(θ))
        for i, (x, y) in enumerate(Xcart[:,:-1]):
            if x < 0:    
                φ[i] = numpy.pi - numpy.arctan2( y, abs(x) )
            else:
                φ[i] = numpy.arctan2( y, x )
        return numpy.concatenate([ [θ], [φ] ],axis=0).T
    
    def get_2D_R_from_Z(self,R,Z,Zshift=0.):
        f =  1 - ( Z - Zshift )**2 / R[2]**2
        return numpy.sqrt(f)*R
    def get_intersection(self, Rleft, Rright, OriginRight=numpy.array([0,0,0]), Z=0.):
        Rleft  = self.get_2D_R_from_Z(Rleft,Z) # determine the radii at a specific heights in 2D
        Rright = self.get_2D_R_from_Z(Rright,Z,OriginRight[-1])
    
        gx     = (Rright[0]/Rleft[0])**2 
        gy     = (Rright[1]/Rleft[1])**2
        xshift = OriginRight[0]
        A  = (gx-gy)
        B  = 2*gy*xshift
        C  = (Rright[0]**2) *(gy-1) - gy*(xshift**2)
    
        Xsolved = numpy.roots([A,B,C])
        Xsolv   = min(Xsolved) if min(Xsolved) < Rleft[0] and min(Xsolved) > -Rleft[0] else Xsolved[1]
    
        Ysolv = Rleft[1]*numpy.sqrt( 1 - (Xsolv/Rleft[0] )**2 )#- (Z/Rleft[1])**2 ) 
        return numpy.array([ [Xsolv,Ysolv,Z],[Xsolv,-Ysolv,Z] ])
    
    
    
    
    
    ##---------------------------##
    ## Post-processing functions ##
    ##---------------------------##
    def visualize(self,*surface):

        H = self.ventricle['left']['outer'].H
        
        # from matplotlib.colors import LightSource
        # from matplotlib import cm
        # light = LightSource(90, 45)
        # illuminated_surface = light.shade(xl, cmap=cm.coolwarm)
        fig = plt.figure(figsize=(8, 6))
        colors = ('green','grey','blue','orange')    
        for j, isurf in enumerate(surface):
            side, surf = isurf.split('-') # split the name
            R = [getattr(self.ventricle[side][surf], i) for i in ('Rx','Ry','Rz')]
            Vertex = self.patchverts[side][surf]
            

            S = self.ventricle[side]['outer'].Origin if side == 'right' else numpy.zeros(3)
        
            
            u  = numpy.linspace(0, 2*numpy.pi, 50)
            vl = numpy.linspace(self.get_θ_from_Z(H,R[2]),  numpy.pi, 50)
            
            xl = R[0] * numpy.outer(numpy.cos(u), numpy.sin(vl)) + S[0] 
            yl = R[1] * numpy.outer(numpy.sin(u), numpy.sin(vl)) + S[1]
            zl = R[2] * numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(vl)) + S[2]
            f  = R[2]*2 # plot limit factor
    

            ax = fig.gca(projection='3d')
            ax.plot(Vertex[:,0], Vertex[:,1], Vertex[:,2], 'rs', label='Vertices')
            #ax.plot_surface(xl, yl, zl, rstride=1, cstride=1, antialiased=False)#, facecolors=illuminated_surface)
            ax.plot_wireframe(xl, yl, zl, rstride=2, cstride=2,color=colors[j])
            ax.view_init(20, -135)
            ax.legend()
            ax.set_xlim(-f, f)
            ax.set_ylim(-f, f)
            ax.set_zlim(-f, f)
        plt.show() 
        return
    
    
    
    
    
#%   
class StoredInput:
    '''Class: StoredInput
        
        Class which stores useful info regarding dimensions. Is able to convert the unit of the stored values. 
        
        Input:
            
            Requires: Rx, Ry, Rz, 
            Optional: 
                TruncH: Boolean which says whether the Input dict conatains a truncated height key 'H'
                Origin: Boolean which says whether the Input contains an origin numpy array key 'Origin'
                
    '''  
    def __init__(self, Input : dict, TruncH : bool = False, Origin : bool = False, Rseptum : bool = False):
        self.Rx = Input['Rx']
        self.Ry = Input['Ry']
        self.Rz = Input['Rz']
        if TruncH: # Truncation height is provided and should be stored
            self.H = Input['H']
        if Origin: # Origin is provided and should be stored
            self.Origin = Input['Origin']
        if Rseptum: # Septum radius is provided and should be stored
            self.Rxs = Input['Rxs']
            
        self.unit = 'm' # Set defaul unit to meter    
        #self.converted = False # Units have not been tampered/converted at this point    
        return 
    def __repr__(self,):
        return "Class: VentricleDataStructure"
    # def convert(self,unit): # Converts the stored variables from an input 'unit = dm,cm,mm,..' to the SI-unit METER
    #     Units = dict(m=1, dm=1e-1, cm=1e-2, mm=1e-3) # unit conversions to meter
    #     assert not self.converted, "Parameters have already been converted ones, multiple conversion are not desired!"
    #     for key, value in vars(self).items():
    #         if key != 'converted':
    #             setattr(self, key, value*Units[unit])
    #     self.converted = True
    #     return
    
    def convert(self,unit): # Converts the stored variables from an input 'unit = dm,cm,mm,..' to the SI-unit METER
        Units = dict(m=1, dm=1e-1, cm=1e-2, mm=1e-3) # unit conversions to meter
        assert unit in Units.keys(), "Unknown unit specified '{}', change to one of the following [m, dm, cm, mm]".format(unit)       
        convfac = Units[unit]/Units[self.unit] # Determine the conversion factor: ( required unit / current unit )
        for key, value in vars(self).items():
            if key != 'unit':
                setattr(self, key, value*convfac)
        self.unit = unit # Change unit type
        return self
    
    
    
    
    
    
    
    
#%    
class MultipatchSurface:
    '''Class:
        
        Extracts surface points for specific input: lv, rv, bi-v 
        
        Input:
            
            Requires Rlv_x, Rlv_y, Rlv_z, R....    
    '''    
    def __init__(self, geom_type, surface, patchverts, patchcon, cps, w, nelems, pnelems, index_lv=None, bnames={}): # in nutils #TODO make nelems and pnelems = ...
        self._geom_type  = geom_type  # geometry type: left- or bi-ventricle
        self._surface    = surface    # surface type: left-outer, left-inner, right-outer, right-inner
        self._patchverts = patchverts # patchvertices of the multipatch surface
        self._patchcon   = patchcon   # patchconnectivity array of the multipatch surface
        self._cps        = cps        # cntrol points of the multipatch surface
        self._w          = w          # weights corresponding the cps
        self._bnames     = self._convert_boundary_names(bnames) # The boundary names of the multipatch surface

        if None in nelems: # If some patches are denoted by None
            nelems = self.expand_nelems(nelems)

        if None in pnelems: # If some patches are denoted by None
            pnelems = self.expand_pnelems(pnelems)

        self._bnelems    = nelems  # elements for each boundary (0,1): 2, ...
        self._pnelems    = pnelems # elements for each patch (u,v)
        
        self._nelems = nelems # Dirty fix for now
        
        
        # Calculate knotvalues and knotmultiplicity (Currently fixed because we have a fixed spline order, quadratic, and we know the number of elements)
        self._knotval, self._knotmult = self.calc_knot_info(nelems)      # knotvalues for each boundary 
        #self._knotmult   =       # knotmultiplicity of the knotvalues for each boundary
        
        self._nrcps      = {int(key[5:]): (val[0]+2)*(val[1]+2) for key, val in self._pnelems.items()}  
        if index_lv:
            self._index_lv = index_lv # array:  indices of the 
            
        # Additional info
        self._input_type = 'nutils'
        return
        
    def __repr__(self,):
        return "Class: MultipatchSurface"
    
    def __str__(self,):
        Text  = "MultipatchSurface class of:\n"
        Text += " -Geometry: {}\n".format(self._geom_type)
        Text += " -Surface : {}\n".format(self._surface)
        return Text
    
    def __getitem__(self, patchID):
        self.assert_patchID(patchID)

        pnelems_patch = self._pnelems[f"patch{patchID}"]
        nu, nv = pnelems_patch # nr of elements in (u,v) directions

        # Create correct nelems dict
        ndims = min( len(self.patchcon(patchID)) // 2, 3) # dimension, curve, surface or solid
        assert ndims == 2, "Something went wrong, only splitting of 2-dimensional  multipatch domains possible"
        bound = []
        connectivity = self.patchcon(patchID).reshape([2]*ndims)
        
        for i in range(2):
            bound += [tuple(connectivity[i,:])]
            bound += [tuple(connectivity[:,i])]

        nelems_patch = {k: v for k, v in self._nelems.items() if k in bound}
        patchSurface = MultipatchSurface(self._geom_type,
                                         f"{self._surface} (patch {patchID})",
                                         self.patchverts(patchID),
                                         self.patchcon(patchID),
                                         self.cps(patchID),
                                         self.weights(patchID),
                                         nelems_patch,
                                         {f"patch{patchID}":pnelems_patch})
        return patchSurface
    
    def __iter__(self):
        self.patch_iter = 0
        return self
    
    def __next__(self):
        if self.patch_iter < self.nrpatches():
            patchSurface = self.__getitem__(self.patch_iter)
            self.patch_iter += 1
            return patchSurface
        else:
            raise StopIteration

    def copy(self,):
        return copy.deepcopy(self)
    
    # Return functions
    def cps(self,patchID=None): # return cps for specific patchID < nr_patches
        '''Return the control points of the given patch.
    
        Parameters
        ----------
        patchID : int :
            The index of the patch. Specify 'None' to obtain all patches.
    
        Returns
        -------
        : numpy.array :
            Array containing the cps
        '''
        if patchID==None:
            return self._cps
        else:
            self.assert_patchID(patchID)
            nrcps = numpy.cumsum( [0] + [j for i, j in self._nrcps.items() ])
            return self._cps[nrcps[patchID]:nrcps[patchID+1]] 
        
    def weights(self,patchID=None):
        '''Return the weights of the given patch.
    
        Parameters
        ----------
        patchID : int :
            The index of the patch. Specify 'None' to obtain all patches.
    
        Returns
        -------
        : numpy.array :
            Array containing the weights
        '''
        if patchID==None:
            return self._w
        else:
            self.assert_patchID(patchID)
            nrw = numpy.cumsum( [0] + [j for i, j in self._nrcps.items() ])
            return self._w[nrw[patchID]:nrw[patchID+1]] 
        
    def surface(self,):
        return self._surface
    def patchverts(self,patchID=None):
        if patchID == None:
            return self._patchverts
        else:
            return self._patchverts[patchID]
    def patchcon(self,patchID=None):
        if patchID == None:
            return self._patchcon
        else:
            return self._patchcon[patchID]
    def bnelems(self,biventr=False):
        if biventr: # Dirty fix for now
            return self._bnelems
        else:
            bnelems = {}
            for key, val in self._knotval.items():
                bnelems[key] = len(val)-1 
            self._bnelems = bnelems 
            return self._bnelems
        #return self._bnelems
    # def pnelems(self,): 
    #     for i in range(self.nrpatches()):           
    #         cps        = self.cps(i).copy
    #         pnelems[i] = 
    #     return self._pnelems
    def knotval(self,boundary=None):
        if boundary == None:
            return self._knotval
        else:
            return self._knotval[boundary]
    def knotmult(self,boundary=None):
        if boundary == None:
            return self._knotmult
        else:
            return self._knotmult[boundary]
        
    def knotvec(self,patchID=None): # Knotvector per patch (u,v)

        #bound = (0,1),(0,2)
        ndims = min( len(self.patchcon(patchID)) // 2, 3) # dimension, curve, surface or solid
        bound = []
        assert ndims == 2, "Only supports 2-dimensional multipatch for knotvec evaluation."
        #if ndims == 3: bound += (0,4)  
        connectivity = self.patchcon(patchID).reshape([2]*ndims)
        bound += [tuple(connectivity[0,:])]
        bound += [tuple(connectivity[:,0])]
        if ndims == 3: bound += [tuple(connectivity[:,:,0])]
        
        knotvec = []
        for ib in bound:
            knotvec += [self._det_knotvec(self.knotval(ib), self.knotmult(ib))]
    
        return knotvec
    

    def expand_nelems(self, nelems):
        
        patches_conn = self.patchcon()
        ndims = min( len(patches_conn[0]) // 2, 3) # dimension, curve, surface or solid
        assert ndims == 2, "Only supports 2-dimensional multipatch for knotvec evaluation."
        bound = []
        for pconn in patches_conn:
            bconn = pconn.reshape([2]*ndims)
            bound = [tuple(i) for i in bconn]+ [tuple(i) for i in bconn.T]
            for b in bound:
                if b not in nelems:
                    nelems[b] = nelems[None]
        nelems.pop(None)            
        return nelems

    def expand_pnelems(self, pnelems):
        for i in range(len(self.patchcon())): # Loop over index of patches
            if f"patch{i}" not in pnelems:
                pnelems[f"patch{i}"] = pnelems[None]
        pnelems.pop(None)                   
        return pnelems
    
    @staticmethod
    def _det_knotvec(knotval, knotmult):
        knotvec = []
        for val, mult in zip(knotval.copy(), knotmult.copy()):
            knotvec += [val]*mult
        return knotvec
    
    @staticmethod
    def calc_knot_info(nelems):
        knotv = {};knotm = {}
        for key, val in nelems.items():
            knotv[key] = [v/val for v in range(val)] + [1] # knotvector (without multiplicity) per boundary
            knotm[key] = [3] + [1]*(val-1) + [3]
        return knotv, knotm
    
    def nrcps(self,patchID=None):
        '''Return the number of cps per patch in dict format (keys are patchID).
        '''
        if patchID == None:
            return self._nrcps
        else:
            return self._nrcps[patchID]
    
    def nrpatches(self,):
        '''Return the number of patches of the multipatch geometry.
        '''
        return self._patchcon.shape[0]
    ## Assert whether the given patchID is allowed
    def assert_patchID(self, patchID):
        assert type(patchID) == int, f"Provided patchID is of the wrong type {type(patchID)}, should be integer." 
        assert patchID < self.nrpatches(), "Provided patchID '{}' exceeds maximum patchID '{}'".format(patchID,self.nrpatches()-1) 
        return    
    
    
    ## Folowing functions are mainly used for modifying the cps, w and knotval/knotmult data
    def dependencies(self,patchID,add_return=False):
        patchcon   = self.patchcon(patchID)
        self.update_patchverts() #patchverts = self._patchverts(patchID)
        IDref    = numpy.array([i for i in range(len(patchcon))])
        boundref = numpy.reshape(numpy.array([0,1,2,3]), [2]*2) # reference surface patch connectivity
        Boundref = numpy.concatenate([ boundref, boundref.T ])
        

        depend = {}
        Depend = {}
        Ncps = numpy.cumsum([0]+[len(self.cps(i)) for i in range(len(self.patchcon()))])
        for ipatch in range(self.nrpatches()):
            boundref = numpy.reshape(numpy.array([0,1,2,3]), [2]*2)

            if ipatch != patchID:
                condition_sub   = [ ( idx in self.patchcon(ipatch) ) for idx in patchcon ] # Mask of sub in main
                condition_main  = [ ( idx in patchcon ) for idx in self.patchcon(ipatch) ] # Mask of main in sub
                if True in condition_main: # If there is an equal index return true
                    key   = (patchID, ipatch) # Keys of the main patch and the patch which it is connected to 
                    #assert self.patchcon(ipatch)[condition_main].tolist() == patchcon[condition_sub].tolist()
                    IDs   = self.patchcon(ipatch)[condition_main].tolist()   
                    #print(self.patchcon(ipatch) )
                    depend[key] = IDs
                    
                    ## Go from the global patchverts combi (IDs) to the local cps indices of that patch
                    IDcps = []
                    for ikey in key:
                        bound    = numpy.reshape( numpy.array([ 0,
                                                                self.bnelems()[tuple(IDs)] + 1,
                                                                self.nrcps(ikey) - self.bnelems()[tuple(IDs)] - 2,
                                                                self.nrcps(ikey) - 1]), [2]*2) 
                        Bound    =  numpy.concatenate([ bound, bound.T ])
                        Boundref = numpy.concatenate([ self.patchcon(ikey).reshape([2]*2), self.patchcon(ikey).reshape([2]*2).T ])
                        
                        # Find corresponding indices
                        mask   = [(b == IDs).all() for b in Boundref ]
                        IDcps += [Bound[mask][0]]
                    
                        # print(Boundref)
                        # print(IDs)
                        # print(mask)
                        # print(Bound)
                        # print(Bound[mask])
                    Depend[key] = IDcps    
                    #print(bound)
        # print(Depend)
        # # find the patchvertices cps and w indices that are matched       
        # Depend = {}
        # #patch_boundaries = numpy.reshape(self.patchcon(ipatch), [2]*2)
        
        
        # for ikey, ival in depend.items():
        #     ID_main=[];ID_sub=[]
        #     for i in ival:
        #         diff_main = numpy.linalg.norm(self.cps(ikey[0]) - self.patchverts(i),axis=1)
        #         ID_main.append(numpy.argwhere( diff_main < 1e-9 )[0][0])
        #         diff_sub  = numpy.linalg.norm(self.cps(ikey[1]) - self.patchverts(i),axis=1)
        #         ID_sub.append(numpy.argwhere( diff_sub < 1e-9 )[0][0])
        #     Depend[ikey]  = [ID_main, ID_sub] # contains the IDs of the complete cps and w arrays that should be equal for the specific patchID   
           
        # Calculate the entire index range between the vertices  
        DependIDs  = {}
        ncps_bound = self.bnelems()
        Ncps = numpy.cumsum([0]+[len(self.cps(i)) for i in range(len(self.patchcon()))]) # nr of cps with increasing patch ID 
        for (CPSkey, CPSval), ( VERkey, VERval ) in zip( Depend.items() , depend.items() ):
            if tuple(VERval) in ncps_bound.keys():
                nelems = ncps_bound[tuple(VERval)]
            elif tuple(VERval[::-1]) in ncps_bound.keys():
                nelems = ncps_bound[tuple(VERval[::-1])]
                
            DependIDs[CPSkey] = [numpy.linspace(CPSval[0][0],CPSval[0][1],nelems+2, dtype=int),
                                 numpy.linspace(CPSval[1][0],CPSval[1][1],nelems+2, dtype=int) + Ncps[CPSkey[1]]]
                
        #print(depend) 
        # print(Depend)
        # print(DependIDs)  
        if add_return:
            return DependIDs, Depend, depend
        else:
            return DependIDs  # Returns dict that gives the dependency of the specified patchID to the remaining ones
                                        # For example, for patchID = 0 it gives:
                                        #    connect = { ( 0, 1 ): [ [ local indices of cps_patchID ], [ corresponding indices of global cps array ] ]
                                        #            ]  }   
    
    def update_cps(self,*patchIDs):
        # update the patches in the specified hierarchical order
        for ipatch in patchIDs: 
            connect = self.dependencies(ipatch)
            for ikey, ival in connect.items():
                self.cps()[ival[1]] =  self.cps(ikey[0])[ival[0]]
        return
    
    def update_w(self,*patchIDs):
        # update the patches in the specified hierarchical order
        for ipatch in patchIDs: 
            connect = self.dependencies(ipatch)
            for ikey, ival in connect.items():
                self.weights()[ival[1]] =  self.weights(ikey[0])[ival[0]]
        return
    
    def update_patchverts(self,):

        if self._geom_type == 'left-ventricle': # Only works for Left ventricle
            patches = 0,2 # Patches which we loop over to get the patch vertices
            bounds = {patches[0]:[ (0,1), (0,2) ], patches[1]:[ (6,7), (6,4) ] } # Boundaries which we will use to determine the number of cps positioned on them
            ind_vert = [[0,1,2,3],[6,7,4,5]]
            
            for i, (ikey, ival) in enumerate( bounds.items() ):
                
                ncps    = ( self.bnelems()[ival[0]] + 2 , self.bnelems()[ival[1]] + 2 )
                cps     = self.cps(ikey).reshape(ncps[0],ncps[1],3)

                self.patchverts()[ind_vert[i]] = numpy.array([ cps[ 0, 0], 
                                                               cps[ 0,-1],
                                                               cps[-1, 0],
                                                               cps[-1,-1] ])
                
            # bounds   = [ (0,1), (0,2) ] # Local boundaries of the patch we will be using
            # ncps_mid = ( self.bnelems()[bounds[0]] , self.bnelems()[bounds[1]] ) # nr of cps in between the boundary cps (quadratic spline)
            # ncps     = (ncps_mid[0]+2)*(ncps_mid[1]+2) # Total nr of cps in that patch
            
            
            # b     = (0,1,2,3)
            # b[1] += ncps_mid
            # b[2] += ncps_mid 
            # for ibound in bounds:
            #     ind_cps_01 = [i for i in range(bounds[0][0], bound[1][1] + ncps_mid[0])]
            #     ind_cps_02  = [i for i in range(bound[0][1], bound[1][1] + ncps_mid[1] + )]
            
            
            # for ikey, ival in bounds.items():
            #     ncps_mid = ( self.bnelems()[ival[0]] , self.bnelems()[ival[1]] ) # nr of cps in between the boundary cps (quadratic spline)
            #     for b, bound in enumerate(ival): 
            #         if b == 1:
            #            ncps  = (ncps_mid[0]+2)*(ncps_mid[1]+2) # Total nr of cps in that patch
            #            ind_cps  = [i for i in range(bound[0], bound[1] + ncps_mid[b] + )]   
            #         ind_cps  = [i for i in range(bound[0],bound[1]+ncps_mid)] 
            #         self.patchverts()[bound] = self.cps(ikey)[indiced_totalcps]
            
            
            
            # for ipatch in patches:
            #     connect_cps, connect_vert, connect_con = self.dependencies(ipatch,add_return=True)
            #     # We are only interested in the connection of patch 0 with 1 and 3 and patch 2 with  and 3
            #     for jpatch in (1,3):
            #         indices_patchvert = connect_con[(ipatch,jpatch)]   # Indices that correspond to specific values in the patchverts array
            #         indiced_totalcps  = connect_vert[(ipatch,jpatch)][0] # The coupled indices of indices_patchvert, which refers to the entire array of the cps
            #         self.patchverts()[indices_patchvert] =  self.cps()[indiced_totalcps] # Ensure the patchvertices coincide with the values inside the total cps array
                
        else:
            raise NotImplementedError("Only works for left-ventricle")
        return
    
    def update(self,*patchIDs): # Update the cps, w and patchverts arrays
        self.update_cps(*patchIDs)
        self.update_w(*patchIDs)
        self.update_patchverts()
        # To be implemented
        # self.update_knotval()
        # self.update_knotmult()
        return
    
    
    def refine(self,nrefine): # Refine entire multipatch surface
        # Convert to Splipy first
        CPS = numpy.array([[0,0,0]])
        W   = numpy.array([0])
        knotval_new  = {}
        knotmult_new = {}
        nelems_new   = {}
        nrcps  = []
        for ipatch in range(self.nrpatches()):

            cps     = self.cps(ipatch)
            w       = self.weights(ipatch)
            knotvec = self.knotvec(ipatch) # <-- This cannot be modified, only knotval() and knotmult()

            #print(self.knotval(ipatch),self.knotmult(ipatch))
            control_net    = numpy.concatenate([ cps*w[:,numpy.newaxis], w[:,numpy.newaxis] ], axis = 1) # Homogenous coordinates
            spline_order   = 2
            basis_u        = BSplineBasis(order=spline_order+1, knots=knotvec[0])
            basis_v        = BSplineBasis(order=spline_order+1, knots=knotvec[1])
            Splipy_surface = Surface(basis_u, basis_v, control_net, rational=True)
            Splipy_surface.refine(nrefine)
            
            # Extract cps and weights
            weight_new  = Splipy_surface[:].reshape(-1,4)[:,-1]
            cps_new     = Splipy_surface[:].reshape(-1,4)[:,:3]/weight_new[:,numpy.newaxis]
            
            CPS = numpy.concatenate([ CPS, cps_new ], axis=0)
            W   = numpy.concatenate([ W,   weight_new ], axis=0)
  
            nrcps += [len(cps_new)]
            

            # Update the knotvector
            knotvec_u = list(Splipy_surface.bases[0])
            knotvec_v = list(Splipy_surface.bases[1])
            
            ndims = min( len(self.patchcon(ipatch)) // 2, 3)
            if ndims == 3: raise NotImplementedError('Refine does not support solids')
            bound = self.patchcon(ipatch).reshape([2]*ndims)
            for i in range(ndims):
                knotval_u, knotmult_u = numpy.unique(knotvec_u, return_counts=True)
                knotval_v, knotmult_v = numpy.unique(knotvec_v, return_counts=True)
                
                knotval_new[tuple(bound[i,:])]  = knotval_u.tolist()
                knotval_new[tuple(bound[:,i])]  = knotval_v.tolist()
                
                knotmult_new[tuple(bound[i,:])] = knotmult_u.tolist()
                knotmult_new[tuple(bound[:,i])] = knotmult_v.tolist()
                
                nelems_new[tuple(bound[i,:])]  = len(knotval_u) - 1
                nelems_new[tuple(bound[:,i])]  = len(knotval_v) - 1

        # Update the knotval and knotmult arrays
        self._w   = W[1:] # remove first one sine it is empty
        self._cps = CPS[1:]
        self._knotmult = knotmult_new
        self._knotval  = knotval_new
        self._nelems   = nelems_new
        for ipatch in range(self.nrpatches()): 
            self._nrcps[ipatch] = nrcps[ipatch]
        return
    
    @staticmethod
    def _convert_boundary_names(boundariesdic):
        boundaries    = {}    
        return {key:','.join([boundaries[key]]+['patch{}-{}'.format(ipatch, bname) for ipatch,bname in val] if key in boundaries else ['patch{}-{}'.format(*v) for v in val]) for key, val in boundariesdic.items()}

    ## UPDATE PATCHVERTS WHEN CPS ARE UPDATED!
    
    # def total_nelems(self,): # Return total nr of elements of the surface
    #     total_nelems = 0
    #     for i in self._nelems.items():
    #         total_nelems += i[0]*i[1]
    #     return total_nelems
    
    # def convert_to_splipy(self,): # converts the cps/w structure to a splipy surface (easier to work with?)
    #     if self.input_type == 'splipy':
    #         print("Already given in SpliPy format")
    #         return
    #     return
    # def convert_to_nutils(self,): # converts the cps/w structure to a nutils input data structure (easier to work with?)
    #     if self.input_type == 'nutils':
    #         print("Already given in Nutils format")
    #         return
    #     return    
    def rotate(self, degr, axis='z', degrees=True):
        R = scipy.spatial.transform.Rotation
        RotM = R.from_euler(axis, degr, degrees=degrees).as_matrix()
        self._patchverts = (self._patchverts @ RotM)   # patchconnectivity array of the multipatch surface
        self._cps        = (self._cps @ RotM)   # cntrol points of the multipatch surface
        return

    def translate(self, x):
        # translate the cps and patchverts by a given vector
        assert type(x) == numpy.ndarray, f"Provided translation vector is of the wrong type {type(x)}, should be a numpy array"
        #assert x.shape
        self._patchverts += x.astype( self._patchverts.dtype )   # patchconnectivity array of the multipatch surface
        self._cps        += x.astype( self._cps.dtype )   # cntrol points of the multipatch surface
        return

    def scale(self, sfactor, direction='x'):
        # translate the cps and patchverts by a given vector
        #assert type(x) == numpy.ndarray, f"Provided translation vector is of the wrong type {type(x)}, should be a numpy array"
        #assert x.shape
        

        # if direction == 'x':
        #     Svec = numpy.array([sfactor, 1, 1 ]) 
        # elif direction == 'y':
        #     Svec = numpy.array([1, sfactor, 1 ])
        # elif direction == 'z':
        #     Svec = numpy.array([1, 1, sfactor ])
        # else:
        #     raise ValueError(f"Unknown direction specified for stretching '{direction}'")

        Svec = numpy.array([sfactor if ax in direction else 1 for ax in "xyz"])
        self._patchverts *= Svec.astype( self._patchverts.dtype )   # patchconnectivity array of the multipatch surface
        self._cps        *= Svec.astype( self._cps.dtype )   # cntrol points of the multipatch surface
        return

    def get_base_center(self,):
        topo, geom = self.get_topo_geom_alt(boundary_names=self._bnames)
        
        ns     = function.Namespace()
        base_coords  = topo.boundary['base'].sample('bezier',10).eval(geom)# Difference with z=0 
        return numpy.mean(base_coords, axis=0)
    

    def get_apex_point(self,):
        # Only for left ventricle multipatch surfaces
        topo, geom = self.get_topo_geom_alt()
        sample_points =  topo.sample('bezier', 40).eval(geom)
        return sample_points[ numpy.argmin( sample_points[:,-1] ) ] # return point with lowest value = apex point

    def get_topo_geom_alt(self, boundary_names: dict = {}):
        # Check if biventricle    
        biventr = False if self._geom_type == 'left-ventricle' else True  
        
        topo, lingeom = mesh.multipatch(patches=self.patchcon(), patchverts=self.patchverts(), nelems=self.bnelems(biventr=biventr))

        if biventr:
            bsplines = topo.basis('spline', degree=2, patchcontinuous=False)
        else:
            bsplines = topo.basis('spline', degree=2, patchcontinuous=False, knotvalues=self.knotval(), knotmultiplicities=self.knotmult())

        weight   = bsplines.dot(self.weights())
        geom     = bsplines.vector(3).dot((self.cps()*self.weights()[:,numpy.newaxis]).ravel())/weight
        
        # Assign boundary names
        if boundary_names:
            topo       = topo.withboundary(**boundary_names)
            
        #print('Value should be zero (or close to): Value=', abs(topo.interfaces['interpatch'].sample('bezier', 5).eval(function.jump(geom))).max()) # New nutils version does not work with this
        print("Number of CPS: {}".format(len(topo.basis('spline', degree=2))))
        print("Number of Elements: {}".format(len(topo.integrate_elementwise(geom, degree=0)) ))
        return topo, geom

    ## A useful function if you want to find a local boundary its connectivity/name
    def print_localbound_names(self,): # function which prints a list of the boundary names and their connectivity. Can be useful  
        for i, j in enumerate(self.patchcon()):
            self.__findboundary(j, i)    
        return   
    def __findboundary(self, patch, patchID=None):
        ndims = 0
        while 2**ndims < len(patch):#patch.shape[1]:
            ndims += 1
        patch_boundaries = numpy.reshape(patch, [2]*ndims)

        if ndims == 1:
            names = ['left','right']
        elif ndims == 2:
            names = ['left','bottom','right','top']
        elif ndims == 3:
            names = ['left','bottom','back','right','top','front']
            
        boundaries = []
        for i in range(2):
            boundaries += [patch_boundaries[i,:]]
            boundaries += [patch_boundaries[:,i]]
            if ndims == 3:
                boundaries += [patch_boundaries[:,:,i]]            
        ## Print results
        if patchID != None:
            print('Patch boundary indices (Patch{})'.format(patchID))
        else:
            print('Patch specific boundary data')
        format_row = "{:>12} : {}"
        for name, indices in zip(names, boundaries):
            ind = [indices if ndims<3 else numpy.concatenate([indices[0], indices[1]], axis=0)]
            print(format_row.format(name, ind[0] ) ) 
        return  
    
    def __conc_conn(self, x): # concatenates the columns
        conn = numpy.ndarray((0,2),int)
        for col in range(x.shape[1]-1):
            conn = numpy.concatenate([ conn, numpy.concatenate([x[:,col,_], x[:,col+1,_]],axis=1)], axis=0)
        return conn
    
    def save_controlnet(self, filename : str = 'filename'):
        '''save_controlnet
        
        Save the control net as a 'filename.vtk' file. Usefull for visualization.

        ''' 
        
        CONN = numpy.ndarray((0,2),int)
        #conf = lambda x : numpy.concatenate([ numpy.concatenate([x[:,0,_], x[:,1,_]],axis=1), numpy.concatenate([x[:,1,_], x[:,2,_]],axis=1)], axis=0)
        assert self.patchcon().shape[1] == 4, "Save control net only works for surfaces"
        id_prev = 0
        for patch in self.patchcon():
            nrcps_u = self.bnelems()[tuple(patch[[0,1]])] + 2
            nrcps_v = self.bnelems()[tuple(patch[[0,2]])] + 2
            id_u, id_v = numpy.meshgrid( range(nrcps_u), range(nrcps_v) )

            id_cps   = numpy.linspace(0,nrcps_u*nrcps_v-1,nrcps_u*nrcps_v).astype(int).reshape(nrcps_u, nrcps_v) + id_prev
            row_conn = numpy.concatenate([self.__conc_conn(id_u)[...,_],   self.__conc_conn(id_v)[...,_]],axis=2)
            col_conn = numpy.concatenate([self.__conc_conn(id_u.T)[...,_], self.__conc_conn(id_v.T)[...,_]],axis=2)

            # row_conn = numpy.concatenate([conf(id_u)[...,_],   conf(id_v)[...,_]],axis=2)
            # col_conn = numpy.concatenate([conf(id_u.T)[...,_], conf(id_v.T)[...,_]],axis=2)


            conn = numpy.concatenate([row_conn, col_conn], axis=0)
            for c in conn:
                c = c.T
                CONN = numpy.concatenate([CONN, id_cps[c[0],c[1]][_]],axis=0)
            
            id_prev += nrcps_u*nrcps_v

        export.vtk( filename, CONN, self.cps())
        return

    def save_vtk(self, filename : str = 'filename', patch = -1, boundary_only: bool = False, nrsamples: int = 15, boundaries: tuple = (), boundary_names: dict = {}):
        '''save_vtk
        
        Construct the Nutils NURBS topology and geometry and export the geometry in a 'filename.vtk' file.
        
        - boundary_only  : Saves all the boundaries of the geometry, not the inside parts
        - boundaries     : Tuple of specific boundary names which are to be saved (are saved in seperate file)
        - boundary_names : Dictionary that renames the standard boundary names in Nutils format
        ''' 
        if type(patch) == int:
            patch = [patch]
          
            
        # # Check if biventricle    
        # biventr = False if self._geom_type == 'left-ventricle' else True  
        
        # topo, lingeom = mesh.multipatch(patches=self.patchcon(), patchverts=self.patchverts(), nelems=self.bnelems(biventr=biventr))

        # if biventr:
        #     bsplines = topo.basis('spline', degree=2, patchcontinuous=False)
        # else:
        #     bsplines = topo.basis('spline', degree=2, patchcontinuous=False, knotvalues=self.knotval(), knotmultiplicities=self.knotmult())

        # weight   = bsplines.dot(self.weights())
        # geom     = bsplines.vector(3).dot((self.cps()*self.weights()[:,numpy.newaxis]).ravel())/weight
        
        # print('Value should be zero (or close to): Value=', abs(topo.interfaces['interpatch'].sample('bezier', 5).eval(function.jump(geom))).max())
        # print("Number of CPS: {}".format(len(topo.basis('spline', degree=2))))
        topo, geom = self.get_topo_geom_alt(boundary_names)
        
        if patch[0] == -1:
            if len(boundaries) != 0:
                for bound in boundaries: # Save specific boundaries
                    bezier     = topo.boundary[bound].sample('bezier', nrsamples)
                    pids       = topo.basis('patch').dot(numpy.arange(len(self.patchcon())))
                    GEOM, PIDS = bezier.eval([geom,pids])
                    Filename   = filename + '_' + bound
                    export.vtk( Filename, bezier.tri, GEOM, PatchID=PIDS)
            else:       
                bezier     = topo.sample('bezier', nrsamples) if not boundary_only else topo.boundary.sample('bezier', nrsamples)
                pids       = topo.basis('patch').dot(numpy.arange(len(self.patchcon())))
                GEOM, PIDS = bezier.eval([geom,pids])
                export.vtk( filename, bezier.tri, GEOM, PatchID=PIDS) 
        else:
            if len(boundaries) != 0:
                NotImplementedError("Saving boundary of specific patch")
            patch_str  = ','.join(['patch{}'.format(ip) for ip in patch])
            bezier     = topo[patch_str].sample('bezier', nrsamples) if not boundary_only else topo[patch_str].boundary.sample('bezier', nrsamples)
            GEOM       = bezier.eval(geom)
            PIDS       = numpy.ones(len(GEOM))
            export.vtk( filename, bezier.tri, GEOM, PatchID=PIDS) 
        
        # if filename!=None:
        #     export.vtk( filename, bezier.tri, GEOM, PatchID=PIDS) 
        
    def save_stl(self, filename: str = 'filename', nrsamples: int = 10, boundaries: tuple = (), boundary_names: dict = {}):
        # TODO Improve        
        topo, geom = self.get_topo_geom_alt(boundary_names)
        

        #bezier  = topo.boundary['base_l,epi,endo_l'].sample('bezier', nrsamples) # Connect
        #bezier  = topo.boundary['base_l,epi,endo_l,base_r,endo_r'].sample('bezier', nrsamples) # Connect
        bezier  = topo.boundary.sample('bezier', nrsamples)
        GEOM    = bezier.eval(geom)                # Coords

        # Re-arrange normals
        bezier_tri = self.get_bezier_tri(bezier)  

        data = numpy.zeros(len(bezier.tri), dtype=stl.mesh.Mesh.dtype)
        stl_mesh = stl.mesh.Mesh(data, remove_empty_areas=False)
        stl_mesh.x[:] = GEOM[:,0][bezier_tri]
        stl_mesh.y[:] = GEOM[:,1][bezier_tri]
        stl_mesh.z[:] = GEOM[:,2][bezier_tri]
        stl_mesh.save(filename + '.stl')
        return

    @staticmethod
    def get_bezier_tri(bezier):
        bezier_tri = bezier.tri.copy()
        if bezier.tri.shape[1] == 3:
            bezier_tri[1::2,1] = bezier.tri[1::2,2] 
            bezier_tri[1::2,2] = bezier.tri[1::2,1]
        return bezier_tri


    def save_mesh_vtk(self,filename: str = 'mesh_lines', nsample : int = 20, nrefine : int = 0):
      # get topo and geometry
      topo, geom = self.get_topo_geom_alt()
      if nrefine: topo = topo.refine(nrefine)
      
      meshrefs  = []
      meshtrans = []
    
      # loop over elements
      for elemtrans, elemref in zip(topo.transforms, topo.references):
        # loop over faces
        for facetrans, faceref in elemref.edges:
          # loop over edges
          for edgetrans, edgeref in faceref.edges:
            meshrefs .append(edgeref)
            meshtrans.append(elemtrans + (facetrans,) + (edgetrans,))
    
      meshrefs  = elementseq.References.from_iter(meshrefs, 1)
      meshtrans = transformseq.PlainTransforms(meshtrans, fromdims=1, todims=topo.ndims)
      space     = topo.space
      mesh_topo = topology.TransformChainsTopology(space=space, references=meshrefs, transforms=meshtrans, opposites=meshtrans)
      
      mbezier = mesh_topo.sample('bezier', nsample)
      mpoints = mbezier.eval(geom)
      
      export.vtk(filename, mbezier.tri, mpoints)
      return 



    
    def save_pickle(self, filename : str = 'filename', boundary : dict = {}):
        '''save_pickle

        Save a new 'filename.pickle' file.
        '''        
        data = [ self._cps, self._w, self._patchverts, self._patchcon, self._bnelems, self.knotval(), self.knotmult()]
        if boundary:
            data += [boundary]
            
        with open(filename + '.pickle', 'wb') as f: 
             pickle.dump(data,f) 
        return
    
    
    def load_pickle(self, filename : str = 'filename'):
        '''load_pickle

        Load an existing 'filename.pickle' file and output its results.
        '''
        with open(filename + '.pickle', 'rb') as f: 
          return pickle.load(f) 

    def save_json(self, filename : str = 'filename', boundary : dict = {}):
        '''save_json

        Saves the geometrical data to a json file.
        '''
        import json

        if filename.split(".")[-1] != "json":
            filename += ".json"

        # Convert np.arrays to lists
        cps = self._cps.tolist()
        w   = self._w.tolist()
        patches    = self._patchcon.tolist()
        patchverts = self._patchverts.tolist()
        
        # Assign everything to one dictionary
        Data = {}
        Data["cps"] = cps
        Data["w"]   = w
        Data["patches"]    = patches
        Data["patchverts"] = patchverts
        Data["nelems"]     = self.__remap_keys(self._bnelems)
        Data["knotval"]    = self.__remap_keys(self.knotval())
        Data["knotmult"]   = self.__remap_keys(self.knotmult())
        if boundary:
            Data["boundaries"] = boundary
        json_object = json.dumps(Data, indent=4)

        with open(filename + ".json", "w") as outfile:
            outfile.write(json_object)
        return

    def load_json(self, filename : str = 'filename', asDict=False): # Return as dictionary True or False
        '''load_json

        Loads the geometrical data from a json file.

        asDict : Bool;  True if user desired the output to be a dictionary, False for individual components
        '''
        import json
        if filename.split(".")[-1] != "json":
            filename += ".json"
        f = open(filename)

        # returns JSON object as 
        # a dictionary
        data = json.load(f)

        # Convert lists to usable np.arrays
        data["cps"] = numpy.array(data["cps"])
        data["w"]   = numpy.array(data["w"])
        data["patchverts"] = numpy.array(data["patchverts"])
        data["patches"] = numpy.array(data["patches"]) 
        
        # Certain values had to be remapped because JSON did not support tuples as keys
        data["nelems"]   = self.__remap_keys_reverse(data["nelems"])
        data["knotval"]  = self.__remap_keys_reverse(data["knotval"])
        data["knotmult"] = self.__remap_keys_reverse(data["knotmult"])

        if asDict:
            return data
        else:
            cps = data["cps"]
            w   = data["w"]
            patchverts = data["patchverts"]
            patches = data["patches"]
            nelems  = data["nelems"]
            knotval = data["knotval"]
            knotmult   = data["knotmult"]
            boundaries = data["boundaries"]
            return cps, w, patchverts, patches, nelems, knotval, knotmult, boundaries
    
    def save_txt(self, filename : str = 'filename', boundary : dict = {}, headerName : str = "Ventricle"):
    
        if filename.split(".")[-1] != "txt":
            filename += ".txt"

        Title   = f"{headerName} multipatch geometry data"
        Data    = [self._cps, self._w, self._patchverts, self._patchcon, self._bnelems, self.knotval(), self.knotmult()]
        
        Headers = ["Control points", "Weights", "Patch vertices", "Patch connectivity", \
                   "Number of elements per boundary", "Knot values per boundary", \
                   "Knot multiplicity per boundary"]

        if boundary:
            Data += [boundary]
            Headers += ["Boundary names"]

        with open(filename, 'w') as f:
            f.write(Title+"\n\n")
            for head, data in zip(Headers,Data):
                f.write(head+"\n")
                if type(data) == numpy.ndarray:
                    lines = "\n".join( [ str(row.tolist()) for row in data ] )
                else:
                    lines = ""
                f.write(lines) 
        
                if type(data) == dict:
                    for key, value in data.items(): 
                        f.write('%s:%s\n' % (key, value))
                f.write("\n\n") 
        return

    @staticmethod
    def load_txt(filename : str = 'filename', returnDict=False):
        import ast

        if filename.split(".")[-1] != "txt":
            filename += ".txt"

        with open(filename, 'r') as f:
            lines = f.readlines()

        Headers = "Control points", "Weights", "Patch vertices", "Patch connectivity", \
                "Number of elements per boundary", "Knot values per boundary", \
                "Knot multiplicity per boundary", "Boundary names" # Ordering is of importance! Should be same as in the save_txt() function
        DictKeys = "control points", "weights", "patch vertices", "patch connectivity", \
                "nelems", "knot values", "knot multiplicity", "boundary names"
        
        # Strip empty spaces and remove '\n'
        lines = [idat.strip("\n") for idat in lines if len(idat.strip("\n")) != 0]   
        catch_line = "View"
        idx = []   
        for i, line in enumerate(lines):
            if line in Headers:
                idx += [i]
        idx += [None]
        loadedData = []
        for k, (i,j) in enumerate(zip(idx[:-1],idx[1:])):
            if k < 4: # We encounter np.arrays
                loadedData += [numpy.array([ast.literal_eval(iline) for iline in lines[(i+1):j]])]  
            # elif k == 3: # We have the weights list (special case)
            #     loadedData += [numpy.array([ast.literal_eval(', '.join(lines[(i+1):j]))])]  
            else: # We encounter dicts
                d = {}
                for b in lines[(i+1):j]:
                    i = b.split(':')
                    if k != len(Headers)-1:
                        d[ast.literal_eval(i[0])] = ast.literal_eval(i[1])
                    else: # Else we have the boundaries dict
                        d[i[0]] = i[1]      
                loadedData += [d]

        if returnDict:
            return {key : value for key, value in zip(DictKeys, loadedData)}
        else:
            return loadedData

    @staticmethod
    def __remap_keys(mapping): # Remapping keys reuired if tuple is a key in the initial dict
        return [{'key':[int(i) for i in k ], 'value': v} for k, v in mapping.items()]    

    @staticmethod
    def __remap_keys_reverse(mapping):
        remapped = {}
        for map in mapping:
            k = map["key"]
            v = map["value"]
            remapped[tuple(k)] = v
        return remapped 
    
    
#%     
class MultipatchSolid(MultipatchSurface):
    '''Class:
        Combines surfaces and solids
        
        Combines the surfaces (inner, outer) of the left-ventricle or right-ventricle. 
        The class is also able to combine the solids of the created left- and right-ventricles
        
        Input:
            
            Requires Rlv_x, Rlv_y, Rlv_z, R....    
    ''' 
    def __init__(self, *surfaces, tnelems = (1, 1)):
        if len(surfaces)==1:
            surfaces = surfaces[0]
        
        # create mpty dictionaries
        self.surf   = dict(left=dict(), right=dict())
        # self.surf['left']['outer'] = surfaces[0]
        # self.surf['left']['inner'] = surfaces[1]
        # self.surfLV = {}
        # self.surfRV = {}

    # TODO Filter the input appropriately with correct asserts
        
        for isurf in surfaces:
            
            # Assign surfaces to corresponding variables (callable in other function)
            if isurf.surface() == 'left-inner':
                self.surf['left']['inner'] = isurf
            elif isurf.surface() == 'left-outer':
                self.surf['left']['outer'] = isurf
            elif isurf.surface() == 'right-inner':
                self.surf['right']['inner'] = isurf
            else:
                self.surf['right']['outer'] = isurf
            
            # if len(surfaces) == 2:
            #     #assert isurf.geom_type == 'left-ventricle', "Incorrect number of surfaces (2), for given geometry type: 'bi-ventricle'"
                
            #     # Assign surfaces to corresponding variables (callable in other function)
            #     if isurf.surface == 'left-inner':
            #         self.surf['left']['inner'] = isurf
            #     else:
            #         self.surf['left']['outer'] = isurf
                    
            # elif len(surfaces) == 4:
            #     assert isurf.geom_type == 'bi-ventricle', "Incorrect number of surfaces (4), for given geometry type: 'left-ventricle'"
            #     # Assign surfaces to corresponding variables (callable in other function)
            #     if isurf.surface == 'left-inner':
            #         self.surf['left']['inner'] = isurf
            #     elif isurf.surface == 'left-outer':
            #         self.surf['left']['outer'] = isurf
            #     elif isurf.surface == 'right-inner':
            #         self.surf['right']['inner'] = isurf
            #     else:
            #         self.surf['right']['outer'] = isurf  
            # else:
            #     raise ValueError("Incorrect number of surfaces, only supports: \n 2 for 'left-ventricle' \n 4 for 'bi-ventricle'")


        ## Construct the solid
        if self.surf['left'] and self.surf['right']:
            #self._combine_surfs('left' , tnelems[0], store=True)
            #self._combine_surfs('right' , tnelems[0], store=True)
            LeftVentricle  = self._combine_surfs('left' , tnelems[0], store=False) # combine left-ventricle surfaces
            RightVentricle = self._combine_surfs('right', tnelems[1], store=False) # combine left-ventricle surfaces
            self.combine_solids(LeftVentricle, RightVentricle, tnelems)
        elif self.surf['left']:
            self._combine_surfs('left', tnelems[0]) # combine left-ventricle surfaces
        else:
            self._combine_surfs('right', tnelems[1]) # combine left-ventricle surfaces      
            
            
            
        #self._combine_rv(tnelems[1]) # combine right-ventricle surfaces
        #self._combine()    # combine left- and right- surfaces
        return
    
    def _combine_surfs(self, solid_type, tnelems, store=True): 
        # solid_type = 'left' or 'right'
        assert self.surf[solid_type]['outer'].bnelems() == self.surf[solid_type]['inner'].bnelems(), "Mismatch in nelems of inner and outer {}-ventricle surface".format(solid_type) # <-- This is not allowed!
 
    
        MaxID       = numpy.max(           self.surf[solid_type]['outer'].patchcon()                        )
        bnelems     = self.__comb_bound(   self.surf[solid_type]['outer'].bnelems().copy() , MaxID, tnelems )
        patchcon    = self.__comb_patches( self.surf[solid_type]['outer'].patchcon()       , MaxID          )
        patchverts  = numpy.concatenate([  self.surf[solid_type]['outer'].patchverts()     , self.surf[solid_type]['inner'].patchverts() ])

        # Knotvalues and multiplicities
        knotv = [v/tnelems for v in range(tnelems)] + [1] # knotvector (without multiplicity) per boundary
        knotm = [3] + [1]*(tnelems-1) + [3]
        knotval     = self.__comb_bound(   self.surf[solid_type]['outer'].knotval().copy() , MaxID, knotv )
        knotmult    = self.__comb_bound(   self.surf[solid_type]['outer'].knotmult().copy(), MaxID, knotm )
        
        # combine cps and weight arrays in correct way
        nrcps_pp = [0] + [j for i, j in self.surf[solid_type]['outer'].nrcps().items() ] # number of cps per surface patch
        nrcps    = numpy.cumsum(nrcps_pp)
        cps      = numpy.zeros((nrcps[-1]*(tnelems+2),3))
        w        = numpy.zeros(nrcps[-1]*(tnelems+2)) # initialize arrays
         

        w_inside   = lambda ipatch: self.__interpolate_values( self.surf[solid_type]['outer'].weights(ipatch), self.surf[solid_type]['inner'].weights(ipatch), tnelems) # determine the weights of the cps inside the wall
        cps_inside = lambda ipatch: self.__interpolate_values( self.surf[solid_type]['outer'].cps(ipatch)    , self.surf[solid_type]['inner'].cps(ipatch)    , tnelems) # determine the cps inside the wall
        #pnelems    = 
        
        for ipatch in range( len( self.surf[solid_type]['outer'].patchcon() ) ) :# loop over each patch
            i_start = nrcps[ipatch]*(tnelems+2)   # start index of ipatch
            i_end   = nrcps[ipatch+1]*(tnelems+2) # end   index of ipatch 
            cps[i_start:i_end, :] = numpy.concatenate([ self.surf[solid_type]['outer'].cps(ipatch)      ] + cps_inside(ipatch) + [ self.surf[solid_type]['inner'].cps(ipatch)     ] ) 
            w[i_start:i_end]      = numpy.concatenate([ self.surf[solid_type]['outer'].weights(ipatch)  ] + w_inside(ipatch)   + [ self.surf[solid_type]['inner'].weights(ipatch) ] ) 
        
        if store:
            self._bnelems    = bnelems
            self._patchverts = patchverts
            self._patchcon   = patchcon
            self._cps   = cps
            self._w     = w
            self._nrcps = {i: j*(tnelems+2) for i, j in enumerate(nrcps_pp[1:])} 
            
            self._knotmult  = knotmult
            self._knotval   = knotval
            self._geom_type = 'left-ventricle'
            return
        else:
            nrcps_dic = {i: j*(tnelems+2) for i, j in enumerate(nrcps_pp[1:])} 
            return dict(solid='{}-ventricle'.format(solid_type), cps=cps, w=w, nrcps=nrcps_dic, patchverts=patchverts, patchcon=patchcon, bnelems=bnelems, knotmult=knotmult, knotval=knotval)
     
    def combine_solids(self, Solid_a, Solid_b, tnelems): # Specifically combining left- and right-ventricles
        # Assign correct names
        LV = Solid_a if 'left-ventricle'  == Solid_a['solid'] else Solid_b
        RV = Solid_a if 'right-ventricle' == Solid_a['solid'] else Solid_b
        
        # 1) Set the patchconnectivity array
        #changeIDX      = {'from': [0,1,6,7], 'to': [6,5,13,12,13]} # Change specific indices to a different value
        #MaxID          = numpy.max( LV['patchcon'] )
        #RV['patchcon'] += MaxID - 1
        Rvpatchcon = numpy.array([ [ 4, 5,32,33, 6, 7,36,37],
                                   [34,35,32,33,38,39,36,37],
                                   [14,15,34,35,12,13,38,39],
                                   [15, 5,35,33,13, 7,39,37]]) # New patchconnectivity of right-ventricle
        # patch_right = numpy.array([[ 4, 5,24,25, 6, 7,28,29],
        #                            [26,27,24,25,30,31,28,29],
        #                            [10,11,26,27, 8, 9,30,31],
        #                            [11, 5,27,25, 9, 7,31,29]])
        
        # Rvpatchcon = numpy.array([ [ 4, 5, 6, 7, 32, 33, 36, 37],
        #                            [34,35,32,33,38,39,36,37],
        #                            #[34,35,38,39, 32, 33, 36, 37],
        #                            [14,15,12,13, 34, 35, 38, 39],
        #                            [15, 5,13, 7, 35, 33, 39, 37]])
                                   #[5,7,33,37,15,13,35,39]])
                                   #[15, 5, 35,33,13,7,39,37]])
                                   #[15, 5,13, 7, 35, 33, 39, 37]]) # New patchconnectivity of right-ventricle
        
                
        
        #self.change_idx(LV['patchcon'], changeID) # change the indices value
        patchcon = numpy.concatenate([ LV['patchcon'], Rvpatchcon], axis=0)
        
        # 2) Set the patchverts array 
        changeIDX      = {'from': [0,1,6,7,8,9,14,15], 'to': [4,5,14,15,6,7,12,13]} # Change specific indices to a different value
                                    # right-ventricle           # left-ventricle
        patchverts = numpy.concatenate([ LV['patchverts'], numpy.delete(RV['patchverts'], changeIDX['from'], axis=0) ], axis=0)
        
        
        
        # 3) Set bnelems array
        # changeIDX      = {'from': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15], 
        #                   'to'  : [4, 5,32,33,34,35,14,15, 6, 7,36,37,38,39,12,13]}
        # RVbnelems = dict()
        # for key, item in RV['bnelems'].items():
        #     newkey    = []
        #     for i, bkey in enumerate(key):
        #         newkey += [changeIDX['to'][bkey]]    
        #     RVbnelems[tuple(newkey)] = item
        # RVbnelems.update(LV['bnelems'])
        
        RVbnelems = self.__comb_solid_dict(LV['bnelems'], RV['bnelems'])
        knotmult  = self.__comb_solid_dict(LV['knotmult'], RV['knotmult'])
        knotval   = self.__comb_solid_dict(LV['knotval'], RV['knotval'])

        # 4) Set the nrcps
        nr_lvpatches = LV['patchcon'].shape[0]
        nrcps = {(key+nr_lvpatches): item for key, item in RV['nrcps'].items()}  # simple update (key+maxPATCHID)
        nrcps.update(LV['nrcps'])        

        # 5) Set the cps and w arrays
        cps = numpy.concatenate([LV['cps'], RV['cps']],axis=0)# Simple concatenate
        w   = numpy.concatenate([LV['w'],   RV['w']]  ,axis=0)# Simple concatenate
        
        nrcps_pp    = [0] + [j for i, j in nrcps.items() ] # number of cps per volume patch
        nrcps_pp_cs = numpy.cumsum(nrcps_pp)

        for patchright, patchleft, bound in zip( (11,13,14), (2,6,9), ( (4,5), (14,15), (15,5) )  ):
            nrcps_surf_right = int( nrcps[patchright] / (tnelems[1]+2)) 
            nrcps_surf_left  = int( nrcps[patchleft] / (tnelems[0]+2))
            nrcps_long       = RVbnelems[bound] + 2
            for i in range(1,tnelems[1]+1):
                start_idx_right = nrcps_pp_cs[patchright] + nrcps_surf_right*i
                end_idx_right   = start_idx_right + nrcps_long
                start_idx_left = nrcps_pp_cs[patchleft] + nrcps_long*i
                end_idx_left   = start_idx_left + nrcps_long
                
                cps[start_idx_right:end_idx_right] = cps[start_idx_left:end_idx_left] 
                w[start_idx_right:end_idx_right]   = w[start_idx_left:end_idx_left]
            
        idx = 1# match the intersection cps and w, use the LV data as main
        
        


        # Store the values
        self._bnelems    = RVbnelems
        self._patchverts = patchverts
        self._patchcon   = patchcon
        self._cps        = cps
        self._w          = w
        self._nrcps      = nrcps 
        self._geom_type  = 'bi-ventricle'    
            
        self._knotmult  = knotmult
        self._knotval   = knotval
        return
        
     
        # Bi-ventricle:
            # Start with right-ventricle surfRV...
            # Compute new knots based on t_nelems in RV
            # Construct the left-ventricle
       # return
    
    # def __combine(self,): # c
    #     ''' Combine the inner and outer surfaces and store the relevant data again.
    #     '''    
    #     self._degree  = self.surfaces['outer'].degree()
        
    #     # We combine the inner and outer surface, by moving from the outer -> inner surface
    #     # One could move from inner -> outer, but we have to chose a direction to move in (changes the code a bit) 
    #     MaxID         = numpy.max(self.surfaces['outer'].patches())
    #     self._nelems  = self.__comb_bound(self.surfaces['outer'].nelems(),MaxID, self._tnelems)
    #     self._patches = self.__comb_patches(self.surfaces['outer'].patches(),MaxID)
    #     self._patchverts = numpy.concatenate([ self.surfaces['outer'].patchverts(), self.surfaces['inner'].patchverts() ])

    #     # Thickness direction knotvectors (only supports quadratic)
    #     knotv_t = [ i/self._tnelems for i in range(0,self._tnelems + 1) ] # knot vector values
    #     knotm_t = [3] + [1]*(len(knotv_t)-2) + [3] # knot multiplicity values
        
    #     self._knotvalues       = self.__comb_bound( self.surfaces['outer'].knotval() ,  MaxID,  knotv_t)
    #     self._knotmultiplicity = self.__comb_bound( self.surfaces['outer'].knotmult(),  MaxID,  knotm_t)
        
    #     # combine cps and weight arrays in correct way
    #     nrcps_pp = [0] + [j for i, j in self.surfaces['outer'].nrcps().items() ] # number of cps per surface patch
    #     nrcps    = numpy.cumsum(nrcps_pp)
    #     cps      = numpy.zeros((nrcps[-1]*(self._tnelems+2),3))
    #     w        = numpy.zeros(nrcps[-1]*(self._tnelems+2)) # initialize arrays
         

    #     w_inside   = lambda ipatch: self.__interpolate_values(self.surfaces['outer'].weights(ipatch), self.surfaces['inner'].weights(ipatch), self._tnelems) # determine the weights of the cps inside the wall
    #     cps_inside = lambda ipatch: self.__interpolate_values(self.surfaces['outer'].cps(ipatch)    , self.surfaces['inner'].cps(ipatch)    , self._tnelems) # determine the cps inside the wall
        
    #     for ipatch in range( len( self.surfaces['outer'].patches() ) ) :# loop over each patch
    #         i_start = nrcps[ipatch]*(self._tnelems+2)
    #         i_end   = nrcps[ipatch+1]*(self._tnelems+2)
    #         cps[i_start:i_end, :] = numpy.concatenate([ self.surfaces['outer'].cps(ipatch)      ] + cps_inside(ipatch) + [ self.surfaces['inner'].cps(ipatch)     ] ) 
    #         w[i_start:i_end]      = numpy.concatenate([ self.surfaces['outer'].weights(ipatch)  ] + w_inside(ipatch)   + [ self.surfaces['inner'].weights(ipatch) ] ) 
    #     self._cps   = cps
    #     self._w     = w
    #     self._nrcps = {i: j*(self._tnelems+2) for i, j in enumerate(nrcps_pp[1:])} 
    #     return

    def __comb_solid_dict(self,LVdict, RVdict):
        '''Return the combined dictionary of the two solids: left ventricle and right ventricle.
    
        Parameters
        ----------
        LVdict  : dict :
            The dictionary that contains patch-boundary specific data (of the left ventricle multipatch solid).
        RVdict  : dict :
            The dictionary that contains patch-boundary specific data (of the right ventricle multipatch solid).
        Returns
        -------
        : dict :
            Updated the dictionary for the entire multipatch solid geometry (biventricle)
        '''         
        changeIDX      = {'from': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15], 
                          'to'  : [4, 5,32,33,34,35,14,15, 6, 7,36,37,38,39,12,13]} # This dict is specific for the multipatch topology used
        NEWdict = dict()
        for key, item in RVdict.items():
            newkey    = []
            for i, bkey in enumerate(key):
                newkey += [changeIDX['to'][bkey]]    
            NEWdict[tuple(newkey)] = item
        NEWdict.update(LVdict)
        return NEWdict.copy()
    
    def __comb_patches(self,patches, MaxID):
        '''Return the combined multipatch connectivity array.
    
        Parameters
        ----------
        patches  : numpy.array :
            The connectivity array of (can be either) the inner or outer surface.
        MaxID  : int :
            Specify the maximum index present inside the patches array.
        Returns
        -------
        : numpy.array :
            New connectivity array of the solid multipatch geometry
        '''        
        return numpy.concatenate([ patches, patches + MaxID + 1 ], axis=1)
    
    def __comb_bound(self,bound_dic, MaxID, add_item):
        '''Return the combined boundary dictionary.
    
        Parameters
        ----------
        bound_dic  : dict :
            The boundary dictionary of (can be either) the inner or outer surface.
        MaxID  : int :
            Specify the maximum index present inside the patches array.
        add_item   : ... :
            Item that is to be added along the boundary in thickness direction
        Returns
        -------
        : dict :
            Updated boundary dictionary for the entire multipatch geometry
        '''         
        for keys in list(bound_dic):
            newkey = numpy.asarray(keys) + MaxID  + 1        # Capture old keys and increment them
            bound_dic[tuple(newkey)] = bound_dic[keys] # Append dict with new key and value
        for i in range(MaxID+1):
            bound_dic[(i,i+MaxID+1)] = add_item
        #bound_dic[None] = add_item # Nr of elements in thickness
        return bound_dic.copy()
    
    def __interpolate_values(self,outer, inner, tnelems):
        '''Return the linearly interpolated values between outer and inner.
    
        Parameters
        ----------
        outer    : numpy.array :
            Outer surface values
        inner    : numpy.array :
            Inner surface values
        tnelems  : int :
            Number of elements in thickness          
        Returns
        -------
        : list of numpy.arrays :
            Interpolated values stored in a list with len(list) = tnelems
        '''   
        # Greville based non-linear interpolation
        knotv    = [0] + [ i/tnelems for i in range(0,tnelems+1) ] + [1]  # knotvector for p=2, with open boundaries
        greville = [0.5*(knotv[i+1]+knotv[i]) for i in range(len(knotv)-1)][1:-1] # greville points, neglect first and last point
        
        inside   = []
        x  = outer
        dx = inner - outer
        for f in greville:
            inside += [ f*dx + x ] 
            
        # Old linear interpolation    
        # inside = []
        # f  = 1 / ( tnelems + 1 )
        # x  = outer
        # dx = inner - outer
        # for i in range(tnelems):
        #     inside += [ (i+1)*f*dx + x ] 
        return inside
    
     ##-------------------------------------------------------------------------------------------------------------------------
    
    
    
    
    def boundary_names(self,): 
        ''' Names that are assigned to the local patch boundaries.

        Returns
        -------
        : dict :
            Dictionary containing the boundary names (keys) and local patch boundary names that are to be renamed (items)
            {inner: 'patch0-left,patch1-left,...', outer: ...}
        '''  
        
        if self.surf['left'] and self.surf['right']: # bi-ventricle
            boundaries    = {}    
            boundariesdic = { 'endo_r'      : [(3,'left'),(4,'left'),(5,'left'),(8,'left'),(11,'right'),(12,'right'),(13,'right'),(14,'right')],
                              'epi'         : [(0,'left'),(1,'left'),(7,'left'),(10,'left'),
                                               (11,'left'),(12,'left'),(13,'left'),(14,'left')],
                              'endo_l'      : [(0,'right'),(1,'right'),(2,'right'),(3,'right'),(4,'right'),
                                               (5,'right'),(6,'right'),(7,'right'),(8,'right'),(9,'right'),(10,'right')],
                              'septum_r'    : [(3,'left'),(4,'left'),(5,'left'),(8,'left')],
                              'endo_r_nosep': [(11,'right'),(12,'right'),(13,'right'),(14,'right')],
                              'base_r'      : [(11,'front'),(12,'front'),(13,'front')],
                              'base_l'      : [(0,'front'),(1,'front'),(2,'front'),(3,'front'),(4,'front'),(5,'front'),(6,'front'),(7,'front')]}                       
            return {key:','.join([boundaries[key]]+['patch{}-{}'.format(ipatch, bname) for ipatch,bname in val] if key in boundaries else ['patch{}-{}'.format(*v) for v in val]) for key, val in boundariesdic.items()}
        elif self.surf['left']: # left-ventricle
            boundaries    = {}    
            boundariesdic = { 'epi':        [(0,'left') ,(1,'left') ,(2,'left') ,(3,'left') ,(4,'left') ],
                              'endo_l':     [(0,'right'),(1,'right'),(2,'right'),(3,'right'),(4,'right')],
                              'base_l':     [(0,'front'),(1,'front'),(2,'front'),(3,'front')]}                       
            return {key:','.join([boundaries[key]]+['patch{}-{}'.format(ipatch, bname) for ipatch,bname in val] if key in boundaries else ['patch{}-{}'.format(*v) for v in val]) for key, val in boundariesdic.items()}

            
            
        # patch_name   = 'patch{}-{}'
        # patch_bound  = ('top','bottom','front','back') # ('-left','-right',) are not reuired as they are interfaces
        # rename_bound = (end, start, inner, outer)      
        # ndim = self.ndims() 
        
        # boundaries = {irename : [] for irename in rename_bound}
        
        # for i, (irename, ibound) in enumerate(zip(rename_bound,patch_bound)): # loop over boundary (re)names and the local boundary names
        #     if i > 1 and ndim != 3: # make distinction between 2D and 3D 
        #         boundaries.pop(irename) # remove the key
        #         continue
        #     else:
        #         for ipatch in range(self.patches().shape[0]): # loop over number of patches
        #                 localname = patch_name.format(ipatch, ibound)    
        #                 boundaries[irename] += [localname] 
        #     boundaries[irename] = ','.join(boundaries[irename]) # join the list of strings with , separator
        # return boundaries

    def get_topo_geom(self,rename=True):   
        ''' Return the nutils topology and geometry
        
        Parameters
        ----------
        rename  : bool :
            Choose if you want to rename some local patch boundaries or not.

        Returns
        -------
        : Topology :
            Nutils topology
        : Geometry :
            Nutils geometry
        '''          
        
        
        if self.surf['left'] and self.surf['right']:
            topo, lingeom = mesh.multipatch(patches=self.patchcon(), patchverts=self.patchverts(), nelems=self._bnelems)
            bsplines = topo.basis('spline', degree=2, patchcontinuous=False)
        else:
            topo, lingeom = mesh.multipatch(patches=self.patchcon(), patchverts=self.patchverts(), nelems=self.bnelems())
            bsplines = topo.basis('spline', degree=2, patchcontinuous=False, knotvalues=self.knotval(), knotmultiplicities=self.knotmult())
        
        weight   = bsplines.dot(self.weights())
        geom     = bsplines.vector(3).dot((self.cps()*self.weights()[:,numpy.newaxis]).ravel())/weight
        
        # Assign boundary names
        if rename:
            boundaries = self.boundary_names()
            topo       = topo.withboundary(**boundaries)
        return topo, geom
    
    def get_volume(self, shift = 0, output=False, wall=True):
        topo, geom = self.get_topo_geom(rename=True)
        
        ns     = function.Namespace()
        # Check if the base is at z=0
        zdiff    = topo.boundary['base_l'].sample('bezier',2).eval(geom)[0,-1] # Difference with z=0 
        if shift == 0 and not numpy.isclose(zdiff,0):
            treelog.info("Detected that the basal plane is not at z=0, shifting it accordingly.")
            ns.H = zdiff
        else:
            ns.H   = shift
            
        ns.Xshift_i = '< 0, 0, H >_i' 
        ns.x   = geom
        
        biventricle = True if self.surf['left'] and self.surf['right'] else False
        if biventricle: # Bi-ventricle
            Volume_rv = topo.boundary['endo_r'].integral('- ( 1 / 3 ) n_i ( x_i - Xshift_i ) d:x'@ ns, degree=8).eval()
            treelog.info('Volume Right-ventricle: {:.2f} [ml]'.format(Volume_rv*1e6))
        Volume_lv = topo.boundary['endo_l'].integral('- ( 1 / 3 ) n_i ( x_i - Xshift_i ) d:x'@ ns, degree=8).eval()
        treelog.info('Volume Left-ventricle:  {:.2f} [ml]'.format(Volume_lv*1e6))
        
        if wall:
            Volume_wall = topo.integral('d:x'@ ns, degree=8).eval()
            treelog.info('Wall volume:  {:.2f} [ml]'.format(Volume_wall*1e6))
            #Volume_lv = topo.boundary['endo_l'].integral('- ( 1 / 3 ) n_i ( x_i - Xshift_i ) d:x'@ ns, degree=8).eval()
        
        
        
        if output:
            if wall:
                if biventricle:
                    return Volume_lv, Volume_rv, Volume_wall
                else:
                    return Volume_lv, Volume_wall
            else:
                if biventricle:
                    return Volume_lv, Volume_rv
                else:
                    return Volume_lv




     
    
    
    
    
    
    
    
    
    
    
    
    
    
    # def get_nelems(self,nelems,MaxID):
    #     for keys in list(nelems):
    #         newkey = numpy.asarray(keys) + MaxID  + 1        # Capture old keys and increment them
    #         nelems[tuple(newkey)] = nelems[keys] # Append dict with new key and value
    #     nelems[None] = 1 # Nr of elements in thickness
    #     return nelems.copy()
    
    # def interpolate_cps_w(self, cpsw1, cpsw2, tnelems):
        
    #     for i in range(tnelems):
    #         cps = i*cpsw1[0] 
        
    #     cps = 1
    #     w = 1
    #     return cps, w
    
    # def surfaces(): # Combine surface inner and outer surfaces of either the left- or right-ventricle
    #     return
    
    # def solids(): # combine solids (left- and right-ventricle)
    #     return
    
    
#%     
if __name__ == '__main__':
    
   ## Left-ventricle geometry 
    mm=1e-3;ml=1e-6  
    C    = 43.0*mm
    H    = 24.8*mm
    Vpap = 4.0 *ml
    Vlv0 = 40.0*ml
    Vw   = 140.*ml
    Vi   = Vlv0 + Vpap
    Ve   = Vw - Vpap

    ξi = 0.3694447442932904 #0.371299?
    ξo = 0.6754874350057226 #0.678357? 
   
    Rleft_inner  = C*numpy.array([numpy.sinh(ξi),numpy.sinh(ξi),numpy.cosh(ξi)])
    Rleft_outer  = C*numpy.array([numpy.sinh(ξo),numpy.sinh(ξo),numpy.cosh(ξo)])

    LV_outer = {'Rx': Rleft_outer[0]   , 'Ry': Rleft_outer[1]  , 'Rz': Rleft_outer[2], 'H': H}
    LV_inner = {'Rx': Rleft_inner[0]   , 'Ry': Rleft_inner[1]  , 'Rz': Rleft_inner[2]} # H can be given, but not required if it is in LV_inner

    leftventricle  = CardiacGeometry('left-ventricle', LVi=LV_inner, LVo=LV_outer) # Constructs the left-ventricle inner and outer surfaces
    surfaceLV_outer = leftventricle.generate_surface('left-outer')
    surfaceLV_inner = leftventricle.generate_surface('left-inner')

    # save_pickle(surfaceLV_outer, filename='LV_left_outer')
    # save_pickle(surfaceLV_inner, filename='LV_left_inner')
    surfaceLV_outer = load_pickle('output/pickle/LV_left_outer')[0]
    surfaceLV_inner = load_pickle('output/pickle/LV_left_inner')[0]

   
    #solid = MultipatchSolid(surfaceLV_outer, surfaceLV_inner)
    
    #solid.surf['left']['inner'].dependencies(4)
    #solid.surf['left']['inner'].update_cps(0)
    #solid.surf['left']['inner'].update_patchverts()
    
    #solid.get_volume(shift = H)
    #solid.save_vtk('output/vtk/LeftventricleGeom') #, boundaries=('base_l','epi','endo_l'), boundary_names=solid.boundary_names()) #boundary_only=True)
    #solid.save_stl('output/vtk/LeftventricleGeom', boundary_names=solid.boundary_names()) 
    #solid.save_mesh_vtk('output/vtk/LeftventricleGeom_mesh')
    #solid.save_pickle('pickle files/LV_GEOMETRY_DATA', boundary=solid.boundary_names())

    # Save mesh refinement
    nelems  = 7
    order   = 0,1,2,3,4
    surfaceLV_outer.refine(nelems)
    surfaceLV_inner.refine(nelems)
    solid   = MultipatchSolid(surfaceLV_outer, surfaceLV_inner, tnelems=(nelems+1,1))
    topo, geom  = solid.get_topo_geom_alt()
    #solid.save_vtk('output/vtk/LeftventricleGeom',boundary_only=True, nrsamples=8) 
    #solid.save_mesh_vtk('output/vtk/LeftventricleGeom_mesh',nsample = 5)    
  
   ## Bi-ventricle geometry  
    # mm=1e-3;ml=1e-6  
    # Folder = 'ECCOMAS geom variations/'
    # Name   = 'LONG'
    # ##  Input
    # C    = 43.0*mm
    # H    = 24.8*mm
    # Vpap = 4.0 *ml
    # Vlv0 = 40.0*ml
    # Vw   = 140.*ml
    # Vi   = Vlv0 + Vpap
    # Ve   = Vw - Vpap

    # ξi = 0.3694447442932904
    # ξo = 0.6754874350057226   

    # dr    = 1.1 #1.#1.1
    # dr2   = 1 #1.1 
    # #Rlong = numpy.array([ 1/dr, 1/dr, dr**2 ])
   
    # Rleft_inner  = C*numpy.array([numpy.sinh(ξi),numpy.sinh(ξi),numpy.cosh(ξi)])
    # Rleft_outer  = C*numpy.array([numpy.sinh(ξo),numpy.sinh(ξo),numpy.cosh(ξo)])
    # Rleft_middl  = 0.5*(Rleft_inner + Rleft_outer)  

    # RV_wall       = 0.005 #0.008 #0.005
    # dZ_right      = 0.005
    # # Rright_inner  = numpy.array([C*numpy.sinh(ξo) + RV_wall,0.7*C*numpy.sinh(ξo)+RV_wall,abs(C*numpy.cosh(ξo)*numpy.cos( 0.85*numpy.pi )) - RV_wall ])
    # # Rright_outer  = Rright_inner + RV_wall
    # Rright_outer  = numpy.array([C*numpy.sinh(ξo) + RV_wall,0.7*C*numpy.sinh(ξo)+RV_wall,abs(Rleft_outer[-1]*numpy.cos( 0.85*numpy.pi )) - dZ_right ])
    # Rright_inner  = Rright_outer - RV_wall
   
    # LV_outer = {'Rx': C*numpy.sinh(ξo)/dr, 'Rxs': 0.8, 'Ry':  C*numpy.sinh(ξo)/dr, 'Rz': C*numpy.cosh(ξo)*dr**2, 'H': H}
    # LV_inner = {'Rx': C*numpy.sinh(ξi)/dr, 'Rxs': 0.8, 'Ry':  C*numpy.sinh(ξi)/dr, 'Rz': C*numpy.cosh(ξi)*dr**2} # H can be given, but not required if it is in LV_inner 
    # RV_inner = {'Rx': Rright_inner[0]/dr2, 'Ry': Rright_inner[1]/dr2, 'Rz': Rright_inner[2]*dr**2}
    # RV_outer = {'Rx': Rright_outer[0]/dr2, 'Ry': Rright_outer[1]/dr2, 'Rz': Rright_outer[2]*dr**2, 'Origin': numpy.array([-Rleft_inner[0],0.,0])} # H can be given, but not required if it is in LV_inner

    # ##-----------------------------------
    # biventricle    = CardiacGeometry('bi-ventricle',   LVi=LV_inner, LVo=LV_outer, RVi=RV_inner, RVo=RV_outer) # Constructs the bi-ventricle inner and outer surfaces    
    
    # surfaceLV_outer = biventricle.generate_surface('left-outer')
    # surfaceLV_inner = biventricle.generate_surface('left-inner')
    # surfaceRV_outer = biventricle.generate_surface('right-outer')
    # surfaceRV_inner = biventricle.generate_surface('right-inner')

    # save_pickle(surfaceLV_outer, filename=Folder+'pickle/BV_left_outer'+'_'+Name)
    # save_pickle(surfaceLV_inner, filename=Folder+'pickle/BV_left_inner'+'_'+Name)
    # save_pickle(surfaceRV_outer, filename=Folder+'pickle/BV_right_outer'+'_'+Name)
    # save_pickle(surfaceRV_inner, filename=Folder+'pickle/BV_right_inner'+'_'+Name)
    # surfaceLV_outer = load_pickle(Folder+'pickle/BV_left_outer'+'_'+Name)[0]
    # surfaceLV_inner = load_pickle(Folder+'pickle/BV_left_inner'+'_'+Name)[0]
    # surfaceRV_outer = load_pickle(Folder+'pickle/BV_right_outer'+'_'+Name)[0]
    # surfaceRV_inner = load_pickle(Folder+'pickle/BV_right_inner'+'_'+Name)[0]
   
    # solid = MultipatchSolid(surfaceLV_outer, surfaceLV_inner, surfaceRV_outer, surfaceRV_inner)
    # solid.get_volume(shift = H)
    # solid.save_vtk(Folder+'vtk/BiventricleGeom'+'_'+Name)#, boundary=True)
    # solid.save_pickle('ECCOMAS geom variations/BV_GEOMETRY_DATA'+'_'+Name, boundary=solid.boundary_names())
   
   ##-----------------------------------
 




  
  

  # long    = 2
  # pnelems_i = {'patch0': [1,long], 'patch1': [1,long], 'patch2': [1,long], 'patch3': [1,long], 'patch4': [long,long]} # General element distribution per patch (should match with multipatches)
  # loftd_i   = {'patch0':  1      , 'patch1':  1      , 'patch2':  0      , 'patch3':  1      , 'patch4':  1}          # lofting direction, 0 or 1
  # lnelems_i = {'patch0':  1      , 'patch1':  1      , 'patch2':  1      , 'patch3':  1      , 'patch4':  1}          # nr of elements in lofting direction, value should be below pnelems
 
  # surfaceLV_outer = leftventricle.generate_surface('left-outer')#, pnelems_i=pnelems_i, loftd_i=loftd_i, lnelems_i=lnelems_i)
  # surfaceLV_inner = leftventricle.generate_surface('left-inner')#, pnelems_i=pnelems_i, loftd_i=loftd_i, lnelems_i=lnelems_i)
   
  # save_pickle(surfaceLV_outer, filename='left_outer')
  # save_pickle(surfaceLV_inner, filename='left_inner')
  
  
  # surfaceLV_outer_s = load_pickle('left_outer')[0]
  # surfaceLV_inner_s = load_pickle('left_inner')[0]
  
  # solid = MultipatchSolid(surfaceLV_outer_s, surfaceLV_inner_s)
  # solid.save_vtk('vtk files/Solid')
  
  
  
  
  
  
  
  
  
  
  
  
     
   # biventricle    = CardiacGeometry('bi-ventricle',   LVi=LV_inner, LVo=LV_outer, RVi=RV_inner, RVo=RV_outer) # Constructs the bi-ventricle inner and outer surfaces
   # biventricle.visualize('left-outer','left-inner','right-outer','right-inner')
   # biventricle.visualize('right-outer')#,'left-inner')
    
    
   # surfaceLV_outer = biventricle.generate_surface('left-outer')
   # surfaceLV_inner = biventricle.generate_surface('left-inner')
   # surfaceRV_outer = biventricle.generate_surface('right-outer')
   # surfaceRV_inner = biventricle.generate_surface('right-inner')
  
   # # # surfaceLV_outer.save_vtk('vtk files/test_outerLV')
   # # # surfaceRV_outer.save_vtk('vtk files/test_outerRV')
   # # # surfaceLV_inner.save_vtk('vtk files/test_innerLV')
   # # # surfaceRV_inner.save_vtk('vtk files/test_innerRV')


   # save_pickle(surfaceLV_outer, filename='BV_left_outer')
   # save_pickle(surfaceLV_inner, filename='BV_left_inner')
   # save_pickle(surfaceRV_outer, filename='BV_right_outer')
   # save_pickle(surfaceRV_inner, filename='BV_right_inner')
  
  
   # surfaceLV_outer = load_pickle('BV_left_outer')[0]
   # surfaceLV_inner = load_pickle('BV_left_inner')[0]
   # surfaceRV_outer = load_pickle('BV_right_outer')[0]
   # surfaceRV_inner = load_pickle('BV_right_inner')[0]
   
  # surfaceLV_outer.save_vtk('vtk files/test_outerLV')
  # surfaceRV_outer.save_vtk('vtk files/test_outerRV')
  # surfaceLV_inner.save_vtk('vtk files/test_innerLV')
  # surfaceRV_inner.save_vtk('vtk files/test_innerRV')


    
  # solid = MultipatchSolid(surfaceLV_outer, surfaceLV_inner, surfaceRV_outer, surfaceRV_inner)
  #solid = MultipatchSolid(surfaceLV_outer, surfaceLV_inner)
  #solid = MultipatchSolid(surfaceRV_outer, surfaceRV_inner)
  #solid.save_vtk('vtk files/BiventricleGeom', patch=6)
  # solid.save_vtk('vtk files/BiventricleGeom')#, boundary=True)
  
  #solid.save_pickle('pickle files/BV_GEOMETRY_DATA_NEW', boundary=solid.boundary_names())
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## Sample only relevant elements
   # Zinterest   = -1
   # topo, geom  = solid.get_topo_geom()
   # topo        = topo.refine(2)
   # bezierverts = topo.sample('bezier', 2)
   # nverts      = 8 # 8 vertices per patch
   # dim         = 3 # spatial dimension
   # X           = bezierverts.eval(geom).reshape(-1,nverts,dim)
   # elemID      = []
   # for ID, elemsi in enumerate(X):
   #     minZ = min(elemsi[:,-1]) 
   #     maxZ = max(elemsi[:,-1])
   #     if Zinterest <= maxZ and Zinterest >= minZ:
   #         elemID += [ID]
     
   # # for i, ID in enumerate(elemID):
   # #     if i == 0:
   # #         topo_sample = topo[[ID]]
   # #     else:
   # #         topo_sample += topo[[ID]]
   # topo_sample = topo[elemID]        
   # bezier     = topo_sample.sample('bezier', 15)
   # Xsample    = bezier.eval(geom)
   # export.vtk('testing', bezier.tri, Xsample)
   
   
  # save_pickle(topo, geom, filename='TopoGeom_Test')          
  
  # bezierER      = topo.boundary['endo_r'].sample('bezier', 15)
  # bezierEL      = topo.boundary['endo_l'].sample('bezier', 15)
  # bezierEPI     = topo.boundary['epi'].sample('bezier', 15) 
  # bezierBAL     = topo.boundary['base_l'].sample('bezier', 15)
  # bezierBAR     = topo.boundary['base_r'].sample('bezier', 15)
  # bezierSEP     = topo.boundary['septum_r'].sample('bezier', 15)
  # bezierSEPr    = topo.boundary['endo_r_nosep'].sample('bezier', 15)

  # export.vtk('BoundER', bezierER.tri, bezierER.eval(geom)) 
  # export.vtk('BoundEL', bezierEL.tri, bezierEL.eval(geom)) 
  # export.vtk('BoundEpic', bezierEPI.tri, bezierEPI.eval(geom)) 
  # export.vtk('BoundBaseL', bezierBAL.tri, bezierBAL.eval(geom)) 
  # export.vtk('BoundBaseR', bezierBAR.tri, bezierBAR.eval(geom)) 
  # export.vtk('BoundSepR', bezierSEP.tri, bezierSEP.eval(geom)) 
  # export.vtk('BoundSepnoR', bezierSEPr.tri, bezierSEPr.eval(geom)) 



  #surfaceLV_inner = leftventricle.generate_surface('left-inner')
  #surfaceLV_inner.save_vtk('vtk files/testinnerLV')
  #surfaceLV_inner.save_pickle('pickle files/testinnerLV')
  
  # store  = StoredInput(LV_outer)
  # store.convert2('mm') 
  #Combine(1,2,3)
  
  ## Possibilities
  #solid = Combine( surf_lv_outer, surf_lv_inner ) # <-- Create left-ventricle
  #solid = Combine( surf_lv_outer, surf_lv_inner, surf_rv_outer, surf_rv_inner ) # <-- Create Biventricle
  
  
  #solid.construct(nelems=2) # Number of elements in thickness
  #solid.construct(nelems={'lv':1,'bv':2}) #
  
  
  
  
  #surfaceLV.cps, surface.w, surface,patchverts, surface.patchconn, surface.nelems
  #
  # Adjust cps, w etc if wanted
  # leftventricle.update() # updates all cps, w and nelems based on adjustments
  
  
  
  # solidLV       = leftventricle.combine.surfaces('left') # Combine surfaces left-ventricle
  # 
  # solidLV       = biventricle.combine.surfaces('left')
  # solidRV       = biventricle.combine.surfaces('right')
  
  # solidBV = solidLV + solidRV