# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 13:52:58 2021

@author: s146407
"""

## Surface class ##
from nutils import mesh, function, solver, export,matrix
import numpy, treelog, re
import matplotlib.pyplot as plt
from splipy import curve_factory, BSplineBasis, Curve, surface_factory
from .loftquadratic import my_loft_function

__matrix__ = matrix.backend('MKL')
#__matrix__ = matrix.backend('scipy')
#__matrix__ = matrix.backend('numpy')
__matrix__.__enter__()

class SurfaceEllipsoid:
    
    def __init__(self, R, patchverts, patch, pnelems={None: (1,1)}, lnelems={None: 1}, loftd=None, S=numpy.zeros(3), bcons=None, bcons_cpsw=None, coons={None: False}):
        """  Construct a Ventricle (left or right) surface based on the input parameters.

        Note: This script does not create an eintire Biventricle, but 
        merely the surface of an arbitrary ellips (subject to certain constrains).
        Required--------------------------------------------------------
        :param R:          Array containing radii in [x,y,z]-directions
        :param patchverts: Array containing all the patchvertices that are used for the mulitpatch geometry
        :param patch:      An n x 4 Connectivity array of the patchvertices, where n is the number of patches 
        
        Optional-------------------------------------------------------
        :param pnelems:    Dictionary containing the total number of element per patch for its parametric coordinates, (u,v) 
                            example: pnelems = {None: (1,1), 0: (1,2), 1: (3,4)}
                            Defaults to a single element patch (1,1)
        :param lnelems:    Dictionary containing the number of elements used to loft a patch. The direction in which these
                            nr. of elements are assigned is given in loftd. The value used in lnelems should he less than 
                            or equal to the the same given in pnelems (total nr. of elements). 
                            Ideally, this should also hold: mod(pnlems,lnelem)=0
        :param loftd:      Dictionary that contains the local direction of lofting (0 or 1) per patch, example: loftd={3: 0, 2: 1, 10: 0}                  
        :param S:          A 1x3 array, containing the physical coordinates of the origin of the ellips 
        :param bcons:      Dictionary containing geometrical constrains specified per line segment. 2 Constrains are available:
                           - Plane: Specify the cutting plane, {(0,1): ('plane', [a,b,c,d])}, where a*x+b*y+c*z=d
                           - Ellips: Specify a second ellips, {(0,1): ('ellips', S, R)}, where S = [x,y,z] (origin), R = [Rx,Ry,Rz] (radii)
        :param bcons_cpsw: Dictionary containing the control points (cps) and weights (w) for a line segment: {(0,1): [cps, w]}
        
        Return---------------------------------------------------------
        :return: Structure containing all patch control points and weight
        :rtype: Dictionary containing Splipy surfaces    
        """
        
        ## Store input
        self.__R          = R           # Radius array: R=[Rx,Ry,Rz]
        self.__S          = S           # Origin of the ellips, default=[0,0,0]
        self.__patchverts = patchverts  # Patch vertices array 
        self.__patch      = patch       # Patch connectivity array, multiple surface patches are possible
        
        
        self.degree       = 2 # the surface is constructed with quadratic splines only
        self.__degree     = 2
        
        ## Change lnelems, pnelems and loftd keys from 'patch0', 'patch1', 'patch2', ... to 0, 1, 2, ..
     
        
        
        self.__lnelems    = dict((int(key[-(len(key)-5):]), value) for (key, value) in lnelems.items()) #lnelems # Lofting elements per patch in the loft-direction, should be read in combination with loftd (the loft direction) 
        #self.__bnelems    = bnelems      # , bnelems={None: 1} Array corresponding each patch, consisting the number of elems in u & v-direction
        self.__pnelems    = dict((int(key[-(len(key)-5):]), value) for (key, value) in pnelems.items()) #pnelems
        self.__coons      = coons 
        
        # Create empty dictionaries if None is specified
        if bcons == None:
            self.__bcons = {} # Empty dict
        else:
            self.__bcons = bcons
        if bcons_cpsw == None:
            self.__bcons_cpsw = {} # Empty dict
        else:
            self.__bcons_cpsw = bcons_cpsw
            
            
        self.run = 0 # Post-processing step initial value    
        ## Calculate curve/surface cps and weights:
        self.get_bound_connectivity()
        
        ## Pre-processing
        loftdN = dict((int(key[-(len(key)-5):]), value) for (key, value) in loftd.items()) 
        self.get_loftdirection(loftdN) # Create default lofting direction (based on best result), can be altered by input loftb (if specified)
        self.generate_nelems_dict()   # generate nelems dict, which agrees with pnelems and bnelems, creates Error otherwise

        ## Processing
        LHSb = self.get_boundary_cpsw()
        self.solve_inner_curves(LHSb)
        #self.construct_surfaces()

        
    def get_bound_connectivity(self,):
        """  Get the (unique) boundary connectivity of the multipatch/patches.
        :store self.nrelemsb:  Dictionary containing number of elements for each local boundary 
                                per patch id: {0: [elem_left,elem_right,elem_bot,elem_top]}
        :store self.boundidx:  Array containing all unique boundary connectivities [[0,1],[2,3],[..]]
        :store self.boundidxp: Array containing the boundary connectivities for each patch [ [ [0,1],[2,3],[..] ], [ [], [], .. ], ..]
        :store self.nbound:    Total number of unique boundaries present
        """
        # Construct boundary connectivity array
        lp            = len(self.__patch)
        boundidxp     = numpy.ones((lp,4,2), dtype=int)
        self.nrelemsb = dict.fromkeys([i for i in range(len(self.__patch))])
        for i, ipatch in enumerate(self.__patch):
            bindx     = numpy.reshape(ipatch, [2]*2)
            boundidxp[i,:] = numpy.array([ bindx[0,:], bindx[1,:], bindx[:,0], bindx[:,1] ]) 
            self.nrelemsb[i] = [1]*4
            # for bj, j in enumerate(boundidxp[i,:]):
            #     if tuple(j) in self.__bnelems:
            #         self.nrelemsb[i][bj] =  self.__bnelems[tuple(j)] 
        boundidx      = numpy.unique(boundidxp.reshape(-1,2), axis=0) # Filter duplicates, make unique

        self.boundidx  = boundidx#[args] # Unique
        self.boundidxp = boundidxp       # Non-unique (per-patch)
        self.nbound    = len(boundidx)   # nr of unique boundaries
        return

    def get_loftdirection(self, loftd, default=False):
        """  Construct dictionary that contains the lofting direction per patches. 
             These directions are either 0 or 1, which refer to the parametric surface directions.
        :param loftd:   Optional dictionary containing the direction of certain patches
        :param default: True or False, returns 'random' (all 0-direction) dictionary if True
        """
        loftdir      = dict.fromkeys([i for i in range(len(self.__patch))])
        self.loftdir = {x: 0 for x in loftdir}
        if default==True: # Return default dictionary ('random')
            return # Just pick the first direction 0

        # # Make distinction for left- and right-ventricle (they have differences)
        # if len(self.__patch) > 5: # In case of left-ventricle
        #     patches = [0]   # List of patches with different direction
        # else: # In case of right-ventricle
        #     patches = [0]
            
        # for i in patches:
        #     self.loftdir[i] = 1
            
        # Change specific patch directions
        if loftd==None:
            return
        else:
            for i in loftd:
                self.loftdir[i] = loftd[i]    
            return


    def assign_nelems(self,boundaryIDX,nelems):
        """  Assign specified number of elements to self.nelems dictionary.
               In case the value to be overwritten is different than the current value or None,
               an exception is given.
        :param boundaryIDX: 1x2 Array containing the oundary connectivity [0,1], [2,3], [0,6] , ..
        :param nelems:      Number of elements to be assigned to this boundary
        """
        if self.nelems[tuple(boundaryIDX)] == None or self.nelems[tuple(boundaryIDX)] == nelems:
            self.nelems[tuple(boundaryIDX)] = nelems
        else:
            raise ValueError('Mismatch in nr of elements between patches, check pnelems.')
        return
    
    def generate_nelems_dict(self,):
        """     Generates the self.nelems dictionary, which is required for the eventual Nutils topology. 
                 It is a dictionary containing the number of element per curve, example: 
                 bnelems = {None: 1, (0,1): 3, (4,5): 2}.
                 
                 The function also asserts the pnelems input, which should agree with the lofting nr of elements
        """         
        ## Check the nr of elements input, raise error if wrong          
        for ip in range(len(self.__patch)):
            
            # Create knot vectors (u and v) direction
            if ip in self.__lnelems:
                m = self.__lnelems[ip] # nr of elements
            else:
                m = 1
            ## Only allow refinement if it is uniformly possible
            if self.__pnelems[ip][0] % ((1-self.loftdir[ip])*(m - 1) + 1)  != 0:# If the residual is 0, this means it is hierarchicaly possible
               raise ValueError('Refinement is not Hierarchical')           
            if self.__pnelems[ip][1] % (self.loftdir[ip]*(m-1) + 1)  != 0:# If the residual is 0, this means it is hierarchicaly possible
               raise ValueError('Refinement is not Hierarchical') 
                                
        ## Continue with the verified input       
        self.nelems  = self.construct_bound_dict()
        #self.nelems = {x: 1 for x in nelems} # All boundaries have as default 1 element  
         
        # Loop over each pnelems patch combi
        for i in self.__pnelems:
            if i == None: # If none, it means use default elements = 1
                for idx in self.boundidxp[i]:
                    self.assign_nelems(idx,1)
            else:
                for j, idx in enumerate(self.boundidxp[i]):
                    if j <= 1:
                        self.assign_nelems(idx,self.__pnelems[i][0])
                    else:
                        self.assign_nelems(idx,self.__pnelems[i][1])   
        return   
        
    def construct_bound_dict(self,):
        """  Construct empty solution dictionary containing both control points 
             and weights (arrays) per unique boundary.
        """
        return dict.fromkeys([tuple(i) for i in self.boundidx])   
    
    def get_boundary_cpsw(self,):
        """  Store the solution (control points and weights) for each unique boundary.
        :return LHSb:  Dictionary containing the solved control points and weights 
                        for each boundary curve: LHSb = { (0,1): (cps,w), (2,3): (cps,w), .. }
        """
        LHSb = self.construct_bound_dict() 

        for j, jb in enumerate(self.boundidx): # loop over unique boundaries (indices)
            
            if tuple(jb) in self.__bcons_cpsw or tuple(jb)[::-1] in self.__bcons_cpsw: # Check whether it has already been determined/stored           
                
                ## Convert to unit sphere, first convert to left ellips, than translate en scal eback to right ventr
                # W   = self.__bcons_cpsw[tuple(jb)][1]
                # CPS = ( self.__bcons_cpsw[tuple(jb)][0].reshape(-1,3)/W[:,numpy.newaxis] ) + self.__S 
                # LHSb[tuple(jb)] = ( (CPS * W[:,numpy.newaxis]).ravel(), W )
                
                LHSb[tuple(jb)] = self.__bcons_cpsw[tuple(jb)]
                
            elif tuple(jb) in self.__bcons or tuple(jb)[::-1] in self.__bcons:# Check whether the boundary has constrains => immediately solve for it
                with treelog.context('Boundary {} (constrained) {}/{}'.format(tuple(jb),j+1,self.nbound)):
                    LHSb[tuple(jb)] = self.get_lhs_boundary(jb, constrain=self.__bcons[tuple(jb)]) 
            
            else: # Else solve for minimum strain constrain
                with treelog.context('Boundary (unconstrained) {}/{}'.format(j+1,self.nbound)): #- len(self.__bcons_cpsw) - len(self.__bcons)
                    lhsbx, lhsbw = self.get_lhs_boundary(jb)
                    LHSb[tuple(jb)] = (lhsbx, lhsbw) 
                             
        self.LHSb = LHSb # Store (can be removed in eventual script)
        return LHSb          
    
    def get_lhs_boundary(self, vertIDX, weight=[1,1], constrain=None, IDX=True):
        """  Solve for each (unique) boundary curve (quadratic, single element) on the unit sphere,
                unless the constrain type is ellips, in that case it is computed for the actual ellipsoid.
        
        Note: Solving on the unit sphere showed faster convergence compared to solving on the actual ellipsoid
                The computed result will be converted to the ellipsoid in a later stage.
                
        :param vertIDX: 1 x 2 Array containing the indices of the patchvertices that form a curve/line
                         OR
                        2 x 3 Array containing the physical coordinates through which the curve passes,
                        if this is the case, IDX should be set to False; IDX=False
        :param contrain: Tuple containing type of constrain (plane or ellips) and corresponding parameters
        :param IDX:      True or False, determines the interpretation of vertIDX (either index or physical coordinates)
            
        :return lhsbx: Array (ravel) containing the control points on the unit sphere 
        :return lhsbw: Array containing the weights
        """
        
        topo, geom = mesh.rectilinear([numpy.linspace(0,1,2)]) # Construct curve geometry
        ns         = self.define_namespace(topo, geom, cons=constrain)
        
        ## Create empty initial arrays & assign vertex constrains
        if ( constrain==None or constrain[0] == 'plane' ) and IDX==True: # If there is no constrain convert to unit sphere coordinates
            X        = self.convert_to_sphere(self.__patchverts[vertIDX])
        elif IDX==True:
            X        = self.__patchverts[vertIDX]  
        elif IDX==False: # Input is vertex coordinates
            X        = vertIDX #- self.__S
        lhsxi    = self.linterp_bound(X, 1).reshape(1,-1)[0] # Initial guess, linearized field
        lhsxi[:3]  = lhsxi[:3]*weight[0]
        lhsxi[-3:] = lhsxi[-3:]*weight[1]
        
        consx    = lhsxi.copy()
        consx[3:-3] = numpy.nan 
        consw     = numpy.ones((self.__degree+1))
        consw[0]  = weight[0]
        consw[-1] = weight[1]
        conswi    = consw.copy()
        
        consw[1:-1] = numpy.nan # Start and end vertices are by defined by the weights, rest nan
        lagr        = conswi.copy()
        lagr[:]     = numpy.nan
        
        ## Solve the remaining cps and w's
        #if constrain==None:
        # lhsi  = numpy.concatenate([lhsxi, conswi, numpy.ones((self.__degree+1))]) # Initial guess, w=1, lagrangian=1
        # #consi = numpy.concatenate([consx, ones, ones]) # Constrain, w=1, lagrangian=1
        # cons  = numpy.concatenate([consx, consw, lagr]) # Note consw is same as conslagrangian
        # #nrgh = topo.integral('( ε_ij ε_ij + λ ellips ) d:u' @ ns, degree=4*self.__degree)
        # #lhsb = solver.optimize('lhs', nrgh, lhs0=lhsi, constrain=consi, droptol=1e-6, tol=1e-6)        

        ## Calculate NURBS cps & w values
        if constrain==None:
            # Constrain arrays
            lhsi  = numpy.concatenate([lhsxi, conswi, numpy.ones((self.__degree+1))]) # Initial guess, w=1, lagrangian=1
            cons  = numpy.concatenate([consx, consw, lagr]) 
            
            nrg  = topo.integral('( 0.5 ε_ij ε_ij + λ sphere ) d:u' @ ns, degree=4*self.__degree)
            #lhs  = solver.optimize('lhs', nrg, lhs0=lhsb, constrain=cons, tol=1e-6)      
            lhs  = solver.optimize('lhs', nrg, lhs0=lhsi, constrain=cons, tol=1e-6)  
            lhsbx, lhsbw = self.split_lhs(lhs, 1)

        elif constrain[0]=='default': # constrain is the actual ellips (not unit sphere)
            lhsi  = numpy.concatenate([lhsxi, conswi, numpy.ones((self.__degree+1))]) # Initial guess, w=1, lagrangian=1
            cons  = numpy.concatenate([consx, consw, lagr]) 
            
            nrg   = topo.integral('( 0.5 ε_ij ε_ij + ( λ ellips )  ) d:u' @ ns, degree=4*self.__degree)    
            lhs   = solver.optimize('lhs', nrg, lhs0=lhsi, constrain=cons, tol=1e-6) 
            lhsbx, lhsbw = self.split_lhs(lhs, 1)
            
        elif constrain[0]=='plane': # Constrain if a plane
            # Constrain arrays
            # lhsi  = numpy.concatenate([lhsxi, conswi, numpy.ones((self.__degree+1)),numpy.ones((self.__degree+1))]) # Initial guess, w=1, lagrangian=1
            # cons  = numpy.concatenate([consx, consw, lagr, lagr]) 
            # #lhsi  = numpy.concatenate([lhsxi, conswi, numpy.ones((self.__degree+1))]) # Initial guess, w=1, lagrangian=1
            # #cons  = numpy.concatenate([consx, consw, lagr]) 
            # ns.geomcons = self.get_geom_constrain(constrain)
            # nrg  = topo.integral('( ε_ij ε_ij + ( λ ellips ) + ( γ geomcons ) ) d:u' @ ns, degree=4*self.__degree)  
            # #nrg  = topo.integral('( 0.5 ellips + ( λ sphere  ) + γ ( geomcons ) ) d:u' @ ns, degree=4*self.__degree)
            # lhs  = solver.optimize('lhs', nrg, lhs0=lhsi, constrain=cons, tol=1e-6) 

            # lhsbx, lhsbw = self.split_lhs(lhs, 1)
            ## Analytically determine the cps location and weight of a quadratic NURBS n a sphere - plane intersection
            plane_str = self.get_geom_constrain(constrain, map2sphere=True) # plane expression in string format
            plane_num = [float(val.replace(" ","")) for val in re.split('x_0 | x_1 | x_2', plane_str)] # convert to useable numbers [a,b,c,d], from a*x+b*y+c*z=d

            n_plane = numpy.array(plane_num[:-1]) # Normal of the plane (unit length)
            δ_plane = -plane_num[-1]               # Distance plane from origin   
            ρ_plane = numpy.sqrt(1-δ_plane**2)    # Radius intersection circle
            e0 = δ_plane*n_plane
            if numpy.isclose(abs(n_plane[-1]),1):
                e1 = numpy.array([ρ_plane,0,0]) # vector inside the plane
                e2 = numpy.array([0,ρ_plane,0]) # second vector inside the plane
            else:
                e1 = ρ_plane/numpy.linalg.norm(n_plane[:-1])*numpy.array([n_plane[1],-n_plane[0],0])
                e2 = numpy.cross(n_plane, e1)

            p1_circle = X[0] - e0
            p2_circle = X[1] - e0
            dp_circle = p2_circle - p1_circle

            sina = numpy.linalg.norm( numpy.cross(p1_circle,dp_circle) )/( numpy.linalg.norm(p1_circle)*numpy.linalg.norm(dp_circle) )

            w   = sina
            cps = ρ_plane/sina * ( p1_circle + p2_circle )/ numpy.linalg.norm( p1_circle + p2_circle ) + e0

            lhsbx = numpy.concatenate([X[0],cps*w,X[1]]).reshape(-1,3)
            lhsbw = numpy.array([1,w,1])
            lhs   = numpy.concatenate([lhsbx.ravel(), lhsbw, numpy.ones((self.__degree+1)),numpy.ones((self.__degree+1))]) # include lagrange multipliers for evaluation step eventhough this is not necessary (old bug fix)   
            # Constrain arrays (does not work anymore with updated nutils... -> search vector does not reduce residual)
            # lhsi  = numpy.concatenate([lhsxi, conswi, numpy.ones((self.__degree+1)),numpy.ones((self.__degree+1))]) # Initial guess, w=1, lagrangian=1
            # cons  = numpy.concatenate([consx, consw, lagr, lagr]) 
            # ns.geomcons = self.get_geom_constrain(constrain, map2sphere=True)
            # nrg  = topo.integral('( 0.5 ε_ij ε_ij + ( λ sphere ) + ( γ geomcons ) ) d:u' @ ns, degree=4*self.__degree)  
            # lhs  = solver.optimize('lhs', nrg, lhs0=lhsi, constrain=cons, tol=1e-6) 

            # lhsbx, lhsbw = self.split_lhs(lhs, 1)
            
            
        elif constrain[0]=='ellips': # Constrain if an ellips
            # Constrain arrays
            lhsi  = numpy.concatenate([lhsxi, conswi,numpy.ones((self.__degree+1))]) # Initial guess, w=1, lagrangian=1
            cons  = numpy.concatenate([consx, consw, lagr]) 
            
            ns.geomcons = self.get_geom_constrain(constrain)
            #nrg  = topo.integral('( ellips^2 + λ geomcons ) d:u' @ ns, degree=4*self.__degree)   
            nrg  = topo.integral('( 0.5 ellips^2 + λ geomcons ) d:u' @ ns, degree=4*self.__degree)   
            lhs  = solver.optimize('lhs', nrg, lhs0=lhsi, constrain=cons, tol=1e-6) 
            lhsbx, lhsbw = self.split_lhs(lhs, 1)

        ## Compute the L2-error and map onto sphere if not yet done
        if constrain==None or constrain[0] == 'plane':
            Error  = topo.integral('( sphere )^2 d:u' @ ns, degree=8*self.__degree)
        else:
            lhsbNH = lhsbx.reshape(-1,3)/lhsbw[:,numpy.newaxis] 
            lhsbxH = self.convert_to_sphere(lhsbNH.ravel())*lhsbw[:,numpy.newaxis] 
            lhsbx  = lhsbxH.ravel()
            #lhsbx  = self.convert_to_sphere(lhsbx).ravel() # Map onto unit sphere because solution is located on ellips
            Error  = topo.integral('( ellips )^2 d:u' @ ns, degree=8*self.__degree)  
            
        treelog.info('The curve L\u00b2-error norm: {}'.format(Error.eval(lhs=lhs)**0.5))
        return lhsbx, lhsbw

    def get_geom_constrain(self, constrain, map2sphere=False):
        """  Convert the geometrical constrains to useable Nutils (namespace) variables.
                
        :param constrain: Array containing the constrain type and corresponding parameters:
                           'plane' : [a,b,c,d], based on a*x+b*x+c*x=d
                           'ellips': [S,R], where S = origin and R = [Rx,Ry,Rz]  
        """        
        if constrain[0] == 'plane':
            constrain_copy = constrain[1].copy()
            if map2sphere: # scale ellips plane to unit sphere plane, will be converted back to ellips at later point
                # Find point on plane
                C_old = constrain_copy[-1]  # rhs constant
                N_old = constrain_copy[:-1] # normal vector
                index = -1 if len(numpy.nonzero(N_old)[0]) > 1 else numpy.nonzero(N_old)[0][0]
                if index == 0: # the x-component of the normal vector is equal to 1, other values are 0
                    X = numpy.array([ ( C_old - N_old[1]*C_old - N_old[2]*C_old ) / N_old[0], 0, 0 ]) 
                elif index == 1: # the y-component of the normal vector is equal to 1, other values are 0
                    X = numpy.array([ 0, ( C_old - N_old[0]*C_old - N_old[2]*C_old ) / N_old[1], 0 ])                 
                else: # If index = 2 or index is something else = None
                    X = numpy.array([ C_old, C_old,  ( C_old - N_old[0]*C_old - N_old[1]*C_old ) / N_old[2] ])
                N_new = ( constrain_copy[:-1]*self.__R )/ numpy.linalg.norm( constrain_copy[:-1]*self.__R ) # rescale the normal vector
                C_new = numpy.dot(N_new, X / self.__R) # determine new rhs constant
                # Specify mapped plane information
                constrain_copy[:-1] = N_new 
                constrain_copy[-1]  = C_new
                
            geomcons = ''
            for i, val in enumerate(constrain_copy):
                sign = '-' if numpy.sign(val) == -1 else '+'
                if i == 0 and sign == '+':
                    sign = ''
                if i != len(constrain_copy)-1:
                    geomcons += sign + ' {} '.format(abs(val)) + 'x_' + str(i) + ' '
                else:
                    sign = '-' if numpy.sign(val) != -1 else '+'
                    geomcons += sign + ' {} '.format(abs(val))  
            return geomcons    
            #return '{} x_0 + {} x_1 + {} x_2 - {}'.format(*constrain[1])
        elif constrain[0] == 'ellips':
            geomcons = ''
            for i, (c,r) in enumerate(zip(constrain[1], constrain[2])):
                if c >= 0:
                    geomcons += (' ( x_' + str(i) + ' - {} )^2 / ( {} )^2 +').format(c,r)
                else:
                    geomcons += (' ( x_' + str(i) + ' + {} )^2 / ( {} )^2 +').format(abs(c),r)
            return geomcons[:-1] + '- 1'
        else:
            raise ValueError('The specified constrain type {}, does not exist'.format(constrain[0]))
        
    
    def split_lhs(self, lhs, nelems):
        """  Split the computed (Nutils) lhs array into a control point and weights array.
                
        :param lhs:    Array containing both cps and weight (concatenated)
        :param nelems: Number of elems, default=1
        """     
        if isinstance(nelems, int):
            indx = (self.__degree+nelems)
        else:
            indx = (self.__degree+nelems[0])*(self.__degree+nelems[1])
        lhsx = lhs[:3*indx]
        lhsb = lhs[3*indx:4*indx]    
        return lhsx, lhsb

    def convert_to_sphere(self, Xell):
        """  Map cartesian points on ellipsoid to unit sphere.
        :param Xell: n x 3 array, containing ellipsoid coordinates to be mapped onto unit sphere
        """ 
        return numpy.dot(Xell.reshape(-1,3) - self.__S ,numpy.diag([1/i for i in self.__R])) 

    # def convert_to_ellips(self, Xsph):
    #     """  Map cartesian points on unit sphere to ellipsoid.
    #     :param Xsph: n x 3 array, containing unit sphere coordinates to be mapped onto the ellipsoid
    #     """ 
    #     return numpy.dot(Xsph.reshape(-1,3) - self.__S ,numpy.diag(self.__R)) + self.__S 


    def solve_inner_curves(self,LHSb): # Solve inner curves and store inner + boundary curves as splipy objects
        """  Solve inner curves of each surface patch and construct a surface from them. Results are stored accordingly
        
        The stored curves and surfaces are all dictionaries containing Splipy objects: self.dic = { patch0: Splipy_object, patch1: Splipy_object, ..}
        :Store self.inner_curves:     Curves to be lofted through, which form the surface: incl. the outer (start and end boundary curves)
        :Store self.boundary_curves:  Boundary curves solved in previous step  
        :Store self.outer_curves:     Start and end curves use for lofting, these are boundary curves defined by the lofting direction
        :Store self.surface:          Actual surface result used for post-processing
                
        :param LHSb: Dictionary containing the control points (cps) and weights (w) of 
                      each boundary curve, which is already solved: LHSb = { (0,1): (cps,w), (2,3): (cps,w), .. }
        """ 
        
        ## First convert and store cps and w arrays as Splipy objects (curves)
        self.inner_curves    = dict.fromkeys([i for i in range(0,len(self.__patch))]) # These are the curves that are used for lofting (incl. 2 patch boundary curves as well)
        self.boundary_curves = dict.fromkeys([i for i in range(0,len(self.__patch))]) # These are the remaing 2 patch curves of the surface which are not used for lofting!
        self.outer_curves    = dict.fromkeys([i for i in range(0,len(self.__patch))])
        self.surface         = dict.fromkeys([i for i in range(0,len(self.__patch))])
        
        basis    = BSplineBasis(order=3, knots=[0,0,0,1,1,1]) # Set default single element quadratic B-spline basis (Splipy)

        for i, patchb in enumerate(self.boundidxp): # Zip boundaries of patch with lofting direction
            self.boundary_curves[i] = []
            self.outer_curves[i]    = []
            self.inner_curves[i]    = []
            
            # Create knot vectors (u and v) direction
            if i in self.__lnelems:
                m = self.__lnelems[i] # nr of elements
            else:
                m = 1 # Default nr of elements
                
            p          = self.__degree+1
            knots      = numpy.concatenate([[0]*p,  numpy.linspace(1/m,1-1/m,m-1) , [1]*(p)  ])

            for j, idx in enumerate(patchb): # patchb is build up as: [ [a,b], [c,d], [a,c], [b,d] ], of an [a,b,c,d]-patch
                                             # Where [a,b] and [c,d] are 0-direction 
                                             # While [a,c] and [b,d] are 1-direction

                cps4D    = numpy.concatenate([LHSb[tuple(idx)][0].reshape(-1,3), LHSb[tuple(idx)][1][:,numpy.newaxis]],axis=1) # CPS + W
                curve    = Curve(basis, controlpoints=cps4D, rational=True) # Create curve
                greville = numpy.array([knots[(i):i+self.__degree].sum()/self.__degree for i in range(1,m+p)])

                if self.loftdir[i] == 0: # if loft direction is in left-right direction              

                    if j <= 1: # First 2 boundaries in patchb are left and right-boundaries (not to be lofted through)
                        curve_bound = curve.clone()    
                        self.boundary_curves[i] += [curve_bound.insert_knot(knots[p:-p])] # Store boundary curves per patch
                    else:
                        curve_outer = curve.clone() # Start and End curves of lofting
                        self.outer_curves[i] += [curve_outer]

                else: # If loft direction is in bot-top direction
                
                    if j > 1: # Last 2 boundaries in patchb are bot and top-boundaries (not to be lofted through)
                        curve_bound = curve.clone()    
                        self.boundary_curves[i] += [curve_bound.insert_knot(knots[p:-p])] # Store boundary curves per patch
                    else:
                        curve_outer = curve.clone() # Start and End curves of lofting
                        self.outer_curves[i] += [curve_outer]            
                
                
            if i in self.__coons: # Create Coons surface if True
                #print(self.outer_curves[i] + self.boundary_curves[i])
                self.surface[i] = surface_factory.edge_curves(self.outer_curves[i] + self.boundary_curves[i]) # Create coons surface

                if self.loftdir[i] == 0:
                    self.surface[i].refine(int(self.__pnelems[i][0]-1),direction=1)
                    self.surface[i].refine(int(self.__pnelems[i][1]-1),direction=0)
                else:
                    self.surface[i].refine(int(self.__pnelems[i][0]-1),direction=0)
                    self.surface[i].refine(int(self.__pnelems[i][1]-1),direction=1)
                
                ## Map all curves and surfaces back to the ellips (from the unit sphere)
                
                self.surface[i].scale(self.__R[0],self.__R[1],self.__R[2]) 
                self.surface[i].translate(self.__S)
                for mc in range(len(self.outer_curves[i])):
                
                    self.outer_curves[i][mc].scale(self.__R[0],self.__R[1],self.__R[2]) # Outer curves are a part of inner, so dont scale it twice
                    self.boundary_curves[i][mc].scale(self.__R[0],self.__R[1],self.__R[2])                  
                    self.outer_curves[i][mc].translate(self.__S) # Translate back
                    self.boundary_curves[i][mc].translate(self.__S)
                    
            # elif i in self.cons_inner: # some inner curves are constrained
            #     # Assume Xstart and Xend 
                
            else: # Create lofted surface instead
                ## Evaluate at greville points and create inner curves
                Xstart = self.boundary_curves[i][0].evaluate(greville[1:-1])
                Xend   = self.boundary_curves[i][1].evaluate(greville[1:-1])
                # Construct paires
                Xinner = numpy.concatenate([ Xstart[:,numpy.newaxis], Xend[:,numpy.newaxis] ],axis=1)
                
                ## Solve the inner curves
                self.inner_curves[i] = [self.outer_curves[i][0]]
                w1 = self.get_boundary_weight(self.boundary_curves[i][0]) # Calculate desired weigths = array
                w2 = self.get_boundary_weight(self.boundary_curves[i][1])
                for t, inner in enumerate(Xinner):
                    #weight   = [self.boundary_curves[i][0][t+1][-1] , self.boundary_curves[i][1][t+1][-1]]   
                    weight   = [w1[t+1],w2[t+1]] # Skip the first weight because it does not belong to inner curves (same for last)
                    with treelog.context('Inner curve (patch{}) {}/{}'.format(i,t+1,len(Xinner))):
                        CPSi, Wi = self.get_lhs_boundary(inner, weight=weight, IDX=False)
                    cps4Di   = numpy.concatenate([CPSi.reshape(-1,3), Wi[:,numpy.newaxis]],axis=1) # CPS + W
    
                    curvei = Curve(basis, controlpoints=cps4Di, rational=True) # Create curve
                    self.inner_curves[i] += [curvei] # Insert knots, such that it corresponds with outer_curves .insert_knot(knots[self.loftdir[i]][p:-p])
                self.inner_curves[i] += [self.outer_curves[i][1]]
                
                ## Create lofted surface through these inner curves
                self.surface[i] = self.loftsurface(self.boundary_curves[i], self.inner_curves[i])    
               
                ## Only allow refinement if it is uniformly possible
                if self.loftdir[i] == 0:
                    ref0 = self.__pnelems[i][0] / ((1-self.loftdir[i])*(m - 1) + 1)
                    self.surface[i].refine(int(ref0-1),direction=1)
                    ref1 = self.__pnelems[i][1] / (self.loftdir[i]*(m-1) + 1) 
                    self.surface[i].refine(int(ref1-1),direction=0)
                else:
                    ref0 = self.__pnelems[i][0] / ((1-self.loftdir[i])*(m - 1) + 1)
                    self.surface[i].refine(int(ref0-1),direction=0)
                    ref1 = self.__pnelems[i][1] / (self.loftdir[i]*(m-1) + 1) 
                    self.surface[i].refine(int(ref1-1),direction=1) 
    
                    
                ## Map all curves and surfaces back to the ellips (from the unit sphere)
                self.surface[i].scale(self.__R)
                self.surface[i].translate(self.__S)
                for mc in range(len(self.inner_curves[i])):
                    
                    self.inner_curves[i][mc].scale(self.__R[0],self.__R[1],self.__R[2]) 
                    self.inner_curves[i][mc].translate(self.__S)
                for mc in range(len(self.boundary_curves[i])):
            
                    #self.outer_curves[i][mc].scale(self.__R[0],self.__R[1],self.__R[2]) # Outer curves are a part of inner, so dont scale it twice
                    self.boundary_curves[i][mc].scale(self.__R[0],self.__R[1],self.__R[2])   
                    self.boundary_curves[i][mc].translate(self.__S)
        return

    def loftsurface(self, boundary_curve, inner_curve):
        """  Create surface by lofting through the inner curves and constraining the free-boundaries by the boundary curves.
        
        The curves are all dictionaries containing splipy objects:  self.dic = { patch0: Splipy_object, patch1: Splipy_object, ..}
        :param inner_curve:     Curves to be lofted through, which form the surface: incl. the outer (start and end boundary curves)
        :param boundary_curve:  Boundary curves solved in previous step, correspond to free boundary  
                
        :return surface: Splipy object that is used to describe a surface
        """ 
        
        surface = my_loft_function(inner_curve) #surface_factory.loft(inner_curve)
        degree = 2
        m = len(inner_curve[0].knots()[0]) - 1 # nr of elements
        ## Constrain, fix free-boundary splines to boundary curves
        for idx, (bot, top) in enumerate(zip(boundary_curve[0][1:-1],boundary_curve[1][1:-1])):
            ib = (idx+1)*(m+degree)
            it = ib + m + degree - 1
            surface[it] += top - surface[it]
            surface[ib] += bot - surface[ib]
        return surface             
    
    def get_boundary_weight(self,curve):
        """  Determine the weights of the start and end cps of the inner curves, such that after lofting (interpolation in 4D space), 
             the respective inner curve its start and end weight corresponds to the boundary curve exactly (incl the cps location.)
             
        :param curve:     Splipy curve representing the inner curve of a patch
              
        :return w: Required weights of the boundary curve, to be assigned to the correct inner curve
        """ 
        basis = curve.bases[0]
        v   = basis.greville() 
        Nv  = basis(v) # Get values of basis functions at greville points
        w   = numpy.asarray([ numpy.dot(Nv,curve[:,-1]) ]).reshape(-1) # multiplyt/dot product the basis function with the weight,
                                                                       # this way we can extract the weight of the curve such that it corresponds 
                                                                       # to a desired weight value upon lofting (curves are interpolated)
        return w
    
    def linterp_bound(self, X, nelems, δ=0):
        """  Linearly interpolate input by number of elements.
                
        :return linearly interpolation between two points in space:
        """ 
        
        n    = self.__degree + nelems # Nr of control points, incl boundary
        subd = n - 2           # Number of subdivisions
        dx   = numpy.diff(X, axis=0)[0]/(subd+1) 
        x    = [i*dx for i in range(1,subd+1)]
        ## Normalized correction
        v = numpy.array([0,0,1])
        n = numpy.cross( dx, v ) / ( numpy.linalg.norm( numpy.cross( dx, v ) ) ) 
        δx = δ*n # Correction, in case a better initial guess is desired
        return numpy.concatenate([X[:1],x + X[0],X[1:]], axis=0)  

    def define_namespace(self, topo, geom, cons=None):
        """  Define the Nutils namespace variables used for computing curves.
        The curves are all dictionaries containing splipy objects:  self.dic = { patch0: Splipy_object, patch1: Splipy_object, ..}
        :param topo:  Nutils topology
        :param geom:  Nutils geometry  
        :param cons:  Constrain type   
            
        :return ns: Namespace
        """ 
        ## Define namespace and Spline basis
        ns      = function.Namespace()
        if cons==None:
            Φ, ω, Λ = function.chain([topo.basis('spline',degree=self.__degree).vector(3), 
                                      topo.basis('spline',degree=self.__degree), 
                                      topo.basis('spline',degree=self.__degree)])
        elif cons[0]=='plane':
            # Φ, ω, Λ = function.chain([topo.basis('spline',degree=self.__degree).vector(3), 
            #                           topo.basis('spline',degree=self.__degree), 
            #                           topo.basis('spline',degree=self.__degree)])
            Φ, ω, Λ, Γ = function.chain([topo.basis('spline',degree=self.__degree).vector(3), 
                                          topo.basis('spline',degree=self.__degree), 
                                          topo.basis('spline',degree=self.__degree),
                                          topo.basis('spline',degree=self.__degree)])
            ns.Γ    = Γ
            ns.γ    = 'Γ_n ?lhs_n'
        else:
            Φ, ω, Λ = function.chain([topo.basis('spline',degree=self.__degree).vector(3), 
                                      topo.basis('spline',degree=self.__degree), 
                                      topo.basis('spline',degree=self.__degree)])
            

            
        ns.Φ    = Φ
        ns.Λ    = Λ
        ns.ω    = ω
        ns.u    = geom
        ns.xh_i = 'Φ_ni ?lhs_n'# B-spline
        ns.w    = 'ω_n ?lhs_n' # Weights
        ns.λ    = 'Λ_n ?lhs_n' # Lagrangian multiplier
        ns.x_i  = 'xh_i / w'   # NURBS
        ns.Rm   = function.diagonalize([1/i for i in self.__R]) # Diagonalized matrix with inverse radii
        ns.R    = function.asarray(list(self.__R))
        ns.S    = function.asarray(list(self.__S))
        #ns.S    = function.asarray([-1,0,0]) 
        ns.sqellips_i = 'Rm_ij ( x_j - S_j )'
        ns.ellips     = 'sqellips_i sqellips_i - 1' # Ellips
        ns.sphere     = 'x_i x_i - 1'               # unit sphere, R=1
        
        ## Define ellipsoid equations
        #ns.ellipsh     = self.__shape #' ( xh_0 / R_0 )^2 + ( xh_1 / R_1 )^2 + ( xh_2 / R_2 )^2 - 1' # homogenous coord. ellips
        #ns.ellips      = ' ( ( x_0 - S_0 ) / R_0 )^2 + ( ( x_1 - S_1 ) / R_1 )^2 + ( ( x_2 - S_2 ) / R_2 )^2 - 1' 
        #ns.ellipsright = self.__shape #' ( ( x_0 + 0.5 ) / R_0 )^2 + ( x_1 / R_1 )^2 + ( x_2 / R_2 )^2 - 1'
    
        #ns.ε  = function.localgradient(ns.x, __ndims=len(geom))# 1 or 2 __ndim
        #ns.εh = function.localgradient(ns.xh, __ndims=len(geom))
        
        ns.ε  = function.grad(ns.x, geom, ndims=len(geom))
        ns.εh = function.grad(ns.x, geom, ndims=len(geom))
        
        return ns    
    
    def visualize(self,*param): # Does only accept dictionaries       
        """  Visualize the results in 3D plots. Input can be of the following form:
        param: - Single Splipy curve
               - Single Splipy surface
               - Array containing multiple Splipy curves
               - DIctionary containing (multiple) Splipy curves or Surfaces
        """ 
            
        ## Parametric values
        u   = numpy.linspace(0,1,20)
        v   = numpy.linspace(0,1,20) # surface.knots(1)[-1]
        
        # Specify ellips
        φ = numpy.linspace(0, 2* numpy.pi, 30)
        θ = numpy.linspace(0, numpy.pi, 30)
        x = self.__R[0] * numpy.outer(numpy.cos(φ), numpy.sin(θ)) + self.__S[0]
        y = self.__R[1] * numpy.outer(numpy.sin(φ), numpy.sin(θ)) + self.__S[1]
        z = self.__R[2] * numpy.outer(numpy.ones(numpy.size(φ)), numpy.cos(θ)) + self.__S[2]
        f = 1.1*self.__R[-1]
        
        ## Plot curve
        fig = plt.figure(figsize=(8, 6))
        ax = fig.gca(projection='3d')     
        
        # loop over input tuple()
        for inp in param:
            # Check if it is a dicitonary
            if isinstance(inp, dict):
                # Loop over the input of the dictionary
                for dic in inp:
                    objectt = inp[dic]
                    try: # it is a surface
                        kn = objectt.knots()
                        surf = objectt.evaluate(u,v)
                        ax.plot_wireframe(surf[:,:,0], surf[:,:,1], surf[:,:,2])                        
                    except AttributeError: # It is a curve
                        for ic in objectt:
                            curve = ic.evaluate(u)
                            ax.plot(curve[:,0], curve[:,1], curve[:,2])
            else: # If it is not a dictionary but a list or something       
             # Check if surface is given 
                try: # it is a surface OR a curve
                    kn   = inp.knots()
                    if len(kn)>1: # It is indeed a surface
                        surf = inp.evaluate(u,v)
                        ax.plot_wireframe(surf[:,:,0], surf[:,:,1], surf[:,:,2])   
                    else: # It is a curve
                        curve = inp.evaluate(u)
                        ax.plot(curve[:,0], curve[:,1], curve[:,2])                         
                except AttributeError: # It is a curve
                    if len(inp) == 1:
                        curve = inp.evaluate(u)
                        ax.plot(curve[:,0], curve[:,1], curve[:,2])                        
                    else:
                        for ic in inp:
                            curve = ic.evaluate(u)
                            ax.plot(curve[:,0], curve[:,1], curve[:,2])               
        ax.plot_wireframe(x, y, z, color='darkgrey')
        #ax.view_init(20, 45)
        #plt.axis('off')
        #plt.grid(b=None)
        ax.set_xlim(-f, f)
        ax.set_ylim(-f, f)
        ax.set_zlim(-f, f)
        plt.show() 
        return
    
    def get_cps_w(self,):
        """  Obtain the control points and weights useable for Nutils purpose (separate concatenated arrays)
                
        :return cps: n x 3 array containing all the control points
        :return w:   1 x n array containing all the weights corresponding the cps
        """ 
        ## This function should be run only ones(!) due to the swapping command used (otherwise should be run twice)
        if self.run == 0:
            self.run += 1
            ## Concatenate all surfaces and weightsin case of mulitple patches
            for i in range(len(self.__patch)):
                if self.loftdir[i] == 0: # Check if loftd corresponds to first patch direction, if so, the cps should be swapped
                    homog_coord = self.surface[i].swap()[:].reshape(-1,4) # Swap directions because of loft curves (not same direction as patch)
                else:
                    homog_coord = self.surface[i][:].reshape(-1,4)
                if i == 0:
                    w   = homog_coord[:,-1]
                    cps = homog_coord[:,:-1]/w[:,numpy.newaxis]
                else:
                    
                    wi   = homog_coord[:,-1]
                    cpsi = homog_coord[:,:-1]/wi[:,numpy.newaxis]
                    w   = numpy.concatenate([w  , wi])
                    cps = numpy.concatenate([cps, cpsi], axis=0)     
            self.cps = cps
            self.w   = w
            return cps, w
        else:
            return self.cps, self.w
    
    def get_topo_geom(self,cps,w,filename=None):
        """  Construct the Nutils NURBS topology and geometry.
                
        :param cps:      n x 3 array containing all the control points
        :param w:        1 x n array containing all the weights corresponding the cps
        :param filename: Name of the file (string) used to save the geometry as .vtk file
        """ 
        topo, lingeom = mesh.multipatch(patches=self.__patch, patchverts=self.__patchverts, nelems=self.nelems)

        bsplines = topo.basis('spline', degree=2, patchcontinuous=False)
        weight   = bsplines.dot(w)
        geom     = bsplines.vector(3).dot((cps*w[:,numpy.newaxis]).ravel())/weight
        
        print('Value should be zero (or close to): Value=', abs(topo.interfaces['interpatch'].sample('bezier', 5).eval(function.jump(geom))).max())
        bezier     = topo.sample('bezier', 15)
        pids       = topo.basis('patch').dot(numpy.arange(len(self.__patch)))
        GEOM, PIDS = bezier.eval([geom,pids])
        
        if filename!=None:
            export.vtk(filename, bezier.tri, GEOM, PatchID=PIDS) 
        return topo, geom
    
    def save(self,filename): # bare in mind that we are mainly interested in saving the actual 3D solid bi-ventricle geometry
        """  Save the result/geometry to a .vtk file with the name: filename.""" 
        cps, w     = self.get_cps_w()
        topo, geom = self.get_topo_geom(cps,w, filename=filename)
        return
