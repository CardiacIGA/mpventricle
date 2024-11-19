# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 13:12:05 2022

@author: s146407
"""


## Adjusted loft function when compared to the standard SpliPy version
## Constrained the loft to only quadratic and not cubic splines (cubis was default for SpliPy)

from splipy import BSplineBasis, Curve, Surface
from splipy.surface_factory import edge_curves
import numpy as np


def my_loft_function(*curves):
    if len(curves) == 1:
        curves = curves[0]

    # clone input, so we don't change those references
    # make sure everything has the same dimension since we need to compute length
    curves = [c.clone().set_dimension(3) for c in curves]
    if len(curves)==2:
        return edge_curves(curves)
    elif len(curves)==3:
        # can't do cubic spline interpolation, so we'll do quadratic
        basis2 = BSplineBasis(3)
        dist  = basis2.greville()
    else:
        x = [c.center() for c in curves]

        # create knot vector from the euclidian length between the curves
        dist = [0]
        for (x1,x0) in zip(x[1:],x[:-1]):
            dist.append(dist[-1] + np.linalg.norm(x1-x0))

        # using "free" boundary condition by setting N'''(u) continuous at second to last and second knot
        knot = [dist[0]]*4 + dist[2:-2] + [dist[-1]]*4
        basis2 = BSplineBasis(4, knot)

        ## Making use of the greville points (assuming start/end of curvesare located here), 2nd order
        knot = [0.]*3 + [1-i*1/(len(curves)-2) for i in range((len(curves)-2)-1,0,-1)] + [1.]*3
        basis2 = BSplineBasis(3, knot)
        dist  = basis2.greville()
        # basis2 = BSplineBasis(3, [0., 0., 0., 0.33333333, 0.66666667,1., 1., 1.])
        # dist  = basis2.greville()
        
        ## Making use of the distance
        # knot = [dist[0]]*3
        # for i,j in zip(dist[1:-2],dist[2:]):
        #     knot+=[i+(j-i)*0.5]
        # knot += [dist[-1]]*3    
        # basis2 = BSplineBasis(3, knot)

        
    n = len(curves)
    for i in range(n):
        for j in range(i+1,n):
            Curve.make_splines_identical(curves[i], curves[j])

    basis1 = curves[0].bases[0]
    m      = basis1.num_functions()
    u      = basis1.greville() # parametric interpolation points
    v      = dist              # parametric interpolation points

    # compute matrices
    Nu     = basis1(u)
    Nv     = basis2(v)
    Nu_inv = np.linalg.inv(Nu)
    Nv_inv = np.linalg.inv(Nv)

    # compute interpolation points in physical space
    x      = np.zeros((m,n, curves[0][0].size))
    for i in range(n):
        x[:,i,:] = Nu @ curves[i].controlpoints

    # solve interpolation problem
    cp = np.tensordot(Nv_inv, x,  axes=(1,1))
    cp = np.tensordot(Nu_inv, cp, axes=(1,1))
    
    # re-order controlpoints so they match up with Surface constructor
    cp = cp.transpose((1, 0, 2))
    cp = cp.reshape(n*m, cp.shape[2])

    return Surface(basis1, basis2, cp, curves[0].rational)