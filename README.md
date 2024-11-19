# NURBS geometry

This project contains the code to construct a template (multipatch) left- or bi-ventricle heart geometry using NURBS. 
Use is made of a multipatch approach, where the user is able to modify the position and numbers of control points (CPS) and weights.
The obtained geometry can be converted to a usable Nutils topology and geometry required for the cardiac IGA model.

## Files

The repository makes use of the following file structure:
The geometry construction module '[vtnurbs](https://github.com/CardiacIGA/multipatch-ventricle/tree/main/vtnurbs)' (ventricle template NURBS), a set of [examples](https://github.com/CardiacIGA/multipatch-ventricle/tree/main/examples) that show the usage of the module, images of two left and bi-ventricle geometry results, and some miscellaneous [utility](https://github.com/CardiacIGA/multipatch-ventricle/tree/main/utils) files.


The ventricle template NURBS ([vtnurbs](https://github.com/CardiacIGA/multipatch-ventricle/tree/main/vtnurbs)) consists of:
- [surface_generator.py](https://github.com/CardiacIGA/multipatch-ventricle/blob/main/vtnurbs/surface_generator.py)            : Calculates the position and weights of a single patch surface on an ellipsoid given specific constrains. It is only limited to quadratic NURBS.
- [cardiac_geometry.py](https://github.com/CardiacIGA/multipatch-ventricle/blob/main/vtnurbs/cardiac_geometry.py)    : Combines the calculated surface information (cps, weights, knotvectors etc.) into a modifiable solid Nutils topology and geometry, 
- [loftquadratic.py](https://github.com/CardiacIGA/multipatch-ventricle/blob/main/vtnurbs/loftquadratic.py)        : Modified Splipy function for lofting quadratic NURBS surfaces. The original Splipy integration elevated the spline degree which is undesirable.

And [example](https://github.com/CardiacIGA/multipatch-ventricle/tree/main/examples) scripts:
- [leftventricle_example.py](https://github.com/CardiacIGA/multipatch-ventricle/blob/main/examples/leftventricle_template.py) : Example script which shows how to generate an idealized left ventricle NURBS geometry,
- [biventricle_example.py](https://github.com/CardiacIGA/multipatch-ventricle/blob/main/examples/biventricle_template.py) : Example script which shows how to generate an idealized bi-ventricle NURBS geometry,
- [geometry_variations_paper.py](https://github.com/CardiacIGA/multipatch-ventricle/blob/main/examples/geometry_variations_paper.py) : Example script which shows how the bi-ventricle variations of the paper by Willems et al. are generated,
- [modify_leftventricle.py](https://github.com/CardiacIGA/multipatch-ventricle/blob/main/examples/modify_leftventricle.py) : Example script in which the left-ventricle its cps and weights are modified, refined and saved.

## Additional information
This repository is part of the COMBAT-VT project (https://combatvt.nl/), and is maintained by R. Willems during his PhD project.

## Citing this module
Please consider citing this module when using it:

Willems, R., Janssens, K. L., Bovendeerd, P. H., Verhoosel, C. V., & van der Sluis, O. (2024). An isogeometric analysis framework for ventricular cardiac mechanics. Computational Mechanics, 73(3), 465-506. DOI: https://doi.org/10.1007/s00466-023-02376-x