# spatial_v2_extended

New algorithms include:
* Methods to compute the Coriolis matrix ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/CoriolisMatrix.m)) [all systems] and Christoffel symbols ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/Christoffel.m)) [Single DoF joints only] ([paper](http://dx.doi.org/10.1115/1.4051169))
* Methods for assessing identifiability ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/identifiability), [paper](https://arxiv.org/abs/1711.03896)) 

New features include:
* Extensions of most algorithms (RNEA, ABA, CRBA, etc.) to address dynamic effects from motor rotors ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/dynamics))
* Updated most algorithms (RNEA, ABA, CRBA, etc.) to handle multi-DoF joints (e.g., spherical or floating base) without needing specialized versions of the algorithms ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/dynamics))
* Regressor calculation algorithms ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/regressor))
* Variety of tools for converting between different representations of orientation ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/orientation_tools))
* Methods for computing the partial derivatives of Inverse Dynamics ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/derivatives/ID_derivatives.m)) 
* Partial compatibility for complex-valued input arguments toward support of complex-step derivative comptuations in unit tests (including the complex step on matrix Lie groups [link](https://ieeexplore-ieee-org.proxy.library.nd.edu/abstract/document/8957301)).
