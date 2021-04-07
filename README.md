# spatial_v2_extended

New algorithms include:
* Methods to compute the Coriolis matrix ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/CoriolisMatrix.m)) [all systems] and Christoffel symbols ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/Christoffel.m)) [Single DoF joints only] ([paper](https://arxiv.org/abs/2010.01033))
* Methods for assessing identifiability ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/identifiability), [paper](https://arxiv.org/abs/1711.03896)) 
* Methods for computing the partial derivatives of Inverse Dynamics ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/derivatives/ID_derivatives.m)) [paper forthcoming]

New features include:
* Extensions of most algorithms (RNEA, ABA, CRBA, etc.) to address dynamic effects from motor rotors ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/dynamics))
* Updated most algorithms (RNEA, ABA, CRBA, etc.) to handle multi-DoF joints (e.g., spherical or floating base) without needing specialized versions of the algorithms ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/dynamics))
* Regressor calculation algorithms ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/regressor))
* Variety of tools for converting between different representations of orientation ([link](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/orientation_tools))
* Partial compatibility for complex-valued input arguments toward support of complex-step derivative comptuations in unit tests.
