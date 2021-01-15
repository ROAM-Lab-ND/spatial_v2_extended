# spatial_v2_extended

Extensions include:
* Extensions of most algorithms (RNEA, ABA, CRBA, etc.) to address dynamic effects from motor rotors [(link)](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/dynamics)
* Updated most algorithms (RNEA, ABA, CRBA, etc.) to handle multi-DoF joints (e.g., spherical or floating base) without needing specialized versions of the algorithms [(link)](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/dynamics)
* New algorithms to compute the Coriolis matrix [(link)](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/CoriolisMatrix.m) [all systems] and Christoffel symbols [(link)](https://github.com/ROAM-Lab-ND/spatial_v2_extended/blob/main/v3/dynamics/Christoffel.m) [Single DoF joints only]
* Regressor calculation algorithms [(link)](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/regressor)
* Algorithms for assessing identifiability [(link)](https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/identifiability)
* Variety of tools for converting between different representations of orientation [https://github.com/ROAM-Lab-ND/spatial_v2_extended/tree/main/v3/orientation_tools]((link))
