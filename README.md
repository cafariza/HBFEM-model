# HBFEMmodel_ElementaryMatrices
Here you will find the MatLab Code to built the HBFEM model of the steel structure given in (Franco et al., 2022). You will obtain
- Dynamic properties (eigen frequencies and mode shapes)
- Top roof displacement from a time history analysis using Newmark

![Modes_num_ana_fem_he](https://user-images.githubusercontent.com/30736872/157275778-f2a449c6-36db-465a-a161-74395c2b4741.svg)

The HBFEM model is a 1D enriched model discretized in finite elements with 3 degrees of freedom at each node. The construction of this numerical formulation is well detailed in (Franco et al., 2022) 
Through a representative example in the framework of linear elasticity, it is shown that this homogenized beam finite element (HBFEM) solution
recovers the analytical results and is close to the finite element solution of the detailed structural model. This new formulation simplifies the implementation of the beam-like model in complex configurations, enables parametric studies, and could be easily used for a wide range of applications in the field of dynamics of solids and structures, other than free vibration analysis, at reduced computational costs.

The implementation of the beam models and the subsequent homogenized beam finite element (HBFEM) model  correctly characterizes the transverse dynamics of real buildings if the structure fulfills specific requirements. Firstly, the analyzed structure must be vertically regular in mass, stiffness and strength to respect the periodicity condition required by the homogenization. Secondly, the building must be tall enough to respect the scale separation condition. Structures of at least five stories (ε ≤ 0.3) fulfill this condition. Finally, the dynamic characterization must be limited to the analysis of the vibration modes with a sufficiently long wavelength. For each vibration type, the maximum number of modes which can be modelled with this approach is approximately N/3, where N is the number of stories. This condition is typically valid in earthquake engineering studies where the lowest frequencies are of interest. In this framework, a lot of the existing reinforced concrete (RC) and steel buildings seem to easily fulfill the stated conditions. 


# How to use the MatLab code

Step 1: 
  Verify you have a MatLab version 2015a or later with the Symbolic Toolbox installed.
  
Step 2:
  Download all the available files in the folder: "HBFEM_MatLabCode"
  
Step 3:
  Place all the donwloaded files in the same folder 

Step 4: 
  Run MAIN_ExampleFrancoetal2022CS.m


# References

Franco, C., Chesnais, C., Semblat, J. F., Giry, C., & Desprez, C. (2022). Finite element formulation of a homogenized beam for reticulated structure dynamics. Computers and Structures, 261–262. https://doi.org/10.1016/j.compstruc.2021.106729

# Authors 


Franco Carolina -  Univ Gustave Eiffel, GERS-SRO and MAST-EMGCU, F-77447 Marne-la-Vallée, France
francoariza@edu.univ-eiffel.fr, carolina.franco@geodynamique.com

Chesnais Céline - Univ Gustave Eiffel, GERS-SRO, F-77454 Marne - la- Vallée, France,
celine.chesnais@univ-eiffel.fr

Semblat Jean-François - IMSIA (UMR9219), CNRS, EDF, CEA, ENSTA Paris, Institut Polytechnique de
Paris, 91120 Palaiseau, France, jean-francois.semblat@ensta-paris.fr

Desprez Cédric - Univ. Lyon, INSA-Lyon, GEOMAS, F-69621, Villeurbanne, France,
cedric.desprez@insa-lyon.fr

Giry Cédric - Université Paris-Saclay, ENS Paris-Saclay, CentraleSupélec, CNRS, LMPS - Laboratoire de
Mécanique Paris-Saclay, 91190, Gif-sur-Yvette, France, cedric.giry@ens-paris-saclay.fr

This work is part of the Ph.D. project entitled: “Multiscale modeling of the seismic response of buildings: Coupling between Homogenization and Multifiber element methods” developed at University Gustave Eiffel, France in the framework of the research initiative 3-COP2017.

# Licence

SPDX-License-Identifier:  EUPL-1.2

# Versions

Creation date: July 2020
Last modified: April 2022


