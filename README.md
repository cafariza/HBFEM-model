# HBFEMmodel
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
  Download all the available files in ## Folder
  
Step 3:
  Run MAIN.m 



# References

Franco, C., Chesnais, C., Semblat, J. F., Giry, C., & Desprez, C. (2022). Finite element formulation of a homogenized beam for reticulated structure dynamics. Computers and Structures, 261–262. https://doi.org/10.1016/j.compstruc.2021.106729

