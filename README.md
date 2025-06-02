# ANSYS-USERELEMENT-PHFLD
The phase field fracture approach is implemented using ANSYS user programmable features. In order to obtain numerical solutions of the coupled system of partial differential equations using the finite element method an additional degree of freedom is introduced for the element. The phase field model is implemented by means of Ansys UEL subroutine which allows for user-defined computation of the element tangent stiffness matrices and the nodal force vectors. We consider isoparametric 2D elements with 3 degrees of freedom per node, i.e. UX,UY and φ, and four integration points. The framework is general, and is supported by addressing several classical 2D boundary value problems as well as the ductile fracture and 3D surface flaws behaviors of particular interest. The framework can be very extended to modern material models like strain gradient plasticity, transgranular and intergranular damage mechanisms under cyclic loading conditions. 


If you using this code for your research purposes please refer following paper:

[A phase field approach implementation in ANSYS and numerical examples](https://doi.org/10.3221/IGF-ESIS.70.08) 

**Some implemented features**
- Elastic, plastic and creep parts of the total strain energy density
- Load as pressure
- The plastic energy dissipation
- Dependence of critical energy release rate from temperature (power law)



In SOURCE folder you can find the APDL script example of 2D cracked plate under tension.  



Parallel version of this library for 3D cases you can find here: https://github.com/DmitryKosov1/Phase-field-fracture-ANSYS
