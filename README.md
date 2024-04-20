# ANSYS-USERELEMENT-PHFLD
Phase field fracture realized by ANSYS custom USERELEMENT subroutine. In order to obtain numerical solutions of the coupled system of partial differential equations using the finite element method an additional degree of freedom is introduced for the element. The phase field model is implemented by means of Ansys UEL subroutine which allows for user-defined computation of the element tangent stiffness matrices and the nodal force vectors. We consider isoparametric 2D quadrilateral elements (linear and quadratic) with 3 degrees of freedom per node, i.e. UX,UY and φ, and four integration points. The framework is general, and is supported by addressing several classical 2D boundary value problems as well as the ductile fracture and 3D surface flaws behaviors of particular interest. The framework can be very easily extended to modern material models like strain gradient plasticity, transgranular and intergranular damage mechanisms under cyclic loading conditions.  





If you using this code for your research problems please refer:

D.Kosov1, A. Tumanov1, V. Shlyannikov "A phase field approach implementation in ANSYS and numerical examples" (2024)


More complex version of this library for 3D problems you can find here: https://github.com/DmitryKosov1/Phase-field-fracture-ANSYS
