            SUBROUTINE obrat(N,A,C)
            ! input ...
            ! a(n,n) - array of coefficients for matrix A
            ! n      - dimension
            ! output ...
            ! c(n,n) - inverse matrix of A
            ! comments ...
            ! the original matrix a(n,n) will be destroyed 
            ! during the calculation
            !===========================================================
            implicit none 
            integer n
            double precision A(n,n), C(n,n), OOO(n,2*n), lkl(n,2*n) 
            double precision ll,tt, gaus(n,2*n)
            integer i, j, pp, kk, kl
                  gaus = 0.0d0
                  do i = 1, N
                      do j = 1, N
                          gaus(i,j) = A(i,j)
                      end do
                  end do
                  do i=1, N
                      gaus(i,i+N) = 1.0d0
                  end do
                  OOO = gaus
                  do i = 1,N
                      lkl = OOO
                      kl = i+1
                      If(OOO(i,i).EQ.0.0d0) then
                          do while (OOO(i,i).EQ.0.0d0)
                              kl = kl
                          do j = 1, 2*N
                              OOO(i,j) = lkl(kl,j)
                              OOO(kl,j) = lkl(i,j)
                          end do
                          If(OOO(i,i).EQ.0.0d0) then
                              OOO = lkl
                          end if    
                          kl = kl +1
                          end do
                      end if    
                      tt = OOO(i,i)
                      do pp = 1, 2*N
                         OOO(i,pp) = OOO(i,pp)/tt 
                      end do   
                      do kk = 1,N-i
                          ll = OOO(kk+i,i)
                          do j = 1,2*N
                              OOO(kk+i,j) = OOO(kk+i,j)-OOO(i,j)*ll
                          end do
                      end do
                  end do   
c          обратный ход
                  do i = -N,-1
                      tt = OOO(-i,-i)
                      do pp = 1, 2*N
                         OOO(-i,pp) = OOO(-i,pp)/tt 
                      end do   
                      do kk = i+1, -1
                          ll = OOO(-kk,-i)
                          do j = 1,2*N
                              OOO(-kk,j) = OOO(-kk,j)-OOO(-i,j)*ll
                          end do
                      end do
                  end do   
                  do i=1,N
                      do j=1,N
                          C(i,j) = OOO(i,j+N)
                      end do
                  end do    
            end subroutine obrat
            
*deck,UserElem     USERDISTRIB                                     jxw
c Copyright ANSYS.  All Rights Reserved.      
   
      subroutine UserElem (elId, matId, keyMtx, lumpm, nDim, nNodes,
     &                     Nodes, nIntPnts, nUsrDof, kEStress, 
     &                     keyAnsMat, keySym, nKeyOpt, KeyOpt,
     &                     temper, temperB, tRef, kTherm, 
     &                     nPress, Press, kPress, nReal, RealConst, 
     &                     nSaveVars, saveVars, xRef, xCur, 
     &                     TotValDofs, IncValDofs, ItrValDofs,
     &                     VelValDofs, AccValDofs,
     &                     kfstps, nlgeom, nrkey, outkey, elPrint, iott,
     &                     keyHisUpd, ldstep, isubst, ieqitr, timval, 
     &                     keyEleErr, keyEleCnv,
     &                     eStiff, eMass, eDamp, eSStiff,
     &                     fExt, fInt, elVol, elMass, elCG, 
     &                     nRsltBsc, RsltBsc, nRsltVar, RsltVar, 
     &                     nElEng, elEnergy)
     &      
*deck,usermat      USERDISTRIB  parallel                                gal
  
      
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERELEM"::UserElem
*deck,usermat      USERDISTRIB  parallel                                gal

 
c*************************************************************************
c
c *** Primary function: General User Element Subroutine
c *** Note:
c       This routine is completed with an example, see more details later.
c
c
c     PROGRAMMER SHOULD NOT CHANGE ANY PURE INPUT ARGUMENTS (marked by ....,in)!
c
c     elId      (int,sc,in)        element number
c     matId     (int,sc,in)        material number of this element
c     keyMtx    (int,ar(10),in)    matrix and load vector form requests
c                                     0 = not requested, 1 = requested
c                                     see below for more details
c     lumpm     (int,sc,in)        mass matrix format
c                                    = 0 no lumped mass matrix
c                                    = 1 lumped mass matrix
c     nDim      (int,sc,in)        number of dimensions of the problem
c                                       (defined on USRELEM command as NDIM)
c                                    = 2 2D
c                                    = 3 3D
c     nNodes    (int,sc,in)        number of nodes of the element
c                                       (defined on USRELEM command as NNODES)
c     Nodes     (int,ar(nNodes),in)node list of this element 
c     nIntPnts  (int,sc,in)        maximum number of integration points
c                                       (defined on USRELEM command as NINTPNTS)
c     nUsrDof   (int,sc,in)        number of DOFs of this element (matrix and 
c                                     load vector size)
c     kEStress  (int,sc,in)        kEStress 
c                                       (defined on USRELEM command as KESTRESS)
c     keyAnsMat (int,sc,in)        key to indicate if ANSYS material
c                                     routine is going to be called
c                                     (defined on USRELEM command as KEYANSMAT)
c                                     = 0, No
c                                     = 1, Yes
c     keySym    (int,sc,in)        key to indicate if element matrices
c                                     is symmetric
c                                       (defined on USRELEM command as KEYSYM)
c                                     = 0, symmetric
c                                     = 1, unsymmetric
c     nKeyOpt   (int,sc,in)        number of element key options able to be
c                                     used in this routine
c     KeyOpt    (int,ar(nKeyOpt),in) values of element key option defined
c                                     by et or keyopt command for the
c                                     user elements, only the first
c                                     nKeyOpt values are passed in and can
c                                     be used to branch the routine for
c                                     different formulations
c     temper    (dp,ar(nNodes),in) nodal temperatures at current time
c     temperB   (dp,ar(nNodes),in) nodal temperatures at the beginning of this
c                                     incremental step (substep)
c     tRef      (dp,sc,in)         reference temperature
c     kTherm    (int,sc,inout)     input:  flag for thermal loading 
c                                      = 1, Temperatures at nodes are different 
c                                      from the reference temperature, 
c                                      thermal loading might be needed.
c                                      = 0, Temperatures at nodes are the same
c                                      as the reference temperature, 
c                                      thermal loading is not needed.
c                                  output:  flag for thermal strains
c     nPress    (int,sc,in)        number of pressure values for this element
c     Press     (dp,ar(nPress),in) applied elemental face load (pressure)
c     kPress    (int,sc,in)        flag for pressure loading 
c                                      = 1, pressure load is applied and 
c                                      equivalent nodal forces should be 
c                                      calculated
c                                      = 0, no pressure loading
c     nReal     (int,sc,in)        number of real constants
c                                       (defined on USRELEM command as NREAL)
c     RealConst (dp,ar(nReal),in)  user defined real constants 
c     nSaveVars (int,sc,in)        number of saved variables
c                                      (defined on USRELEM command as NSAVEVARS)
c     saveVars  (dp,ar(nSaveVars),inout) user saved variables
c     xRef      (dp,ar(nDim,nNodes),in)
c                                  nodal coordinates in initial configuration
c     xCur      (dp,ar(nDim,nNodes),in)
c                                  nodal coordinates in current configuration
c     TotValDofs (dp,ar(nUsrDof),in) total values of DOFs (displacements) 
c                                     from time = 0
c     IncValDofs (dp,ar(nUsrDof),in) incremental values of DOFs (displacements) 
c                                     for the current step
c     ItrValDofs (dp,ar(nUsrDof),in) iterative values of DOFs (displacements)
c                                     for the current iteration
c                                     (normally needed for debug only)
c     VelValDofs (dp,ar(nUsrDof),in) first time derivatives of DOFs 
c                                             (velocities) (normally not needed)
c     AccValDofs (dp,ar(nUsrDof),in) second time derivatives of DOFs 
c                                          (accelerations) (normally not needed)
c     kfstps    (int,sc,in)        key for the first iteration of first 
c                                     substep of the first load step
c                                     = 1 yes
c                                     = 0 no
c     nlgeom    (int,sc,in)        large deformation key [from nlgeom command]
c                                     = 0 NLGEOM,OFF
c                                     = 1 NLGEOM, ON
c     nrkey     (int,sc,in)        key to indicate a newton-raphson
c                                     (incremental) procedure
c                                     = 0 No
c                                     = 1 Yes
c     outkey    (int,sc,in)        key to indicate if any element output is
c                                     to be placed on the print file or the 
c                                     result file
c                                     = 0 No
c                                     = 1 Yes
c     elPrint   (int,sc,in)        key to indicate if any element output is 
c                                     to be placed on the print file
c                                     = 0 No
c                                     = 1 Yes
c     iott      (int,sc,in)        print output file unit number
c     keyHisUpd (int,sc,in)        key to indicate if history-dependent
c                                    variables need to be updated, like
c                                    equivalent plastic strain, back stress
c                                    etc. since the iteration is already
c                                    converged
c                                     = 0 not converged, don't need to update
c                                         history dependent variables
c                                     = 1 yes, converged, need to update
c                                         history dependent variables
c
c     --- The following 7 variable group can usually be ignored.
c     --- The variables are used for debug, timing, and convergence control.
c     ldstep    (int,sc,in)        current load step number
c     isubst    (int,sc,in)        current substep number
c     ieqitr    (int,sc,in)        current equilibium iteration  number
c     timval    (int,sc,in)        current time value
c     keyEleErr (int,sc,inout)     key to indicate if there is any element 
c                                     formulation error, like negative Jacobian.
c                                     The error could be caused by too
c                                     large incremental step, illegal model.
c                                     = 0 no error (preset value before calling)
c                                     = 1 some error happens. ANSYS will
c                                     decide to stop the analysis or cutback
c                                     the substep (bi-section) based on other
c                                     user input and information at higher
c                                     level.
c     keyEleCnv (int,sc,inout)     key to flag if this element satisfies
c                                     the user defined element convergence
c                                     criterion. 
c                                     = 1, yes, the criterion is satisfied
c                                       or don't have any criterion at all
c                                       it is preset value before calling
c                                     = 0, no, the element doesn't satisfy
c                                       element convergence criterion. If
c                                       this is the case, the iteration will
c                                       not converge even when both force
c                                       and displacement converge 
c       ---- end of 7 variable group -----
c
c                                                                  requested if
c     eStiff(dp,ar(nUsrDof,nUsrDof),inout) stiffness matrix         keyMtx(1)=1
c     eMass (dp,ar(nUsrDof,nUsrDof),inout) mass matrix              keyMtx(2)=1
c     eDamp (dp,ar(nUsrDof,nUsrDof),inout) damping matrix           keyMtx(3)=1
c     eSStiff(dp,ar(nUsrDof,nUsrDof),inout)stress stiffness matrix  keyMtx(4)=1
c     fExt      (dp,ar(nUsrDof),out)       applied f vector         keyMtx(5)=1
c     fInt      (dp,ar(nUsrDof),out)       internal force vector    keyMtx(6)=1

c     elVol     (dp,sc,out)        element volume
c     elMass    (dp,sc,out)        element mass
c     elCG      (dp,ar(3),out)     element centroid coordinates in current
c                                     configuration
c     nRsltBsc  (dp,sc,in)         number of basic elemental results saved in
c                                   result files 
c     RsltBsc   (dp,ar(nRsltBsc),out) basic elemental results 
c                                       (see EXPLANATION below)
c     nRsltVar  (int,sc,in)        number of elemental results saved in 
c                                     result file as non-summable miscellaneous
c                                     variables 
c                                       (defined on USRELEM command as NRSLTVAR)
c     RsltVar   (dp,ar(nRsltVar),out) variables to saved in result files as
c                                      non-summable miscellaneous variables 
c                                      requested when outkey = 1
c
c     nElEng    (int,sc,in)        number of energies (fixed at 3)
c
c     elEnergy  (dp,ar(nElEng),out) elemental energy
c                                     elEnergy(1), element strain energy
c                                     elEnergy(2), element plastic strain energy
c                                     elEnergy(3), element creep strain energy
c
c     EXPLANATION OF RsltBsc
c     
c       Basic element results are saved and total number of result 
c     quantities is nRsltBsc, where:
c            nRsltBsc = (7+7)* number of corner nodes at one element.
c       To process the quantites by post processing properly, the results 
c     must be in the following order:
c       1.) Stresses: Sx Sy Sz Sxy Syz Sxz Seqv at all corner points,
c     followed by:
c       2.) Strains : Ex Ey Ez Exy Eyz Exz Eeqv at all corner points
c     where Seqv and Eeqv = equivalent stress and strain respectively
c
c
************************************************************************
c

#include "impcom.inc"
c
      EXTERNAL         ElemGetMat  ! this routine may be user-programmed

      INTEGER          elId, matId, keyMtx(10), lumpm,nDim, nNodes,
     &                 Nodes(nNodes), nIntPnts, nUsrDof, kEStress, 
     &                 keyAnsMat, keySym, nKeyOpt, KeyOpt(nKeyOpt),
     *                 kTherm, nPress, kPress, nReal, nSaveVars, 
     &                 kfstps, nlgeom, nrkey, outkey, jdim, 
     &                 elPrint, iott, keyHisUpd, l, inod,
     &                 ldstep, isubst, ieqitr, keyEleErr, keyEleCnv,
     &                 nRsltBsc, nRsltVar, nElEng, gggggg, intPnttt


      DOUBLE PRECISION temper(nNodes), temperB(nNodes), tRef, 
     &                 Press(nPress), RealConst(nReal),
     &                 saveVars(nSaveVars), klk(30),
     &                 xRef(nDim,nNodes), xCur(nDim,nNodes),
     &                 TotValDofs(nUsrDof), IncValDofs(nUsrDof), 
     &                 ItrValDofs(nUsrDof), VelValDofs(nUsrDof),
     &                 AccValDofs(nUsrDof),      timval,
     &                 eStiff(nUsrDof,nUsrDof), eMass(nUsrDof,nUsrDof), 
     &                 eDamp(nUsrDof,nUsrDof), eSStiff(nUsrDof,nUsrDof), 
     &                 fExt(nUsrDof), fInt(nUsrDof), 
     &                 elVol, elMass, elCG(3),
     &                 RsltBsc(nRsltBsc), RsltVar(nRsltVar), 
     &                 elEnergy(nElEng)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c *** CODE EXAMPLE ***
c
c --- The element code is only to show how to use the routine to create user 
c     elements.  Two element types are coded.   Only the stiffness matrix, mass 
c     matrix and internal load vector are shown.
c
c       When KeyOpt(1) = 0, it is a structural 2D plane strain element 
c                                          with 4 nodes and 4 integration points
c       When KeyOpt(1) = 1, it is a structural 3D solid elements 
c                                         with 20 nodes and 8 integration points
c       No advanced element technology is employed,
c                         and they are only coded for geometric linear analysis.
c
c --- This demonstration code only shows how to create eStiff, eMass, and fInt.
c --- Other matrices and/or vectors should be created similarly. 
c --- This is coded for good readability.
c
c --- Included decks, functions, variables defined for the example

#include "locknm.inc"
      EXTERNAL         vzero, vmove, vmult, vdot, vidot,
     &                 maxv, matxb, matba, maat, matsym, getMatProp,
     &                 erhandler, equivStrain, ElemJac, ElemMass,
     &                 ElemRsltNode, ElemShpFn, pplock, ppunlock

      DOUBLE PRECISION vdot, vidot

      INTEGER          nUsrDof2, intPnt, iNode, nTens, flgSingular,
     &                 k1, k2, k3, nComp, iDim, iDim1, iComp,numU,
     &                 nNodesCorner, nDirect, kThermIP, i, j, kFlagF,
     &                 kFlagS, typemodel, linear, drivingforce,
     &                 func_name, lll, llll   
      DOUBLE PRECISION BMat(nDim*2,nUsrDof), Ex, nu, density, G, workDb,
     &                 con1, con2, cMat(nDim*2,nDim*2), shIsoC(nNodes),
     &                 shIso(nNodes), shDerIso(nDim,nNodes), wtIP(1),
     &                 workArr(360), elJac(nDim*nDim), detJac, dperr(2),
     &                 shDerEl(nDim,nNodes), dVol, Strain(nDim*2), 
     &                 Stress(nDim*2), wStrain(48), wStress(48),
     &                 nStrain(48), nStress(48), sigm, tem, prop(3),
     &                 IncStrain(nDim*2),  defG(3,3), 
     &                 defG0(3,3), xCurIP(nDim), TemperIP, 
     &                 TemperIPB, StressTh(nDim*2), MatProp(5),
     &                 StrainPl(nDim*2), StrainCr(nDim*2), 
     &                 StrainTh(nDim*2), StrainSw, StressBk(nDim*2),
     &                 MatRotGlb(3,3), wStrainTh(48), wStrainPl(48),
     &                 wStrainCr(48), eMassb(nNodes,nNodes), EnergyD(3),
     &                 phik(nNodes), phi, dNdx(nDim,nNodes) , phin, H,               !
     &                 du(nDim*nNodes), Gc, xlc, xk, dstran(nDim*2,1),             !
     &                 Strain1(nDim*2), Stress1(nDim*2), Hn, Psi,                    !     
     &                 Stress2(nDim*2,1), dN(nNodes,1), w, dw, ddw, cw,                             !
     &                 b(nDim*2,nUsrDof-nNodes), rhs(nUsrdof), trE,                      !
     &                 amatrx(nUsrDof,nUsrDof),eStiff1(nUsrDof,nUsrDof),             !   
     &                 u(nUsrDof-nNodes), pl, pln, plast, Psi1, sedd,
     &                 alph, alphn, alphBn, Fdeg, Ac, alphT, seddn,                ! 
     &                 alphB, coordx, pi, bulk, Hmin, e, d, gNum, dgNum, 
     &                 ddgNum, gDen, dgDen, ddgDen, a, StrainEl(nDim*2), 
     &                 gf, dgf, sigc, ddgf, mm, EdevS, Edev(3), phinn,
     &                 psip, psin, eg, trEp2, trEn2, trEp, trEn, ggf
      
      
      CHARACTER*4      label(3)

c --- temporary debug key
      INTEGER debug, ix
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c --- B E G I N   E X E C U T A B L E  C O D I N G
c
c --- initialization
      Ex = RealConst(1) 
      nu = RealConst(2)
      xlc= RealConst(3)
      Gc=  RealConst(4)
      xk=  RealConst(5)
      kFlagF = 0
      kFlagS = 0
      typemodel = 3
      linear = 1
      drivingforce = 4
      Ac = 0
      nTens = nDim*2
      nComp = nDim*nDim
      nDirect = 3     
      amatrx = 0.d0
      rhs = 0.d0
      numU = nUsrDof/nNodes
      sigc = 1000
      pi = 3.14159
      bulk = Ex/3/(1-2*nu)
      lll = 3
      llll = 4
      
      IF (numU.EQ.3.d0) THEN
          KeyOpt(1) = 0.0d0
      ELSE
          KeyOpt(1) = 1.0d0
      END IF    
      
      if (kFlagF.eq.1) then
        alphT=Gc/(2.d0*6.d0*xlc)
      else
        alphT=1.d10
      endif   
      
      keyMtx = 0.d0
      keyMtx(1) = 1.d0
      keyMtx(6) = 1.d0
      nUsrDof2 = nUsrDof*nUsrDof
      CALL vzero (BMat(1,1),nUsrDof*nTens)
      IF (keyMtx(1).EQ.1) CALL vzero (eStiff(1,1),nUsrDof2)
      IF (keyMtx(2).EQ.1) CALL vzero (eMass(1,1) ,nUsrDof2)
      IF (keyMtx(5).EQ.1) CALL vzero (fExt(1)    ,nUsrDof)
      IF (keyMtx(6).EQ.1) CALL vzero (fInt(1)    ,nUsrDof)
      nlgeom = 0
      IF (nlgeom.EQ.0) THEN
         DO iDim = 1, 3
            DO iDim1 = 1, 3
               defG0(iDim, iDim1) = 0.0D0
            END DO
            defG0(iDim, iDim) = 1.0D0
         END DO
         CALL vmove (defG0(1,1),defG(1,1),9)
      ELSE
c        Nonlinear logic not defined here
      END IF

c --- start integration loop

      DO intPnt = 1, nIntPnts

c --- obtain shape functions and derivatives of shape functions

         IF (KeyOpt(1).EQ.0) THEN
            CALL ElemShpFn (1, intPnt, 1, shIso(1), nNodes)
            CALL ElemShpFn (1, intPnt, 2, shDerIso(1,1), nUsrDof)
            CALL ElemShpFn (1, intPnt, 3, wtIP(1), 1)
         ELSE
            CALL ElemShpFn (2, intPnt, 1, shIso(1), nNodes)
            CALL ElemShpFn (2, intPnt, 2, shDerIso(1,1), nUsrDof)
            CALL ElemShpFn (2, intPnt, 3, wtIP(1), 1)
         END IF
         
         dN = 0.d0
         do i=1,nNodes
           dN(i,1) = shIso(i)
         end do  

c --- coordinates at integration points

         DO iDim = 1, nDim
            xCurIP(iDim) =  vidot(shIso(1), 1, xCur(iDim,1), nDim,
     &                 nNodes)
         END DO

c --- temperatures at integration points

         TemperIP = vdot(shIso(1), temper(1), nNodes)
         TemperIPB = vdot(shIso(1), temperB(1), nNodes)
         

c --- derivatives of shape functions

         CALL vzero (workArr(1), nComp)
         iComp = 1
         DO iDim = 1, nDim
            DO iDim1 = 1, nDim
               DO iNode = 1, nNodes
                  workArr(iComp) = workArr(iComp) 
     &                + shDerIso(iDim1,iNode)*xCur(iDim,iNode)
               END DO
               iComp = iComp + 1
            END DO
         END DO
         CALL ElemJac (workArr(1), elJac(1), nDim, detJac, 
     &                   flgSingular)
         IF (flgSingular.LE.0) THEN
            dperr(1) = detJac
            dperr(2) = elId
            CALL erhandler ('UserElem', 1100, 3, 'Negative element
     &                       Jacobian value %I at element %I. This
     &                       is due to wrong element order or bad  
     &                       mesh.',dperr(1),' ') 
            GOTO 990
         END IF
         DO iNode = 1, nNodes
            IF (KeyOpt(1).EQ.0) THEN
               shDerEl(1,iNode) =  elJac(1)*shDerIso(1,iNode)
     &                           + elJac(3)*shDerIso(2,iNode)
               shDerEl(2,iNode) =  elJac(2)*shDerIso(1,iNode)
     &                           + elJac(4)*shDerIso(2,iNode)
            ELSE
               shDerEl(1,iNode) =  elJac(1)*shDerIso(1,iNode)
     &                           + elJac(4)*shDerIso(2,iNode)
     &                           + elJac(7)*shDerIso(3,iNode)
               shDerEl(2,iNode) =  elJac(2)*shDerIso(1,iNode)
     &                           + elJac(5)*shDerIso(2,iNode)
     &                           + elJac(8)*shDerIso(3,iNode)
               shDerEl(3,iNode) =  elJac(3)*shDerIso(1,iNode)
     &                           + elJac(6)*shDerIso(2,iNode)
     &                           + elJac(9)*shDerIso(3,iNode)
            END IF
         END DO
         
         dNdx = 0.d0
         do i=1,nDim
             do j=1,nNodes
                 dNdx(i,j) = shDerEl(i,j) 
             end do
         end do    
         dVol = detJac*wtIP(1)

c --- create B matrix
         BMat = 0.d0 
         k1 = 1
         DO iNode = 1,nNodes
            k2  = k1 + 1
            BMat(1,k1) =  shDerEl(1,iNode)
            BMat(2,k2) =  shDerEl(2,iNode)
            BMat(4,k1) =  shDerEl(2,iNode)
            BMat(4,k2) =  shDerEl(1,iNode)
            IF (KeyOpt(1).EQ.1) THEN
               k3 = k2 + 1
               BMat(3,k3) =  shDerEl(3,iNode)
               BMat(5,k2) =  shDerEl(3,iNode)
               BMat(5,k3) =  shDerEl(2,iNode)
               BMat(6,k3) =  shDerEl(1,iNode)
               BMat(6,k1) =  shDerEl(3,iNode)
            END IF         
            k1 = k1 + nDim
         END DO
         
          b = 0.d0
          k1 = 1
          DO iNode = 1,nNodes
            k2  = k1 + 1
            b(1,k1) =  shDerEl(1,iNode)
            b(2,k2) =  shDerEl(2,iNode)
            b(4,k1) =  shDerEl(2,iNode)
            b(4,k2) =  shDerEl(1,iNode)
            IF (KeyOpt(1).EQ.1) THEN
               k3 = k2 + 1
               b(3,k3) =  shDerEl(3,iNode)
               b(5,k2) =  shDerEl(3,iNode)
               b(5,k3) =  shDerEl(2,iNode)
               b(6,k3) =  shDerEl(1,iNode)
               b(6,k1) =  shDerEl(3,iNode)
            END IF         
            k1 = k1 + nDim
          END DO
          
!         calculate phi
          
          DO iNode = 1,nNodes
              phik(iNode) = TotValDofs(numU*iNode)
          end do    
          phi = 0.d0
          do inod=1,nNodes
              phi=phi+shIso(inod)*phik(inod)
          end do
          if (phi.gt.1.d0) then
              phi=0.999d0
          end if
          if (TotValDofs(2).ne.0.d0) then
              cMat(4,4) = G
          end if
          
!         calculate dstran
          
          du = 0.d0
          do i=1,nUsrDof
               IF (MOD(i,numU).EQ.0.d0) THEN
                  cMat(4,4) = G
                  cMat(4,4) = G
              ELSE 
                  du(i-INT(i/numU)) = ItrValDofs(i)
                  u(i-INT(i/numU)) = TotValDofs(i)
              END IF    
          end do 
c          dstran=0.d0
c          dstran= matmul(b,du)
          
          do i=1,ntens
              Stress1(i) = saveVars(20*(intPnt-1) + i)
              Strain1(i) = saveVars(20*(intPnt-1) + i + 6)
          end do    
          phin = saveVars(20*(intPnt-1) + 13)
          Hn = saveVars(20*(intPnt-1) + 14)
          alphn = saveVars(20*(intPnt-1) + 15)
          alphBn = saveVars(20*(intPnt-1) + 16) 
          sedd = saveVars(20*(intPnt-1) + 18)
          seddn = sedd  
          phinn = phin
          
          coordx=0.d0
          if (phi.ge.0.95d0) then
            do i=1,nNodes
                coordx=coordx+dN(i,1)*xCur(1,i)
            enddo
            if (coordx.gt.Ac) then
                Ac=coordx
            endif
          endif 

c --- calculate strains and stress

         CALL maxv (b(1,1), u(1), IncStrain(1), nTens, 
     &              nUsrDof-nNodes)
         
c         CALL straininc(ntens,ndim,nNodes,nUsrDof,b,u,defG0)
         
         keyAnsMat= 1.0d0
         IF (keyAnsMat.EQ.1) THEN
c           ---- Use standard ANSYS material (METHOD 1)
c          ---- USERMAT is called from here
            CALL ElemGetMat (elId, matId, nDim, nTens, nDirect,
     &                         intPnt, xCurIP(1), TemperIP,
     &                         TemperIPB, kThermIP, IncStrain(1),
     &                         defG0(1,1), defG(1,1),
     &                         cMat(1,1), MatProp(1), Stress(1), 
     &                         Strain(1), StressTh(1), StrainTh(1),
     &                         StrainPl(1), StrainCr(1), 
     &                         StressBk(1), StrainSw, EnergyD(1),
     &                         MatRotGlb(1,1))
            if (kThermIP .eq. 1) kTherm = 1
c            call get_ElmData ('SVAR', elId,intPnt,30, klk)
                        
         ELSE
c          ---- Make up your own material (METHOD 2)
            IF (nlgeom.EQ.0) CALL vmove (IncStrain(1), Strain(1),
     &                                   nTens)
            Stress(1) = con1*Strain(1) + con2*(Strain(2)+Strain(3))
            Stress(2) = con1*Strain(2) + con2*(Strain(3)+Strain(1))
            Stress(3) = con1*Strain(3) + con2*(Strain(1)+Strain(2))
            Stress(4) = G*Strain(4)
            IF (KeyOpt(1).EQ.1) THEN
               Stress(5) = G*Strain(5)
               Stress(6) = G*Strain(6)
            END IF
         END IF 
         
         if(kflagS.eq.0) then
             phin=phi
         end if 

******Выбор параметров модели****************************************
         if (typemodel .EQ. 1) then
             gf = (1.d0-phin)**2 + xk
             ggf = (1.d0-phinn)**2 + xk
             dgf = 2*(phin -1)
             ddgf = 2         
             w = phin**2
             dw = 2*phin
             ddw = 2
             cw = 0.5
c             Hmin = 0
             Hmin = 3.d0*Gc/(16.d0*xlc)
         else if (typemodel .EQ. 2) then
             gf = (1.d0-phin)**2 + xk
             dgf = 2*(phin -1)
             ddgf = 2         
             w = phin
             dw = 1
             ddw = 0
             cw = 2/3
             Hmin = 0 
         else
            mm=4.d0*Ex*Gc/(xlc*sigc**2) 
            w=2.d0*phi-phi**2
            dw=2.d0-2.d0*phi
            ddw=-2.d0
            cw = 4.d0/3.d0
            if (linear .EQ. 1) then
                e=-0.5d0
                d=2.d0
            else
                e=2.d0**(5.d0/3.d0)-3.d0
                d=2.5d0
            end if
            gNum=(1.d0-phi)**d
            dgNum=-d*(1.d0-phi)**(d-1.d0);
            ddgNum=d*(d-1.d0)*(1.d0-phi)**(d-2.d0)
            gDen=gNum+mm*phi+mm*e*phi**2.d0
            dgDen=dgNum+mm+2.d0*mm*e*phi
            ddgDen=ddgNum+2.d0*mm*e

            gf=gNum/gDen + xk
            dgf=(dgNum*gDen-gNum*dgDen)/(gDen**2.d0)
            ddgf=((ddgNum*gDen-gNum*ddgDen)*gDen-2.d0*
     1 (dgNum*gDen-gNum*dgDen)*dgDen)/(gDen**3.d0)
            
            Hmin = 0.5d0*sigc**2/Ex
         end if
         
c********Расчёт движущей силы разрушения*********************************
         psip = 0
         psin = 0
         StrainEl=Strain-StrainPl       
         trE=StrainEl(1)+StrainEl(2)+StrainEl(3)
         trEp=0.5d0*(trE+abs(trE))
         trEn=0.5d0*(trE-abs(trE))
         Edev(1:3)=StrainEl(1:3)-trE/3.d0
         EdevS=Edev(1)**2+Edev(2)**2+Edev(3)**2 
         eg = Ex/(1.d0+nu)/2.d0
         if (drivingforce.eq.2) then ! Amor et al. 
             psip=0.5d0*bulk*trEp**2+eg*EdevS
             psin=0.5d0*bulk*trEn**2
         elseif (drivingforce.eq.3) then ! Miehe et al.
             trEp2=0.d0
             trEn2=0.d0
             do i=1,3
                 trEp2=trEp2+(StrainEl(i)+abs(StrainEl(i)))**2.d0/4.d0
                 trEn2=trEn2+(StrainEl(i)-abs(StrainEl(i)))**2.d0/4.d0
             end do
             psip=nu*eg/(1d0-2d0*nu)*trEp**2d0+eg*trEp2
             psin=nu*eg/(1d0-2d0*nu)*trEn**2d0+eg*trEn2
         elseif  (drivingforce.eq.1) then! no split
             do i=1,ntens
c                 psip=psip+Stress(i)*StrainEl(i)*0.5d0/gf
                  psip = EnergyD(1)
             end do
         else 
            psip=0.5d0*max(Stress(1),Stress(2),Stress(3))**2/Ex
         endif            
    
         H=max(Hmin,psip,Hn)
         
         if(kflagS.eq.0) then
             H=H
         else
             H=Hn   
         end if
c         if (H.lt.Hn) H=Hn
c*************************************************************************
          
c********Функция деградации усталость*************************************  
         alph=H*((1-phin)**2+xk)

         if (alph.ge.alphn) then
             alphB = alphBn+abs(alph-alphn)
         else
             alphB=alphBn
         endif            
          
          if (alphB.lt.alphT) then
              Fdeg= 1.d0
          else
              Fdeg=(2.d0*alphT/(alphB+alphT))**2.d0
          endif      
c*************************************************************************
          
c*********Запись данных***************************************************  
          DO i=1,ntens
              saveVars(20*(intPnt-1) + i) = Stress(i)
              saveVars(20*(intPnt-1) + i + 6) = Strain(i)
          END DO    
          saveVars(20*(intPnt-1) + 13) = phi
          saveVars(20*(intPnt-1) + 14) = H    
          saveVars(20*(intPnt-1) + 15) = alph 
          saveVars(20*(intPnt-1) + 16) = alphB
          saveVars(20*(intPnt-1) + 17) = Ac
          saveVars(20*(intPnt-1) + 18) = psip
          saveVars(20*(intPnt-1) + 19) = H + EnergyD(2) 
c*************************************************************************
          
c*********Расчёт матрицы узловых нагрузок и матрицы жёсткости элемента          
           IF (keyMtx(1).EQ.1) THEN
              amatrx(1:ndim*nNodes,1:ndim*nNodes)=
     1      amatrx(1:ndim*nNodes,1:ndim*nNodes)+
     1      dvol*(ggf*
     1      matmul(matmul(transpose(b),cMat),b)) 
            
          
              amatrx((ndim*nNodes+1):nUsrDof,(ndim*nNodes+1):nUsrDof)=
     1  amatrx((ndim*nNodes+1):nUsrDof,(ndim*nNodes+1):nUsrDof)
     1 +dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc*Fdeg/2.d0/cw
     2 +matmul(dN,transpose(dN))*(Gc/xlc*Fdeg/4.d0/cw*ddw
     2  +ddgf*(H+EnergyD(2))))
              
          END IF  
          
          IF (keyMtx(6).EQ.1) THEN 
              rhs(1:ndim*nNodes)=rhs(1:ndim*nNodes)+
     1 dvol*(ggf*matmul(transpose(b),Stress)) 
              
              rhs((ndim*nNodes+1):nUsrDof)=
     1 +rhs((ndim*nNodes+1):nUsrDof)
     1 +dvol*(matmul(transpose(dNdx),matmul(dNdx,phik(1:nNodes)))*Fdeg
     2 *Gc*xlc/2.d0/cw+dN(1:nNodes,1)*(Gc/xlc*Fdeg/4.d0/cw*dw
     2  +dgf*(H+EnergyD(2))))
            
          END IF  
            

          
c*************************************************************************

            

c --- create stiffness matrix

c         IF (keyMtx(1).EQ.1) CALL matba (BMat(1,1), cMat(1,1), 
c     &                            eStiff(1,1), nTens, nTens, nUsrDof, 
c     &                            nTens, nUsrDof, workArr(1), dVol)

c --- prepare to create mass matrix

c         IF (keyMtx(2).EQ.1) THEN
c            IF (density.NE.0.0d0) THEN
c               workDb = density*dVol
c               CALL maat (shIso(1),eMassb(1,1),nNodes,nNodes,workDb)
c            ENDIF
c         ENDIF

c --- create external force vector

c         IF (keyMtx(5).EQ.1 .AND. kThermIP.EQ.1 .AND. outkey.EQ.0) THEN
c            CALL vmult (StressTh(1), workArr(1), nTens, dVol)
c            CALL matxb (BMat(1,1), workArr(1), fExt(1), nTens, nTens,
c     &                  nUsrDof, nUsrDof, 1, -nTens)
c         END IF

c --- create internal force vector

c         IF (keyMtx(6).EQ.1) THEN
c            CALL vmult (Stress(1), workArr(1), nTens, dVol)
c            CALL matxb (BMat(1,1), workArr(1), fInt(1), nTens, nTens,
c     &                  nUsrDof, nUsrDof, 1, -nTens)
c         END IF

c --- calculate other element quantities

!         elVol = elVol+dVol
!         elMass = elMass+dVol*density
!         IF (keyAnsMat.EQ.0) elEnergy(1) = elEnergy(1)
!     &                 + 0.5d0*dVol*vdot(Strain(1), Stress(1), nTens)
!         k1 = (intPnt-1)*nTens+1
!c         CALL vmove (Stress(1), saveVars(k1), nTens)
!         IF (outkey.EQ.1) THEN
!            CALL vmove (Strain(1), wStrain(k1), nTens)
!            CALL vmove (Stress(1), wStress(k1), nTens)
!            IF (keyAnsMat.EQ.1) THEN
!               CALL vmove (StrainTh(1), wStrainTh(k1), nTens)
!               CALL vmove (StrainPl(1), wStrainPl(k1), nTens)
!               CALL vmove (StrainCr(1), wStrainCr(k1), nTens)
!            END IF
!            IF (debug.EQ.1) THEN
!               write (*,3010) intPnt, (Strain(ix),ix=1,nTens)
!               write (*,3020) (Stress(ix),ix=1,nTens)
!               write (*,3030) (StrainPl(ix),ix=1,nTens)
! 3010          FORMAT (/1x, 'intPnt=',i2, 'Strain=',6(e15.8,2x))
! 3020          FORMAT (1x, 8x, 'Stress=',6(e15.8,2x))
! 3030          FORMAT (1x, 8x, 'StrainPl=',6(e15.8,2x))
!            END IF
!         END IF
      END DO
!      
          do i=1,nUsrDof
               IF (MOD(i,numU).EQ.0.d0) THEN
                  fint(i) = rhs(ndim*nNodes + INT(i/numU))
              ELSE 
                  fint(i) = rhs(i-INT(i/numU))
              END IF    
          end do
          
          do i=1,nUsrDof
               IF (MOD(i,numU).EQ.0.d0) THEN
                   do j=1,nUsrDof
                       IF (MOD(j,numU).EQ.0.d0) THEN
                          eStiff(i,j) = 
     &    amatrx(ndim*nNodes + INT(i/numU),ndim*nNodes + INT(j/numU))
                       ELSE
                          eStiff(i,j) = 0.d0
                       END IF    
                   end do    
              ELSE 
                   do j=1,nUsrDof
                       IF (MOD(j,numU).EQ.0.d0) THEN
                          eStiff(i,j) = 0.d0 
                       ELSE
                          eStiff(i,j) = 
     &                        amatrx(i-INT(i/numU),j-INT(j/numU))     
                       END IF    
                   end do   
              END IF    
          end do
      

          
          
          
!      
!
!c --- symmetricize eStiff
!
!c      IF (keyMtx(1).EQ.1) CALL matsym (eStiff(1,1), nUsrDof, nUsrDof)
!
!c --- create mass matrix
!
!c      IF (keyMtx(2).EQ.1 .AND. density.NE.0.0d0) THEN
!c         CALL ElemMass (eMassb(1,1), nNodes, nDim, nUsrDof, eMass(1,1))
!c      ENDIF
!
!c --- calculate strains and stresses at nodes by extrapolating and output
!c       to result files
!
!      IF (outkey.EQ.1) THEN
!         IF (KeyOpt(1).EQ.0) THEN
!            nNodesCorner = nNodes
!         ELSE
!            nNodesCorner = 8
!         END IF
!         CALL ElemRsltNode (KeyOpt(1), nrkey, nTens, wStress(1), 
!     &                        wStrain(1), nIntPnts, nStress(1), 
!     &                        nStrain(1), nNodesCorner)
!
!c --- only calculate basic result variables when it is allowed
!c        and necessary
!
!         IF (nRsltBsc.GT.0) THEN
!            DO iNode = 1, nNodesCorner
!               k1 = (iNode-1)*nTens + 1
!               k2 = (iNode-1)*7 + 1
!               CALL vmove (nStress(k1), RsltBsc(k2), nTens)
!               sigm = (nStress(k1)+nStress(k1+1)+nStress(k1+2))/3.0d0
!               IF (KeyOpt(1).EQ.0) THEN
!                  RsltBsc(k2+4) = 0.0d0
!                  RsltBsc(k2+5) = 0.0d0
!                  RsltBsc(k2+6) = SQRT(1.5d0*
!     &                    ( (nStress(k1)-sigm)*(nStress(k1)-sigm)
!     &                    + (nStress(k1+1)-sigm)*(nStress(k1+1)-sigm)
!     &                    + (nStress(k1+2)-sigm)*(nStress(k1+2)-sigm)
!     &                    + 2.0d0*nStress(k1+3)*nStress(k1+3)))
!               ELSE
!                  RsltBsc(k2+6) = SQRT(1.5d0*
!     &                    ( (nStress(k1)-sigm)*(nStress(k1)-sigm)
!     &                    + (nStress(k1+1)-sigm)*(nStress(k1+1)-sigm)
!     &                    + (nStress(k1+2)-sigm)*(nStress(k1+2)-sigm)
!     &                    + 2.0d0*(nStress(k1+3)*nStress(k1+3)
!     &                    +        nStress(k1+4)*nStress(k1+4)
!     &                    +        nStress(k1+5)*nStress(k1+5))))
!               END IF
!
!               k2 = (nNodesCorner+iNode-1)*7 + 1
!               CALL vmove (nStrain(k1), RsltBsc(k2), nTens)
!               IF (KeyOpt(1).EQ.0) THEN
!                  RsltBsc(k2+4) = 0.0d0
!                  RsltBsc(k2+5) = 0.0d0
!               END IF
!               CALL equivStrain (nu, nStrain(k1), nTens, 
!     &                            RsltBsc(k2+6))
!            END DO    
!         END IF
!         k1 = nNodesCorner*nTens
!c         CALL vmove (nStrain(1), RsltVar(1), k1)
!c         CALL vmove (nStress(1), RsltVar(k1+1), k1)
!
!         IF (elPrint .EQ. 1) THEN
!c --- print out the results in OUT file (requested by the OUTPR command)
!            CALL pplock (LOCKOT)
!            WRITE (iott,2000) elId
! 2000       FORMAT (/1x, 'Material Point output for element',I8)
!            WRITE (iott,2100)
! 2100       FORMAT(/4x, 'Intg.Pt. "S"     Stresses')
!            DO intPnt = 1, nIntPnts
!               k1 = (intPnt-1)*nTens
!               WRITE (iott, 2110) intPnt, (wStress(k1+k2),k2=1,nTens)
! 2110          FORMAT (4x,I4, 4x, 6(E12.5,1x))
!            END DO
!            WRITE (iott,2200)
! 2200       FORMAT(/4x, 'Intg.Pt. "EPTO"     Strains')
!            DO intPnt = 1, nIntPnts
!               k1 = (intPnt-1)*nTens
!               WRITE (iott, 2110) intPnt, (wStrain(k1+k2),k2=1,nTens)
!            END DO
!            IF (keyAnsMat.EQ.1) THEN
!               workDb = 0.0d0
!               DO k1 = 1, nIntPnts*nTens
!                  workDb = workDb + ABS(wStrainPl(k1))
!               END DO
!               IF (workDb.GT.0.0d0) THEN
!                  WRITE (iott,2300)
! 2300             FORMAT(/4x, 'Intg.Pt. "EPPL"     Strains')
!                  DO intPnt = 1, nIntPnts
!                     k1 = (intPnt-1)*nTens
!                     WRITE (iott, 2110) intPnt, (wStrainPl(k1+k2),
!     &                                  k2=1,nTens)
!                  END DO
!               END IF
!               workDb = 0.0d0
!               DO k1 = 1, nIntPnts*nTens
!                  workDb = workDb + ABS(wStrainCr(k1))
!               END DO
!               IF (workDb.GT.0.0d0) THEN
!                  WRITE (iott,2400)
! 2400             FORMAT(/4x, 'Intg.Pt. "EPCR"     Strains')
!                  DO intPnt = 1, nIntPnts
!                     k1 = (intPnt-1)*nTens
!                     WRITE (iott, 2110) intPnt, (wStrainCr(k1+k2),
!     &                                  k2=1,nTens)
!                  END DO
!               END IF
!               workDb = 0.0d0
!               DO k1 = 1, nIntPnts*nTens
!                  workDb = workDb + ABS(wStrainTh(k1))
!               END DO
!               IF (workDb.GT.1.0d-12) THEN
!                  WRITE (iott,2500)
! 2500             FORMAT(/4x, 'Intg.Pt. "EPTH"     Strains')
!                  DO intPnt = 1, nIntPnts
!                     k1 = (intPnt-1)*nTens
!                     WRITE (iott, 2110) intPnt, (wStrainTh(k1+k2),
!     &                                  k2=1,nTens)
!                  END DO
!               END IF
!               write (iott,3000)
! 3000          format(2/)
!            END IF
!            CALL ppunlock (LOCKOT)
!         END IF
!      END IF
     
 990  CONTINUE
      RsltVar = saveVars 
      RETURN
      END

       subroutine straininc(ntens,ndim,nNodes,nUsrDof,bmat,utmp,xx1)
c
c     Notation:  
c       dstran(i)  incremental strain component 
c       note:      i = 1   xx direction 
c                    = 2   yy direction 
c                    = 3   zz direction
c                    = 4   xy direction
c    u() - displacement
c   du() - increment of displacement in the last inc.
c   
#include "impcom.inc"

      DOUBLE PRECISION bmat(ndim,nnodes),
     & xx1(3,3),utmp(nUsrDof-nNodes)
      
      INTEGER k1,i, ntens,ndim,nNodes,nUsrDof, nodi
      
      DOUBLE PRECISION dNidx, dNidy

      ! set xx1 to Identity matrix
      xx1=0.d0
      do k1=1,3
       xx1(k1,k1)=1.d0       
      end do

c************************************
c    Compute incremental strains
c************************************
c
      do nodi=1,nNodes

       dNidx=bmat(1,nodi)
       dNidy=bmat(2,nodi)


c     deformation gradient
       xx1(1,1)=xx1(1,1)+dNidx*utmp(1)
       xx1(1,2)=xx1(1,2)+dNidy*utmp(1)
       xx1(2,1)=xx1(2,1)+dNidx*utmp(2)
       xx1(2,2)=xx1(2,2)+dNidy*utmp(2)
c

      end do
c
      return
      end