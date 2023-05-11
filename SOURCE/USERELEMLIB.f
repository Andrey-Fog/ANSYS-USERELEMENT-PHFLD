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
      
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"USERELEM"::UserElem

#include "impcom.inc"
c
      EXTERNAL         ElemGetMat  ! this routine may be user-programmed

      INTEGER          elId, matId, keyMtx(10), lumpm,nDim, nNodes,
     &                 Nodes(nNodes), nIntPnts, nUsrDof, kEStress, 
     &                 keyAnsMat, keySym, nKeyOpt, KeyOpt(nKeyOpt),
     *                 kTherm, nPress, kPress, nReal, nSaveVars, 
     &                 kfstps, nlgeom, nrkey, outkey, jdim, 
     &                 elPrint, iott, keyHisUpd, l, inod,nTens,
     &                 ldstep, isubst, ieqitr, keyEleErr, keyEleCnv,
     &                 nRsltBsc, nRsltVar, nElEng


      DOUBLE PRECISION temper(nNodes), temperB(nNodes), tRef, 
     &                 Press(nPress), RealConst(nReal),
     &                 saveVars(nSaveVars), klk(10),
     &                 xRef(nDim,nNodes), xCur(nDim,nNodes),
     &                 TotValDofs(nUsrDof), IncValDofs(nUsrDof), 
     &                 ItrValDofs(nUsrDof), VelValDofs(nUsrDof),
     &                 AccValDofs(nUsrDof), timval,
     &                 eStiff(nUsrDof,nUsrDof), eMass(nUsrDof,nUsrDof), 
     &                 eDamp(nUsrDof,nUsrDof), eSStiff(nUsrDof,nUsrDof), 
     &                 fExt(nUsrDof), fInt(nUsrDof), 
     &                 elVol, elMass, elCG(3),
     &                 RsltBsc(nRsltBsc), RsltVar(nRsltVar), 
     &                 elEnergy(nElEng)
c
#include "locknm.inc"


      EXTERNAL         vzero, vmove, vdot, vidot,
     &                 ElemShpFn

      DOUBLE PRECISION vdot, vidot

      INTEGER          nUsrDof2, intPnt, iNode,
     &                 k1, k2, k3, nComp, iDim, iDim1, iComp, i, j, k,
     &                 nDirect, kThermIP
      DOUBLE PRECISION cMat(ndim*2,ndim*2),shIsoC(nNodes),
     &                 shIso(nNodes),dVol, Strain(ndim*2), 
     &                 Stress(ndim*2), prop(5),Bmat(nDim*2,8),
     &                 IncStrain(ndim*2),defG(3,3),coord24(2,4),
     &                 defG0(3,3), xCurIP(ndim), TemperIP,
     &                 TemperIPB, StressTh(ndim*2), MatProp(5),
     &                 StrainPl(ndim*2), StrainCr(ndim*2),u(8,1), 
     &                 StrainTh(ndim*2), StrainSw, StressBk(ndim*2),
     &                 MatRotGlb(3,3), EnergyD(3), 
     &                 xjac(ndim,ndim), xjaci(ndim,ndim),
     &                 dN(nNodes,1),dNdz(ndim,nNodes),eg2,elam,
     &                 dNdx(ndim,nNodes), djac, gaussCoord, phi,
     &                 Stress1(ndim*2), Strain1(ndim*2), H, phin,
     &                 Hn, Psi, phik(4),xx1(3,3),xx1Old(3,3),
     &                 du(8,1), stress2(ndim*2,1),
     &                 rhs(nUsrDof),dstran(4,1),
     &                 amatrx(nUsrDof,nUsrDof),stran(4,1) 
      
c --- Real constants      
      DOUBLE PRECISION Ex, nu, xlc, Gc,xk
c --- Flags
      INTEGER ELTYPE
c --- temporary debug key
      INTEGER debug, ix, debugELM 
      
      DOUBLE PRECISION IPCOORD(nDim,nIntPnts),s,t,gg,hh 
      
      
      data  coord24 /-1.d0, -1.d0,
     2                1.d0, -1.d0,
     3                1.d0,  1.d0,
     4                -1.d0, 1.d0/ 
      
      parameter (gaussCoord=0.577350269d0)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Ex = RealConst(1) 
      nu = RealConst(2)
      xlc= RealConst(3)
      Gc=  RealConst(4)
      xk=  RealConst(5)
      ELTYPE=  RealConst(6)
c      ELID near the crack tip 578, and 611 is far away
      debugELM = 578
      IF (elID.eq.debugELM) THEN
          debugELM=0.d0
      END IF
            nTens = nDim*2
      nComp = nDim*nDim
      nDirect = 3
      nlgeom = 0
      amatrx = 0.d0
      rhs = 0.d0
      eStiff = 0.d0
      fInt = 0.d0
      nUsrDof2 = nUsrDof*nUsrDof
      CALL vzero (Bmat(1,1),nUsrDof*nTens)
      IF (keyMtx(1).EQ.1) CALL vzero (eStiff(1,1),nUsrDof2)
      IF (keyMtx(2).EQ.1) CALL vzero (eMass(1,1) ,nUsrDof2)
      IF (keyMtx(3).EQ.1) CALL vzero (eDamp(1,1) ,nUsrDof2)
      IF (keyMtx(5).EQ.1) CALL vzero (fExt(1)    ,nUsrDof)
      IF (keyMtx(6).EQ.1) CALL vzero (fInt(1)    ,nUsrDof)
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
      elVol  = 0.d0
      elMass = 0.d0
      Strain = 0.d0
      Stress = 0.d0
      CALL vzero (elEnergy(1), nElEng)
      
c --- \\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////   
c ---  \\\\\\\\Start loop on material integration points/////       
c ---   \\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////       
      DO intPnt = 1, nIntPnts
c --- temperatures at integration points

         TemperIP = vdot(shIso(1), temper(1), nNodes)
         TemperIPB = vdot(shIso(1), temperB(1), nNodes)
         
          gg=coord24(1,intPnt)*gaussCoord
          hh=coord24(2,intPnt)*gaussCoord
          
          IF (ELTYPE.EQ.0) THEN
          !shape functions  KOSOV    
              gg=coord24(1,intPnt)*gaussCoord
              hh=coord24(2,intPnt)*gaussCoord
        !     
              dN(1,1)=(1.d0-gg)*(1.d0-hh)/4.d0
              dN(2,1)=(1.d0+gg)*(1.d0-hh)/4.d0
              dN(3,1)=(1.d0+gg)*(1.d0+hh)/4.d0
              dN(4,1)=(1.d0-gg)*(1.d0+hh)/4.d0
          
c --- derivatives of shape functions
        !     derivative d(Ni)/d(gg)
              dNdz(1,1)=-(1.d0-hh)/4.d0
              dNdz(1,2)=(1.d0-hh)/4.d0
              dNdz(1,3)=(1.d0+hh)/4.d0
              dNdz(1,4)=-(1.d0+hh)/4.d0

        !     derivative d(Ni)/d(h)
              dNdz(2,1)=-(1.d0-gg)/4.d0
              dNdz(2,2)=-(1.d0+gg)/4.d0
              dNdz(2,3)=(1.d0+gg)/4.d0
              dNdz(2,4)=(1.d0-gg)/4.d0


              xjac=0.d0

          do iNode=1,nNodes
            do idim=1,ndim
              do jdim=1,ndim
                xjac(jdim,idim)=xjac(jdim,idim)+
     1        dNdz(jdim,iNode)*xCur(idim,iNode)      
              end do
            end do 
          end do
    !
              djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
              if (djac.gt.0.d0) then ! jacobian is positive - o.k.
               xjaci(1,1)=xjac(2,2)/djac
               xjaci(2,2)=xjac(1,1)/djac
               xjaci(1,2)=-xjac(1,2)/djac
               xjaci(2,1)=-xjac(2,1)/djac
              endif
  
              dNdx=matmul(xjaci,dNdz)
          ELSEIF (ELTYPE.EQ.1) THEN
c --- create B matrix TUMANOV
         if (intPnt==1) then
          s=-0.577350269189626
          t=-0.577350269189626
         elseif (intPnt==2) then 
          s=0.577350269189626
          t=-0.577350269189626
         elseif (intPnt==3) then
          s=0.577350269189626
          t=0.577350269189626  
         elseif (intPnt==4) then      
          s=-0.577350269189626
          t=0.577350269189626   
         end if
         
c     adopted linear shape functions 

         DO I =1,4
          dN(I,1)=(1.d0+coord24(1,intPnt)*gaussCoord*s)*
     &             (1.d0+coord24(2,intPnt)*gaussCoord*t)/4.d0
         END do 
            dNdz(1,1)=0.083333333333333401314*t-0.1443375672974065 
            dNdz(2,1)=0.083333333333333401314*s-0.1443375672974065
             
            dNdz(1,2)=-0.083333333333333401314*t+0.1443375672974065
            dNdz(2,2)=-0.083333333333333401314*s-0.1443375672974065
             
            dNdz(1,3)=0.083333333333333401314*t+0.1443375672974065
            dNdz(2,3)=0.083333333333333401314*s+0.1443375672974065
             
            dNdz(1,4)=-0.083333333333333401314*t-0.1443375672974065
            dNdz(2,4)=-0.083333333333333401314*s+0.1443375672974065
             
            xjac(1,1)= dNdz(1,1)*xCur(1,1)+dNdz(1,2)*xCur(1,2)
     1 +dNdz(1,3)*xCur(1,3)+dNdz(1,4)*xCur(1,4)
         
            xjac(1,2)= dNdz(1,1)*xCur(2,1)+dNdz(1,2)*xCur(2,2)
     1 +dNdz(1,3)*xCur(2,3)+dNdz(1,4)*xCur(2,4)
         
            xjac(2,1)= dNdz(2,1)*xCur(1,1)+dNdz(2,2)*xCur(1,2)
     1 +dNdz(2,3)*xCur(1,3)+dNdz(2,4)*xCur(1,4)
          
            xjac(2,2)= dNdz(2,1)*xCur(2,1)+dNdz(2,2)*xCur(2,2)
     1 +dNdz(2,3)*xCur(2,3)+dNdz(2,4)*xCur(2,4)
         
            djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1) 
          
            xjaci(1,1)=xjac(2,2)/djac 
            xjaci(1,2)=-xjac(1,2)/djac  
            xjaci(2,1)=-xjac(2,1)/djac   
            xjaci(2,2)=xjac(1,1)/djac
            !dN/dx
            dNdx(1,1)=xjaci(1,1)*dNdz(1,1)+xjaci(1,2)*dNdz(2,1) 
            dNdx(1,2)=xjaci(1,1)*dNdz(1,2)+xjaci(1,2)*dNdz(2,2) 
            dNdx(1,3)=xjaci(1,1)*dNdz(1,3)+xjaci(1,2)*dNdz(2,3) 
            dNdx(1,4)=xjaci(1,1)*dNdz(1,4)+xjaci(1,2)*dNdz(2,4) 
            !dN/dy
            dNdx(2,1)=xjaci(2,1)*dNdz(1,1)+xjaci(2,2)*dNdz(2,1) 
            dNdx(2,2)=xjaci(2,1)*dNdz(1,2)+xjaci(2,2)*dNdz(2,2)  
            dNdx(2,3)=xjaci(2,1)*dNdz(1,3)+xjaci(2,2)*dNdz(2,3) 
            dNdx(2,4)=xjaci(2,1)*dNdz(1,4)+xjaci(2,2)*dNdz(2,4) 
            
          END  IF
          dVol = djac*1.0d0
    !     form B-matrix
          Bmat=0.d0
          do iNode=1,nNodes
              Bmat(1,2*iNode-1)=dNdx(1,iNode)
              Bmat(2,2*iNode)=dNdx(2,iNode)
              Bmat(4,2*iNode-1)=dNdx(2,iNode)
              Bmat(4,2*iNode)=dNdx(1,iNode)
          end do
!     compute from nodal values	          
          phik(1) = TotValDofs(3)
          phik(2) = TotValDofs(6)
          phik(3) = TotValDofs(9)
          phik(4) = TotValDofs(12)
          phi = 0.d0
          do inod=1,nNodes
              phi=phi+dN(inod,1)*phik(inod)
          end do

          du = 0.d0
          u = 0.d0
          if (phi.gt.1.d0) phi=1
          do i=1,2
              du(i,1) = ItrValDofs(i)
              u(i,1) = IncValDofs(i)
          end do
          do i=4,5
              du(i-1,1) = ItrValDofs(i)
              u(i-1,1) = IncValDofs(i)
          end do  
          do i=6,7
              du(i-1,1) = ItrValDofs(i+1)
              u(i-1,1) = IncValDofs(i+1)
          end do  
          do i=8,9
              du(i-1,1) = ItrValDofs(i+2)
              u(i-1,1) = IncValDofs(i+2)
          end do
          
 !     compute the increment of strain and recover history variables          
          dstran = 0.d0
          dstran= matmul(Bmat,du)
          stran = matmul(Bmat,u)
          
          do i=1,ntens
              Stress1(i) = saveVars(10*(intPnt-1) + i)
              Strain1(i) = saveVars(10*(intPnt-1) + i + 4)
          end do    
          phin = saveVars(10*(intPnt-1) + 9)
          Hn = saveVars(10*(intPnt-1) + 10)
          Psi=0.d0
          do k1=1,ntens
              Psi=Psi+Stress1(k1)*Strain1(k1)*0.5d0
          end do
c --- Current IP coords          
       xCurIP=0.d0
       do k1=1,nnodes
        do k2=1,ndim
         xCurIP(k2)=xCurIP(k2)+dN(k1,1)*xRef(k2,k1)
        end do
       end do
c     calculate incremental strains from nodal values
       call straininc(ntens,ndim,nnodes,nUsrDof,dNdx,du,dstran,u,defG0,
     * xx1Old)
       IncStrain(1:4)=dstran(1:4,1)
c --- Get stresses and material stiffnes matrix   
       CALL ElemGetMat (elId, matId, nDim, nTens, nDirect,
     &                         intPnt, xCurIP(1), TemperIP,
     &                         TemperIPB, kThermIP, IncStrain(1),
     &                         defG0(1,1), defG(1,1),
     &                         cMat(1,1), MatProp(1), Stress(1), 
     &                         Strain(1), StressTh(1), StrainTh(1),
     &                         StrainPl(1), StrainCr(1), 
     &                         StressBk(1), StrainSw, EnergyD(1),
     &                         MatRotGlb(1,1))
c --- ************************************************************************ 
c ---   ************************ ElemGetMat description ********************
c --- ************************************************************************
c     input arguments
c     ===============
c     elId        (int,sc)           element number
c     matId       (int,sc)           material number of this element
c     nDim        (int,sc)           number of dimensions of the problem
c                                    = 2 2D
c                                    = 3 3D
c     nTens       (int,sc)           number of stress/strain components
c     nDirect     (int,sc)           number of stress/strain direct 
c                                      components
c     intPnt      (int,sc)           current integration point number
c     xCurIP      (dp,ar(3))         coordinates of integration point
c     TemperIP    (dp,sc)            integration point  temperatures at 
c                                      current time
c     TemperIPB   (dp,sc)            integration point  temperatures at 
c                                      the end of last incremetal step
c     IncStrain   (dp,ar(nTens))     strain for the current substep step when
c                                       nlgeom = on
c                                    total strain when nlgeom = off
c     defG0       (dp,ar(3x3))       deformation gradient tensor at the end
c                                       of last incremental step 
c
c     input output arguments         input desc     / output desc
c     ======================         ==========       ===========
c     defG        (dp, ar(3x3))      deformation gradient tensor at current
c                                      time, updated for thickness change in
c                                      plane stress when nlgeom=on
c     kTherm      (int,sc)           flag for thermal loading 
c                                      input as:
c                                      = 0 if temp = tref
c                                      = 1 if temp .ne. tref
c                                      gets reset to 0
c                                                   if ALPX, ALPY, and ALPZ = 0
c                                     
c     output arguments
c     ================
c     cMat        (nTens*nTens)      material Jacobian matrix
c     MatProp     (dp,ar(5))         commonly used materail properties
c                                    MatProp(1),Gxz: shear modulus
c                                    MatProp(2),Gyz: shear modulus
c                                    MatProp(3),Gxy: shear modulus
c                                    MatProp(4), density
c                                    MatProp(5), nuxy
c     Stress      (dp,ar(nTens))     total stress
c     Strain      (dp,ar(nTens))     total strain
c     StressTh    (dp,ar(nTens))     thermal stress
c     StrainTh    (dp,ar(nTens))     thermal strain
c     StrainPl    (dp,ar(nTens))     plastic strain
c     StrainCr    (dp,ar(nTens))     creep strain
c     StressBk    (dp,ar(nTens))     back stress for kinematic hardening
c     StrainSw    (dp,sc)            isotropic swelling strain 
c                                        (swelling capability not available yet)
c     EnergyD      (dp,ar(3))        energy density 
c                                    EnergyD(1) elastic energy density
c                                    EnergyD(2) plastic energy density  
c                                    EnergyD(3) creep energy density
c     MatRotGlb   (dp,ar(3,3))       The rotation matrix from global 
c                                     to  material coordinate system
c --- ************************************************************************
            kTherm = 0
!c
        do i=1,nTens
                Strain(i) = Strain1(i) + Strain(i)
                Stress(i) = Stress1(i) + Stress(i)
        end do
        if (Psi.gt.Hn) then
                  H=Psi
              else
                  H=Hn
        endif
        do i=1,ntens
            saveVars(10*(intPnt-1) + i) = Stress(i) 
            saveVars(10*(intPnt-1) + i + 4) = Strain(i)
        end do    
        saveVars(10*(intPnt-1) + 9) = phi
        saveVars(10*(intPnt-1) + 10) = H
        
        amatrx(1:8,1:8)=amatrx(1:8,1:8)+
     1      dvol*(((1.d0-phin)**2+xk)*
     1      matmul(matmul(transpose(Bmat),cMat),Bmat)) 
            
        rhs(1:8)=rhs(1:8)+
     1 dvol*(matmul(transpose(Bmat),Stress)*((1.d0-phin)**2+xk))       
            
        amatrx(9:12,9:12)=amatrx(9:12,9:12)+
     1    dvol*(matmul(transpose(dNdx),dNdx)*Gc*xlc+
     2    matmul(dN,transpose(dN))*(Gc/xlc+2.d0*H))
           
        rhs(9:12)=rhs(9:12)+
     1    dvol*(matmul(transpose(dNdx),matmul(dNdx,phik(1:4)))
     2    *Gc*xlc+dN(1:4,1)*((Gc/xlc+2.d0*H)*phi-2.d0*H)) 



    
        do i=1,ntens
              saveVars(40 + 10*(intPnt-1) + i) = 
     1        Stress(i)*((1.d0-phin)**2+xk) 
              saveVars(40 + 10*(intPnt-1) + i + 4) = Strain(i)
        end do

      end do    ! end loop on material integration points 
c ---   ////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\   
c ---  ////////end loop on material integration points\\\\\\\       
c --- //////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\ 
      do i=1,4
          fint(i*3) = rhs(8+i)
      end do
      do i=1,2
            fint(i) = rhs(i)
      end do
      do i=4,5
            fint(i) = rhs(i-1)
      end do  
      do i=7,8
            fint(i) = rhs(i-2)
      end do  
      do i=10,11
            fint(i) = rhs(i-3)
      end do   
      
      do k1=1,4
          do i=1,12
              estiff(3*k1,i) = 0.d0
          end do    
      end do
      do k1=1,4
              estiff(3*k1,3) = amatrx(8+k1,9)
              estiff(3*k1,6) = amatrx(8+k1,10)
              estiff(3*k1,9) = amatrx(8+k1,11)
              estiff(3*k1,12) = amatrx(8+k1,12)
      end do    
      do k1=1,2
              estiff(k1,1) = amatrx(k1,1)
              estiff(k1,2) = amatrx(k1,2)
              estiff(k1,3) = 0.d0
              estiff(k1,4) = amatrx(k1,3)
              estiff(k1,5) = amatrx(k1,4)
              estiff(k1,6) = 0.d0
              estiff(k1,7) = amatrx(k1,5)
              estiff(k1,8) = amatrx(k1,6)
              estiff(k1,9) = 0.d0
              estiff(k1,10) = amatrx(k1,7)
              estiff(k1,11) = amatrx(k1,8)
              estiff(k1,12) = 0.d0
      end do
      do k1=4,5
              estiff(k1,1) = amatrx(k1-1,1)
              estiff(k1,2) = amatrx(k1-1,2)
              estiff(k1,3) = 0.d0
              estiff(k1,4) = amatrx(k1-1,3)
              estiff(k1,5) = amatrx(k1-1,4)
              estiff(k1,6) = 0.d0
              estiff(k1,7) = amatrx(k1-1,5)
              estiff(k1,8) = amatrx(k1-1,6)
              estiff(k1,9) = 0.d0
              estiff(k1,10) = amatrx(k1-1,7)
              estiff(k1,11) = amatrx(k1-1,8)
              estiff(k1,12) = 0.d0
      end do 
      do k1=7,8
              estiff(k1,1) = amatrx(k1-2,1)
              estiff(k1,2) = amatrx(k1-2,2)
              estiff(k1,3) = 0.d0
              estiff(k1,4) = amatrx(k1-2,3)
              estiff(k1,5) = amatrx(k1-2,4)
              estiff(k1,6) = 0.d0
              estiff(k1,7) = amatrx(k1-2,5)
              estiff(k1,8) = amatrx(k1-2,6)
              estiff(k1,9) = 0.d0
              estiff(k1,10) = amatrx(k1-2,7)
              estiff(k1,11) = amatrx(k1-2,8)
              estiff(k1,12) = 0.d0
      end do
      do k1=10,11
              estiff(k1,1) = amatrx(k1-3,1)
              estiff(k1,2) = amatrx(k1-3,2)
              estiff(k1,3) = 0.d0
              estiff(k1,4) = amatrx(k1-3,3)
              estiff(k1,5) = amatrx(k1-3,4)
              estiff(k1,6) = 0.d0
              estiff(k1,7) = amatrx(k1-3,5)
              estiff(k1,8) = amatrx(k1-3,6)
              estiff(k1,9) = 0.d0
              estiff(k1,10) = amatrx(k1-3,7)
              estiff(k1,11) = amatrx(k1-3,8)
              estiff(k1,12) = 0.d0
      end do            

     
 990  CONTINUE
      ! do intPnt=1,4

      !enddo
      RsltVar = saveVars 
c      
      RsltBsc=0
      do intPnt=1,4
       do j=1,4
         RsltBsc(14*(intPnt-1) + j )=  saveVars(40+10*(intPnt-1) + j) 
         RsltBsc(14*(intPnt-1) + j +7) = 
     &                             saveVars(40+10*(intPnt-1) + j + 4) 
       enddo
      end do    
           
      return
      
      end
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
      subroutine straininc(ntens,ndim,nnode,mlvarx,bmat,du,dstran,u,xx1,
     1 xx1Old)
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

      DOUBLE PRECISION dstran(ntens),bmat(ndim,nnode),
     & du(mlvarx,*),xdu(3),
     & xx1(3,3),u(mlvarx,*),utmp(3),utmpOld(3),xx1Old(3,3)
      INTEGER k1,i, nnode, nodi, incr_row, mlvarx,
     & ntens,ndim
      DOUBLE PRECISION dNidx, dNidy
      dstran=0.d0
      ! set xx1 to Identity matrix
      xx1=0.d0
      xx1Old=0.d0
      do k1=1,3
       xx1(k1,k1)=1.d0
       xx1Old(k1,k1)=1.d0         
      end do

c************************************
c    Compute incremental strains
c************************************
c
      do nodi=1,nnode
           
       incr_row=(nodi-1)*ndim
       do i=1,ndim
        xdu(i)=du(i+incr_row,1)
        utmp(i)=u(i+incr_row,1)
        utmpOld(i)=utmp(i)-xdu(i)
       end do

       dNidx=bmat(1,nodi)
       dNidy=bmat(2,nodi)

       dstran(1)=dstran(1)+dNidx*xdu(1)
       dstran(2)=dstran(2)+dNidy*xdu(2)
       dstran(4)=dstran(4)+dNidy*xdu(1)+dNidx*xdu(2)  

c     deformation gradient
       xx1(1,1)=xx1(1,1)+dNidx*utmp(1)
       xx1(1,2)=xx1(1,2)+dNidy*utmp(1)
       xx1(2,1)=xx1(2,1)+dNidx*utmp(2)
       xx1(2,2)=xx1(2,2)+dNidy*utmp(2)
c
       xx1Old(1,1)=xx1Old(1,1)+dNidx*utmpOld(1)
       xx1Old(1,2)=xx1Old(1,2)+dNidy*utmpOld(1)
       xx1Old(2,1)=xx1Old(2,1)+dNidx*utmpOld(2)
       xx1Old(2,2)=xx1Old(2,2)+dNidy*utmpOld(2)

      end do
c
      return
      end
