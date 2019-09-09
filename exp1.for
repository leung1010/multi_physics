!***********************************************************************

        module global

        ! This module is used to transfer SDV's from the UEL
        !  to the UVARM so that SDV's can be visualized on a
        !  dummy mesh
        !
        !  globalSdv(X,Y,Z)
        !   X - element pointer
        !   Y - integration point pointer
        !   Z - SDV pointer
        !
        !  numElem
        !   Total number of elements in the real mesh, the dummy
        !   mesh needs to have the same number of elements, and 
        !   the dummy mesh needs to have the same number of integ
        !   points.  You must set that parameter value here.
        !
        !  ElemOffset
        !   Offset between element numbers on the real mesh and
        !    dummy mesh.  That is set in the input file, and 
        !    that value must be set here the same.
  
        integer numElem,ElemOffset,err
  
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Set the number of UEL elements used here
        parameter(numElem=352)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Set the offset here for UVARM plotting, must match input file!
        parameter(ElemOffset=1000)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
        real*8, allocatable :: globalSdv(:,:,:)
  
        end module global

***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.
     
      use global
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

      uvar(1) = globalSdv(noel-ElemOffset,npt,1)
c      for example
c      uvar(2) = globalSdv(noel-ElemOffset,npt,2)
c      uvar(3) = globalSdv(noel-ElemOffset,npt,3)
c      uvar(4) = globalSdv(noel-ElemOffset,npt,4)

      RETURN
      END SUBROUTINE UVARM
*************************************************************************
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1  PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2  KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3  NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4  PERIOD)

        use global

C
        IMPLICIT NONE
C
C       VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
C
        REAL(8) :: RHS,AMATRX,SVARS,ENERGY
C
C       VARIABLES PASSED INTO UEL
C
        REAL(8) :: PROPS,COORDS,U,DU,V,A,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
        INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,
     1  KSTEP,KINC,JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,
     2  MDLOAD,JPROPS,NJPROP
C
        DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1 SVARS(*),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2 U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3 PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4 DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5 JPROPS(*)

        INTEGER NINTPT, NDIM, lenJobName, lenOutDir
        character*256 jobName,outDir,fileName
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        PARAMETER(NINTPT = 9) ! number of integration points !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
        !============================debug===========================
          logical :: FirstCall=.true.
          integer :: dummyVar
          if (FirstCall==.true.) then
              write(*,*) 'please input an integer'
              read(*,*) dummyVar
              FirstCall = .false.  
          end if
          dummyVar = 1234 
        !==========================debug============================
          


        !!---------------------------------------------------------
        !  Preform the initial checks
        !
        !
        ! Open the debug/error message file
        !
        call getJobName(jobName,lenJobName)
        call getOutDir(outDir,lenOutDir)
        fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
        open(unit=80,file=fileName,status='unknown')
        !
        !  check the procedure type, this should be a coupled
        !  tempreture displacement or pore pressure displacement
        !  which are any of the following (64,65,72,73)
        !
        IF((LFLAGS(1).EQ.64).OR.(LFLAGS(1).EQ.65).OR.
     +  (LFLAGS(1).EQ.72).OR.(LFLAGS(1).EQ.73)) THEN
            !
            ! any is OK
            !
        ELSE
            WRITE(*,*) 'ABAQUS does not have the right procedure'
            WRITE(*,*) 'go back and chekc the procedure type'
            WRITE(*,*) 'lflags(1)=', lflags(1)
        ENDIF
        !
        !  do nothing if a "dummy" step
        !
        IF(DTIME.EQ.0.0) RETURN
        !
        !  Done with the initial checks
        !!----------------------------------------------------------


        !!----------------------------------------------------------
        ! Call the special element to perform the analysis
        !
        IF(JTYPE.EQ.1) THEN
            !
            !  this is a plane strain analysis
            !
            NDIM = 2
            CALL CHEMOMECHAN8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     +               NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,
     +               A,JTYPE,TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,
     +               JDLTYP,ADLMAG,PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,
     +               MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD,
     +               NDIM,NINTPT)
        ELSE
            !
            !  emmm...I think it should take some time to solve the problem
            !
            WRITE(*,*) 'Element type not supported, jtype = ', jtype
        ENDIF
        !
        !  Done with this element, RHS and AMATRX already returned
        !  as output from the spcific element subroutine called
        !!-------------------------------------------------------------


        RETURN
        END SUBROUTINE UEL






*************************************************************************************
        SUBROUTINE CHEMOMECHAN8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     1   NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,
     2   DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3   NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4   PERIOD,
     5   NDIM,NINTPT)
!        PROGRAM CHEMOMECHAN8
!
        use global
        IMPLICIT NONE
!
!       VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
!
        REAL(8) :: RHS,AMATRX,SVARS,ENERGY
!
!       VARIABLES PASSED INTO UEL
!
        REAL(8) :: PROPS,COORDS,U,DU,V,A,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
        INTEGER :: NDOFEL,NRHS,MCRD,NNODE,JTYPE,
     1  KSTEP,KINC,JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,
     2  MDLOAD,JPROPS,NJPROP,NSVARS,NPROPS
        

!
        DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1  SVARS(*),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2  U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3  PARAMS(*),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4  DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5  JPROPS(*)
        


        REAL*8 Unew(NNODE,2), DUnew(NNODE,2), Uold(NNODE,2)
        REAL*8 Cnew(NNODE), DCnew(NNODE), Cold(NNODE), C_tau, C_t, DCDT
        REAL*8 Kuu(2*NNODE,2*NNODE), Kuc(2*NNODE,NNODE)
        REAL*8 Kcc(NNODE,NNODE), Ru(2*NNODE,1), Rc(NNODE,1)
        REAL*8 GPTX(NINTPT,1), GPTY(NINTPT,2), N(NNODE), w(NINTPT)
        REAL*8 DNST(NNODE,2), DNXY(NNODE,2)
        REAL*8 U_tau, U_t, DUDT, E, MU, Cmat(3,3), Bmatu(3,2*NNODE)
        REAL*8 Bmatm(1,2*NNODE), Uall(2*NNODE), Rgas, Tem, Cini
        REAL*8 Uallold(2*NNODE), DUDTmat(2*NNODE), DCDTmat(NNODE)
        REAL*8 Bmatc(2,NNODE), Phi, PHImat(2,2), detmapJ, C, Le
        REAL*8 S(3,1), BmUDT, SHAF(1,8), Ec(2,1), Sphi(2,1)
        REAL*8 Kcu(NNODE,2*NNODE), SIGMA(3)
        REAL*8 BuC(2*NNODE,3), BcPhi(NNODE,2) !CHANGE


!===========================  DEBUG  ==================================

!===================================================================



        INTEGER i, j, k, kk, NINTPT, NDIM, INTPT
        INTEGER nlSdv, ngSdv, stat, jj
        

        REAL*8 ZERO, ONE, TWO, THREE, HALF, PI, THIRD
        PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0,
     +   HALF=1.D0/2.D0, THIRD=1.D0/3.D0, PI=3.141592653D0)

     
!===========================  DEBUG  ==================================

!===================================================================
        ! Get element parameters
        !
        nlSdv = JPROPS(1) ! number of local sdv's per integ point
        ngSdv = JPROPS(2) ! number of global sdv's per integ point
        ! Allocate memory for the globalSdv's used for viewing
        !  results on the dummy mesh
        !
        if(.not.allocated(globalSdv)) then
            !
            ! allocate memory for the globalSdv's
            !
            ! numElem needs to be set in the MODULE
            ! nInt needs to be set in the UEL
            !
            stat=0        
c           allocate(globalSdv(numElem,nInt,ngSdv))
c           deallocate(globalSdv)
            allocate(globalSdv(numElem,NINTPT,ngSdv),stat=err)

            if(stat.ne.0) then

                write(*,*) '///////////////////////////////////////////'
                write(*,*) 'error when allocating globalSdv'
                write(*,*) '//////////////////////////////////////////'
                write(*,*) '   stat=',stat
                write(*,*) '  ngSdv=',ngSdv
                write(*,*) '   NINTPT=',NINTPT
                write(*,*) 'numElem=',numElem
                write(*,*) '  nNode=',nNode
                write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
                write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
                write(*,*) '///////////////////////////////////////////'
            endif
            write(*,*) '-----------------------------------------------'
            write(*,*) '----------- globalSDV ALLOCATED ---------------'
            write(*,*) '-----------------------------------------------'
            write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS ---------'
            write(*,*) '---------- numElem=',numElem
            write(*,*) '---------- UPE4 ELEMENTS ----------------------'
            write(*,*) '-----------------------------------------------'
            write(*,*) '---------- YOU PUT NUMBER OF POINTS -----------'
            write(*,*) '---------- NINTPT =',NINTPT
            write(*,*) '-----------------------------------------------'
            write(*,*) '---------- YOU PUT NUMBER OF SDVs -------------'
            write(*,*) '---------- ngSdv=',ngSdv
            write(*,*) '-----------------------------------------------'
        endif

        E = PROPS(1)
        MU = PROPS(2)
        Rgas = PROPS(3)
        Tem = PROPS(4)
        Cini = PROPS(5)
        Phi = PROPS(6)


        !  Initial the residual and tangent matrices to zero
        !
        Kuu = ZERO
        Kuc = ZERO
        Kcu = ZERO
        Kcc = ZERO
        Ru = ZERO
        Rc = ZERO

        !  Obtain the nodal displacement and concentration
        !
        k = 0
        DO i = 1, NNODE
            DO j = 1, NDIM
                k = k + 1
                Unew(i,j) = U(k)
                DUnew(i,j) = DU(k,1)
                Uold(i,j) = Unew(i,j) - DUnew(i,j)
            ENDDO
            k = k + 1
            Cnew(i) = U(k)
            DCnew(i) = DU(k,1)
            Cold(i) = Cnew(i) - DCnew(i)
        ENDDO

        
        Uall = ZERO
        i = 0
        DO i = 1, NNODE
            Uall(2*i-1) = Unew(i,1)
            Uall(2*i) = Unew(i,2)
            Uallold(2*i-1) = Uold(i,1)
            Uallold(2*i) = Uold(i,2)
        ENDDO


        ! Impose any time-stepping changes on the increments of
        !  concentration or displacement if you wish
        !
        !  concentration increment
        !
        i = 0
        DO i=1,NNODE
            IF(DABS(DCnew(i)).GT.1.D6) THEN
            PNEWDT = 0.5
            !RETURN
            ENDIF
        ENDDO
        !
        ! displacement increment, based on element diagonal
        !
        Le = DSQRT(((COORDS(1,1)-COORDS(1,4))**TWO) + 
     +     ((COORDS(2,1)-COORDS(2,4))**TWO))
        !
        i = 0
        j = 0
        DO i=1,NNODE
            DO j=1,NDIM
                IF(DABS(DUnew(i,j)).GT.10.0*Le) THEN
                PNEWDT = 0.5
                !RETURN
                ENDIF
            ENDDO
        ENDDO

        
        !!------------------------------------------------------------------------
        !
        !  Begin the loop over integration points
        !
        !  Obtain the integration point local coordinates and weights
        !
        IF(NNODE.EQ.8) THEN
            IF(NINTPT.EQ.9) THEN
                CALL XINT2D9PT(GPTX,GPTY,w,NINTPT) ! 9-pt integration
            ELSEIF(NINTPT.EQ.4) THEN
                CALL XINT2D4PT(GPTX,GPTY,w,NINTPT) ! 4-pt integration
            ELSE
                WRITE(*,*) 'Invalid number of integration points, 
     +                      NINTPT',NINTPT
            ENDIF

        ELSE
            WRITE(*,*) 'Incorrect number of nodes: NNODE.EQ.8'
        ENDIF

        !  Loop over the integration point
        !
        jj = 0
        INTPT = 0 !CHANGE
        DO INTPT = 1, NINTPT

            ! Obtain state variables from previous increment
            !
            if((kinc.le.1).and.(kstep.eq.1)) then
                !
                ! this is the first increment, of the first step
                !  give initial conditions
                !
                C_t  = Cini
                !
            else
                !
                ! this is not the first increment, read old values
                !
                C_t  = svars(1+jj)
                !
            endif

            !  Obtain the shape function
            !
            IF(NNODE.EQ.8) THEN
                CALL CALCSHAPE2DBILINEAR(NINTPT,GPTX,GPTY,INTPT,N,DNST)
            ELSE
                WRITE(*,*) 'Incorrect number of nodes: NNODE.EQ.8'
            ENDIF

            !  Map shape functions from local to global reference coordinate system
            !
            CALL MAPSHAPE2D(NNODE,DNST,COORDS,DNXY,detmapJ)

            !  Obtain the concentration at this INTPT at the beginning and 
            !  end of the increment
            U_tau = ZERO
            U_t = ZERO
            C_tau = ZERO
            !C_t = ZERO
            DCDTmat = ZERO
            k = 0
            DO k = 1, NNODE
                C_tau = C_tau + Cnew(k)*N(k)
                C_t = C_t + Cold(k)*N(k)
                !DCDTmat(k) = (Cnew(k) - Cold(k))/DTIME
                DCDTmat(k) = DCnew(k)/DTIME
            ENDDO

            do k = 1, NNODE
                C_tau = C_t + DCnew(k)*N(k)
            enddo

            DCDT = (C_tau - C_t)/DTIME
            
            k = 0
            DUDTmat = ZERO
            DO k = 1, 2*NNODE
                DUDTmat(k) = (Uall(k) - Uallold(k))/DTIME
            ENDDO

            ! Save the state variables at this integ point
            !  at the end of the increment
            !
            svars(1+jj) = C_tau
            jj = jj + nlSdv ! setup for the next intPt


            ! Save the state variables at this integ point in the
            !  global array used for plotting field output
            !
            globalSdv(JELEM,INTPT,1) = C_tau   ! polymer volume fraction




            C = E*(1-MU)/((1+MU)*(1-2*MU))
            Cmat = ZERO
            Cmat(1,1) = 1.D0*C
            Cmat(1,2) = MU/(1-MU)*C
            Cmat(2,1) = MU/(1-MU)*C
            Cmat(2,2) = 1.D0*C
            Cmat(3,3) = (1-2*MU)/(2*(1-MU))*C

            PHImat = ZERO
            PHImat(1,1) = Phi
            PHImat(2,2) = Phi


            !  Compute/update the displacement residual vector
            !
            Bmatu = ZERO
            kk = 0
            DO kk = 1, NNODE
                Bmatu(1,1+NDIM*(kk-1)) = DNXY(kk,1)
                Bmatu(2,2+NDIM*(kk-1)) = DNXY(kk,2)
                Bmatu(3,1+NDIM*(kk-1)) = DNXY(kk,2)
                Bmatu(3,2+NDIM*(kk-1)) = DNXY(kk,1)
            ENDDO

            Bmatm = ZERO
            kk = 0
            DO kk = 1, NNODE
                Bmatm(1,1+NDIM*(kk-1)) = DNXY(kk,1)
                Bmatm(1,2+NDIM*(kk-1)) = DNXY(kk,2)
            ENDDO

            SIGMA = MATMUL(MATMUL(Cmat,Bmatu),Uall)
            i = 0
            DO i = 1,3
                S(i,1) = SIGMA(i)
            ENDDO
    

            Ru = Ru +  (
     +                  -MATMUL(TRANSPOSE(Bmatu),S)+
     +                      Rgas*Tem*TRANSPOSE(Bmatm)*C_tau
     +                                  )*detmapJ*w(INTPT)

            !  Compute/update the concentration residual vector
            !
            Bmatc = ZERO
            kk = 0
            DO kk = 1, NNODE
                Bmatc(1,kk) = DNXY(kk,1)
                Bmatc(2,kk) = DNXY(kk,2)
            ENDDO

            BmUDT = 0.D0
            i = 0
            DO i = 1,16
                BmUDT = BmUDT + Bmatm(1,i)*DUDTmat(i)
            ENDDO
    
            SHAF = 0.D0
            i = 0
            DO i = 1, 8
                SHAF(1,i) = N(i)
            ENDDO
        
            Ec = 0.D0
            i = 0
            j = 0
            DO i = 1,2
                DO j = 1,8
                    Ec(i,1) = Ec(i,1) + Bmatc(i,j)*Cnew(j)
                ENDDO
            ENDDO
    
            Sphi = MATMUL(Phimat,Ec) 

!            Rc = Rc + (
!     +                  Cini*TRANSPOSE(SHAF)*BmUDT-TRANSPOSE(SHAF)*DCDT-
!     +                  MATMUL(TRANSPOSE(Bmatc),Sphi)
!     +                                          )*detmapJ*w(INTPT)

            Rc = Rc + (
     +                  DTIME*MATMUL(TRANSPOSE(Bmatc),Sphi)
     +                                          )*detmapJ*w(INTPT)
            !  Compute/update the displacement tangent matrix
            !
            BuC = MATMUL(TRANSPOSE(Bmatu),Cmat)
            Kuu = Kuu + (
     +                   MATMUL(BuC,Bmatu)
     +                                        )*detmapJ*w(INTPT)


            !  Compute/update the concentration tangent matrix
            !
            BcPhi = MATMUL(TRANSPOSE(Bmatc),PHImat)
            Kcc = Kcc + (
     +                   -MATMUL(TRANSPOSE(SHAF),SHAF)-
     +                   MATMUL(BcPhi,Bmatc)*DTIME
     +                                       )*detmapJ*w(INTPT)

            !  Compute/update the displacement-concentration tangent matrix
            !
            Kuc = Kuc + (
     +                  -Rgas*Tem*MATMUL(TRANSPOSE(Bmatm),SHAF)
     +                        )*detmapJ*w(INTPT)

            !  Compute/update the concentration-displacement tangent matrix
            !
            Kcu = Kcu + (
     +                    Cini*MATMUL(TRANSPOSE(SHAF),Bmatm)
     +                            )*detmapJ*w(INTPT)

            !
            !  End loop over integration points
            !!------------------------------------------------------------------


            !!------------------------------------------------------------------
            !  Return ABAQUS the RHS vector and the AMATRX stiffiness matrix
            !
            CALL ASSEMBLEELEMENT(NDIM,NNODE,NDOFEL,
     +                           Ru,Rc,Kuu,Kcc,Kuc,Kcu,
     +                           RHS,AMATRX)
            !
            !  End return the RHS and AMATRX
            !!------------------------------------------------------------------

        ENDDO
!        WRITE(*,*) RHS
        RETURN
        END SUBROUTINE CHEMOMECHAN8
!        END PROGRAM CHEMOMECHAN8

************************************************************************************

        SUBROUTINE XINT2D4PT(Gptx,Gpty,w,nIntPt)
        !
        ! This subroutine will get the integration point locations
        !  and corresponding gauss quadrature weights for 2D elements
        !  using 4 gauss points for integration
        !
        !  xi(nIntPt,2): xi,eta coordinates for the integration pts
        !  w(nIntPt):    corresponding integration weights
        !
        IMPLICIT NONE
        !
        INTEGER nIntPt,nDim
        !
        REAL*8 xi(4,2), w(4), Gptx(4,1), Gpty(4,1)


        ! Initialize
        !
        w = 0.D0
        Gptx = 0.D0
        Gpty = 0.D0


        ! Number of Gauss points
        !
        nIntPt = 4


        ! Gauss weights
        !
        w(1) = 1.d0
        w(2) = 1.d0
        w(3) = 1.d0
        w(4) = 1.d0
    

        ! Gauss pt locations in master element
        !
        Gptx(1,1) = -dsqrt(1.d0/3.d0)
        Gpty(1,1) = -dsqrt(1.d0/3.d0)
        Gptx(2,1) = dsqrt(1.d0/3.d0)
        Gpty(2,1) = -dsqrt(1.d0/3.d0)
        Gptx(3,1) = -dsqrt(1.d0/3.d0)
        Gpty(3,1) = dsqrt(1.d0/3.d0)
        Gptx(4,1) = dsqrt(1.d0/3.d0)
        Gpty(4,1) = dsqrt(1.d0/3.d0)


        RETURN
        END SUBROUTINE XINT2D4PT


************************************************************************************

        SUBROUTINE XINT2D9PT(Gptx,Gpty,w,nIntPt)
        !
        ! This subroutine will get the integration point locations
        !  and corresponding gauss quadrature weights for 2D elements
        !  using 4 gauss points for integration
        !
        !  xi(nIntPt,2): xi,eta coordinates for the integration pts
        !  w(nIntPt):    corresponding integration weights
        !
        IMPLICIT NONE
        !
        INTEGER nIntPt,nDim
        !
        REAL*8 Gpty(9,1), Gptx(9,1), w(9)        !!!!!!HERE
        REAL*8 Wa, Wb, Wc
        PARAMETER(Wa=5.D0/9.D0, Wb=8.D0/9.D0, Wc=5.D0/9.D0)


        ! Initialize
        !
        w = 0.D0
        Gptx = 0.D0
        Gpty = 0.D0
        


        ! Number of Gauss points is 9
        !
        
        

        ! Gauss weights
        !
        CALL GSWT(w)


        ! Gauss pt locations in master element
        !
        CALL GetGptFI(Gptx,Gpty)


        RETURN
        END SUBROUTINE XINT2D9PT
*****************************************************************************
        SUBROUTINE GSWT(GWEI) 
        IMPLICIT NONE        
        REAL*8 GWEI(9), GWE(3), FIVE, EIGHT, NINE
        INTEGER I, J, NUMGP
        !
        PARAMETER(FIVE=5.D0,EIGHT=8.D0,NINE=9.D0) 
    !
    ! GWEI: the weights of Gaussian points 
    !
        GWE(1) = FIVE/NINE 
        GWE(2) = EIGHT/NINE 
        GWE(3) = FIVE/NINE 
        I = 0
        J = 0
        DO I=1,3 
            DO J=1,3 
                NUMGP=(I-1)*3+J 
                GWEI(NUMGP)=GWE(I)*GWE(J)
            END DO
        END DO 
        RETURN 
        END SUBROUTINE GSWT
*********************************************************************************************
        SUBROUTINE GetGptFI(Gptx,Gpty)
        IMPLICIT NONE
        REAL*8 AR(3), Gptx(9), Gpty(9), R
        REAL*8 ZERO, NEGONE, ONE, THREE, FIVE
        INTEGER NUMGP, I, J
    !
        PARAMETER(ZERO = 0.D0, NEGONE = -1.D0, ONE = 1.D0,  
     +              THREE = 3.D0, FIVE = 5.D0)
    !
    
        R = SQRT(THREE/FIVE)
    
    ! the sign of Gauss points
    !
        AR(1) = NEGONE
        AR(2) = ZERO
        AR(3) = ONE


    ! Gptx: x coordinates of Gaussian points 
    ! Gpty: y coordinates of Gaussian points
    !
        I = 0
        J = 0
        DO I=1,3 
            DO J=1,3
                NUMGP=(I-1)*3+J 
                Gptx(NUMGP)=AR(I)*R 
                Gpty(NUMGP)=AR(J)*R
            END DO
        END DO 

        RETURN
        END SUBROUTINE GetGptFI



*********************************************************************************************      
        SUBROUTINE MAPSHAPE2D(nNode,DNST,COORDS,DNXY,detmapJ)
        !
        ! Map derivatives of shape fns from xi-eta-zeta domain
        !  to x-y-z domain.
        !
        IMPLICIT NONE
        !
        INTEGER i,j,k,nNode
        !
        REAL*8 DNST(nNode,2),DNXY(nNode,2),COORDS(2,nNode),mapJ(2,2),
     +          mapJ_inv(2,2),detmapJ
        !
        REAL*8 ZERO,ONE,TWO,HALF,FOURTH,EIGHTH
        PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,HALF=0.5D0,FOURTH=0.25D0,
     +            EIGHTH=1.d0/8.D0)


        ! Calculate the mapping Jacobian matrix:
        !
        mapJ = ZERO
        i = 0
        j = 0
        DO i=1,2
            DO j=1,2
                DO k=1,nNode
                    mapJ(i,j) = mapJ(i,j) + DNST(k,i)*COORDS(j,k)
                ENDDO
            ENDDO
        ENDDO


      ! Calculate the inverse and the determinant of Jacobian
      !
      CALL matInv2D(mapJ,mapJ_inv,detmapJ)
      


      ! Calculate first derivatives wrt x, y, z
      !
      DNXY = transpose(matmul(mapJ_inv,transpose(DNST)))
      

      RETURN
      END SUBROUTINE MAPSHAPE2D



*************************************************************************************************  
      SUBROUTINE matInv2D(A,A_inv,det_A)
        !
        ! Returns A_inv, the inverse, and det_A, the determinant
        ! Note that the det is of the original matrix, not the
        ! inverse
        !
          IMPLICIT NONE
        !
          REAL*8 A(2,2),A_inv(2,2),det_A,det_A_inv
    
        
        
          det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
          
          IF (det_A .le. 0.d0) then
            WRITE(*,*) 'WARNING: subroutine matInv2D:'
            WRITE(*,*) 'WARNING: det of mat=',det_A
            RETURN
          ENDIF
              
          det_A_inv = 1.D0/det_A
            
          A_inv(1,1) =  det_A_inv*A(2,2)
          A_inv(1,2) = -det_A_inv*A(1,2)
          A_inv(2,1) = -det_A_inv*A(2,1)
          A_inv(2,2) =  det_A_inv*A(1,1)
    
    
          RETURN
          END SUBROUTINE matInv2D


**********************************************************************************************
        SUBROUTINE ASSEMBLEELEMENT(NDIM,NNODE,NDOFEL,
     +              Ru,Rc,Kuu,Kcc,Kuc,Kcu,
     +              RHS,AMATRX)

        !
        !  Subroutine to assemble the local elements residual and tangent
        !
        IMPLICIT NONE

        INTEGER i, j, k, A11, A12, B11, B12, NDIM, NNODE, NDOFEL,NDOFN

        REAL*8 Ru(NDIM*NNODE,1), Rc(NNODE,1), Kuu(NDIM*NNODE,NDIM*NNODE)
        REAL*8 Kcc(NNODE,NNODE), Kuc(NDIM*NNODE,NNODE), RHS(NDOFEL,1)
        REAL*8 Kcu(NNODE,NDIM*NNODE), AMATRX(NDOFEL,NDOFEL)

        !  Total number of degrees of freedom per node 
        !
        NDOFN = NDOFEL/NNODE

        !  initial
        !
        RHS(:,1) = 0.D0
        AMATRX = 0.D0

        !  Assemble the element level residual
        !
        i = 0
        DO i = 1, NNODE
            A11 = NDOFN*(i-1) + 1
            A12 = NDIM*(i-1) + 1
            !
            !  displacement
            !
            RHS(A11,1) = Ru(A12,1)
            RHS(A11+1,1) = Ru(A12+1,1)
            !
            !  concentration
            !
            RHS(A11+2,1) = Rc(i,1)
        ENDDO
        !
        !  Assemble the element level tangent matrix
        !
        i = 0
        j = 0
        DO i = 1, NNODE
            DO j = 1, NNODE
                A11 = NDOFN*(i-1) + 1
                A12 = NDIM*(i-1) + 1
                B11 = NDOFN*(j-1) + 1
                B12 = NDIM*(j-1) + 1
                !
                !  displacement
                !
                AMATRX(A11,B11) = Kuu(A12,B12)
                AMATRX(A11+1,B11) = Kuu(A12+1,B12)
                AMATRX(A11,B11+1) = Kuu(A12,B12+1)
                AMATRX(A11+1,B11+1) = Kuu(A12+1,B12+1)
                !
                !  concentration 
                !
                AMATRX(A11+2,B11+2) = Kcc(i,j)
                !
                !  displacement - concentration
                !
                AMATRX(A11,B11+2) = Kuc(A12,j)
                AMATRX(A11+1,B11+2) = Kuc(A12+1,j)
                !
                !  concentration - displacement
                !
                AMATRX(A11+2,B11) = Kuu(i,B12)
                AMATRX(A11+2,B11+1) = Kuu(i,B12+1)

            ENDDO
        ENDDO

        RETURN
        END SUBROUTINE ASSEMBLEELEMENT


**********************************************************************************************   
        SUBROUTINE CALCSHAPE2DBILINEAR(nIntPt,GPTX_int,GPTY_int,
     +                intpt,N,DNST)
        !
        ! Calculate the shape functions and their derivatives at the
        ! given integration point in the master element


        ! Calculate the shape functions and their derivatives at the
        ! given integration point in the master element
        !
        !         8                eta
        !   1-----------3          |
        !   |           |          |
        !   |           |          |
        !  5|           |7         |
        !   |           |          |
        !   |           |          O--------- xi
        !   2-----------4        origin at center
        !         6
        !
        ! sh(i) = shape function of node i at the intpt.
        ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
!
        IMPLICIT NONE
!
        INTEGER intpt,nDim,nIntPt
!
        REAL*8 GPTX_int(nIntPt,1),GPTY_int(nIntPt,1),N(8),DNST(8,2),s,t
!
        REAL*8 ZERO, ONE, FOURTH, SECOND, FOUR, TWO
        PARAMETER(ZERO=0.D0,ONE=1.D0,FOURTH=1.D0/4.D0,SECOND=1.D0/2.D0,
     +  FOUR=4.D0,TWO=2.D0)
            
      
        ! Location in the master element
        !
        s = GPTX_int(intpt,1)
        t = GPTY_int(intpt,1)
            
            
        ! The shape functions
        !    
        N(1) = FOURTH*(ONE-s)*(ONE-t)*(-s-t-ONE)
        N(2) = SECOND*(ONE-t)*(ONE+s)*(ONE-s)
        N(3) = FOURTH*(ONE+s)*(ONE-t)*(s-t-ONE)
        N(4) = SECOND*(ONE+s)*(ONE+t)*(ONE-t)
        N(5) = FOURTH*(ONE+s)*(ONE+t)*(s+t-ONE)
        N(6) = SECOND*(ONE+t)*(ONE+s)*(ONE-s)
        N(7) = FOURTH*(ONE-s)*(ONE+t)*(-s+t-ONE)
        N(8) = SECOND*(ONE-s)*(ONE+t)*(ONE-t)
        
        
        


        ! The first derivatives 13572468
        !
        DNST(1,1) = - (s/FOUR - FOURTH)*(t - ONE) - ((t - ONE)*
     +                (s + t + ONE))/FOUR
        DNST(5,1) = (t/TWO - SECOND)*(s - ONE) + 
     +                (t/TWO - SECOND)*(s + ONE)
        DNST(2,1) = ((t - ONE)*(t - s + ONE))/FOUR - 
     +              (s/FOUR + ONE/FOUR)*(t- ONE)
        DNST(6,1) = -((t - ONE)*(t + ONE))/TWO
        DNST(3,1) = (s/FOUR + FOURTH)*(t + ONE) + ((t + ONE)*
     +                (s + t - ONE))/FOUR
        DNST(7,1) = - (t/TWO + SECOND)*(s - ONE) - 
     +                (t/TWO + SECOND)*(s + ONE)
        DNST(4,1) = (s/FOUR - FOURTH)*(t + ONE) + 
     +              ((t + ONE)*(s - t + ONE))/FOUR
        DNST(8,1) = ((t - ONE)*(t + ONE))/TWO
        
        DNST(1,2) = - (s/FOUR - FOURTH)*(t - ONE) - 
     +         (s/FOUR - FOURTH)*(s + t + ONE)
        DNST(5,2) = ((s - ONE)*(s + ONE))/TWO
        DNST(2,2) = (s/FOUR + FOURTH)*(t - s + ONE) + 
     +                (s/FOUR + FOURTH)*(t - ONE)
        DNST(6,2) = - (s/TWO + SECOND)*(t - ONE) - 
     +                (s/TWO + SECOND)*(t + ONE)
        DNST(3,2) = (s/FOUR + FOURTH)*(t + ONE) +
     +             (s/FOUR + FOURTH)*(s + t - ONE)
        DNST(7,2) = -((s - ONE)*(s + ONE))/TWO
        DNST(4,2) = (s/FOUR - FOURTH)*(s - t + ONE) - 
     +                      (s/FOUR - FOURTH)*(t + ONE)
        DNST(8,2) = (s/TWO - SECOND)*(t - ONE) + 
     +                (s/TWO - SECOND)*(t + ONE)

        RETURN
        END SUBROUTINE CALCSHAPE2DBILINEAR