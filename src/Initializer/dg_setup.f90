!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Verena Hermann (hermann AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/hermann)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2009-2017, SeisSol Group
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice,
!!    this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions and the following disclaimer in the documentation
!!    and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its
!!    contributors may be used to endorse or promote products derived from this
!!    software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!! @section DESCRIPTION
!=============================================================================!
!!                                                                           !!
!!  ADER Discontinuous Galerkin Module                                       !!
!!                                                                           !!
!!---------------------------------------------------------------------------!!
!!                                                                           !!
!!  Equations                :  Linear Hyperbolic Equations                  !!
!!  Number of dimensions     :  3D                                           !!
!!  Element types supported  :  tetrahedrons                                 !!
!!  Spatial accuracy         :  arbitrary                                    !!
!!  Temporal accuracy        :  arbitrary                                    !!
!!                                                                           !!
!!---------------------------------------------------------------------------!!

! define preprocessor variables
#include <Initializer/preProcessorMacros.fpp>

MODULE dg_setup_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  INTERFACE closeGalerkin3D_us
     MODULE PROCEDURE closeGalerkin3D_us
  END INTERFACE
  INTERFACE iniGalerkin3D_us_level1_new
     MODULE PROCEDURE iniGalerkin3D_us_level1_new
  END INTERFACE
  INTERFACE iniGalerkin3D_us_level2_new
     MODULE PROCEDURE iniGalerkin3D_us_level2_new
  END INTERFACE
  INTERFACE iniGalerkin3D_us_intern_new
     MODULE PROCEDURE iniGalerkin3D_us_intern_new
  END INTERFACE
  INTERFACE icGalerkin3D_us_new
     MODULE PROCEDURE icGalerkin3D_us_new
  END INTERFACE

!  INTERFACE NonConformingGPEvaluation3D
!     MODULE PROCEDURE NonConformingGPEvaluation3D
!  END INTERFACE
  ! Private procedures and functions
  INTERFACE BuildSpecialDGGeometry3D_new
     MODULE PROCEDURE BuildSpecialDGGeometry3D_new
  END INTERFACE

  !
  !---------------------------------------------------------------------------!
  PUBLIC  ::                               &
       closeGalerkin3D_us,                 &
       iniGalerkin3D_us_level1_new,        &
       iniGalerkin3D_us_level2_new,        &
       iniGalerkin3D_us_intern_new,        &
       icGalerkin3D_us_new
!      NonConformingGPEvaluation3D
  !---------------------------------------------------------------------------!

CONTAINS

  !===========================================================================!
  !!                                                                         !!
  !! closeGalerkin deallocates the arrays used by RKDG method                !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE closeGalerkin3D_us(EQN,MESH,DISC,IO)
    !-------------------------------------------------------------------------!

    USE COMMON_operators_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)         :: EQN
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tDiscretization)    :: DISC
    TYPE(tInputOutput)       :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                  :: i,j,k,iElem,iSide
    INTEGER                  :: LocElemType
    INTEGER                  :: iLocalNeighborSide
    INTEGER                  :: iLocalNeighborVrtx
    !-------------------------------------------------------------------------!
    INTENT(IN)               :: IO
    INTENT(INOUT)            :: DISC
    !-------------------------------------------------------------------------!
    !
    IF(.NOT.DISC%Galerkin%init) THEN
       logError(*) 'closeGalerkin: SeisSol Interface not initialized!!'
       call exit(134)
    ENDIF
    !
    logInfo(*) 'Enter closeGalerkin...'
    !
    DEALLOCATE( DISC%Galerkin%geoNormals         )
    DEALLOCATE( DISC%Galerkin%geoTangent1        )
    DEALLOCATE( DISC%Galerkin%geoTangent2        )
    DEALLOCATE( DISC%Galerkin%geoSurfaces        )
    !
    IF (ASSOCIATED( DISC%Galerkin%Faculty)) DEALLOCATE(DISC%Galerkin%Faculty)
    IF (ASSOCIATED( DISC%Galerkin%TimeGaussP)) DEALLOCATE(DISC%Galerkin%TimeGaussP)
    IF (ASSOCIATED( DISC%Galerkin%TimeGaussW)) DEALLOCATE(DISC%Galerkin%TimeGaussW)
    ! Interoperability with C needs continuous arrays in memory.
    deallocate(disc%galerkin%dgvar)

    IF (ASSOCIATED( DISC%Galerkin%intGaussP_Hex)) DEALLOCATE(DISC%Galerkin%intGaussP_Hex)
    IF (ASSOCIATED( DISC%Galerkin%intGaussW_Hex)) DEALLOCATE(DISC%Galerkin%intGaussW_Hex)
    IF (ASSOCIATED( DISC%Galerkin%bndGaussP_Hex)) DEALLOCATE(DISC%Galerkin%bndGaussP_Hex)
    IF (ASSOCIATED( DISC%Galerkin%bndGaussW_Hex)) DEALLOCATE(DISC%Galerkin%bndGaussW_Hex)


    IF (ASSOCIATED( DISC%Galerkin%intGaussP_Tet)) DEALLOCATE(DISC%Galerkin%intGaussP_Tet)
    IF (ASSOCIATED( DISC%Galerkin%intGaussW_Tet)) DEALLOCATE(DISC%Galerkin%intGaussW_Tet)
    IF (ASSOCIATED( DISC%Galerkin%bndGaussP_Tet)) DEALLOCATE(DISC%Galerkin%bndGaussP_Tet)
    IF (ASSOCIATED( DISC%Galerkin%bndGaussW_Tet)) DEALLOCATE(DISC%Galerkin%bndGaussW_Tet)
    !
    DISC%Galerkin%init = .FALSE.
    !
    logInfo(*) 'closeGalerkin successful '
    !
  END SUBROUTINE closeGalerkin3D_us

!
!   unified solver uses routine_new
!

  !===========================================================================!
  !!                                                                         !!
  !! iniGalerkin3D_us_level1 initializes all essential DG variables that     !!
  !! do NOT depend on the mesh                                               !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE iniGalerkin3D_us_level1_new(OptionalFields,EQN,DISC,MESH,BND,IC,SOURCE,MPI,IO)

    USE DGBasis_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tUnstructOptionalFields)   :: OptionalFields
    TYPE(tEquations)                :: EQN
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tDiscretization)           :: DISC
    TYPE(tBoundary)                 :: BND
    TYPE(tInitialCondition)         :: IC
    TYPE(tSource)                   :: SOURCE
    TYPE(tMPI)                      :: MPI
    TYPE(tInputOutput)              :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: i, j, iElem, iDirac, iRicker
    REAL                            :: x,y,z
    REAL                            :: pi, earth_rad
    REAL                            :: e1(3), e2(3), e3(3)
    ! ------------------------------------------------------------------------!

    DISC%Galerkin%init = .TRUE.

    DISC%Galerkin%init     = .TRUE.
    DISC%Galerkin%nDegFr   = (DISC%Galerkin%nPoly+1)*(DISC%Galerkin%nPoly+2)*(DISC%Galerkin%nPoly+3)/6   !**3

    DISC%Galerkin%nPolyRec  = DISC%Galerkin%nPoly
    DISC%Galerkin%nDegFrRec = DISC%Galerkin%nDegFr


    DISC%Galerkin%nDegFrST = DISC%Galerkin%nDegFrRec*(DISC%Galerkin%nPolyRec+1)

    logInfo(*) 'Interface SEISSOL successful '

  END SUBROUTINE iniGalerkin3D_us_level1_new

  !===========================================================================!
  !!                                                                         !!
  !! iniGalerkin3D_us_level2 initializes all essential DG variables that     !!
  !! DO depend on the mesh                                                   !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE iniGalerkin3D_us_level2_new(OptionalFields,EQN,DISC,MESH,BND,IC,SOURCE,MPI,IO)

    USE DGBasis_mod
!    USE DGSponge_mod ! not yet done for hybrids/unified version
#ifdef PARALLEL
    USE MPIExchangeValues_mod
#endif
    use iso_c_binding, only: c_loc, c_null_char, c_bool
    use f_ftoc_bind_interoperability
    use calc_deltaT_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tUnstructOptionalFields)   :: OptionalFields
    TYPE(tEquations)                :: EQN
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tDiscretization)           :: DISC
    TYPE(tBoundary)                 :: BND
    TYPE(tInitialCondition)         :: IC
    TYPE(tSource)                   :: SOURCE
    TYPE(tMPI)                      :: MPI
    TYPE(tInputOutput)              :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: i, j, k, l, iElem, iDirac, iRicker
    INTEGER                         :: iDRFace
    INTEGER                         :: iCurElem
    INTEGER                         :: nDGWorkVar
    INTEGER                         :: allocstat
    INTEGER                         :: iIntGP,iTimeGP,NestSize,BlockSize,iFace
    INTEGER                         :: iDegFr_xi,iDegFr_tau,iDegFr_tau2
    REAL                            :: x,y,z
    REAL                            :: pi_const, earth_rad
    REAL                            :: e1(3), e2(3), e3(3)
    REAL                            :: xiGP,etaGP,zetaGP,tauGP
    REAL                            :: phi_i,psi_i,phi_j,psi_j
    REAL                            :: phi_xi_i(3), phi_xi_j(3)
    REAL                            :: psi_tau_i, psi_tau_j,psi_tau_grad_i2
    REAL                            :: TimeStiff(DISC%Galerkin%nPoly+1,DISC%Galerkin%nPoly+1)
    REAL                            :: TimeMass(DISC%Galerkin%nPoly+1,DISC%Galerkin%nPoly+1)
    REAL                            :: TimeF0(DISC%Galerkin%nPoly+1)
    INTEGER, POINTER                :: MPI_Dirac_Element(:,:)
    INTEGER                         :: iError,iCPU
    real                            :: l_timeStepWidth
    real                            :: l_loads(3), l_scalings(3), l_cuts(2), l_timeScalings(2), l_gts
    integer                         :: iObject, iSide, iNeighbor, MPIIndex
    real, target                    :: materialVal(EQN%nBackgroundVar)
    logical(kind=c_bool)                        :: enableFreeSurfaceIntegration
    
    ! ------------------------------------------------------------------------!
    !
    CALL iniGalerkin3D_us_intern_new(EQN, DISC, MESH, BND, IC, SOURCE, OptionalFields, IO, MPI)
    !
    CALL BuildSpecialDGGeometry3D_new(OptionalFields%BackgroundValue,EQN,MESH,DISC,BND,MPI,IO)

    ! we have the material parameters and insphere diameters, let's get some time step widths out of this

    ! Set used initial conditions.
    call c_interoperability_setInitialConditionType(trim(IC%cICType) // c_null_char)
    if (IC%cICType == "Travelling") then
      call c_interoperability_setTravellingWaveInformation(IC%origin, IC%kVec, IC%ampField)
    endif

    ! malloc fortran arrays
    call ini_calc_deltaT( eqn%eqType,     &
                          optionalFields, &
                          eqn,            &
                          mesh,           &
                          io,             &
                          mpi                )

    ! get the time step width for every tet
    call cfl_step( optionalFields, &
                   eqn,            &
                   mesh,           &
                   disc,           &
                   io,             &
                   mpi                )

    ! get gts time step width
    l_gts = minval( optionalFields%dt_convectiv(:) )
    if (l_gts .le. 0.0) then
      logError(*) 'Invalid timestep width'
      call MPI_ABORT(MPI%commWorld, 134)
    endif

#ifdef PERIODIC_LTS_SCALING
    ! compute total load per half-sapce
    !         _____________________
    !        /         /     /    /|
    !       /         /     /    / |
    !      /         /     /    /  |
    !     /         /     /    /   |
    !    /         /     /    /    |
    !   /_________/_____/____/     |
    !  |         |     |     |     |
    !  |    T    |  S  |  F  |     |
    !  |    h    |  e  |  i  |     |
    !  |    i    |  c  |  r  |     |_ _ _ _ _ _ _ _ _ _
    !  |    r    |  o  |  s  |    /                    /
    !  |    d    |  n  |  t  |   /                    /
    !  |         |  d  |     |  /       Symmetry     /
    !  |         |     |     | /                    /
    !  |_________|_____|_____|/ _ _ _ _ _ _ _ _ _ _/
    !  |s^2/2 + s|  s  |s^2/2|     |     |         |
    !  |         |     |     |     |     |         |
    ! -50     -cut2  -cut1   0    cut1  cut2      50

    ! derive scaling of time step widths
    l_timeScalings(1) = PERIODIC_LTS_SCALING * PERIODIC_LTS_SCALING
    l_timeScalings(2) = PERIODIC_LTS_SCALING

    ! derive scalings
    l_scalings(1) = l_timeScalings(2) / 2
    l_scalings(2) = 1
    l_scalings(3) = l_timeScalings(1) + l_timeScalings(2)

    ! derive cuts
    l_cuts(1) = 50.0d0 * ( l_scalings(1)                 ) / ( l_scalings(1) + l_scalings(2) + l_scalings(3) )
    l_cuts(2) = 50.0d0 * ( l_scalings(1) + l_scalings(2) ) / ( l_scalings(1) + l_scalings(2) + l_scalings(3) )
#endif

    ! propagate the time step width to time manager
    do iElem = 1, mesh%nElem
      l_timeStepWidth = optionalFields%dt_convectiv(iElem)

#ifdef PERIODIC_LTS_SCALING
      ! perform the scaling of the time step width
      if(     mesh%elem%xybary(1,iElem) > -l_cuts(2) .and.  mesh%elem%xybary(1,iElem) < -l_cuts(1) .or. \
              mesh%elem%xybary(1,iElem) >  l_cuts(1) .and.  mesh%elem%xybary(1,iElem) <  l_cuts(2) )    \
      then
        l_timeStepWidth = l_gts / ( l_timeScalings(2) - 1.0E-2 )
      elseif( mesh%elem%xybary(1,iElem) > -l_cuts(1) .and.  mesh%elem%xybary(1,iElem) <  0.d0 .or.      \
              mesh%elem%xybary(1,iElem) >  0.d0      .and.  mesh%elem%xybary(1,iElem) <  l_cuts(1) )    \
      then
        l_timeStepWidth = l_gts / ( l_timeScalings(1) + 1.0E-2 )
      else
        l_timeStepWidth = l_gts
      endif
#endif

      call c_interoperability_setTimeStepWidth( i_meshId        = iElem,          &
                                                i_timeStepWidth = l_timeStepWidth &
                                              )
    enddo

    enableFreeSurfaceIntegration = (io%surfaceOutput > 0)
    ! put the clusters under control of the time manager
    call c_interoperability_initializeClusteredLts(&
            i_clustering = disc%galerkin%clusteredLts, &
            i_enableFreeSurfaceIntegration = enableFreeSurfaceIntegration, &
            usePlasticity = logical(EQN%Plasticity == 1, 1))

    !
    SELECT CASE(DISC%Galerkin%DGMethod)
    CASE(3)
        nDGWorkVar = EQN%nVar+EQN%nAneFuncperMech
    CASE DEFAULT
        nDGWorkVar = EQN%nVarTotal
    END SELECT
    !
    ! Allocation of arrays
    ALLOCATE(                                                                                   &
         DISC%Galerkin%dgvar( DISC%Galerkin%nDegFr,EQN%nVarTotal,MESH%nElem,DISC%Galerkin%nRK), &
         STAT = allocstat                                                                       )
    IF(allocStat .NE. 0) THEN
       logError(*) 'could not allocate all variables!'
       call MPI_ABORT(MPI%commWorld, 134)
    END IF
    !
    IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        ALLOCATE( DISC%Galerkin%DGTaylor(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly,MESH%nElem), &
                  STAT = allocstat )
        IF(allocStat .NE. 0) THEN
           logError(*) 'could not allocate DISC%Galerkin%DGTayl.'
           call MPI_ABORT(MPI%commWorld, 134)
        END IF
    ENDIF

#ifdef PARALLEL
    IF(MPI%nCPU.GT.1) THEN                                          !
                                                                    !
    ! We comunicate background fields                               !
    ! Needed for discontinuous Riemann problems                     !
    !                                                               !
    CALL MPIExchangeBackground(DISC           = DISC,             & !
                               EQN            = EQN,              & !
                               BND            = BND,              & !
                               MESH           = MESH,             & !
                               IO             = IO,               & !
                               OptionalFields = OptionalFields,   & !
                               MPI            = MPI               ) !
                                                                    !
    ENDIF                                                           !
#endif
    !
    ! Initialize sparse star matrices
    !
  do iElem = 1, MESH%nElem
    iSide = 0

    materialVal = OptionalFields%BackgroundValue(iElem,:)
    call c_interoperability_setMaterial( i_elem = iElem,                                         \
                                         i_side = iSide,                                         \
                                         i_materialVal = materialVal,                            \
                                         i_numMaterialVals = EQN%nBackgroundVar                  )

    do iSide = 1,4
      IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
          iObject         = MESH%ELEM%BoundaryToObject(iSide,iElem)
          MPIIndex        = MESH%ELEM%MPINumber(iSide,iElem)
          materialVal = BND%ObjMPI(iObject)%NeighborBackground(1:EQN%nBackgroundVar,MPIIndex) ! rho,mu,lambda
      ELSE
          SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
          CASE(0)
              iNeighbor       = MESH%ELEM%SideNeighbor(iSide,iElem)
              materialVal = OptionalFields%BackgroundValue(iNeighbor,:)
          CASE DEFAULT ! For boundary conditions take inside material
              materialVal = OptionalFields%BackgroundValue(iElem,:)
          END SELECT
      ENDIF
      call c_interoperability_setMaterial( i_elem = iElem,                        \
                                           i_side = iSide,                        \
                                           i_materialVal = materialVal,       \
                                           i_numMaterialVals = EQN%nBackgroundVar )

    enddo
  enddo

#ifdef USE_MPI
  ! synchronize redundant cell data
  logInfo0(*) 'Synchronizing copy cell material data.';
  call c_interoperability_synchronizeCellLocalData(logical(EQN%Plasticity == 1, 1))
#endif

    call c_interoperability_initializeMemoryLayout(&
            clustering = disc%galerkin%clusteredLts, &
            enableFreeSurfaceIntegration = enableFreeSurfaceIntegration, &
            usePlasticity = logical(EQN%Plasticity == 1, 1))

  ! Initialize source terms
  select case(SOURCE%Type)
    case(0)
      ! No source terms
      ! Do nothing
    case(42)
      call c_interoperability_setupNRFPointSources(trim(SOURCE%NRFFileName) // c_null_char)
    case(50)
      call c_interoperability_setupFSRMPointSources( momentTensor           = SOURCE%RP%MomentTensor,      &
                                                     solidVelocityComponent = SOURCE%RP%SolidVelocityComponent, &
                                                     pressureComponent      = SOURCE%RP%PressureComponent, &
                                                     fluidVelocityComponent = SOURCE%RP%FluidVelocityComponent, &
                                                     numberOfSources   = SOURCE%RP%nSbfs(1),          &
                                                     centres           = SOURCE%RP%SpacePosition,     &
                                                     strikes           = SOURCE%RP%Strks,             &
                                                     dips              = SOURCE%RP%Dips,              &
                                                     rakes             = SOURCE%RP%Rake,              &
                                                     onsets            = SOURCE%RP%Tonset,            &
                                                     areas             = SOURCE%RP%Area,              &
                                                     timestep          = SOURCE%RP%t_samp,            &
                                                     numberOfSamples   = SOURCE%RP%nsteps,            &
                                                     timeHistories     = SOURCE%RP%TimeHist       )
    case default
      logError(*) 'Generated Kernels: Unsupported source type: ', SOURCE%Type
      call MPI_ABORT(MPI%commWorld, 134)
  end select

  if (DISC%Galerkin%FluxMethod .ne. 0) then
    logError(*) 'Generated kernels currently supports Godunov fluxes only.'
    call MPI_ABORT(MPI%commWorld, 134)
  endif

  call c_interoperability_initializeEasiBoundaries(trim(EQN%BoundaryFileName) // c_null_char)
  call c_interoperability_initializeGravitationalAcceleration(EQN%GravitationalAcceleration)

  logInfo0(*) 'Initializing element local matrices.'
  call c_interoperability_initializeCellLocalMatrices(logical(EQN%Plasticity == 1, 1))

  IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        ALLOCATE( DISC%LocalIteration(MESH%nElem) )
        ALLOCATE( DISC%LocalTime(MESH%nElem)      )
        ALLOCATE( DISC%LocalDt(MESH%nElem)        )
        DISC%LocalIteration(:)  = 0.
        DISC%LocalTime(:)       = 0.
    ENDIF
    !
    IF(EQN%DR.EQ.1) THEN

      ALLOCATE(DISC%DynRup%SlipRate1(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%SlipRate2(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%Slip(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%Slip1(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%Slip2(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%TracXY(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%TracXZ(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%Mu(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%PeakSR(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%rupture_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%dynStress_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      
      ! TODO: Transpose StateVar
      ALLOCATE(DISC%DynRup%StateVar(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      !

      ! Initialize w/ first-touch
      !$omp parallel do schedule(static)
      DO i=1,MESH%fault%nSide
          DISC%DynRup%SlipRate1(:,i) = EQN%IniSlipRate1
          DISC%DynRup%SlipRate2(:,i) = EQN%IniSlipRate2
          DISC%DynRup%Slip(:,i) = 0.0
          DISC%DynRup%Slip1(:,i) = 0.0
          DISC%DynRup%Slip2(:,i) = 0.0
          DISC%DynRup%TracXY(:,i) = 0.0
          DISC%DynRup%TracXZ(:,i) = 0.0
          DISC%DynRup%StateVar(:,i) = EQN%IniStateVar(:,i)
          DISC%DynRup%Mu(:,i) = EQN%IniMu(:,i)
          DISC%DynRup%PeakSR(:,i) = 0.0
          DISC%DynRup%rupture_time(:,i) = 0.0
          DISC%DynRup%dynStress_time(:,i) = 0.0
      END DO

      allocate(disc%DynRup%output_Mu(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_Strength(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_Slip(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_Slip1(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_Slip2(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_rupture_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_PeakSR(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_dynStress_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))      
      allocate(disc%DynRup%output_StateVar(DISC%Galerkin%nBndGP,MESH%Fault%nSide))

      ! Initialize w/ first-touch
      !$omp parallel do schedule(static)
      DO i=1,MESH%fault%nSide
          disc%DynRup%output_Mu(:,i) = 0.0
          disc%DynRup%output_Strength(:,i) = 0.0
          disc%DynRup%output_Slip(:,i) = 0.0
          disc%DynRup%output_Slip1(:,i) = 0.0
          disc%DynRup%output_Slip2(:,i) = 0.0
          disc%DynRup%output_rupture_time(:,i) = 0.0
          disc%DynRup%output_PeakSR(:,i) = 0.0
          disc%DynRup%output_dynStress_time(:,i) = 0.0
          disc%DynRup%output_StateVar(:,i) = 0.0
      END DO

    else
        ! Allocate dummy arrays to avoid debug errors
        allocate(DISC%DynRup%SlipRate1(0,0), &
            DISC%DynRup%SlipRate2(0,0),      &
            DISC%DynRup%Slip(0,0),           &
            DISC%DynRup%Slip1(0,0),          &
            DISC%DynRup%Slip2(0,0),          &
            DISC%DynRup%Mu(0,0),             &
            DISC%DynRup%StateVar(0,0),       &
            DISC%DynRup%PeakSR(0,0),         &
            DISC%DynRup%Strength(0,0),       &
            DISC%DynRup%rupture_time(0,0),   &
            DISC%DynRup%dynStress_time(0,0)  )
        allocate(DISC%DynRup%output_Mu(0,0),      &
            DISC%DynRup%output_StateVar(0,0),     &
            DISC%DynRup%output_Strength(0,0),     &
            DISC%DynRup%output_Slip(0,0),         &
            DISC%DynRup%output_Slip1(0,0),        &
            DISC%DynRup%output_Slip2(0,0),        &
            DISC%DynRup%output_rupture_time(0,0), &
            DISC%DynRup%output_PeakSR(0,0),       &
            DISC%DynRup%output_dynStress_time(0,0))
    ENDIF
    !
    IF(DISC%Galerkin%CKMethod.EQ.1) THEN ! not yet done for hybrids
        print*,' ERROR in SUBROUTINE iniGalerkin3D_us_level2_new'
        PRINT*,' DISC%Galerkin%CKMethod.EQ.1 not implemented'
        call MPI_ABORT(MPI%commWorld, 134)
        !
    ENDIF
  END SUBROUTINE iniGalerkin3D_us_level2_new


 !===========================================================================!
  !!                                                                         !!
  !!  iniGalerkin3D_us_intern initializes the private data                   !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE iniGalerkin3D_us_intern_new(EQN, DISC, MESH, BND, IC, SOURCE, OptionalFields, IO, MPI)
    !-------------------------------------------------------------------------!

    USE DGBasis_mod
    USE COMMON_operators_mod
    USE QuadPoints_mod
    USE JacobiNormal_mod

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)         :: EQN
    TYPE(tDiscretization)    :: DISC
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tBoundary)          :: BND
    TYPE(tInitialCondition)  :: IC
    TYPE(tSource)            :: SOURCE
    TYPE(tUnstructOptionalFields)   :: OptionalFields
    TYPE(tInputOutput)       :: IO
    TYPE(tMPI)               :: MPI
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER :: allocstat                                  ! Allocation status !
    INTEGER :: stat                                       ! IO status         !
    INTEGER :: iDegFr                                     ! Loop counter      !
    INTEGER :: iElem, iSide, iBndGP
    INTEGER :: iNeighbor, iLocalNeighborSide, iNeighborSide, iNeighborVertex
    INTEGER :: iLocalNeighborVrtx
    INTEGER :: nDegFr, MaxDegFr, nDGWorkVar
    INTEGER :: i,j,k,l,m,r,r1,r2,r3                       ! Loop counter      !
    INTEGER :: iPoly, iXi, iEta, iZeta, iIntGP            ! Loop counter      !
    REAL    :: xi, eta, zeta, tau, chi, tau1, chi1
    REAL    :: xiS, etaS, zetaS, xP, yP, zP
    REAL    :: phi, gradphixieta(3)
    REAL                            :: rho,mu,lambda,c0,ce(6,6)
    REAL                            :: x(MESH%nVertexMax)
    REAL                            :: y(MESH%nVertexMax)
    REAL                            :: z(MESH%nVertexMax)
    REAL                            :: JacobiT(EQN%Dimension,EQN%Dimension)
    REAL                            :: auxMatrix(EQN%nVar,EQN%nVar)
    REAL                            :: A_Star(EQN%nVar,EQN%nVar)
    REAL                            :: B_Star(EQN%nVar,EQN%nVar)
    REAL                            :: C_Star(EQN%nVar,EQN%nVar)
    REAL                            :: A(EQN%nVar,EQN%nVar)
    REAL                            :: B(EQN%nVar,EQN%nVar)
    REAL                            :: C(EQN%nVar,EQN%nVar)
    REAL                            :: nx, ny, nz
    REAL                            :: sx, sy, sz
    REAL                            :: tx, ty, tz
    REAL                            :: locA(EQN%nVar,EQN%nVar), locabsA(EQN%nVar,EQN%nVar)
    REAL                            :: T(EQN%nVar,EQN%nVar)
    REAL                            :: iT(EQN%nVar,EQN%nVar)
    !
    INTEGER                         :: nTens3GP, indx
    REAL,POINTER                    :: Tens3BaseFunc(:,:)
    REAL,POINTER                    :: Tens3BaseGrad(:,:,:)
    REAL,POINTER                    :: Tens3GaussP(:,:)
    REAL,POINTER                    :: Tens3GaussW(:)
    INTEGER                         :: nTFMGaussP, ngll
    REAL,POINTER                    :: TFMGaussP(:,:)
    REAL,POINTER                    :: TFMGaussW(:)
    REAL                            :: phi_m, phi_l, phi_k
    REAL                            :: phigrad(EQN%Dimension)
    REAL                            :: chiGP, tauGP
    REAL                            :: VAND(DISC%Galerkin%nPoly+1,DISC%Galerkin%nPoly+1),Temp(DISC%Galerkin%nPoly+1,1)
    REAL                            :: grad(EQN%Dimension,EQN%Dimension)
    !
    LOGICAL                         :: configexist
    CHARACTER(LEN=600)              :: FileName_Tet, FileName_Hex, FileName_Time
    ! CHARACTER(LEN=600)              :: TimeFile
    !
    ! Dynamic Rupture variables
    INTEGER                         :: iFace
    INTEGER                         :: iDegFr2
    INTEGER                         :: MPIIndex, iObject
    REAL                            :: xGP,yGP,zGP
    REAL                            :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
    REAL                            :: x_host(MESH%GlobalVrtxType), y_host(MESH%GlobalVrtxType), z_host(MESH%GlobalVrtxType)
    REAL                            :: phi_iDegFr,phi_iDegFr2
    !-------------------------------------------------------------------------!

    IF(.NOT.DISC%Galerkin%init) THEN
       logError(*) 'iniGalerkin: SeisSol Interface not initialized!!'
       call MPI_ABORT(MPI%commWorld, 134)
    ENDIF

    IF(MESH%nElem_Tet.EQ.0 .AND. MESH%nElem_Hex.EQ.0) THEN
       logError(*) 'Quadraturefree ADER-DG is only implemented for tetrahedral and hexahedral.'
       call MPI_ABORT(MPI%commWorld, 134)
    ENDIF

    ! Reading polynomial coefficients and mass matrices

!===================================================================================!
!--------------------------------Tetrahedral Elements-------------------------------!
!===================================================================================!

    IF(MESH%nElem_Tet .GT. 0)THEN
        MaxDegFr = (DISC%Galerkin%nPolyRec+1)*(DISC%Galerkin%nPolyRec+2)*(DISC%Galerkin%nPolyRec+3)/6
        ALLOCATE(DISC%Galerkin%BndGaussP_Tet(EQN%Dimension-1,DISC%Galerkin%nBndGP),                                      &
                 DISC%Galerkin%BndGaussW_Tet(Disc%Galerkin%nBndGP),                                                    &
                 STAT = allocstat                                                                                      )
        IF(allocStat .NE. 0) THEN
           logError(*) 'could not allocate all variables!'
           call MPI_ABORT(MPI%commWorld, 134)
        END IF
    ENDIF ! Tets

  !==============================================================================!
  !------------------------Source Terms------------------------------------------!
  !==============================================================================!

   !
   ! Allocate Gausspoints for ADER-DG Time Integration of source terms
   !
#ifndef NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
   DISC%Galerkin%nTimeGP = DISC%Galerkin%nPoly+1
#else
   disc%galerkin%nTimeGp = NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
#endif
   ALLOCATE(DISC%Galerkin%TimeGaussP(DISC%Galerkin%nTimeGP),                          &
            DISC%Galerkin%TimeGaussW(DISC%Galerkin%nTimeGP),                          &
            DISC%Galerkin%dtPowerFactor(0:DISC%Galerkin%nPoly,DISC%Galerkin%nTimeGP), &
            DISC%Galerkin%Faculty(0:CONVERGENCE_ORDER+5)  )
   !
   ! Precalculate the Faculty of i because otherwise this calculation takes a loooong time...
   !
   DISC%Galerkin%Faculty(0) = 1.
   DO i = 1, CONVERGENCE_ORDER+5
       DISC%Galerkin%Faculty(i) = DISC%Galerkin%Faculty(i-1)*REAL(i)
   ENDDO

!===================================================================================!
!---------------------------Quadrature-free ADER DG---------------------------------!
!===================================================================================!

    DISC%Galerkin%nRK = 1

    ! Attention: Don't change Nr of GP here since some routine depend on these numbers
    DISC%Galerkin%nIntGP = (DISC%Galerkin%nPoly + 2)**3

    SELECT CASE(DISC%Galerkin%DGMethod)
    CASE(3)
        nDGWorkVar = EQN%nVar+EQN%nAneFuncperMech
    CASE DEFAULT
        nDGWorkVar = EQN%nVarTotal
    END SELECT

    IF(MESH%nElem_Tet.GT.0) THEN

#ifdef USE_DR_CELLAVERAGE
        call CellCentresOfSubdivision(DISC%Galerkin%nPoly + 1, DISC%Galerkin%BndGaussP_Tet)
        DISC%Galerkin%BndGaussW_Tet = 1.e99 ! blow up solution if used
#else
        ! Compute and store surface gaussian integration points
        CALL TriangleQuadraturePoints(                         &
                 nIntGP     = DISC%Galerkin%nBndGP,            &
                 IntGaussP  = DISC%Galerkin%BndGaussP_Tet,     &
                 IntGaussW  = DISC%Galerkin%BndGaussW_Tet,     &
                 M          = DISC%Galerkin%nPoly+2,           &
                 IO         = IO,                              &
                 quiet      = .TRUE.,                          &
                 MPI        = MPI                              )
#endif

#ifndef NDEBUG
        ! assert contant material parameters per element
        if ( disc%galerkin%nDegFrMat .ne. 1 ) then
          logError(*) 'iniGalerkin3D_us_intern_new, disc%galerkin%nDegFrMat not equal 1.', disc%galerkin%nDegFrMat
          call MPI_ABORT(MPI%commWorld, 134)
        endif

        ! assert 4 sides for tetrahedrons
        if ( mesh%nSides_tet .ne. 4 ) then
          logError(*) 'iniGalerkin3D_us_intern_new, mesh%nSides_tet not equal 4.', mesh%nSides_tet
          call MPI_ABORT(MPI%commWorld, 134)
        endif

        ! assert 3 vertices for triangles
         if ( mesh%nVertices_tri .ne. 3 ) then
          logError(*) 'iniGalerkin3D_us_intern_new, mesh%nVertices_tri not equal 3.', mesh%nVertices_tri
          call MPI_ABORT(MPI%commWorld, 134)
        endif
#endif
    ENDIF ! Tetras
    !
    logInfo(*) 'iniGalerkin successful '
    !
  END SUBROUTINE iniGalerkin3D_us_intern_new


  !===========================================================================!
  !!                                                                         !!
  !!  icGalerkin3D_us initializes the degrees of freedom at time t=0.0       !!
  !!  by L2 projection                                                       !!
  !!                                                                         !!
  !===========================================================================!


  SUBROUTINE icGalerkin3D_us_new(EQN, DISC, MESH, IC, SOURCE, IO, MPI)
    !-------------------------------------------------------------------------!
    use iso_c_binding, only: c_loc
    use f_ftoc_bind_interoperability
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)         :: EQN
    TYPE(tDiscretization)    :: DISC
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tInitialCondition)  :: IC
    TYPE(tSource)            :: SOURCE
    TYPE(tInputOutput)       :: IO
    TYPE(tMPI)               :: MPI
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER :: iElem                                                          ! Element number
    INTEGER :: iIntGP                                                         ! Index of internal Gausspoint
    INTEGER :: iDegFr                                                         ! Degree of freedom
    INTEGER :: iVar                                                           ! Variable number
    INTEGER :: iVert                                                          ! Vertex counter
    INTEGER :: iPoly, nIntGP, nDegFr, eType
    INTEGER :: LocPoly, LocDegFr                                              ! Variables for p-adaptivity
    REAL    :: xi, eta, zeta                                                  ! Reference coordinates
    REAL    :: xGP, yGP, zGP                                                  ! Physical coordinates
    REAL    :: x(MESH%nVertexMax)                                             ! Element vertices in physical coordinates system
    REAL    :: y(MESH%nVertexMax)                                             ! Element vertices in physical coordinates system
    REAL    :: z(MESH%nVertexMax)                                             ! Element vertices in physical coordinates system
    REAL    :: iniGP_Plast(6)                                                 ! Initial stress loading for plastic calculations
    REAL    :: phi                                                            ! Value of the base function at GP      !
    REAL    :: Einv, v                                                        ! Inverse of Young's modulus, Poisson ratio v
    !
    ! temporary degrees of freedom
    real    :: l_initialLoading( NUMBER_OF_BASIS_FUNCTIONS, 6 )
    REAL    :: oneRankedShaped_iniloading(NUMBER_OF_BASIS_FUNCTIONS*6)        ! l_iniloading to one rank array  (allows removing warning we running with plasticity))
    real    :: l_plasticParameters(2)
    !-------------------------------------------------------------------------!
    !
    IF(.NOT.DISC%Galerkin%init) THEN
       logError(*) 'icGalerkin: SeisSol Interface not initialized!!'
       call MPI_ABORT(MPI%commWorld, 134)
    ENDIF
    !
    ALLOCATE(EQN%Energy(3,1:MESH%nElem))
    EQN%Energy = 0.


    IF(EQN%Plasticity.EQ.1) THEN
      ALLOCATE(DISC%Galerkin%DOFStress(DISC%Galerkin%nDegFr,6,MESH%nElem), DISC%Galerkin%pstrain(7, MESH%nElem),&
               DISC%Galerkin%PlasticParameters(4,1:MESH%nElem), DISC%Galerkin%Strain_Matrix(6,6))
        !Initialization
        DISC%Galerkin%DOFStress = 0.
        DISC%Galerkin%pstrain = 0.
        DISC%Galerkin%PlasticParameters = 0.
        DISC%Galerkin%Strain_Matrix = 0.

        !Initialize the stress-strain relation matrix (mu and lambda should be element dependent)
        Einv = (EQN%lambda+EQN%mu)/(EQN%mu*(3*EQN%lambda+2*EQN%mu))!Inv of the Young's modulus
        v = EQN%lambda/(2*(EQN%lambda+EQN%mu)) !Poisson's ratio

        DISC%Galerkin%Strain_Matrix(1,1) = Einv
        DISC%Galerkin%Strain_Matrix(2,2) = Einv
        DISC%Galerkin%Strain_Matrix(3,3) = Einv
        DISC%Galerkin%Strain_Matrix(4,4) = 1/(2*EQN%mu)
        DISC%Galerkin%Strain_Matrix(5,5) = 1/(2*EQN%mu)
        DISC%Galerkin%Strain_Matrix(6,6) = 1/(2*EQN%mu)
        DISC%Galerkin%Strain_Matrix(1,2) = -v*Einv
        DISC%Galerkin%Strain_Matrix(1,3) = -v*Einv
        DISC%Galerkin%Strain_Matrix(2,1) = -v*Einv
        DISC%Galerkin%Strain_Matrix(2,3) = -v*Einv
        DISC%Galerkin%Strain_Matrix(3,1) = -v*Einv
        DISC%Galerkin%Strain_Matrix(3,2) = -v*Einv
    ENDIF


    logInfo0(*) 'DG initial condition projection... '

    call c_interoperability_projectInitialField()
    !
    iPoly  = DISC%Galerkin%nPoly
    nIntGP = DISC%Galerkin%nIntGP
    nDegFr = DISC%Galerkin%nDegFr

    !$omp parallel do schedule(static) shared(eqn, disc, mesh, ic, source, io, iPoly, nIntGp, nDegFr) private(iElem, iIntGP, iDegFr, iVar, iVert, eType, locPoly, locDegFr, xi, eta, zeta, xGp, yGp, zGp, x, y, z, phi, l_initialLoading,oneRankedShaped_iniloading, l_plasticParameters, iniGP_plast) default(none)
    DO iElem = 1,MESH%nElem
        l_initialLoading=0
        l_plasticParameters=0
        
        IF(EQN%Plasticity == 1) THEN !high-order points approach
           !elementwise assignement of the initial loading
           l_initialLoading(1,1:6) = EQN%IniStress(1:6,iElem)

           ! initialize the element dependent plastic parameters
           l_plasticParameters(1) = EQN%PlastCo(iElem) !element-dependent plastic cohesion
           l_plasticParameters(2) = EQN%BulkFriction(iElem) !element-dependent bulk friction
        
           ! initialize loading in C
           oneRankedShaped_iniloading = pack( l_initialLoading, .true. )
           call c_interoperability_setInitialLoading( i_meshId = iElem, \
                                                      i_initialLoading = oneRankedShaped_iniloading)

           !initialize parameters in C
           call c_interoperability_setPlasticParameters( i_meshId            = iElem, \
                                                         i_plasticParameters = l_plasticParameters )

       END IF
    ENDDO ! iElem

    IF (EQN%Plasticity == 1) THEN
      call c_interoperability_setTv( tv = EQN%Tv )
      ! TODO: redundant (see iniGalerkin3D_us_level2_new) call to ensure correct intitial loading in copy layers.
      call c_interoperability_synchronizeCellLocalData(logical(.true., 1));

  END IF

    logInfo0(*) 'DG initial condition projection done. '

  END SUBROUTINE icGalerkin3D_us_new

  SUBROUTINE BuildSpecialDGGeometry3D_new(MaterialVal,EQN,MESH,DISC,BND,MPI,IO)
    USE iso_c_binding, only: c_loc, c_null_char, c_bool
    USE common_operators_mod
    USE DGbasis_mod
    USE ini_faultoutput_mod
    USE f_ftoc_bind_interoperability
#ifdef HDF
    USE hdf_faultoutput_mod
#endif
    use, intrinsic :: iso_c_binding

    !-------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------!
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    ! Argument list declaration
    TYPE(tEquations)         :: EQN
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tDiscretization)    :: DISC
    TYPE(tBoundary)          :: BND
    TYPE(tMPI)               :: MPI
    TYPE(tInputOutput)       :: IO
    REAL                     :: MaterialVal(MESH%nElem,EQN%nBackgroundVar)          !
    !-------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER :: iElem, iSide, iDomain            ! Loop counter                !
    REAL    :: sidevec(3,2)                     ! Tmp vertex connection vector!
    INTEGER :: VertexSide_Tet(MESH%nSides_Tet,MESH%nVertices_Tri)
    INTEGER :: VertexSide_Hex(MESH%nSides_Hex,MESH%nVertices_Quad)              ! # sides = 6, # vertices per side = 4
    INTEGER :: allocstat                        ! Status of allocation        !
    REAL    :: Length                           ! Length of tangent vector    !
    REAL    :: BaryVec(EQN%Dimension,MESH%nSideMax), Dist(MESH%nSideMax)
    REAL    :: minv
    INTEGER :: minl(1)
    INTEGER :: i, j, iLayer, iZone
    REAL    :: nx, ny, nz, sx, sy, sz, tx, ty, tz, rho, amax, a
    REAL    :: c(6,6), Voigt_rot(6,6), TT(6,6), T(6,6)
    REAL    :: coefficients(4), coefficients2(5)
    REAL    :: Re_solution(3), Im_solution(3)
    REAL    :: K_F, K_S, K_Mean, MM, Poro, Alpha(6)  !Porous parameters to obtain undrained c_ij parameters
    REAL    :: rho_F, rho_S, nu, Kappa(3), Tor(3)
    REAL    :: Rho1, Rho2, Beta1, Beta2, solution2(4)
    REAL, POINTER :: zone_minh(:), zone_maxh(:), zone_deltah(:), zone_deltap(:)
    COMPLEX :: solution(3)
    INTEGER :: nDOF,TotDOF, PoroFlux
    REAL    :: elementWaveSpeeds(4)
    !
    INTEGER :: iErr,iPoly,iVrtx
    INTEGER :: nLocPolyElem(0:100), TempInt(MESH%nSideMax)
    REAL    :: min_h, max_h, deltah, deltap
    INTEGER :: LocPoly
    REAL    :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType), Tmp(3,3),D(2,2,2,3)
    REAL    :: xi, eta, zeta, xP, yP, zP
    INTEGER :: ngll, k, Fix1(2), Fix2(2)
    CHARACTER(LEN=200) :: Filename
    !
    ! variables for fault output
    REAL, POINTER :: S_inc(:)
    REAL, ALLOCATABLE :: chi_vector(:), tau_vector(:)
    REAL    :: S_tmp, chi, tau, phi1, phi2
    integer :: hasDR
    INTEGER :: iNeighbor, iNeighborSide, NeigBndGP, l, iFault, iDegFr, iP, iBndGP, iPlusElem
    !
    !-------------------------------------------------------------------------!
    INTENT(IN)                :: MaterialVal, EQN, IO
    INTENT(INOUT)             :: DISC, BND, MESH
    !-------------------------------------------------------------------------!
    !                                                                         !
    !
    ! The unit tetrahedron has the following 4 local vertices:                !
    !                                                                         !
    ! 1 = (0,0,0)                                                             !
    ! 2 = (1,0,0)                                                             !
    ! 3 = (0,1,0)                                                             !
    ! 4 = (0,0,1)                                                             !
    !                                                                         !
    VertexSide_Tet(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of side I        !
    VertexSide_Tet(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of side II       !
    VertexSide_Tet(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of side III      !
    VertexSide_Tet(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of side IV       !
    ! The unit hexahedron has the following 6 local vertices:                 !
    !                                                                         !
    ! 1 = (0,0,0)                                                             !
    ! 2 = (1,0,0)                                                             !
    ! 3 = (0,1,0)                                                             !
    ! 4 = (1,1,0)                                                             !
    ! 5 = (0,0,1)                                                             !
    ! 6 = (1,0,1)                                                             !
    ! 7 = (0,1,1)                                                             !
    ! 8 = (1,1,1)                                                             !
    !                                                                         !
    !         7           8         1          1,2,6,5
    !         *-----------*         2          2,4,8,6
    !        /|          /|         3          4,3,7,8
    !       / |         / |         4          3,1,5,7
    !     5*-----------*6 |         5          2,1,3,4
    !      |  |        |  |         6          5,6,8,7
    !      |  |        |  |
    !      | 3*--------|--*4        Each face is characterized by the sum of it local vertex numbers,
    !      | /         | /          i.e. side1=14, side2=20, side3=22, side4=16, side5=10, side6=26
    !      |/          |/
    !     1*-----------*2
    !
    VertexSide_Hex(1,:) =  (/ 1,2,6,5 /)   ! Local hex. vertices of side I        !
    VertexSide_Hex(2,:) =  (/ 2,4,8,6 /)   ! Local hex. vertices of side II       !
    VertexSide_Hex(3,:) =  (/ 4,3,7,8 /)   ! Local hex. vertices of side III      !
    VertexSide_Hex(4,:) =  (/ 3,1,5,7 /)   ! Local hex. vertices of side IV       !
    VertexSide_Hex(5,:) =  (/ 2,1,3,4 /)   ! Local hex. vertices of side V        !
    VertexSide_Hex(6,:) =  (/ 5,6,8,7 /)   ! Local hex. vertices of side VI       !
    !
    ALLOCATE( DISC%Galerkin%geoNormals(EQN%Dimension,MESH%nSideMax,MESH%nElem),  &
              DISC%Galerkin%geoTangent1(EQN%Dimension,MESH%nSideMax,MESH%nElem),  &
              DISC%Galerkin%geoTangent2(EQN%Dimension,MESH%nSideMax,MESH%nElem),  &
              DISC%Galerkin%geoSurfaces(MESH%nSideMax,MESH%nElem),  &
              STAT=allocstat )
    IF (allocStat .NE. 0) THEN
       logError(*) 'Interface SeisSol: could not allocate all variables!'
       call MPI_ABORT(MPI%commWorld, 134)
    END IF

    ! Calculating boundary surfaces (3D)

    DO iElem=1,MESH%nElem
       DO iSide=1,MESH%LocalElemType(iElem)
          SELECT CASE(MESH%LocalElemType(iElem))
          CASE(4)
              ! Boundary side vector pointing in chi-direction
              sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,2),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem))
              ! Boundary side vector pointing in tau-direction
              sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,3),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem))
              ! Normal vector computed by cross product
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  sidevec(:,1).x.sidevec(:,2)
              ! Triangle surface = 0.5 * cross_product
              DISC%Galerkin%geoSurfaces(iSide,iElem)  =  0.5*SQRT(                                     &
                                                         DISC%Galerkin%geoNormals(1,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(2,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(3,iSide,iElem)**2    )
              ! Normalize normal vector to length 1
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  0.5*DISC%Galerkin%geoNormals(:,iSide,iElem) / &
                                                         DISC%Galerkin%geoSurfaces(iSide,iElem)
              !
              ! Compute MinDistBarySide :
              ! 1. Compute vector connecting barycenter of tetrahedron and local point 1 of the side
              BaryVec(:,iSide) = MESH%ELEM%xyBary(:,iElem) - MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem))
              ! 2. Compute scalar product of previous vector and normal vector (already normalized to 1)
              !    and take the absolute value.
              Dist(iSide) = ABS(DOT_PRODUCT(BaryVec(:,iSide),DISC%Galerkin%geoNormals(:,iSide,iElem)))
          CASE(6)
              ! Boundary side vector pointing in chi-direction
              sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,2),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,1),iElem))
              ! Boundary side vector pointing in tau-direction
              sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,3),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,1),iElem))
              ! Normal vector computed by cross product
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  sidevec(:,1).x.sidevec(:,2)
              ! Triangle's surface = 0.5 * cross_product
              DISC%Galerkin%geoSurfaces(iSide,iElem)  =  0.5*SQRT(                                     &
                                                         DISC%Galerkin%geoNormals(1,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(2,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(3,iSide,iElem)**2    )

              ! Boundary side vector pointing in chi-direction
              sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,4),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,3),iElem))
              ! Boundary side vector pointing in tau-direction
              sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,2),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,3),iElem))
              ! Normal vector computed by cross product
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  sidevec(:,1).x.sidevec(:,2)
              ! Second triangle's surface = 0.5 * cross_product
              DISC%Galerkin%geoSurfaces(iSide,iElem)  =  DISC%Galerkin%geoSurfaces(iSide,iElem) +      &
                                                         0.5*SQRT(                                     &
                                                         DISC%Galerkin%geoNormals(1,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(2,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(3,iSide,iElem)**2    )

              ! Normalize normal vector to length 1
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  DISC%Galerkin%geoNormals(:,iSide,iElem) / &
                                                         SQRT(SUM(DISC%Galerkin%geoNormals(:,iSide,iElem)**2))
              !
              ! Compute MinDistBarySide :
              ! 1. Compute vector connecting barycenter of tetrahedron and local point 1 of the side
              BaryVec(:,iSide) = MESH%ELEM%xyBary(:,iElem) - MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,1),iElem))
              ! 2. Compute scalar product of previous vector and normal vector (already normalized to 1)
              !    and take the absolute value.
              Dist(iSide) = ABS(DOT_PRODUCT(BaryVec(:,iSide),DISC%Galerkin%geoNormals(:,iSide,iElem)))

              !
              CONTINUE
              !
          END SELECT
          ! Compute vector inside the triangle's plane for the rotation matrix
          DISC%Galerkin%geoTangent1(:,iSide,iElem) = sidevec(:,1)
          ! Normalize to 1
          Length = SQRT(DISC%Galerkin%geoTangent1(1,iSide,iElem)**2 + &
                        DISC%Galerkin%geoTangent1(2,iSide,iElem)**2 + &
                        DISC%Galerkin%geoTangent1(3,iSide,iElem)**2   )
          DISC%Galerkin%geoTangent1(:,iSide,iElem) = DISC%Galerkin%geoTangent1(:,iSide,iElem)/Length
          ! Compute second vector in the plane, orthogonal to the normal and tangent 1 vectors
          ! using the crossproduct
          DISC%Galerkin%geoTangent2(:,iSide,iElem) = DISC%Galerkin%geoNormals( :,iSide,iElem) .x. &
                                                     DISC%Galerkin%geoTangent1(:,iSide,iElem)
          ! Check: Length must be 1.
          Length = DISC%Galerkin%geoTangent2(1,iSide,iElem)**2 + &
                   DISC%Galerkin%geoTangent2(2,iSide,iElem)**2 + &
                   DISC%Galerkin%geoTangent2(3,iSide,iElem)**2
      ENDDO
      !
      SELECT CASE(MESH%LocalElemType(iElem))
      CASE(6)
          ! Minimum distance of Barycenter to a face
          MESH%ELEM%MinDistBarySide(iElem) = MINVAL(Dist(:))
      CASE(4)
          ! Insphere radius of tetrahedron (is larger the the previous minimum distance !!!)
          MESH%ELEM%MinDistBarySide(iElem) = 3*MESH%ELEM%Volume(iElem)/SUM(DISC%Galerkin%geoSurfaces(1:4,iElem))
      END SELECT
      !
    ENDDO
    !
    ! Report mesh quality
    !
    minl = MINLOC(MESH%ELEM%Volume(:))
    minv = MINVAL(MESH%ELEM%Volume(:))

    logInfo(*) 'Smallest volume found in tetraedron number : ', minl(1)
    logInfo(*) 'Smallest volume is                         : ', minv

    minl = MINLOC(MESH%ELEM%MinDistBarySide(:))
    minv = MINVAL(MESH%ELEM%MinDistBarySide(:))

    logInfo(*) 'Smallest insphere found in tetraedron number : ', minl(1)
    logInfo(*) 'Smallest insphere radius is                  : ', minv

    IF(minv.LE.1e-15) THEN
        logError(*) 'Mesh contains a singular tetrahedron with radius ', minv
        logError(*) 'Element number and position : ', minl(1), MESH%ELEM%xyBary(:,minl(1))
        call MPI_ABORT(MPI%commWorld, 134)
    ENDIF
    DISC%DynRup%DynRup_out_elementwise%DR_pick_output = .FALSE.
    DISC%DynRup%DynRup_out_elementwise%nDR_pick       = 0
    !
    !
    !
    ! Initialize fault rupture output
    ! only in case Dynamic rupture is turned on, and for + elements assigned to the fault
    IF(EQN%DR.EQ.1 .AND. DISC%DynRup%DR_output) THEN
        ! Case 3
        ! output at certain positions specified in the *.dyn file
        IF(DISC%DynRup%OutputPointType.EQ.3) THEN
            !
            DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output = .TRUE.
            DISC%DynRup%DynRup_out_atPickpoint%nDR_pick       = DISC%DynRup%DynRup_out_atPickpoint%nOutPoints
            !
            ! test if fault pickpoints are on the fault (within a tolerance) and find corresponding "+"-element (iElem)

#ifdef HDF
            CALL ini_fault_receiver_hdf(EQN, MESH, DISC, IO, MPI)
!#else
#endif
            CALL ini_fault_receiver(EQN,MESH,BND,DISC,IO,MPI)

        ! Case 4
        ! for full fault output without pickpoints
        ELSEIF(DISC%DynRup%OutputPointType.EQ.4) THEN
            !
            DISC%DynRup%DynRup_out_elementwise%DR_pick_output = .TRUE.
            DISC%DynRup%DynRup_out_elementwise%nDR_pick       = 0
            CALL ini_fault_subsampled(EQN,MESH,BND,DISC,IO,MPI)
        ! Case 5
        ! for full fault output and pickpoints
        ELSEIF(DISC%DynRup%OutputPointType.EQ.5) THEN
            !
            DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output = .TRUE.
            DISC%DynRup%DynRup_out_atPickpoint%nDR_pick       = DISC%DynRup%DynRup_out_atPickpoint%nOutPoints
            !
            ! test if fault pickpoints are on the fault (within a tolerance) and find corresponding "+"-element (iElem)
            CALL ini_fault_receiver(EQN,MESH,BND,DISC,IO,MPI)
            !
            !
            DISC%DynRup%DynRup_out_elementwise%DR_pick_output = .TRUE.
            DISC%DynRup%DynRup_out_elementwise%nDR_pick       = 0
            CALL ini_fault_subsampled(EQN,MESH,BND,DISC,IO,MPI)
        ENDIF ! DISC%DynRup%OutputPointType
    ENDIF ! end initialize fault output
    !
    !
    !
    !
    ! Allocate rest of MPI communication structure
    logInfo(*) 'Allocation of remaining MPI communication structure '
    logInfo(*) '  General info: ', BND%NoMPIDomains,DISC%Galerkin%nDegFr,EQN%nVar
    DO iDomain = 1, BND%NoMPIDomains
        logInfo(*) 'Bnd elements for domain ', iDomain, ' : ',  BND%ObjMPI(iDomain)%nElem
        IF(DISC%Galerkin%DGMethod.EQ.3) THEN
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborDOF(DISC%Galerkin%nDegFrST,EQN%nVarTotal,BND%ObjMPI(iDomain)%nElem) )
        ELSE
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborDOF(DISC%Galerkin%nDegFrRec,EQN%nVarTotal,BND%ObjMPI(iDomain)%nElem) )
        ENDIF
        ALLOCATE( BND%ObjMPI(iDomain)%NeighborBackground(EQN%nBackgroundVar,BND%ObjMPI(iDomain)%nElem)     )
        BND%ObjMPI(iDomain)%Init = .FALSE.
        IF(DISC%Galerkin%DGMethod.EQ.3) THEN
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborDuDt(DISC%Galerkin%nDegFr,EQN%nVar+EQN%nAneFuncperMech,BND%ObjMPI(iDomain)%nElem) )
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborTime( BND%ObjMPI(iDomain)%nElem) )
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborDt(   BND%ObjMPI(iDomain)%nElem) )
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborUpdate(BND%ObjMPI(iDomain)%nElem))
            BND%ObjMPI(iDomain)%NeighborDuDt(:,:,:) = 0.
            BND%ObjMPI(iDomain)%NeighborTime(:)     = -1e10
            BND%ObjMPI(iDomain)%NeighborDt(:)       = -2e10
            BND%ObjMPI(iDomain)%NeighborUpdate(:)   = -1
        ENDIF
    ENDDO ! iDomain
    !
    IF(DISC%Galerkin%ZoneOrderFlag.EQ.1) THEN
        ALLOCATE( zone_minh(MESH%nZones), zone_maxh(MESH%nZones), zone_deltah(MESH%nZones), zone_deltap(MESH%nZones) )
        DO iZone = 1, MESH%nZones
            min_h = MINVAL( MESH%ELEM%Volume(:), MASK = MESH%ELEM%Reference(0,:).EQ.iZone )
            max_h = MAXVAL( MESH%ELEM%Volume(:), MASK = MESH%ELEM%Reference(0,:).EQ.iZone )
#ifdef PARALLEL
            CALL MPI_ALLREDUCE(min_h,zone_minh(iZone),1,MPI%MPI_AUTO_REAL,MPI_MIN,MPI%commWorld,iErr)
            CALL MPI_ALLREDUCE(max_h,zone_maxh(iZone),1,MPI%MPI_AUTO_REAL,MPI_MAX,MPI%commWorld,iErr)
#else
            zone_minh(iZone) = min_h
            zone_maxh(iZone) = max_h
#endif
            zone_deltah(iZone) = zone_maxh(iZone) - zone_minh(iZone)
            zone_deltap(iZone) = DISC%Galerkin%ZoneMaxPoly(iZone) - DISC%Galerkin%ZoneMinPoly(iZone)
        ENDDO
    ELSE
        min_h = MINVAL( MESH%ELEM%MinDistBarySide(:) )
        max_h = MAXVAL( MESH%ELEM%MinDistBarySide(:) )
#ifdef PARALLEL
        CALL MPI_ALLREDUCE(min_h,MESH%min_h,1,MPI%MPI_AUTO_REAL,MPI_MIN,MPI%commWorld,iErr)
        CALL MPI_ALLREDUCE(max_h,MESH%max_h,1,MPI%MPI_AUTO_REAL,MPI_MAX,MPI%commWorld,iErr)
#else
        MESH%min_h = min_h
        MESH%max_h = max_h
#endif
        deltah = MESH%max_h - MESH%min_h
        deltap = DISC%Galerkin%nPoly - DISC%Galerkin%nMinPoly
    ENDIF


    ! Allocate objects for local p-adaptivity (not yet done for hybrid meshes)
    IF(DISC%Galerkin%pAdaptivity.GE.1) THEN
        ALLOCATE( DISC%Galerkin%LocPoly(MESH%nElem) )
        !
        nLocPolyElem(:) = 0
        ! Distribute local polynomial degree
        !
        DO iElem = 1, MESH%nElem
            IF(DISC%Galerkin%ZoneOrderFlag.EQ.1)THEN
                !
                !  Distribute local polynomial degree locally by zone
                !  ----------------------------------------------------
                iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
                IF(DISC%Galerkin%ZoneMinPoly(iLayer).GT.-1)THEN
                     SELECT CASE(DISC%Galerkin%ZoneOrder(iLayer))
                     CASE(1)
                         LocPoly = NINT( DISC%Galerkin%ZoneMinPoly(iLayer) + zone_deltap(iLayer)* &
                                         ((MESH%ELEM%Volume(iElem)-zone_minh(iLayer))/zone_deltah(iLayer))**DISC%Galerkin%ZonePower(iLayer) )
                     CASE(-1)
                         LocPoly = NINT( DISC%Galerkin%ZoneMaxPoly(iLayer) - zone_deltap(iLayer)* &
                                         ((MESH%ELEM%Volume(iElem)-zone_minh(iLayer))/zone_deltah(iLayer))**DISC%Galerkin%ZonePower(iLayer) )
                     END SELECT
                ELSE
                   PRINT *, ' ERROR: local order must not be less or equal to zero! ', iLayer
                   call MPI_ABORT(MPI%commWorld, 134)
                ENDIF
            ELSE
                !
                !  Distribute local polynomial degree globally by size
                !  ----------------------------------------------------
                LocPoly = NINT( DISC%Galerkin%nMinPoly + deltap*((MESH%ELEM%MinDistBarySide(iElem)-MESH%min_h)/deltah)**(1) )
            ENDIF
            !
            DISC%Galerkin%LocPoly(iElem) = LocPoly
            nLocPolyElem( LocPoly ) = nLocPolyElem( LocPoly ) + 1
            !
        ENDDO

        TotDOF = 0
        DO iPoly = DISC%Galerkin%nMinPoly, DISC%Galerkin%nPoly
            logInfo('(a,i3,a,i6)') 'Nr of elements with degree ', iPoly, ' : ', nLocPolyElem(iPoly)
            nDOF = nLocPolyElem(iPoly)*(iPoly+1)*(iPoly+2)*(iPoly+3)/6
            logInfo(*) 'Nr of DOF for this p-zone:   ', nDOF
            TotDOF = TotDOF + nDOF
            IF(nLocPolyElem(iPoly).GT.0) THEN
                WRITE(FileName,'(a10,i2.2,a4)') 'MeshZone-p',iPoly,'.dat'
                OPEN(UNIT=999,FILE=TRIM(FileName),RECL=300)
                WRITE(999,*) ' TITLE = "Mesh Subzone for p-Adaptivity" '
                WRITE(999,*) ' VARIABLES = "x" "y" "z" "N" '
                WRITE(999,*) ' ZONE N=  ', nLocPolyElem(iPoly)*4, &
                                  ' E=  ', nLocPolyElem(iPoly),   &
                                  ' F=FEPOINT  ET=TETRAHEDRON '
                DO iElem = 1, MESH%nElem
                   IF(DISC%Galerkin%LocPoly(iElem).EQ.iPoly) THEN
                      DO iVrtx = 1, MESH%GlobalElemType
                          WRITE(999,*) MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(iVrtx,iElem)), REAL(iPoly)
                      ENDDO
                   ENDIF
                ENDDO
                i = 0
                DO iElem = 1, MESH%nElem
                   IF(DISC%Galerkin%LocPoly(iElem).EQ.iPoly) THEN
                      DO iVrtx = 1, MESH%GlobalElemType
                        i = i + 1
                        TempInt(iVrtx) = i
                      ENDDO
                      WRITE(999,*) TempInt(:)
                   ENDIF
                ENDDO
                CLOSE(999)
            ENDIF
        ENDDO
        logInfo(*) '======================================================= '
        logInfo(*) 'Total number of degrees of freedom in domain:   ', TotDOF
        logInfo(*) '======================================================= '
    ELSE
        ALLOCATE( DISC%Galerkin%LocPoly(1) )
    ENDIF

    IF(DISC%Galerkin%ZoneOrderFlag.EQ.1) THEN
        DEALLOCATE( zone_minh, zone_maxh, zone_deltah, zone_deltap )
    ENDIF

    CONTINUE
    !
  END SUBROUTINE BuildSpecialDGGeometry3D_new

END MODULE dg_setup_mod
