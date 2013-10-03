!=============================================================================
!
! $Id: smv2chem1.h,v 1.1.1.1 2008-02-12 16:06:36 trayanov Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Mark Jacobson, UCLA/Stanford)
!   jrt@llnl.gov
!
! FILE
!   Smv2Chem1_mod.F90
!
! DESCRIPTION
!   This include file contains the data for SmvgearII for the
!   variables that need to be saved between calls.
!
!   Note that during analysis for adding openmp to the Physproc routine,
!   smv2chem.h was split into smv2chem1.h & smv2chem2.h, to make the handling
!   of data cleaner and clearer.
!=============================================================================

module Smv2Chem1_mod

   implicit none
   private

   public :: Smv2Chem1_type

   type Smv2Chem1_type

      ! if 1, then reorder grid-cells by stiffness
      integer :: reorderGridCellsStiffness
      ! identifies spc # of water vapor
      integer :: speciesNumOfWaterVapor
      integer :: imgas
      ! identifies spc # of nitrogen gas
      integer :: speciesNumOfNitrogen
      ! identifies spc # of oxygen   gas
      integer :: speciesNumOfOxygen
      ! intended # of grid-cells in a grid-block
      integer :: intendedNumGridCellsInBlock
      ! logical unit number to write to when pr_smv2 is true
      integer :: unitNumberPrSmv2
      ! identifies gas chemistry type (1..NCSGAS)
      integer :: gasChemistryType
      ! jphotrat  : tbd
      integer, allocatable :: jphotrat (:) !(ICS)
      ! # of kinetic rxns (non-photo)
      integer, allocatable :: numKineticRxns   (:) !(ICS)
      ! ntloopncs : tbd
      integer, allocatable :: ntloopncs(:) !(ICS)
      ! # of active + inactive gases
      integer, allocatable :: numActiveAndInactiveGases   (:) !(ICS)
      ! original spc # of each new jnew spc
      integer, allocatable :: inewold  (:,:) !(MXGSAER, ICS)
      !  npphotrat : tbd
      integer, allocatable :: npphotrat (:,:) !(IPHOT,   ICS)
      ! fraction time step is decreased in Smvgear if convergence
      real*8  :: fractionDecreaseConvergence
      ! max time step for night all chem (s)
      real*8  :: maxTimeStepNight
   end type Smv2Chem1_type

      contains

end module Smv2Chem1_mod

