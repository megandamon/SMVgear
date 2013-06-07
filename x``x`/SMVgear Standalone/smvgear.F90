!=============================================================================
!
! $Id: smvgear.F90,v 1.1.1.1 2008-02-12 16:06:36 trayanov Exp $
!
! CODE DEVELOPER
!   Original code from Mark Z. Jacobson ((C) COPYRIGHT, 1993).
!   LLNL modifications:  John Tannahill
!                        jrt@llnl.gov
!
! FILE
!   smvgear.F
!
! ROUTINES
!   Smvgear
!   doMassBalanceAccounting
!   initializeFirstTimeStep
!   doSubStep
!   correctorStep
!   advanceTimeStep
!   initConcentrationArray
!   storeInitConcAndDerivatives
!   computeErrorFromCorrected1stDeriv
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Smvgear
!
! DESCRIPTION
!   This routine is the driver for the Smvgear (Sparse Matrix Vector Gear code)
!   chemistry solver.  It uses a Gear-type integrator that solves first order
!   ordinary differential equations with initial value boundary conditions.
!   Smvgear differs from an original Gear code in that it uses sparse matrix
!   and vectorization techniques to improve its computational speed.
!
!   This version is Smvgear II, 9/96.  It has been modified to include
!   grid-cell reordering prior to each time interval and different chemistry
!   for different atmospheric regions.  The purpose of the reordering is to
!   equations together.  This reordering can save signifcant computer time
!   (e.g., speed the code by a factor of two or more), depending on the
!   variation in stiffness throughout the grid-domain.  When the stiffness is
!   the same throughout the grid-domain (e.g., if all concentrations and rates
!   are the same), then reordering is unnecessary and will not speed solutions.
!
!   This version includes a variable absolute error tolerance.  The absolute
!   tolerance is recalculated every few Gear time steps.  This version also
!   contains different sets of chemistry for different regions of the
!   atmosphere.  Thus, urban, free tropospheric, and stratospheric chemistry
!   can be solved during the same model run.
!
!   References =>
!   ----------
!
!     Jacobson M. Z. (1997) Improvement in Smvgear II through Absolute
!     Error Tolerance Control; in submission.
!
!     Jacobson M. Z. (1995) Computation of Global Photochemistry with Smvgear
!     II, Atmos. Environ., 29a, 2541-2546.
!
!     Jacobson M. Z. (1994) Developing, Coupling, and Applying a Gas, Aerosol,
!     Transport, and Radiation Model to Studying Urban and Regional Air
!     Pollution, PhD thesis, University of California, Los Angeles.
!
!     Jacobson M. Z. and Turco R. P. (1994) Smvgear: A Sparse Matrix,
!     Vectorized Gear Code for Atmospheric Models, Atmos. Environ. 28a,
!     273-284.
!
!     The origins of the Gear integrator used in Smvgear are found in:
!       Gear C. W. (1971) Numerical Initial Value Problems in Ordinary
!       Differential Equations, Prentice-Hall, NJ, pp. 158-166.
!
!     Finally, in subroutine Smvgear, the following ideas originated from
!     Lsodes, the Livermore solver for ordinary differential with sparse
!     matrices (Hindmarsh A. C. and Sherman A. H.):
!       (a) predicting the first time-step;
!       (b) determining corrector convergence differently than in Gear's
!           original code (goc);
!       (c) determining error differently than in goc;
!       (d) summing up the pascal matrix differently than in goc.
!
!     References for the 1987 Lsodes version include:
!
!       Sherman A. H. and Hindmarsh A. C. (1980) Gears: A Package for the
!       Solution of Sparse, Stiff Ordinary Differential Equations, Lawrence
!       Livermore National Laboratory Report UCRL-84102.
!
!       Hindmarsh A. C. (1983) Odepack, A Systematized Collection of ODE
!       Solvers, Scientific Computing, R.S. Stepleman et. al., eds.,
!       North-Holland, Amsterdam, pp. 55-74.
!
! ARGUMENTS
!   do_qqjk_inchem   : if pr_qqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   do_semiss_inchem : do surface emissions inside the chemistry solver, or
!                      outside?
!   pr_qqjk  : should the periodic qqjk output file be written?
!   pr_smv2  : should the SmvgearII     output file be written
!              (non-parallel mode only)?
!   ifsun    : identifies whether sun is up (=1) or down (=2)
!   ilat     : # of latitudes
!   ilong    : # of longitudes
!   ivert    : # of vertical layers
!   reorder   : 1 => reorder grid-cells and blocks for chemistry
!              2 => solve chemistry
!   numZones   : # of zones (ilong * ilat * ivert)
!   jlooplo  : low ntloop grid-cell - 1 in a grid-block
!   numGridCellsInBlock   : # of grid-cells in a grid-block
!   lunsmv   : logical unit number to write to when pr_smv2 is true
!   numActiveReactants    : # of active rxns
!   nfdh2    : nfdh3 + # of rxns with two   active reactants
!   nfdh3    :         # of rxns with three active reactants
!   nfdl1    : nfdh2 + 1
!   nfdl2    : nfdh3 + 1
!   nfdrep   : nfdh3 + # of rxns with two active reactants that are not
!              followed by a rxn with the same reactants
!   nfdrep1  : nfdrep + 1 ! MRD: No longer needed. Remove it.
!   timeStepDecreaseFraction  : fraction time step is decreased in Smvgear if convergence
!              test fails
!   maxTimeStepNight  : max time step for night all chem (s)
!   pr_nc_period     : NetCDF output period
!   tdt      : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   jreorder : gives original grid-cell from re-ordered grid-cell
!   jphotrat : tbd
!   origJnewSpcNumber  : original spc # of each new jnew spc
!   airDensity   : density of air (molec/cm^3)
!   corig    : original gas-phase concs used to restart Smvgear if a failure
!              occurs (molec/cm^3)
!   pratk1   : tbd
!   surfaceEmissions    : surface emissions (units?)
!   amountAddedToEachSpecies    : amount added to each spc at each grid-cell, for mass balance
!              accounting (# cm^-3 for gas chemistry (?))
!   nfdh1    : nfdh2 + # of rxns with one   active reactant
!   stiffness   : measure of stiffness/nearness to convergence of each block
!              sum ydot/y for all species (MRD per Kareem Sorathia)
!   valuesDecomposedMatrix      : array holding values of decomposed matrix
!   cnew     : stores cnewDerivatives (y (estimated)) (molec/cm^3)
!   gloss    : value of first derivatives on output from velocity; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!   vdiag    : 1 / current diagonal term of the decomposed matrix
!   rateConstants    : rate constants
!-----------------------------------------------------------------------------

      subroutine Smvgear  &
     &  (do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2, ifsun,  &
     &   ilat, ilong, ivert, reorder, numZones, jlooplo, numGridCellsInBlock, lunsmv,  &
     &   numActiveReactants, gasChemistryType, nfdh2, nfdh3, nfdl1, nfdl2, nfdrep, nfdrep1,  &
     &   timeStepDecreaseFraction, maxTimeStepNight, pr_nc_period, tdt, do_cell_chem, irma, irmb,  &
     &   irmc, jreorder, jphotrat, ntspec, origJnewSpcNumber, airDensity, corig,  &
     &   pratk1, surfaceEmissions, amountAddedToEachSpecies, nfdh1, stiffness, valuesDecomposedMatrix, cnew, gloss, vdiag,  &
     &   rateConstants, &
     &   yda, qqkda, qqjda, qkgmi, qjgmi, &
     &   CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
     &   num_qjo, num_qks, num_qjs, num_active, prDiag)

      use GmiPrintError_mod, only : GmiPrintError
      use GmiMechanism_mod
      use GmiManager_mod
      use GmiSparseMatrix_mod
      use Smv2Chem2_mod

      implicit none

#     include "smv2chem_par.h"

      ! Argument declarations.
      integer, intent(in) :: CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2
      integer, intent(in) :: num_qjo, num_qks, num_qjs, num_active
      real*8 , intent(in   ) :: qjgmi(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qjo)
      real*8 , intent(in   ) :: qkgmi(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qks)
      real*8 , intent(inout) :: qqjda(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qjs)
      real*8 , intent(inout) :: qqkda(CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_qks)
      real*8 , intent(inout) :: yda  (CTMi1:CTMi2, CTMju1:CTMj2, CTMk1:CTMk2, num_active)
      logical, intent(in)  :: do_qqjk_inchem
      logical, intent(in)  :: do_semiss_inchem
      logical, intent(in)  :: pr_qqjk
      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: ifsun
      integer, intent(in)  :: ilat, ilong, ivert
      integer, intent(in)  :: reorder
      integer, intent(in)  :: numZones
      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: numGridCellsInBlock
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: numActiveReactants
      integer, intent(in)  :: gasChemistryType
      integer, intent(in)  :: nfdh2,  nfdh3
      integer, intent(in)  :: nfdl1,  nfdl2
      integer, intent(in)  :: nfdrep, nfdrep1
      real*8,  intent(in)  :: timeStepDecreaseFraction
      real*8,  intent(in)  :: maxTimeStepNight
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      logical, intent(in)  :: do_cell_chem(ilong, ilat, ivert)
		!K: irm[a,b,c] can be removed as they are now in GenChem and don't vary with block
      integer, intent(in)  :: irma    (NMTRATE)
      integer, intent(in)  :: irmb    (NMTRATE)
      integer, intent(in)  :: irmc    (NMTRATE)
      integer, intent(in)  :: jreorder(numZones)
      integer, intent(in)  :: jphotrat(ICS)
      integer, intent(in)  :: ntspec  (ICS)
      integer, intent(in)  :: origJnewSpcNumber (MXGSAER, ICS)
      real*8,  intent(in)  :: airDensity  (numGridCellsInBlock)
      real*8,  intent(in)  :: corig   (KBLOOP, MXGSAER)
      real*8,  intent(in)  :: pratk1  (KBLOOP, IPHOT)
      real*8,  intent(in)  :: surfaceEmissions   (ilat*ilong, IGAS)
      real*8,  intent(inout) :: stiffness(numZones)
		!K: Why is valuesDecomposedMatrix passed in/out?  Seems silly
		!K: Same question for everything but cnew
      real*8,  intent(inout) :: valuesDecomposedMatrix   (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: cnew  (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: gloss (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: amountAddedToEachSpecies (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: vdiag (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: rateConstants (KBLOOP, NMTRATE)
      integer, intent(out) :: nfdh1
      logical, intent(in)  :: prDiag



      ! Variable declarations.
      integer :: diurnalGasChemType  ! gasChemistryType       => for daytime   gas chemistry
                       ! gasChemistryType + ICS => for nighttime gas chemistry
      integer :: ibcb(IGAS)
      integer :: concAboveAbtolCount(KBLOOP, 5) ! counts # of concs above abtol(i), i = 1..
      integer :: numFinalMatrixPositions ! excluding diagonal
      integer :: numLoopCycles
      integer :: evaluatePredictor
      integer :: kloop ! loops over cell block


      real*8  :: r1delt
      real*8  :: xtimestep ! for diagnostics
      real*8, parameter  :: MAX_REL_CHANGE = 0.3d0

      real*8  :: dely  (KBLOOP)
      real*8  :: absoluteErrTolerance (KBLOOP) !(molec/cm^-3 for gases)
      real*8  :: accumulatedErrorStorage  (KBLOOP, MXGSAER) ! stores value of accumulatedError when numSuccessStepsBeforeReTest = 1
      real*8  :: explic (KBLOOP, MXGSAER)

      ! an array of length num1stOEqnsSolve*(MAXORD+1) that carries the
      ! derivatives of cnew, scaled by currentTimeStep^j/factorial(j), where j is
      ! the jth derivative; j varies from 1 to orderOfIntegrationMethod; e.g., cnewDerivatives(jspc,2)
      ! stores currentTimeStep*y' (estimated)
      real*8  :: cnewDerivatives  (KBLOOP, MXGSAER*7)

      type (Mechanism) :: mechanismObject
      type (Manager) :: managerObject
      type (SparseMatrix) :: sparseMatrixObject ! not used yet

#     include "setkin_ibcb.h"




      call initializeMechanism (mechanismObject, numGridCellsInBlock, irma, &
                              &  irmb, irmc, nfdh2, nfdh3, nfdrep, rateConstants)
      call resetGear (managerObject, diurnalGasChemType, gasChemistryType, ifsun, maxTimeStepNight)




      numLoopCycles = -1
      do while (managerObject%timeRemainingInChemInterval > MINIMUM_TIME_STEP)
         numLoopCycles = numLoopCycles + 1

         call calculateNewRmsError (managerObject, numGridCellsInBlock, dely, managerObject%correctorIterations)

         if (numLoopCycles == 0) then ! Restarting with new cell block.

           call startTimeInterval (managerObject, gasChemistryType)
           call initConcentrationArray (numGridCellsInBlock, cnew, corig, managerObject)
           call initializeFirstTimeStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
           call doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)

         else

            ! MRD: beware of misleading comments regarding convergence
            if (managerObject%dcon > CONVERGENCE_THRESHOLD) then ! NON-CONVERGENCE

              ! If nonconvergence after one step, re-evaluate first derivative with new values of cnew.
              if (managerObject%correctorIterations == 1) then

                  call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

                  !---!
                  cycle
                  !---!

              ! If the predictor matrix is more than one step old, update it and cycle
              else if (evaluatePredictor == DO_NOT_EVAL_PREDICTOR) then ! This path doesn't seem to be taken by current testing

                  managerObject%numFailOldJacobian = managerObject%numFailOldJacobian + 1
                  evaluatePredictor = EVAL_PREDICTOR
                  call correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
                  call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

                  !---!
                  cycle
                  !---!

              end if

              ! This path doesn't seem to be taken by current testing

              ! If the Predictor is current, then reduce the time step, reset the accumulated derivatives to their values before the failed step, and retry with the smaller step.
              evaluatePredictor     = EVAL_PREDICTOR
              call updateAfterNonConvTightenLimits (managerObject, 2.0d0, managerObject%told, timeStepDecreaseFraction)
              call resetCnewDerivatives (managerObject, cnewDerivatives, numGridCellsInBlock)
              ! these three routines are called in succession several times below
              call doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
              call correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
              call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

              !---!
              cycle
              !---!

            end if  !NON-CONVERGENCE



            ! The corrector iteration CONVERGED.
            evaluatePredictor = DO_NOT_EVAL_PREDICTOR

            if (managerObject%correctorIterations > 1) call testAccumulatedError (managerObject, numGridCellsInBlock, dely)


            ! The accumulated error test failed.
            if (managerObject%der2max > managerObject%enqq) then

               call updateAfterAccumErrorTestFails (managerObject)
               call resetCnewDerivatives (managerObject, cnewDerivatives, numGridCellsInBlock)

               if (managerObject%numFailuresAfterVelocity <= NUM_FAILURES_SMALL) then

                  managerObject%ifsuccess = 0
                  managerObject%timeStepRatioHigherOrder   = 0.0d0

                  call estimateTimeStepRatio (managerObject, numGridCellsInBlock, dely, cnewDerivatives)




                  ! This path isn't being executed by testing data
                  ! If the last step was successful and timeStepRatio is small, keep the
                  ! urrent step and order, and allow three successful steps before
                  ! re-checking the time step and order.
                  if ((managerObject%timeStepRatio < 1.1d0) .and. (managerObject%ifsuccess == 1)) then

                    ! --------------------------------------------------------------------!
                    ! The code here is a duplicate of a block below
                    managerObject%numSuccessStepsBeforeReTest = 3

                    call doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
                    call correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
                    call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

                    !---!
                    cycle
                    !---!

                  ! If the maximum time step ratio is that of one order lower than
                  ! the current order, decrease the order.  Do not minimize timeStepRatio
                  ! to <= 1, when ifsuccess = 0 since this is less efficient.
                  else if (managerObject%timeStepRatio == managerObject%timeStepRatioLowerOrder) then
                    managerObject%orderOfIntegrationMethod = managerObject%orderOfIntegrationMethod - 1

                  else if (managerObject%timeStepRatio == managerObject%timeStepRatioHigherOrder) then
                     call increaseOrderAndAddDerivativeTerm (managerObject, cnewDerivatives, numGridCellsInBlock)
                  end if

                  ! If the last two steps have failed, re-set numSuccessStepsBeforeReTest to the current order + 1.
                  ! (Do not minimize timeStepRatio if managerObject%numFailuresAfterVelocity >= 2 since tests show
                  ! that this merely leads to additional computations.)
                  managerObject%numSuccessStepsBeforeReTest = managerObject%orderOfIntegrationMethod + 1

                  call doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
                  call correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
                  call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

                  ! End duplicate code block
                  ! --------------------------------------------------------------------!

                  !---!
                  cycle
                  !---!




               ! This path is not being executed with testing data
               ! if the first attempts fail, retry the step at timeStepDecreaseFraction
               else if (managerObject%numFailuresAfterVelocity <= NUM_FAILURES_LARGE) then

                  managerObject%ifsuccess = 0
                  managerObject%timeStepRatio = timeStepDecreaseFraction

                  call doSubStep(ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
                  call correctorStep(valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
                  call advanceTimeStep(ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

                  !---!
                  cycle
                  !---!

               ! this path isn't being executed by testing data (we need bad input data)
               else ! numFailuresAfterVelocity >= NUM_FAILURES_LARGE

                  call resetTermsBeforeStartingOver (managerObject, cnew, cnewDerivatives, &
                                                & numGridCellsInBlock, lunsmv, pr_smv2)

                  if (managerObject%numExcessiveFailures == LIMIT_EXCESSIVE_FAILURES) then
                     if (pr_smv2) Write(*,*) "Smvgear:  Stopping because of excessive errors."
                     call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
                  end if


                  call initializeFirstTimeStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
                  call doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
                  call correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
                  call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

                  !---!
                  cycle
                  !---!

               end if



            else ! The accumulated error test did not fail


              ! After a successful step, do diagnostics
              if (pr_qqjk .and. do_qqjk_inchem) then
                xtimestep = managerObject%elapsedTimeInChemInterval - managerObject%told

                call Do_Smv2_Diag  &

                 &      (jlooplo, numGridCellsInBlock, pr_nc_period, tdt, managerObject%told, do_cell_chem,  &
                 &       jreorder, origJnewSpcNumber, airDensity, cnew, xtimestep, &
                 &       yda, qqkda, qqjda, qkgmi, qjgmi, &
                 &       ilong, ilat, ivert, numZones, &
                 &       CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
                 &       num_qjo, num_qks, num_qjs, num_active)
              end if

              ! update the concentration and all derivatives (reset told, set ifsuccess = 1, increment numSuccessTdt)
              call updateAndResetAfterSucessfulStep (managerObject, cnewDerivatives, numGridCellsInBlock)
              call doMassBalanceAccounting (managerObject%integrationOrderCoeff1, managerObject%num1stOEqnsSolve, numGridCellsInBlock, &
                                             amountAddedToEachSpecies, managerObject%accumulatedError, explic, cnewDerivatives)

              ! Exit smvgear if a time interval has been completed.
              managerObject%timeRemainingInChemInterval = managerObject%chemTimeInterval - managerObject%elapsedTimeInChemInterval
              if (managerObject%timeRemainingInChemInterval <= 1.0d-06) then

               !----!
               return
               !----!

              endif



              ! numSuccessStepsBeforeReTest counts the number of successful steps before re-testing the
              ! if numSuccessStepsBeforeReTest > 1, decrease it and go on to the next time step with
              ! the current step-size and order;
              if (managerObject%numSuccessStepsBeforeReTest > 1) then
                call storeAccumErrorAndSetTimeStepRatio (managerObject, accumulatedErrorStorage, numGridCellsInBlock)
                call doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
                call correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
                call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

                !---!
                cycle
                !---!

              end if

            end if

            call estimateTimeStepRatioOneOrderHigher (managerObject, MAXORD, numGridCellsInBlock, &
                        dely, accumulatedErrorStorage)

            call estimateTimeStepRatio (managerObject, numGridCellsInBlock, dely, cnewDerivatives)


            ! If the last step was successful and timeStepRatio is small, keep the
            ! current step and order, and allow three successful steps before
            ! re-checking the time step and order.
            if ((managerObject%timeStepRatio < 1.1d0) .and. (managerObject%ifsuccess == 1)) then

              ! --------------------------------------------------------------------!
              ! The code here is a duplicate of a block above
              managerObject%numSuccessStepsBeforeReTest = 3

              call doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
              call correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
              call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

              !---!
              cycle
              !---!

            ! If the maximum time step ratio is that of one order lower than
            ! the current order, decrease the order.  Do not minimize timeStepRatio
            ! to <= 1, when ifsuccess = 0 since this is less efficient.
            else if (managerObject%timeStepRatio == managerObject%timeStepRatioLowerOrder) then
              managerObject%orderOfIntegrationMethod = managerObject%orderOfIntegrationMethod - 1

            else if (managerObject%timeStepRatio == managerObject%timeStepRatioHigherOrder) then
               call increaseOrderAndAddDerivativeTerm (managerObject, cnewDerivatives, numGridCellsInBlock)
            end if

            ! If the last two steps have failed, re-set numSuccessStepsBeforeReTest to the current
            ! order + 1.  Do not minimize timeStepRatio if managerObject%numFailuresAfterVelocity >= 2 since tests show
            ! that this merely leads to additional computations.
            managerObject%numSuccessStepsBeforeReTest = managerObject%orderOfIntegrationMethod + 1

            call doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
            call correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
            call advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)


            ! End duplicate code block - MRD: there may need to be a cycle here!?
            ! --------------------------------------------------------------------!
         end if



      end do

      return

   end subroutine Smvgear



      subroutine initializeFirstTimeStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
         use GmiManager_mod
         use GmiMechanism_mod
         use GmiSparseMatrix_mod
         use Smv2Chem2_mod
         implicit none

#     include "smv2chem_par.h"
         integer, intent(in) :: ilat
         integer, intent(in) :: ilong
         integer, intent(in) :: numZones
         integer, intent(in) :: origJnewSpcNumber(MXGSAER, ICS)
         integer, intent(in) :: reorder
         integer, intent(in) :: jlooplo
         integer, intent(in) :: jphotrat(ICS)
         integer, intent(in) :: jreorder(numZones)
         integer, intent(in) :: numGridCellsInBlock
         integer, intent(in) :: lunsmv
         integer :: concAboveAbtolCount(KBLOOP, 5)
         integer :: evaluatePredictor
         integer, intent(in) :: gasChemistryType
         integer :: diurnalGasChemType
         integer, intent(out) :: nfdh1
         integer, intent(in) :: ntspec(ICS)
         integer, intent(in) :: numActiveReactants
         integer :: numFinalMatrixPositions
         real*8 :: absoluteErrTolerance(KBLOOP)
         real*8 :: valuesDecomposedMatrix(KBLOOP, 0:MXARRAY)
         real*8 :: cnew(KBLOOP, MXGSAER)
         real*8 :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(in) :: corig(KBLOOP, MXGSAER)
         real*8 :: dely(KBLOOP)
         real*8 :: stiffness(numZones)
         real*8 :: explic(KBLOOP, MXGSAER)
         real*8 :: gloss(KBLOOP, MXGSAER)
         real*8 :: MAX_REL_CHANGE
         real*8, intent(in) :: pratk1(KBLOOP, IPHOT)
         real*8 :: r1delt
         real*8 :: vdiag(KBLOOP, MXGSAER)
         real*8, intent(in) :: surfaceEmissions(ilat*ilong, IGAS)
         logical, intent(in) :: do_semiss_inchem
         logical, intent(in) :: prDiag
         logical, intent(in) :: pr_smv2
         type(Mechanism) :: mechanismObject
         type(Manager) :: managerObject

         integer :: kloop

         call resetBeforeUpdate (managerObject)

         !!DIR$ INLINE
         call updatePhotoDissRates  (mechanismObject, numGridCellsInBlock, numActiveReactants, gasChemistryType, diurnalGasChemType, jphotrat, pratk1)
         !!DIR$ NOINLINE

         call velocity(mechanismObject, managerObject%num1stOEqnsSolve, diurnalGasChemType, cnew, gloss, nfdh1)
         managerObject%numCallsVelocity = managerObject%numCallsVelocity + 1

         call setBoundaryConditions (mechanismObject, numZones, jreorder, jlooplo, ilat, &
               ilong, ntspec, gasChemistryType, origJnewSpcNumber, do_semiss_inchem, gloss, surfaceEmissions)

         call determineInitialAbTol (managerObject, cnew, concAboveAbtolCount, reorder, numGridCellsInBlock, gasChemistryType, absoluteErrTolerance)

         call calculateErrorTolerances (managerObject, numGridCellsInBlock, cnew, gloss, dely)

         if (reorder /= SOLVE_CHEMISTRY) then
            do kloop = 1, numGridCellsInBlock
               stiffness(jlooplo+kloop) = dely(kloop)
            end do
            return
         end if

         call calcInitialTimeStepSize (managerObject, numGridCellsInBlock, dely, gasChemistryType)

         call setInitialOrder (managerObject, evaluatePredictor)

         call storeInitConcAndDerivatives(managerObject%num1stOEqnsSolve, numGridCellsInBlock, &
                              cnewDerivatives, cnew, managerObject%currentTimeStep, gloss)

      end subroutine initializeFirstTimeStep


      subroutine doSubStep (ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)
         use GmiManager_mod
         use GmiMechanism_mod
         use GmiSparseMatrix_mod
         use Smv2Chem2_mod
         implicit none

#     include "smv2chem_par.h"
         integer, intent(in) :: ilat
         integer, intent(in) :: ilong
         integer, intent(in) :: numZones
         integer, intent(in) :: origJnewSpcNumber(MXGSAER, ICS)
         integer, intent(in) :: reorder
         integer, intent(in) :: jlooplo
         integer, intent(in) :: jphotrat(ICS)
         integer, intent(in) :: jreorder(numZones)
         integer :: kloop
         integer, intent(in) :: numGridCellsInBlock
         integer, intent(in) :: lunsmv
         integer :: evaluatePredictor
         integer :: concAboveAbtolCount(KBLOOP, 5)
         integer, intent(in) :: gasChemistryType
         integer :: diurnalGasChemType
         integer, intent(out) :: nfdh1
         integer, intent(in) :: ntspec(ICS)
         integer, intent(in) :: numActiveReactants
         integer :: numFinalMatrixPositions
         real*8 :: absoluteErrTolerance(KBLOOP)
         real*8 :: valuesDecomposedMatrix(KBLOOP, 0:MXARRAY)
         real*8 :: cnew(KBLOOP, MXGSAER)
         real*8 :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(in) :: corig(KBLOOP, MXGSAER)
         real*8 :: dely(KBLOOP)
         real*8 :: stiffness(numZones)
         real*8 :: explic(KBLOOP, MXGSAER)
         real*8 :: gloss(KBLOOP, MXGSAER)
         real*8 :: MAX_REL_CHANGE
         real*8, intent(in) :: pratk1(KBLOOP, IPHOT)
         real*8 :: r1delt
         real*8 :: vdiag(KBLOOP, MXGSAER)
         real*8, intent(in) :: surfaceEmissions(ilat*ilong, IGAS)
         logical, intent(in) :: do_semiss_inchem
         logical, intent(in) :: prDiag
         logical, intent(in) :: pr_smv2
         type(Manager) :: managerObject
         type(Mechanism) :: mechanismObject

         do while (.true.) ! could be do while managerObject%currentTimeStep < HMIN

            if (managerObject%orderOfIntegrationMethod /= managerObject%oldOrderOfIntegrationMethod) then
               call updateCoefficients (managerObject)
            endif

            call calculateTimeStep (managerObject, evaluatePredictor, MAX_REL_CHANGE)

            if (managerObject%currentTimeStep < HMIN) then !HMIN is a minimum acceptable time step

               call tightenErrorTolerance (managerObject, pr_smv2, lunsmv, gasChemistryType)
               call startTimeInterval (managerObject, gasChemistryType)
               call initConcentrationArray(numGridCellsInBlock, cnew, corig, managerObject)

               call initializeFirstTimeStep(ilat, ilong, numZones, absoluteErrTolerance, valuesDecomposedMatrix, &
                        & cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, stiffness, &
                        & evaluatePredictor, explic, gloss, origJnewSpcNumber, reorder, jlooplo, jphotrat, jreorder, &
                        & kloop, numGridCellsInBlock, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, gasChemistryType, &
                        & diurnalGasChemType, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, surfaceEmissions)

            else ! order of old integration method is same as current order

               if (managerObject%timeStepRatio /= 1.0d0) then
                  call scaleDerivatives (managerObject, numGridCellsInBlock, cnewDerivatives)
               end if

               if (managerObject%ifsuccess == 1) then
                  managerObject%rdelmax = 10.0d0

                  if (Mod (managerObject%numSuccessTdt, 3) == 2) then
                     call calcNewAbsoluteErrorTolerance (managerObject, cnew, concAboveAbtolCount, &
                        numGridCellsInBlock, absoluteErrTolerance, gasChemistryType)
                  end if

                  call updateChold (managerObject, numGridCellsInBlock, cnew, absoluteErrTolerance)
               end if

               call predictConcAndDerivatives (managerObject, cnewDerivatives, explic, numGridCellsInBlock, prDiag)
               call correctorStep(valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
               call advanceTimeStep(ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)

               exit

            end if
         end do

      end subroutine doSubStep


      subroutine correctorStep (valuesDecomposedMatrix, cnew, cnewDerivatives, evaluatePredictor, numGridCellsInBlock, managerObject, mechanismObject, diurnalGasChemType, numFinalMatrixPositions, r1delt, vdiag)
         use GmiManager_mod
         use GmiMechanism_mod
         use GmiSparseMatrix_mod
         use Smv2Chem2_mod
         implicit none

#     include "smv2chem_par.h"

         real*8 :: valuesDecomposedMatrix(KBLOOP, 0:MXARRAY)
         real*8 :: cnew(KBLOOP, MXGSAER)
         real*8 :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer, intent(inout) :: evaluatePredictor
         integer, intent(in) :: numGridCellsInBlock
         integer :: diurnalGasChemType
         integer :: numFinalMatrixPositions
         type(Manager) :: managerObject
         type(Mechanism) :: mechanismObject

         real*8 :: r1delt
         real*8 :: vdiag(KBLOOP, MXGSAER)

         call initCorrector (managerObject, numGridCellsInBlock, cnew, cnewDerivatives)

         ! Re-evaluate predictor matrix before starting the corrector iteration.
         if (evaluatePredictor == EVAL_PREDICTOR) then

            r1delt = -managerObject%integrationOrderCoeff1 * managerObject%currentTimeStep
            numFinalMatrixPositions  = sparseMatrixDimension(diurnalGasChemType) - managerObject%num1stOEqnsSolve

            call calculatePredictor (numFinalMatrixPositions, sparseMatrixDimension(diurnalGasChemType), &
                  numGridCellsInBlock, cnew, npdhi(diurnalGasChemType), npdlo(diurnalGasChemType), r1delt, valuesDecomposedMatrix, &
                  mechanismObject%rateConstants)

            managerObject%numCallsPredict = managerObject%numCallsPredict + 1
            evaluatePredictor  = PREDICTOR_JUST_CALLED

         ! K: Consider un-inlining this.
         !!DIR$   INLINE
            call LU_Decomp  (managerObject%num1stOEqnsSolve, numGridCellsInBlock, diurnalGasChemType, valuesDecomposedMatrix, vdiag)
         !!DIR$   NOINLINE

            call setConvergenceTerms (managerObject, MBETWEEN)

         end if

      end subroutine correctorStep


      subroutine advanceTimeStep (ilat, ilong, numZones, valuesDecomposedMatrix, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, origJnewSpcNumber, jlooplo, jreorder, numGridCellsInBlock, managerObject, mechanismObject, gasChemistryType, diurnalGasChemType, nfdh1, ntspec, vdiag, surfaceEmissions)
         use GmiManager_mod
         use GmiMechanism_mod
         use GmiSparseMatrix_mod
         implicit none

#     include "smv2chem_par.h"
         integer, intent(in) :: ilat
         integer, intent(in) :: ilong
         integer, intent(in) :: numZones
         real*8, intent(inout) :: valuesDecomposedMatrix(KBLOOP, 0:MXARRAY)
         real*8, intent(inout) :: cnew(KBLOOP, MXGSAER)
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(inout) :: dely(KBLOOP)
         logical, intent(in) :: do_semiss_inchem
         real*8 :: gloss(KBLOOP, MXGSAER)
         integer, intent(in) :: origJnewSpcNumber(MXGSAER, ICS)
         integer, intent(in) :: jlooplo
         integer, intent(in) :: jreorder(numZones)
         integer, intent(in) :: numGridCellsInBlock
         type(Manager) :: managerObject
         type(Mechanism) :: mechanismObject
         integer, intent(in) :: gasChemistryType
         integer :: diurnalGasChemType
         integer, intent(out) :: nfdh1
         integer, intent(in) :: ntspec(ICS)
         real*8 :: vdiag(KBLOOP, MXGSAER)
         real*8, intent(in) :: surfaceEmissions(ilat*ilong, IGAS)

         ! Evaluate the first derivative using corrected values of cnew.
         call velocity (mechanismObject, managerObject%num1stOEqnsSolve, diurnalGasChemType, cnew, gloss, nfdh1)
         managerObject%numCallsVelocity = managerObject%numCallsVelocity + 1

         call setBoundaryConditions (mechanismObject, numZones, jreorder, jlooplo, ilat, &
                  ilong, ntspec, gasChemistryType, origJnewSpcNumber, do_semiss_inchem, gloss, surfaceEmissions)

         call computeErrorFromCorrected1stDeriv (managerObject%num1stOEqnsSolve, numGridCellsInBlock, &
                  gloss, managerObject%currentTimeStep, cnewDerivatives, managerObject%accumulatedError)

         call Backsub (managerObject%num1stOEqnsSolve, numGridCellsInBlock, diurnalGasChemType, valuesDecomposedMatrix, vdiag, gloss)

         call sumAccumulatedError (managerObject, cnew, cnewDerivatives, dely, gloss, numGridCellsInBlock)

      end subroutine advanceTimeStep











      subroutine doMassBalanceAccounting (integrationOrderCoeff, num1stOEqnsSolve, numGridCellsInBlock, amountAddedToEachSpecies, &
                                       &  accumulatedError, explic, cnewDerivatives)

         implicit none
#     include "smv2chem_par.h"

         real*8, intent(in) :: integrationOrderCoeff
         integer, intent(in) :: num1stOEqnsSolve
         integer, intent(in) :: numGridCellsInBlock
         real*8, intent(inout) :: amountAddedToEachSpecies(KBLOOP, MXGSAER)
         real*8, intent(in)  :: accumulatedError (KBLOOP, MXGSAER) ! on a successful return; accumulatedError(kloop,i) contains the estimated one step local error in cnew
         real*8, intent(in) :: explic(KBLOOP, MXGSAER)
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)


         integer :: i, kloop
         real*8 :: dtasn1

         if (integrationOrderCoeff == 1.0d0) then
            do i = 1, num1stOEqnsSolve
              do kloop = 1, numGridCellsInBlock
                amountAddedToEachSpecies(kloop,i) = amountAddedToEachSpecies(kloop,i) + accumulatedError(kloop,i) + explic(kloop,i)
                cnewDerivatives(kloop,i) = cnewDerivatives(kloop,i) + accumulatedError(kloop,i)
              end do
            end do

          else

            do i = 1, num1stOEqnsSolve
              do kloop = 1, numGridCellsInBlock
                dtasn1         = integrationOrderCoeff * accumulatedError(kloop,i)
                amountAddedToEachSpecies(kloop,i) = amountAddedToEachSpecies(kloop,i) + dtasn1 + explic(kloop,i)
                cnewDerivatives (kloop,i) = cnewDerivatives (kloop,i) + dtasn1
              end do
            end do

          end if

      end subroutine doMassBalanceAccounting







!-----------------------------------------------------------------------------
!
! ROUTINE
!   computeErrorFromCorrected1stDeriv
! DESCRIPTION
! In the case of the chord method, compute error (gloss) from the
! corrected calculation of the first derivative.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
      subroutine computeErrorFromCorrected1stDeriv (num1stOEqnsSolve, numGridCellsInBlock, gloss, &
                  & currentTimeStep, cnewDerivatives, accumulatedError)

         implicit none
#     include "smv2chem_par.h"

         integer, intent(in) :: num1stOEqnsSolve, numGridCellsInBlock
         real*8, intent(inout) :: gloss (KBLOOP, MXGSAER)
         real*8, intent(in) :: currentTimeStep
         real*8, intent(in) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(in)  :: accumulatedError (KBLOOP, MXGSAER)!

         integer jspc, j, kloop

         do jspc = 1, num1stOEqnsSolve
           j = jspc + num1stOEqnsSolve
           do kloop = 1, numGridCellsInBlock
             gloss(kloop,jspc) = (currentTimeStep * gloss(kloop,jspc)) -  &
                  (cnewDerivatives(kloop,j) + accumulatedError(kloop,jspc))
           end do
         end do

      end subroutine computeErrorFromCorrected1stDeriv



!-----------------------------------------------------------------------------
!
! ROUTINE
!   storeInitConcAndDerivatives
! DESCRIPTION
! Store initial concentration and first derivatives x time step.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
      subroutine storeInitConcAndDerivatives(num1stOEqnsSolve, numGridCellsInBlock, cnewDerivatives, cnew, currentTimeStep, gloss)

         implicit none
#     include "smv2chem_par.h"

         integer, intent(in) :: num1stOEqnsSolve
         integer, intent(in) :: numGridCellsInBlock
         real*8, intent(out) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(in) :: cnew(KBLOOP, MXGSAER)
         real*8, intent(in) :: currentTimeStep
         real*8, intent(in) :: gloss(KBLOOP, MXGSAER)

         ! local variables
         integer :: kloop, jspc, j

         do jspc = 1, num1stOEqnsSolve
            j = jspc + num1stOEqnsSolve

           do kloop = 1, numGridCellsInBlock
             cnewDerivatives(kloop,jspc) = cnew(kloop,jspc)
             cnewDerivatives(kloop,j)    = currentTimeStep * gloss(kloop,jspc)
           end do

         end do
      end subroutine storeInitConcAndDerivatives




      subroutine initConcentrationArray(numGridCellsInBlock, concentrationsNew, concentrationsOld, managerObject)

         use GmiManager_mod
         implicit none

#     include "smv2chem_par.h"

         integer, intent(in) :: numGridCellsInBlock
         real*8, intent(out) :: concentrationsNew(KBLOOP, MXGSAER)
         real*8,  intent(in)  :: concentrationsOld(KBLOOP, MXGSAER)
         type (Manager) :: managerObject

         integer :: jnew, kloop

         do jnew = 1, managerObject%num1stOEqnsSolve
           do kloop = 1, numGridCellsInBlock
             concentrationsNew(kloop, jnew) = concentrationsOld(kloop, jnew) ! why save this?
           end do
         end do

      end subroutine





