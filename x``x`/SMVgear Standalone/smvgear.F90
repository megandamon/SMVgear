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
!   Backsub
!   Smvgear
!   Decomp
!   Update
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
!   ireord   : 1 => reorder grid-cells and blocks for chemistry
!              2 => solve chemistry
!   itloop   : # of zones (ilong * ilat * ivert)
!   jlooplo  : low ntloop grid-cell - 1 in a grid-block
!   ktloop   : # of grid-cells in a grid-block
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
!   hmaxnit  : max time step for night all chem (s)
!   pr_nc_period     : NetCDF output period
!   tdt      : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   jreorder : gives original grid-cell from re-ordered grid-cell
!   jphotrat : tbd
!   inewold  : original spc # of each new jnew spc
!   denair   : density of air (molec/cm^3)
!   corig    : original gas-phase concs used to restart Smvgear if a failure
!              occurs (molec/cm^3)
!   pratk1   : tbd
!   yemis    : surface emissions (units?)
!   amountAddedToEachSpecies    : amount added to each spc at each grid-cell, for mass balance
!              accounting (# cm^-3 for gas chemistry (?))
!   nfdh1    : nfdh2 + # of rxns with one   active reactant
!   errmx2   : measure of stiffness/nearness to convergence of each block
!              sum ydot/y for all species (MRD per Kareem Sorathia)
!   cc2      : array holding values of decomposed matrix
!   cnew     : stores cnewDerivatives (y (estimated)) (molec/cm^3)
!   gloss    : value of first derivatives on output from velocity; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!   vdiag    : 1 / current diagonal term of the decomposed matrix
!   rrate    : rate constants
!   trate    : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?)) !REMOVED!
!
!-----------------------------------------------------------------------------

      subroutine Smvgear  &
     &  (do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2, ifsun,  &
     &   ilat, ilong, ivert, ireord, itloop, jlooplo, ktloop, lunsmv,  &
     &   numActiveReactants, ncs, nfdh2, nfdh3, nfdl1, nfdl2, nfdrep, nfdrep1,  &
     &   timeStepDecreaseFraction, hmaxnit, pr_nc_period, tdt, do_cell_chem, irma, irmb,  &
     &   irmc, jreorder, jphotrat, ntspec, inewold, denair, corig,  &
     &   pratk1, yemis, amountAddedToEachSpecies, nfdh1, errmx2, cc2, cnew, gloss, vdiag,  &
     &   rrate, &
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
      integer, intent(in)  :: ireord
      integer, intent(in)  :: itloop
      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: numActiveReactants
      integer, intent(in)  :: ncs
      integer, intent(in)  :: nfdh2,  nfdh3
      integer, intent(in)  :: nfdl1,  nfdl2
      integer, intent(in)  :: nfdrep, nfdrep1
      real*8,  intent(in)  :: timeStepDecreaseFraction
      real*8,  intent(in)  :: hmaxnit
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      logical, intent(in)  :: do_cell_chem(ilong, ilat, ivert)
		!K: irm[a,b,c] can be removed as they are now in GenChem and don't vary with block
      integer, intent(in)  :: irma    (NMTRATE)
      integer, intent(in)  :: irmb    (NMTRATE)
      integer, intent(in)  :: irmc    (NMTRATE)
      integer, intent(in)  :: jreorder(itloop)
      integer, intent(in)  :: jphotrat(ICS)
      integer, intent(in)  :: ntspec  (ICS)
      integer, intent(in)  :: inewold (MXGSAER, ICS)
      real*8,  intent(in)  :: denair  (ktloop)
      real*8,  intent(in)  :: corig   (KBLOOP, MXGSAER)
      real*8,  intent(in)  :: pratk1  (KBLOOP, IPHOT)
      real*8,  intent(in)  :: yemis   (ilat*ilong, IGAS)
      real*8,  intent(inout) :: errmx2(itloop)
		!K: Why is cc2 passed in/out?  Seems silly
		!K: Same question for everything but cnew
      real*8,  intent(inout) :: cc2   (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: cnew  (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: gloss (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: amountAddedToEachSpecies (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: vdiag (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: rrate (KBLOOP, NMTRATE)
      integer, intent(out) :: nfdh1
      logical, intent(in)  :: prDiag

      ! Variable declarations.
      integer :: evaluatePredictor
      integer :: kloop ! loops over cell block
      integer :: ncsp  ! ncs       => for daytime   gas chemistry
                       ! ncs + ICS => for nighttime gas chemistry
      integer :: ibcb(IGAS)
      ! counts # of concs above abtol(i), i = 1..
      integer :: concAboveAbtolCount(KBLOOP, 5)
      real*8  :: r1delt
      real*8  :: xtimestep

      real*8, parameter  :: MAX_REL_CHANGE = 0.3d0

      real*8  :: dely  (KBLOOP)
      real*8  :: absoluteErrTolerance (KBLOOP) !(molec/cm^-3 for gases)
      real*8  :: accumulatedErrorStorage  (KBLOOP, MXGSAER) ! stores value of accumulatedError when idoub = 1
      real*8  :: explic(KBLOOP, MXGSAER)
      ! an array of length num1stOEqnsSolve*(MAXORD+1) that carries the
      ! derivatives of cnew, scaled by currentTimeStep^j/factorial(j), where j is
      ! the jth derivative; j varies from 1 to orderOfIntegrationMethod; e.g., cnewDerivatives(jspc,2)
      ! stores currentTimeStep*y' (estimated)
      real*8  :: cnewDerivatives  (KBLOOP, MXGSAER*7)

      type (Mechanism_type) :: mechanismObject
      type (Manager_type) :: managerObject
      integer :: numFinalMatrixPositions ! excluding diagonal

#     include "setkin_ibcb.h"

      call initializeMechanism (mechanismObject, ktloop, irma, &
                              &  irmb, irmc, nfdh2, nfdh3, nfdrep, rrate)
      call resetGear (managerObject, ncsp, ncs, ifsun, hmaxnit)


   do while (.true.)


      if (managerObject%dcon == 0.0d0) then

        call calculateNewRmsError (managerObject, ktloop, dely, managerObject%correctorIterations)
        call startTimeInterval (managerObject, ncs)
        call initConcentrationArray(ktloop, cnew, corig, managerObject)

         ! Restarting with new cell block.
         call old150CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
         call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)

         cycle
      else

         call calculateNewRmsError (managerObject, ktloop, dely, managerObject%correctorIterations)

         ! MRD: the comments below may be misleading.
         if (managerObject%dcon > 1.0d0) then ! NON-CONVERGENCE

           ! If nonconvergence after one step, re-evaluate first derivative with new values of cnew.
           if (managerObject%correctorIterations == 1) then

             ! re-eval first derivative (calls velocity)
             call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)
             cycle

            ! If the Jacobian (predictor?) matrix is more than one step old, update it,
            !  and try convergence again.
           else if (evaluatePredictor == DO_NOT_EVAL_PREDICTOR) then ! This path doesn't seem to be taken by current testing

             managerObject%numFailOldJacobian = managerObject%numFailOldJacobian + 1
             evaluatePredictor = EVAL_PREDICTOR

             ! calls corrector loop / re-evaluates predictor
            call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
            call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

            cycle
           end if

            ! This path doesn't seem to be taken by current testing
            ! if the Jacobian is current, then reduce the time step,
            ! reset the accumulated derivatives to their values before the failed step,
            ! and retry with the smaller step.
            evaluatePredictor     = EVAL_PREDICTOR
            call updateAfterNonConvTightenLimits (managerObject, 2.0d0, managerObject%told, timeStepDecreaseFraction)
            call resetCnewDerivatives(managerObject, cnewDerivatives, ktloop)

           !go to 200 ! tighten limit, then re-evaluate the predictor, then call velocity
           call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
           call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
           call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

           cycle
         end if

         ! The corrector iteration CONVERGED.
         evaluatePredictor = DO_NOT_EVAL_PREDICTOR
         if (managerObject%correctorIterations > 1) then
            call testAccumulatedError (managerObject, ktloop, dely)
         end if

         ! The accumulated error test failed.
         if (managerObject%der2max > managerObject%enqq) then

            ! In all cases, reset the derivatives to their values before the last time step.
            call updateAfterAccumErrorTestFails (managerObject)
            call resetCnewDerivatives (managerObject, cnewDerivatives, ktloop)

            ! MRD: magic numbers
            ! re-estimate a time step at the same or one lower order and retry the step;
            if (managerObject%numFailuresAfterVelocity <= 6) then

               managerObject%ifsuccess = 0
               managerObject%timeStepRatioHigherOrder   = 0.0d0

               call estimateTimeStepRatio (managerObject, ktloop, dely, cnewDerivatives)

               ! this path isn't being executed by testing data
               !     If the last step was successful and timeStepRatio is small, keep the
               !     current step and order, and allow three successful steps before
               !     re-checking the time step and order.
               if ((managerObject%timeStepRatio < 1.1d0) .and. (managerObject%ifsuccess == 1)) then

                 managerObject%idoub = 3

                 call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
                 call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
                 call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

                  cycle

               ! If the maximum time step ratio is that of one order lower than
               ! the current order, decrease the order.  Do not minimize timeStepRatio
               ! to <= 1, when ifsuccess = 0 since this is less efficient.
               else if (managerObject%timeStepRatio == managerObject%timeStepRatioLowerOrder) then
                 managerObject%orderOfIntegrationMethod = managerObject%orderOfIntegrationMethod - 1

               else if (managerObject%timeStepRatio == managerObject%timeStepRatioHigherOrder) then
                  call increaseOrderAndAddDerivativeTerm (managerObject, cnewDerivatives, ktloop)
               end if

               !     If the last two steps have failed, re-set idoub to the current
               !     order + 1.  Do not minimize timeStepRatio if managerObject%numFailuresAfterVelocity >= 2 since tests show
               !     that this merely leads to additional computations.
               managerObject%idoub = managerObject%orderOfIntegrationMethod + 1

               call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
               call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
               call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

               cycle

            ! this path is not being executed with testing data
            ! if the first attempts fail, retry the step at timeStepDecreaseFraction
            else if (managerObject%numFailuresAfterVelocity <= 20) then

               managerObject%ifsuccess = 0
               managerObject%timeStepRatio     = timeStepDecreaseFraction

               call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
               call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
               call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

               cycle

            ! this path isn't being executed by testing data
            else

               call resetTermsBeforeStartingOver (managerObject, cnew, cnewDerivatives, &
                                             & ktloop, lunsmv, pr_smv2)

               if (managerObject%numExcessiveFailures == LIMIT_EXCESSIVE_FAILURES) then
                  if (pr_smv2) Write(*,*) "Smvgear:  Stopping because of excessive errors."
                  call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
               end if

               !GO TO 150
               call old150CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
               call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
               call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
               call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

               cycle
            end if

         ! The accumulated error test did not fail
         else

         !       After a successful step, update the concentration and all
         !       derivatives, reset told, set ifsuccess = 1, increment numSuccessTdt,
           if (pr_qqjk .and. do_qqjk_inchem) then
             xtimestep = managerObject%elapsedTimeInChemInterval - managerObject%told

             call Do_Smv2_Diag  &

        &      (jlooplo, ktloop, pr_nc_period, tdt, managerObject%told, do_cell_chem,  &
        &       jreorder, inewold, denair, cnew, xtimestep, &
        &       yda, qqkda, qqjda, qkgmi, qjgmi, &
        &       ilong, ilat, ivert, itloop, &
        &       CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
        &       num_qjo, num_qks, num_qjs, num_active)
           end if

            call updateAndResetAfterSucessfulStep (managerObject, cnewDerivatives, ktloop)
            call doMassBalanceAccounting (managerObject%integrationOrderCoeff1, managerObject%num1stOEqnsSolve, ktloop, &
                                          amountAddedToEachSpecies, managerObject%accumulatedError, explic, cnewDerivatives)

            !       Exit smvgear if a time interval has been completed.
           managerObject%timeRemainingInChemInterval = managerObject%chemTimeInterval - managerObject%elapsedTimeInChemInterval
           if (managerObject%timeRemainingInChemInterval <= 1.0d-06) return


           !       idoub counts the number of successful steps before re-testing the
           !       step-size and order: if idoub > 1, decrease idoub and go on to the next time step with
           !       the current step-size and order;
           if (managerObject%idoub > 1) then
             call storeAccumErrorAndSetTimeStepRatio(managerObject, accumulatedErrorStorage, ktloop)
             call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
             call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
             call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

             cycle
           end if

         end if

         !     Test whether to change the step-size and order.
         !     Determine the time step at (a) one order lower than, (b) the same
         !     order as, and (c) one order higher than the current order.  In the
         !     case of multiple grid-cells in a grid-block, find the minimum
         !     step size among all the cells for each of the orders.  Then, in
         !     all cases, choose the longest time step among the three steps
         !     paired with orders, and choose the order allowing this longest
         !     step.
         call estimateTimeStepRatioOneOrderHigher (managerObject, MAXORD, ktloop, &
                     dely, accumulatedErrorStorage)

         call estimateTimeStepRatio (managerObject, ktloop, dely, cnewDerivatives)

         !     If the last step was successful and timeStepRatio is small, keep the
         !     current step and order, and allow three successful steps before
         !     re-checking the time step and order.
         if ((managerObject%timeStepRatio < 1.1d0) .and. (managerObject%ifsuccess == 1)) then

           managerObject%idoub = 3

           call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
           call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
           call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

            cycle

         ! If the maximum time step ratio is that of one order lower than
         ! the current order, decrease the order.  Do not minimize timeStepRatio
         ! to <= 1, when ifsuccess = 0 since this is less efficient.
         else if (managerObject%timeStepRatio == managerObject%timeStepRatioLowerOrder) then
           managerObject%orderOfIntegrationMethod = managerObject%orderOfIntegrationMethod - 1

         else if (managerObject%timeStepRatio == managerObject%timeStepRatioHigherOrder) then
            call increaseOrderAndAddDerivativeTerm (managerObject, cnewDerivatives, ktloop)
         end if

         !     If the last two steps have failed, re-set idoub to the current
         !     order + 1.  Do not minimize timeStepRatio if managerObject%numFailuresAfterVelocity >= 2 since tests show
         !     that this merely leads to additional computations.
         managerObject%idoub = managerObject%orderOfIntegrationMethod + 1

         call old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
         call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
         call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)

      end if
end do
      return

   end subroutine Smvgear

      subroutine old150CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
         use GmiManager_mod
         use GmiMechanism_mod
         use GmiSparseMatrix_mod
         use Smv2Chem2_mod
         implicit none

#     include "smv2chem_par.h"
         integer, intent(in) :: ilat
         integer, intent(in) :: ilong
         integer, intent(in) :: itloop
         real*8 :: absoluteErrTolerance(KBLOOP)
         real*8 :: cc2(KBLOOP, 0:MXARRAY)
         real*8 :: cnew(KBLOOP, MXGSAER)
         real*8 :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer :: concAboveAbtolCount(KBLOOP, 5)
         real*8, intent(in) :: corig(KBLOOP, MXGSAER)
         real*8 :: dely(KBLOOP)
         logical, intent(in) :: do_semiss_inchem
         real*8 :: errmx2(itloop)
         integer :: evaluatePredictor
         real*8 :: explic(KBLOOP, MXGSAER)
         real*8 :: gloss(KBLOOP, MXGSAER)
         integer, intent(in) :: inewold(MXGSAER, ICS)
         integer, intent(in) :: ireord
         integer, intent(in) :: jlooplo
         integer, intent(in) :: jphotrat(ICS)
         integer, intent(in) :: jreorder(itloop)
         integer :: kloop
         integer, intent(in) :: ktloop
         integer, intent(in) :: lunsmv
         type(manager_type) :: managerObject
         real*8 :: MAX_REL_CHANGE
         type(mechanism_type) :: mechanismObject
         integer, intent(in) :: ncs
         integer :: ncsp
         integer, intent(out) :: nfdh1
         integer, intent(in) :: ntspec(ICS)
         integer, intent(in) :: numActiveReactants
         integer :: numFinalMatrixPositions
         logical, intent(in) :: pr_smv2
         real*8, intent(in) :: pratk1(KBLOOP, IPHOT)
         logical, intent(in) :: prDiag
         real*8 :: r1delt
         real*8 :: vdiag(KBLOOP, MXGSAER)
         real*8, intent(in) :: yemis(ilat*ilong, IGAS)
         call resetBeforeUpdate (managerObject)

         !!DIR$ INLINE
         call updatePhotoDissRates  (mechanismObject, ktloop, numActiveReactants, ncs, ncsp, jphotrat, pratk1)
         !!DIR$ NOINLINE

         call velocity(mechanismObject, managerObject%num1stOEqnsSolve, ncsp, cnew, gloss, nfdh1)
         managerObject%numCallsVelocity = managerObject%numCallsVelocity + 1
         call setBoundaryConditions (mechanismObject, itloop, jreorder, jlooplo, ilat, &
               ilong, ntspec, ncs, inewold, do_semiss_inchem, gloss, yemis)

         do kloop = 1, ktloop
            dely(kloop) = 0.0d0
         end do

         call determineInitialAbTol(managerObject, cnew, concAboveAbtolCount, ireord, ktloop, ncs, absoluteErrTolerance)
         call calculateErrorTolerances (managerObject, ktloop, cnew, gloss, dely)

         if (ireord /= SOLVE_CHEMISTRY) then
            do kloop = 1, ktloop
               errmx2(jlooplo+kloop) = dely(kloop)
            end do
            return
         end if

         call calcInitialTimeStepSize (managerObject, ktloop, dely, ncs)
         call setInitialOrder (managerObject, evaluatePredictor)
         call storeInitConcAndDerivatives(managerObject%num1stOEqnsSolve, ktloop, &
                              cnewDerivatives, cnew, managerObject%currentTimeStep, gloss)
      end subroutine


      subroutine old200CodeBlock(ilat, ilong, itloop, absoluteErrTolerance, cc2, cnew, cnewDerivatives, concAboveAbtolCount, corig, dely, do_semiss_inchem, errmx2, evaluatePredictor, explic, gloss, inewold, ireord, jlooplo, jphotrat, jreorder, kloop, ktloop, lunsmv, managerObject, MAX_REL_CHANGE, mechanismObject, ncs, ncsp, nfdh1, ntspec, numActiveReactants, numFinalMatrixPositions, pr_smv2, pratk1, prDiag, r1delt, vdiag, yemis)
         use GmiManager_mod
         use GmiMechanism_mod
         use GmiSparseMatrix_mod
         use Smv2Chem2_mod
         implicit none

#     include "smv2chem_par.h"
         integer, intent(in) :: ilat
         integer, intent(in) :: ilong
         integer, intent(in) :: itloop
         real*8 :: absoluteErrTolerance(KBLOOP)
         real*8 :: cc2(KBLOOP, 0:MXARRAY)
         real*8 :: cnew(KBLOOP, MXGSAER)
         real*8 :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer :: concAboveAbtolCount(KBLOOP, 5)
         real*8, intent(in) :: corig(KBLOOP, MXGSAER)
         real*8 :: dely(KBLOOP)
         logical, intent(in) :: do_semiss_inchem
         real*8 :: errmx2(itloop)
         integer :: evaluatePredictor
         real*8 :: explic(KBLOOP, MXGSAER)
         real*8 :: gloss(KBLOOP, MXGSAER)
         integer, intent(in) :: inewold(MXGSAER, ICS)
         integer, intent(in) :: ireord
         integer, intent(in) :: jlooplo
         integer, intent(in) :: jphotrat(ICS)
         integer, intent(in) :: jreorder(itloop)
         integer :: kloop
         integer, intent(in) :: ktloop
         integer, intent(in) :: lunsmv
         type(manager_type) :: managerObject
         real*8 :: MAX_REL_CHANGE
         type(mechanism_type) :: mechanismObject
         integer, intent(in) :: ncs
         integer :: ncsp
         integer, intent(out) :: nfdh1
         integer, intent(in) :: ntspec(ICS)
         integer, intent(in) :: numActiveReactants
         integer :: numFinalMatrixPositions
         logical, intent(in) :: pr_smv2
         real*8, intent(in) :: pratk1(KBLOOP, IPHOT)
         logical, intent(in) :: prDiag
         real*8 :: r1delt
         real*8 :: vdiag(KBLOOP, MXGSAER)
         real*8, intent(in) :: yemis(ilat*ilong, IGAS)

         do while (.true.)

               if (managerObject%orderOfIntegrationMethod /= managerObject%oldOrderOfIntegrationMethod) then
                  call updateCoefficients (managerObject)
               endif


               call calculateTimeStep (managerObject, evaluatePredictor, MAX_REL_CHANGE)
               if (managerObject%currentTimeStep < HMIN) then !HMIN is a minimum acceptanbletime step
                 call tightenErrorTolerance (managerObject, pr_smv2, lunsmv, ncs)

                  !     Start time interval or re-enter after total failure.
                  call startTimeInterval (managerObject, ncs)
                  call initConcentrationArray(ktloop, cnew, corig, managerObject)

                  !     Re-enter here if total failure or if restarting with new cell block.
                 call resetBeforeUpdate (managerObject)

         !!DIR$ INLINE
                  call updatePhotoDissRates  (mechanismObject, ktloop, numActiveReactants, ncs, ncsp, jphotrat, pratk1)
         !!DIR$ NOINLINE

                  call velocity(mechanismObject, managerObject%num1stOEqnsSolve, ncsp, cnew, gloss, nfdh1)
                  managerObject%numCallsVelocity = managerObject%numCallsVelocity + 1
                  call setBoundaryConditions (mechanismObject, itloop, jreorder, jlooplo, ilat, &
                        ilong, ntspec, ncs, inewold, do_semiss_inchem, gloss, yemis)

                  do kloop = 1, ktloop
                     dely(kloop) = 0.0d0
                  end do

                  call determineInitialAbTol(managerObject, cnew, concAboveAbtolCount, ireord, ktloop, ncs, absoluteErrTolerance)
                  call calculateErrorTolerances (managerObject, ktloop, cnew, gloss, dely)

                  if (ireord /= SOLVE_CHEMISTRY) then
                     do kloop = 1, ktloop
                        errmx2(jlooplo+kloop) = dely(kloop)
                     end do
                     return
                  end if

                  call calcInitialTimeStepSize (managerObject, ktloop, dely, ncs)
                  call setInitialOrder (managerObject, evaluatePredictor)
                  call storeInitConcAndDerivatives(managerObject%num1stOEqnsSolve, ktloop, &
                                    cnewDerivatives, cnew, managerObject%currentTimeStep, gloss)

               else
                  if (managerObject%timeStepRatio /= 1.0d0) then
                     call scaleDerivatives (managerObject, ktloop, cnewDerivatives)
                  end if

                  if (managerObject%ifsuccess == 1) then
                     managerObject%rdelmax = 10.0d0

                     if (Mod (managerObject%numSuccessTdt, 3) == 2) then
                        call calcNewAbsoluteErrorTolerance (managerObject, cnew, concAboveAbtolCount, &
                           ktloop, absoluteErrTolerance, ncs)
                     end if

                     call updateChold (managerObject, ktloop, cnew, absoluteErrTolerance)
                  end if

                  call predictConcAndDerivatives (managerObject, cnewDerivatives, explic, ktloop, prDiag)

                  call old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
                  call old300Block(ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)
                  exit

            end if
         end do
      end subroutine old200CodeBlock


      subroutine old250Block(cc2, cnew, cnewDerivatives, evaluatePredictor, ktloop, managerObject, mechanismObject, ncsp, numFinalMatrixPositions, r1delt, vdiag)
         use GmiManager_mod
         use GmiMechanism_mod
         use GmiSparseMatrix_mod
         use Smv2Chem2_mod
         implicit none

#     include "smv2chem_par.h"

         real*8 :: cc2(KBLOOP, 0:MXARRAY)
         real*8 :: cnew(KBLOOP, MXGSAER)
         real*8 :: cnewDerivatives(KBLOOP, MXGSAER*7)

         integer :: evaluatePredictor
         integer, intent(in) :: ktloop
         type(manager_type) :: managerObject
         type(mechanism_type) :: mechanismObject
         integer :: ncsp
         integer :: numFinalMatrixPositions

         real*8 :: r1delt
         real*8 :: vdiag(KBLOOP, MXGSAER)

         call initCorrector (managerObject, ktloop, cnew, cnewDerivatives)

         ! Re-evaluate predictor matrix before starting the corrector iteration.
         if (evaluatePredictor == EVAL_PREDICTOR) then

            r1delt = -managerObject%integrationOrderCoeff1 * managerObject%currentTimeStep
            numFinalMatrixPositions  = sparseMatrixDimension(ncsp) - managerObject%num1stOEqnsSolve

            call calculatePredictor (numFinalMatrixPositions, sparseMatrixDimension(ncsp), &
                  ktloop, cnew, npdhi(ncsp), npdlo(ncsp), r1delt, cc2, &
                  mechanismObject%rateConstants)

            managerObject%numCallsPredict = managerObject%numCallsPredict + 1
            evaluatePredictor  = PREDICTOR_JUST_CALLED

         ! K: Consider un-inlining this.
         !!DIR$   INLINE
           call Decomp  (managerObject%num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag)
         !!DIR$   NOINLINE

           call setConvergenceTerms (managerObject, MBETWEEN)

         end if
      end subroutine old250Block


      subroutine old300Block (ilat, ilong, itloop, cc2, cnew, cnewDerivatives, dely, do_semiss_inchem, gloss, inewold, jlooplo, jreorder, ktloop, managerObject, mechanismObject, ncs, ncsp, nfdh1, ntspec, vdiag, yemis)
         use GmiManager_mod
         use GmiMechanism_mod
         implicit none

#     include "smv2chem_par.h"
         integer, intent(in) :: ilat
         integer, intent(in) :: ilong
         integer, intent(in) :: itloop
         real*8, intent(inout) :: cc2(KBLOOP, 0:MXARRAY)
         real*8, intent(inout) :: cnew(KBLOOP, MXGSAER)
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(inout) :: dely(KBLOOP)
         logical, intent(in) :: do_semiss_inchem
         real*8 :: gloss(KBLOOP, MXGSAER)
         integer, intent(in) :: inewold(MXGSAER, ICS)
         integer, intent(in) :: jlooplo
         integer, intent(in) :: jreorder(itloop)
         integer, intent(in) :: ktloop
         type(manager_type) :: managerObject
         type(mechanism_type) :: mechanismObject
         integer, intent(in) :: ncs
         integer :: ncsp
         integer, intent(out) :: nfdh1
         integer, intent(in) :: ntspec(ICS)
         real*8 :: vdiag(KBLOOP, MXGSAER)
         real*8, intent(in) :: yemis(ilat*ilong, IGAS)

         ! Evaluate the first derivative using corrected values of cnew.
         call velocity (mechanismObject, managerObject%num1stOEqnsSolve, ncsp, cnew, gloss, nfdh1)
         managerObject%numCallsVelocity = managerObject%numCallsVelocity + 1

         call setBoundaryConditions (mechanismObject, itloop, jreorder, jlooplo, ilat, &
                  ilong, ntspec, ncs, inewold, do_semiss_inchem, gloss, yemis)

         call computeErrorFromCorrected1stDeriv (managerObject%num1stOEqnsSolve, ktloop, &
                  gloss, managerObject%currentTimeStep, cnewDerivatives, managerObject%accumulatedError)

         call Backsub (managerObject%num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag, gloss)

         call sumAccumulatedError (managerObject, cnew, cnewDerivatives, dely, gloss, ktloop)
      end subroutine old300Block











      subroutine doMassBalanceAccounting (integrationOrderCoeff, num1stOEqnsSolve, ktloop, amountAddedToEachSpecies, &
                                       &  accumulatedError, explic, cnewDerivatives)

         implicit none
#     include "smv2chem_par.h"

         real*8, intent(in) :: integrationOrderCoeff
         integer, intent(in) :: num1stOEqnsSolve
         integer, intent(in) :: ktloop
         real*8, intent(inout) :: amountAddedToEachSpecies(KBLOOP, MXGSAER)
         real*8, intent(in)  :: accumulatedError (KBLOOP, MXGSAER) ! on a successful return; accumulatedError(kloop,i) contains the estimated one step local error in cnew
         real*8, intent(in) :: explic(KBLOOP, MXGSAER)
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)


         integer :: i, kloop
         real*8 :: dtasn1

         if (integrationOrderCoeff == 1.0d0) then
            do i = 1, num1stOEqnsSolve
              do kloop = 1, ktloop
                amountAddedToEachSpecies(kloop,i) = amountAddedToEachSpecies(kloop,i) + accumulatedError(kloop,i) + explic(kloop,i)
                cnewDerivatives(kloop,i) = cnewDerivatives(kloop,i) + accumulatedError(kloop,i)
              end do
            end do

          else

            do i = 1, num1stOEqnsSolve
              do kloop = 1, ktloop
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
      subroutine computeErrorFromCorrected1stDeriv (num1stOEqnsSolve, ktloop, gloss, &
                  & currentTimeStep, cnewDerivatives, accumulatedError)

         implicit none
#     include "smv2chem_par.h"

         integer, intent(in) :: num1stOEqnsSolve, ktloop
         real*8, intent(inout) :: gloss (KBLOOP, MXGSAER)
         real*8, intent(in) :: currentTimeStep
         real*8, intent(in) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(in)  :: accumulatedError (KBLOOP, MXGSAER)!

         integer jspc, j, kloop

         do jspc = 1, num1stOEqnsSolve
           j = jspc + num1stOEqnsSolve
           do kloop = 1, ktloop
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
      subroutine storeInitConcAndDerivatives(num1stOEqnsSolve, ktloop, cnewDerivatives, cnew, currentTimeStep, gloss)

         implicit none
#     include "smv2chem_par.h"

         integer, intent(in) :: num1stOEqnsSolve
         integer, intent(in) :: ktloop
         real*8, intent(out) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         real*8, intent(in) :: cnew(KBLOOP, MXGSAER)
         real*8, intent(in) :: currentTimeStep
         real*8, intent(in) :: gloss(KBLOOP, MXGSAER)

         ! local variables
         integer :: kloop, jspc, j

         do jspc = 1, num1stOEqnsSolve
            j = jspc + num1stOEqnsSolve

           do kloop = 1, ktloop
             cnewDerivatives(kloop,jspc) = cnew(kloop,jspc)
             cnewDerivatives(kloop,j)    = currentTimeStep * gloss(kloop,jspc)
           end do

         end do
      end subroutine storeInitConcAndDerivatives




      subroutine initConcentrationArray(ktloop, concentrationsNew, concentrationsOld, managerObject)

         use GmiManager_mod
         implicit none

#     include "smv2chem_par.h"

         integer, intent(in) :: ktloop
         real*8, intent(out) :: concentrationsNew(KBLOOP, MXGSAER)
         real*8,  intent(in)  :: concentrationsOld(KBLOOP, MXGSAER)
         type (Manager_type) :: managerObject

         integer :: jnew, kloop

         do jnew = 1, managerObject%num1stOEqnsSolve
           do kloop = 1, ktloop
             concentrationsNew(kloop, jnew) = concentrationsOld(kloop, jnew) ! why save this?
           end do
         end do

      end subroutine





!-----------------------------------------------------------------------------
!
! ROUTINE
!   Backsub
!
! DESCRIPTION
!   This routine performs back-substitutions on the decomposed matrix.  It
!   solves the linear set of equations Ax = B FOR x, the correction vector,
!   where "A" is the L-U decompostion of the original matrix =>
!
!     P = I - H x Bo x J
!
!   I = identity matrix, H = time step, Bo = a coefficient that depends on
!   the order of the integration method, and J is the matrix of partial
!   derivatives.  B is sent from Smvgear as a corrected value of the first
!   derivatives of the ordinary differential equations.  Decomp solved for
!   "A", the decomposed matrix.  See Press, et. al. (1992), Numerical
!   Recipes, Cambridge University Press, for a better description of the
!   back-substitution process.
!
!   This back-substitution process uses sparse matrix techniques,
!   vectorizes around the grid-cell dimension, and uses no partial
!   pivoting.  Tests by Sherman & Hindmarsh (1980), Lawrence Livermore
!   Livermore Laboratory, Rep. UCRL-84102, and by us have confirmed that
!   the removal of partial pivoting has little effect on results.
!
!   Backsub loop # 1 =>
!     First, adjust right side of Ax = B using lower triangular matrix.
!     Sum 1,2,3,4, or 5 terms at a time to improve vectorization.
!
! ARGUMENTS
!   num1stOEqnsSolve : # of first-order eqns to solve, = # of spc = order of original
!            matrix; num1stOEqnsSolve has a different value for day and night, and for
!            gas- and aqueous-phase chemistry;
!            # spc with prod or loss terms in Smvgear (?)
!   ktloop : # of grid-cells in a grid-block
!   ncsp   : ncs       => for daytime   gas chemistry
!            ncs + ICS => for nighttime gas chemistry
!   cc2    : array holding values of decomposed matrix
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!   gloss  : first derivative = sum of prod minus loss rates for a spc
!
!-----------------------------------------------------------------------------

      subroutine Backsub  &
     &  (num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag, gloss)

      use Smv2Chem2_mod
      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: num1stOEqnsSolve
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: ncsp
      real*8,  intent(in)  :: cc2  (KBLOOP, 0:MXARRAY)
      real*8,  intent(in)  :: vdiag(KBLOOP, MXGSAER)

      real*8,  intent(inout) :: gloss(KBLOOP, MXGSAER)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, ij
      integer :: ij0, ij1, ij2, ij3, ij4

      integer :: j0, j1, j2, j3, j4

      integer :: k, kc, kzt
      integer :: kh1, kh2, kh3, kh4, kh5
      integer :: kl1, kl2, kl3, kl4, kl5

      integer :: mc, mzt
      integer :: mh1, mh2, mh3, mh4, mh5
      integer :: ml1, ml2, ml3, ml4, ml5


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Backsub called.'


      ij = 1


!     ==========================================
      KZTLOOP: do kzt = kztlo(ncsp), kzthi(ncsp)
!     ==========================================

        i = ikztot(kzt)

        kl5 = kbl5(kzt)
        kh5 = kbh5(kzt)
        kl4 = kbl4(kzt)
        kh4 = kbh4(kzt)
        kl3 = kbl3(kzt)
        kh3 = kbh3(kzt)
        kl2 = kbl2(kzt)
        kh2 = kbh2(kzt)
        kl1 = kbl1(kzt)
        kh1 = kbh1(kzt)

!       -- Sum 5 terms at a time. --

        do kc = kl5, kh5

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij3 = ij + 3
          ij4 = ij + 4
          ij  = ij + 5

          j0  = kzeroa(kc)
          j1  = kzerob(kc)
          j2  = kzeroc(kc)
          j3  = kzerod(kc)
          j4  = kzeroe(kc)

			 !K: Hot loop
          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2)) -  &
     &        (cc2(k,ij3) * gloss(k,j3)) -  &
     &        (cc2(k,ij4) * gloss(k,j4))
          end do

        end do

!       -- Sum 4 terms at a time. --

        do kc = kl4, kh4

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij3 = ij + 3
          ij  = ij + 4

          j0  = kzeroa(kc)
          j1  = kzerob(kc)
          j2  = kzeroc(kc)
          j3  = kzerod(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2)) -  &
     &        (cc2(k,ij3) * gloss(k,j3))
          end do

        end do

!       -- Sum 3 terms at a time. --

        do kc = kl3, kh3

          ij0 = ij
          ij1 = ij + 1
          ij2 = ij + 2
          ij  = ij + 3

          j0  = kzeroa(kc)
          j1  = kzerob(kc)
          j2  = kzeroc(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1)) -  &
     &        (cc2(k,ij2) * gloss(k,j2))
          end do

        end do

!       -- Sum 2 terms at a time. --

        do kc = kl2, kh2

          ij0 = ij
          ij1 = ij + 1
          ij  = ij + 2

          j0  = kzeroa(kc)
          j1  = kzerob(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0)) -  &
     &        (cc2(k,ij1) * gloss(k,j1))
          end do

        end do

!       -- Sum 1 term at a time. --

        do kc = kl1, kh1

          ij0 = ij
          ij  = ij + 1

          j0  = kzeroa(kc)

          do k = 1, ktloop
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (cc2(k,ij0) * gloss(k,j0))
          end do

        end do

!     ==============
      end do KZTLOOP
!     ==============


!     ---------------------------------------------------------------
!     Backsub loop # 2.
!
!     Backsubstite with upper triangular matrix to find solution.
!     Again, sum up several terms at a time to improve vectorization.
!     ---------------------------------------------------------------

!     ===========================
      ILOOP: do i = num1stOEqnsSolve, 1, -1
!     ===========================

        mzt = imztot(i,ncsp)

!       ===================
        MZTIF: if (mzt > 0) then
!       ===================

          ml5 = mbl5(mzt)
          mh5 = mbh5(mzt)
          ml4 = mbl4(mzt)
          mh4 = mbh4(mzt)
          ml3 = mbl3(mzt)
          mh3 = mbh3(mzt)
          ml2 = mbl2(mzt)
          mh2 = mbh2(mzt)
          ml1 = mbl1(mzt)
          mh1 = mbh1(mzt)

!         -- Sum 5 terms at a time. --

          do mc = ml5, mh5

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij3 = ij + 3
            ij4 = ij + 4
            ij  = ij + 5

            j0  = mzeroa(mc)
            j1  = mzerob(mc)
            j2  = mzeroc(mc)
            j3  = mzerod(mc)
            j4  = mzeroe(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2)) -  &
     &          (cc2(k,ij3) * gloss(k,j3)) -  &
     &          (cc2(k,ij4) * gloss(k,j4))
            end do

          end do

!         -- Sum 4 terms at a time. --

          do mc = ml4, mh4

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij3 = ij + 3
            ij  = ij + 4

            j0  = mzeroa(mc)
            j1  = mzerob(mc)
            j2  = mzeroc(mc)
            j3  = mzerod(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2)) -  &
     &          (cc2(k,ij3) * gloss(k,j3))
            end do

          end do

!         -- Sum 3 terms at a time. --

          do mc = ml3, mh3

            ij0 = ij
            ij1 = ij + 1
            ij2 = ij + 2
            ij  = ij + 3

            j0  = mzeroa(mc)
            j1  = mzerob(mc)
            j2  = mzeroc(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1)) -  &
     &          (cc2(k,ij2) * gloss(k,j2))
            end do

          end do

!         -- Sum 2 terms at a time. --

          do mc = ml2, mh2

            ij0 = ij
            ij1 = ij + 1
            ij  = ij + 2

            j0  = mzeroa(mc)
            j1  = mzerob(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0)) -  &
     &          (cc2(k,ij1) * gloss(k,j1))
            end do

          end do

!         -- Sum 1 term at a time. --

          do mc = ml1, mh1

            ij0 = ij
            ij  = ij + 1

            j0  = mzeroa(mc)

            do k = 1, ktloop
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (cc2(k,ij0) * gloss(k,j0))
            end do

          end do

!       ============
        end if MZTIF
!       ============

!       -- Adjust gloss with diagonal element. --

        do k = 1, ktloop
          gloss(k,i) = gloss(k,i) * vdiag(k,i)
        end do

!     ============
      end do ILOOP
!     ============


      return

      end subroutine Backsub


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Decomp
!
! DESCRIPTION
!   This routine decomposes the sparse matrix "P" into the matrix "A" in
!   order to solve the linear set of equations Ax = B for x, which is a
!   correction vector.  Ax = B is solved in Backsub, the original matrix
!   "P" is =>
!
!     P = I - H x Bo x J
!
!   where I = identity matrix, H = time step, Bo = a coefficient that
!   depends on the order of the integration method, and J is the matrix of
!   partial derivatives.  See Press, et. al. (1992), Numerical Recipes,
!   Cambridge University Press, for a better description of the L-U
!   decompostion process.
!
!   This L-U decompostion process uses sparse matrix techniques, vectorizes
!   around the grid-cell dimension, and uses no partial pivoting.  Tests by
!   Sherman & Hindmarsh (1980), Lawrence Livermore National Laboratory,
!   Rep. UCRL-84102, and by us have confirmed that the removal of partial
!   pivoting has little effect on results.
!
! ARGUMENTS
!   num1stOEqnsSolve : # of first-order eqns to solve, = # of spc = order of original
!            matrix; num1stOEqnsSolve has a different value for day and night, and for
!            gas- and aqueous-phase chemistry;
!            # spc with prod or loss terms in Smvgear (?)
!   ktloop : # of grid-cells in a grid-block
!   ncsp   : ncs       => for daytime   gas chemistry
!            ncs + ICS => for nighttime gas chemistry
!   cc2    : array of sparseMatrixDimension units holding values of each matrix
!            position actually used; originally,
!            cc2 = P = I - currentTimeStep * coeffsForIntegrationOrder(orderOfIntegrationMethod,1) * partial_derivatives;
!            however, cc2 is decomposed here
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!
!-----------------------------------------------------------------------------
! MRD: LU Decomp  - should go into sparseMatrix module
      subroutine Decomp  &
     &  (num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag)

      use Smv2Chem2_mod
      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: num1stOEqnsSolve
      integer, intent(in)  :: ktloop
      integer, intent(in)  :: ncsp

      real*8,  intent(inout) :: cc2  (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: vdiag(KBLOOP, MXGSAER)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: iar, ic
      integer :: ih1, ih2, ih3, ih4, ih5
      integer :: ij, ija, ijt
      integer :: ik0, ik1, ik2, ik3, ik4
      integer :: il1, il2, il3, il4, il5
      integer :: j, jc, jh, jl
      integer :: k
      integer :: kj0, kj1, kj2, kj3, kj4


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Decomp called.'


!     -----------------------------------------------------------
!     First loop of L-U decompostion.
!
!     Sum 1,2,3,4, OR 5 terms at a time to improve vectorization.
!     -----------------------------------------------------------
!     =======================
      JLOOP: do j = 1, num1stOEqnsSolve !num species with reaction, we think
!     =======================

!       ==============================================
        IJTLOOP: do ijt = ijtlo(j,ncsp), ijthi(j,ncsp)
!       ==============================================

         !MRD: all things with 5 terms
         ! should be part of sparse matrix type
          ij  = ijval(ijt)
          il5 = idl5 (ijt)
          ih5 = idh5 (ijt)
          il4 = idl4 (ijt)
          ih4 = idh4 (ijt)
          il3 = idl3 (ijt)
          ih3 = idh3 (ijt)
          il2 = idl2 (ijt)
          ih2 = idh2 (ijt)
          il1 = idl1 (ijt)
          ih1 = idh1 (ijt)

!         -- Sum 5 terms at a time. --
          ! MRD: does this unrolling really help...?
			 !K: Unclear, I'd rather remove it for clarity if nothing else
			 !K: But this involves unifying ikdeca,etc

          ! should be part of sparse matrix type
          do ic = il5, ih5

            ik0 = ikdeca(ic)
            ik1 = ikdecb(ic)
            ik2 = ikdecc(ic)
            ik3 = ikdecd(ic)
            ik4 = ikdece(ic)

            kj0 = kjdeca(ic)
            kj1 = kjdecb(ic)
            kj2 = kjdecc(ic)
            kj3 = kjdecd(ic)
            kj4 = kjdece(ic)
!K: Hot loop
            do k = 1, ktloop
              cc2(k,ij) =  & !ij is nth location of this matrix
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2)) -  &
     &          (cc2(k,ik3) * cc2(k,kj3)) -  &
     &          (cc2(k,ik4) * cc2(k,kj4))
            end do

          end do

!         -- Sum 4 terms at a time. --

          do ic = il4, ih4

            ik0 = ikdeca(ic)
            ik1 = ikdecb(ic)
            ik2 = ikdecc(ic)
            ik3 = ikdecd(ic)

            kj0 = kjdeca(ic)
            kj1 = kjdecb(ic)
            kj2 = kjdecc(ic)
            kj3 = kjdecd(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2)) -  &
     &          (cc2(k,ik3) * cc2(k,kj3))
            end do

          end do

!         -- Sum 3 terms at a time. --

          do ic = il3, ih3

            ik0 = ikdeca(ic)
            ik1 = ikdecb(ic)
            ik2 = ikdecc(ic)

            kj0 = kjdeca(ic)
            kj1 = kjdecb(ic)
            kj2 = kjdecc(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1)) -  &
     &          (cc2(k,ik2) * cc2(k,kj2))
            end do

          end do

!         -- Sum 2 terms at a time. --

          do ic = il2, ih2

            ik0 = ikdeca(ic)
            ik1 = ikdecb(ic)

            kj0 = kjdeca(ic)
            kj1 = kjdecb(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0)) -  &
     &          (cc2(k,ik1) * cc2(k,kj1))
            end do

          end do

!         -- Sum 1 term  at a time. --

          do ic = il1, ih1

            ik0 = ikdeca(ic)

            kj0 = kjdeca(ic)

            do k = 1, ktloop
              cc2(k,ij) =  &
     &          cc2(k,ij) -  &
     &          (cc2(k,ik0) * cc2(k,kj0))
            end do

          end do

!       ==============
        end do IJTLOOP
!       ==============

        iar = diagonalTermDecomp(j,ncsp)

        do k = 1, ktloop
          vdiag(k,j) = 1.0d0 / cc2(k,iar)
        end do

!       ----------------------------
!       Second loop of decompostion.
!       ----------------------------

        jl = jloz1(j,ncsp)
        jh = jhiz1(j,ncsp)

        do jc = jl, jh

          ija = jzeroa(jc)

          do k = 1, ktloop
            cc2(k,ija) = cc2(k,ija) * vdiag(k,j)
          end do

        end do

!     ============
      end do JLOOP
!     ============


      return

      end subroutine Decomp




