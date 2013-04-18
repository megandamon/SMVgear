
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
!   Backsub
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
!   group cells with stiff equations together and those with non-stiff
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
!   nallr    : # of active rxns
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!   nfdh2    : nfdh3 + # of rxns with two   active reactants
!   nfdh3    :         # of rxns with three active reactants
!   nfdl1    : nfdh2 + 1
!   nfdl2    : nfdh3 + 1
!   nfdrep   : nfdh3 + # of rxns with two active reactants that are not
!              followed by a rxn with the same reactants
!   nfdrep1  : nfdrep + 1 ! MRD: No longer needed. Remove it.
!   fracdec  : fraction time step is decreased in Smvgear if convergence
!              test fails
!   hmaxnit  : max time step for night all chem (s)
!   pr_nc_period     : NetCDF output period
!   tdt      : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   irma,b,c : spc # of each reactant; locates reordered active spc #s
!   jreorder : gives original grid-cell from re-ordered grid-cell
!   jphotrat : tbd
!   ntspec   : # of active + inactive gases
!   inewold  : original spc # of each new jnew spc
!   denair   : density of air (molec/cm^3)
!   corig    : original gas-phase concs used to restart Smvgear if a failure
!              occurs (molec/cm^3)
!   pratk1   : tbd
!   yemis    : surface emissions (units?)
!   smvdm    : amount added to each spc at each grid-cell, for mass balance
!              accounting (# cm^-3 for gas chemistry (?))
!   nfdh1    : nfdh2 + # of rxns with one   active reactant
!   errmx2   : tbd
!   cc2      : array holding values of decomposed matrix
!   cnew     : stores conc (y (estimated)) (molec/cm^3)
!   gloss    : value of first derivatives on output from velocity; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!   vdiag    : 1 / current diagonal term of the decomposed matrix
!   rrate    : rate constants
!   trate    : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!   urate    : term of Jacobian (J) = partial derivative
!
!-----------------------------------------------------------------------------

      subroutine Smvgear  &
     &  (do_qqjk_inchem, do_semiss_inchem, pr_qqjk, pr_smv2, ifsun,  &
     &   ilat, ilong, ivert, ireord, itloop, jlooplo, ktloop, lunsmv,  &
     &   nallr, ncs, nfdh2, nfdh3, nfdl1, nfdl2, nfdrep, nfdrep1,  &
     &   fracdec, hmaxnit, pr_nc_period, tdt, do_cell_chem, irma, irmb,  &
     &   irmc, jreorder, jphotrat, ntspec, inewold, denair, corig,  &
     &   pratk1, yemis, smvdm, nfdh1, errmx2, cc2, cnew, gloss, vdiag,  &
     &   rrate, trate, urate, &
     &   yda, qqkda, qqjda, qkgmi, qjgmi, &
     &   CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
     &   num_qjo, num_qks, num_qjs, num_active)

      use GmiPrintError_mod, only : GmiPrintError
      use GmiMechanism_mod
      use GmiSparseMatrix_mod


      implicit none

#     include "smv2chem_par.h"
#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

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
      integer, intent(in)  :: nallr
      integer, intent(in)  :: ncs
      integer, intent(in)  :: nfdh2,  nfdh3
      integer, intent(in)  :: nfdl1,  nfdl2
      integer, intent(in)  :: nfdrep, nfdrep1
      real*8,  intent(in)  :: fracdec
      real*8,  intent(in)  :: hmaxnit
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      logical, intent(in)  :: do_cell_chem(ilong, ilat, ivert)
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
      real*8,  intent(inout) :: cc2   (KBLOOP, 0:MXARRAY)
      real*8,  intent(inout) :: cnew  (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: gloss (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: smvdm (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: vdiag (KBLOOP, MXGSAER)
      real*8,  intent(inout) :: rrate (KBLOOP, NMTRATE)
      real*8,  intent(inout) :: trate (KBLOOP, NMTRATE*2)
      real*8,  intent(inout) :: urate (KBLOOP, NMTRATE, 3)

      integer, intent(out) :: nfdh1


!     ----------------------
!     Variable declarations.
!     ----------------------

!     ------------------------------------------------------------------------
!     idoub     : records # of steps since the last change in step size or
!                 order; it must be at least kstep = nqq+1 before doubling is
!                 allowed
!
!     numFailOldJacobian     : # of times corrector failed to converge while the Jacobian
!                                was old
!     numFailErrorTest     : # of times accumulated error test failed
!     numFailAfterPredict     : # of times correcter failed to converge after Pderiv was
!                 called
!
!     ifsuccess : identifies whether step is successful (=1) or not (=0)
!     num1stOEqnsSolve    : # of first-order eqns to solve, = # of spc = order of
!                 original matrix; num1stOEqnsSolve has a different value for day and
!                 night and for gas- and aqueous-phase chemistry;
!                 # spc with prod or loss terms in Smvgear (?)
!     jeval     :  1 => call Pderiv the next time through the corrector steps;
!                  0 => last step successful and do not need to call Pderiv;
!                 -1 => Pderiv just called, and do not need to call again
!                  until jeval switched to 1
!     jrestar   : counts # of times Smvgear starts over at order 1 because of
!                 excessive failures
!     kstep     : nqq + 1
!
!     numCallsPredict   : total # of times predictor is called
!     numCallsVelocity   : total # of times velocity is called
!
!     nqqisc    : nqq * num1stOEqnsSolve
!     nqqold    : value of nqq during last time step
!     nqq       : order of integration method; varies between 1 and MAXORD
!     nslp      : last time step # during which Pderiv was called
!     numSuccessTdt    : total # of successful time steps taken
!     ------------------------------------------------------------------------

      integer :: i, j, k
      integer :: i1, i2
      integer :: idoub
      integer :: numFailOldJacobian, jfail, numFailErrorTest, numFailAfterPredict
      integer :: ifsuccess
      integer :: num1stOEqnsSolve
      integer :: jb
      integer :: jeval
      integer :: jg1
      integer :: jgas
      integer :: jnew
      integer :: jrestar
      integer :: jspc
      integer :: k1, k2, k3, k4, k5
      integer :: kloop
      integer :: kstep, kstepisc
      integer :: l3
      integer :: nact
      integer :: ncsp  ! ncs       => for daytime   gas chemistry
                       ! ncs + ICS => for nighttime gas chemistry
      integer :: numCallsPredict, numCallsVelocity
      integer :: nqisc, nqqisc, nqqold
      integer :: nqq
      integer :: nslp
      integer :: numSuccessTdt
      integer :: numErrTolDecreases

      integer :: ibcb(IGAS)

!     ------------------------------------------------
!     kgrp : counts # of concs above abtol(i), i = 1..
!     ------------------------------------------------

      integer :: kgrp(KBLOOP, 5)

!     ------------------------------------------------------------------------
!     asn1      : value of aset(nqq,1)
!     delt      : current time step (s)
!     drate     : parameter which is used to determine whether convergence
!                 has occurred
!     edwn      : pertst^2*order for one order lower  than current order
!     enqq      : pertst^2*order for current order
!     eup       : pertst^2*order for one order higher than current order
!     maxTimeStep      : max time step at a given time (s)
!     hratio    : relative change in delt*aset(1) each change in step or order
!                 when Abs(hratio-1) > MAX_REL_CHANGE, reset jeval = 1 to call Pderiv
!     MAX_REL_CHANGE     : max relative change in delt*aset(1) before Pderiv is called
!     order     : floating point value of num1stOEqnsSolve, the order of # of ODEs
!     rdelmax   : max factor by which delt can be increased in a single step;
!                 as in Lsodes, set it to 1d4 initially to compensate for the
!                 small initial delt, but then set it to 10 after successful
!                 steps and to 2 after unsuccessful steps
!     rdelt     : factor (time step ratio) by which delt is increased or
!                 decreased
!     rdeltdn   : time step ratio at one order lower  than current order
!     rdeltsm   : time step ratio at current order
!     rdeltup   : time step ratio at one order higher than current order
!     rmsrat    : ratio of current to previous rms scaled error; if this
!                 ratio decreases, then convergence is occuring
!     timremain : remaining time in an chem interval (s)
!     chemTimeInterval : total chem time interval; same as chemintv (s)
!     told      : stores last value of xelaps in case current step fails
!     xelaps    : elapsed time in chem interval (s)
!     failureFraction      : = 1 originially, but is decreased if excessive failures
!                 occur in order to reduce absolute error tolerance
!     ------------------------------------------------------------------------

      real*8  :: abtoler1
      real*8  :: asn1, asnqqj
      real*8  :: cnewylow
      real*8  :: cnw
      real*8  :: conp1, conp2, conp3
      real*8  :: consmult
      real*8  :: dcon
      real*8  :: delt, delt1
      real*8  :: der1max, der2max, der3max
      real*8  :: drate
      real*8  :: dtasn1
      real*8  :: edwn, enqq, eup
      real*8  :: initialError, initialError_inv
      real*8  :: errmax_ncs_inv, errymax
      real*8  :: maxTimeStep
      real*8  :: hmtim
      real*8  :: hratio
      real*8  :: iabove
      real*8  :: order, order_inv
      real*8  :: r1delt, rdelmax, rdelt, rdelta
      real*8  :: rdeltdn, rdeltsm, rdeltup
      real*8  :: real_kstep
      real*8  :: reltol1, reltol2, reltol3
      real*8  :: rmsError, rmsErrorPrevious, rmsrat, rmstop
      real*8  :: timremain, chemTimeInterval
      real*8  :: told
      real*8  :: xelaps, xtimestep
      real*8  :: failureFraction

      real*8, parameter  :: MAX_REL_CHANGE = 0.3d0

!     -------------------------------------------------------------------------
!     dely   : tbd
!     yabst  : absolute error tolerance (molec/cm^-3 for gases)
!
!     cest   : stores value of dtlos when idoub = 1
!     chold  : 1 / (reltol * cnew + abtol); multiply chold by local errors in
!              different error tests
!     dtlos  : an array of length num1stOEqnsSolve, used for the accumulated corrections;
!              on a successful return; dtlos(kloop,i) contains the estimated
!              one step local error in cnew
!     explic : tbd
!     conc   : an array of length num1stOEqnsSolve*(MAXORD+1) that carries the
!              derivatives of cnew, scaled by delt^j/factorial(j), where j is
!              the jth derivative; j varies from 1 to nqq; e.g., conc(jspc,2)
!              stores delt*y' (estimated)
!     -------------------------------------------------------------------------

      real*8  :: dely  (KBLOOP)
      real*8  :: yabst (KBLOOP)

      real*8  :: cest  (KBLOOP, MXGSAER)
      real*8  :: chold (KBLOOP, MXGSAER)
      real*8  :: dtlos (KBLOOP, MXGSAER)
      real*8  :: explic(KBLOOP, MXGSAER)

      real*8  :: conc  (KBLOOP, MXGSAER*7)

      type (Mechanism_type) :: mechanismObject
      integer :: nondiag     ! # of final matrix positions, excluding diagonal

      mechanismObject%numGridCellsInBlock = ktloop
      mechanismObject%speciesNumberA = irma
      mechanismObject%speciesNumberB = irmb
      mechanismObject%speciesNumberC = irmc
      mechanismObject%numRxns2 = nfdh2
      mechanismObject%numRxns3 = nfdh3
      mechanismObject%numRxns3Drep = nfdrep

! routine start initializeGear
!     ----------------
!     Begin execution.
!     ----------------

!    Write (6,*) 'Smvgear called.'

      nact = nnact ! MRD: why is this here? removing it affects results!

!     =======================
#     include "setkin_ibcb.h"
!     =======================

      !MRD: local variables initializations
      print*, "initializing things"
      numFailOldJacobian     = 0
      jfail     = 0
      numFailErrorTest     = 0
      numFailAfterPredict     = 0
      numCallsPredict   = 0
      numSuccessTdt    = 0
      numCallsVelocity   = 0
      numErrTolDecreases  = 0

      ! MRD: local variables
      ! MAX_REL_CHANGE moved to parameter above
      rmsError    = 1.0d0

      ! MRD: derived from common block
      num1stOEqnsSolve    = ischang(ncs)
      order     = num1stOEqnsSolve
      order_inv = 1.0d0 / num1stOEqnsSolve

      ! MRD: derived from common bloc
      chemTimeInterval = timeintv(ncs)

      ! ifsun is an argument to Smvgear
      ! ncs is argument to Smvgear
      ! ICS is from common block
      ncsp      = (ifsun - 1) * ICS + ncs

      maxTimeStep = hmaxnit
      if (ifsun == 1) maxTimeStep = hmaxday(ncs)

      failureFraction   = 1.0d0
      iabove = order * 0.4d0

      initialError     = Min (errmax(ncs), 1.0d-03)
      initialError_inv = 1.0d0 / initialError

      errmax_ncs_inv = 1.0d0 / errmax(ncs)
! routine end initializeGear


! routine start startTimeInterval must call: restartTimeInterval,
   !calcInitialTimeStepSize, updateLimitTighten, and itself (if/else)
! 100 calls 150
!     ========
 100  continue
!     ========
!     ----------------------------------------------------
!     Start time interval or re-enter after total failure.
!     ----------------------------------------------------
      print*, "in 100"


      idoub     = 2
      nslp      = MBETWEEN
      jrestar   = 0
      xelaps    = 0.0d0
      told      = 0.0d0
      timremain = chemTimeInterval

      reltol1   = failureFraction * initialError_inv

      reltol2   = failureFraction * errmax_ncs_inv

      reltol3   = errmax_ncs_inv

      abtoler1  = abtol(6,ncs) * reltol1

!     -------------------------------
!     Initialize concentration array.
!     -------------------------------

      do jnew = 1, num1stOEqnsSolve
        do kloop = 1, ktloop
          cnew(kloop, jnew) = corig(kloop, jnew)
        end do
      end do

      ! then call 150
! routine end startTimeInterval


! routine start restartTimeInterval
!     --------------------------------------------------------------------
!     Re-enter here if total failure or if restarting with new cell block.
!     --------------------------------------------------------------------
! 150 resets some stuff, then calls update
!     ========
 150  continue
!     ========

      print*, "in 150"


      hratio    = 0.0d0
      asn1      = 1.0d0
      ifsuccess = 1
      rdelmax   = 1.0d4

!     ---------------------
!     Initialize photrates.
!     ---------------------

!DIR$ INLINE
! MRD: Update will be moved outside the solver
! It will exist in the GMI driver
! Some GMI chem wrapper will initalize this, and then called the solver
! Move this routine into the mechanism

print*, "Calling update"

!     ===========
      call Update  &
!     ===========
     &  (ktloop, nallr, ncs, ncsp, jphotrat, pratk1, rrate, trate)
!DIR$ NOINLINE


      mechanismObject%rateConstants = rrate
      mechanismObject%numActiveReactants = nallr

   print*, "calling velocity"
      call velocity (mechanismObject, num1stOEqnsSolve, ncsp, cnew, gloss, trate, nfdh1)

      numCallsVelocity = numCallsVelocity + 1

      ! MRD: can this be removed?
      mechanismObject%rateConstants = rrate


!     ----------------------------------------------------------------
!     Zero first derviatives in surface zones for species with fixed
!     concentration boundary conditions (LLNL addition, PSC, 5/14/99).
!     ----------------------------------------------------------------
! MRD: these go into a combo of mechanism and sparse matrix
! but probably in mechanism
print*, " Zero first derviatives"
      do kloop = 1, ktloop

        if (jreorder(jlooplo+kloop) <= (ilat*ilong)) then

          do jspc = 1, ntspec(ncs)

            jgas = inewold(jspc,1)

            if (ibcb(jgas) == 1) then
              gloss(kloop,jspc) = 0.0d0
            else if (do_semiss_inchem) then
              gloss(kloop,jspc) =  &
     &          gloss(kloop,jspc) +  &
     &          yemis(jreorder(jlooplo+kloop),jgas)
            end if

          end do

        end if

      end do

! MRD: Take the reordering and error/tolerance calculations and keep them in the solver for now
!     -------------------------------------------
!     Determine initial absolute error tolerance.
!     -------------------------------------------

      do kloop = 1, ktloop
        dely(kloop) = 0.0d0
      end do

   print*, "determining if ireordif: ", ireord
!     ==========================
      IREORDIF: if (ireord /= 1) then
!     ==========================

        do k = 1, 5
          do kloop = 1, ktloop
            kgrp(kloop,k) = 0
          end do
        end do

        do jspc = 1, num1stOEqnsSolve
          do kloop = 1, ktloop

            cnw = cnew(kloop,jspc)

            if (cnw > abtol(1,ncs)) then
              kgrp(kloop,1) = kgrp(kloop,1) + 1
            else if (cnw > abtol(2,ncs)) then
              kgrp(kloop,2) = kgrp(kloop,2) + 1
            else if (cnw > abtol(3,ncs)) then
              kgrp(kloop,3) = kgrp(kloop,3) + 1
            else if (cnw > abtol(4,ncs)) then
              kgrp(kloop,4) = kgrp(kloop,4) + 1
            else if (cnw > abtol(5,ncs)) then
              kgrp(kloop,5) = kgrp(kloop,5) + 1
            end if

          end do
        end do

        do kloop = 1, ktloop

          k1 = kgrp(kloop,1)
          k2 = kgrp(kloop,2) + k1
          k3 = kgrp(kloop,3) + k2
          k4 = kgrp(kloop,4) + k3
          k5 = kgrp(kloop,5) + k4

          if (k1 > iabove) then
            yabst(kloop) = abtol(1,ncs) ! MRD: these yabst should be passed in
          else if (k2 > iabove) then    ! does the driver pass them in?
            yabst(kloop) = abtol(2,ncs) ! or does the mechanism specify them
          else if (k3 > iabove) then    ! tabled for now
            yabst(kloop) = abtol(3,ncs)
          else if (k4 > iabove) then
            yabst(kloop) = abtol(4,ncs)
          else if (k5 > iabove) then
            yabst(kloop) = abtol(5,ncs)
          else
            yabst(kloop) = abtol(6,ncs)
          end if

        end do

!c
        do kloop = 1, ktloop
          do jspc = 1, num1stOEqnsSolve
            cnewylow    = cnew (kloop,jspc) + (yabst(kloop) * reltol1)
            errymax     = gloss(kloop,jspc) / cnewylow
            dely(kloop) = dely (kloop) + (errymax * errymax) ! this is an error (not a tolerance)
          end do
        end do


!     ====
      else
!     ====
   print*, "I am in else because ireord is not equal to 1", ireord
!       ------------------------------------------------------
!       Use lowest absolute error tolerance when reordering.
!       If reordering, set errmx2 then return to Physproc.
!
!       abtoler1 = failureFraction * abtol(6,ncs) / Min (errmax, 1.0d-03)
!       ------------------------------------------------------
! MRD: it is conceviable that there could be different norms
! MRD: or error criteria
! MRD: this is part of the management system
!c
        do kloop = 1, ktloop
          do jspc = 1, num1stOEqnsSolve
            errymax     = gloss(kloop,jspc) /  &
     &                    (cnew(kloop,jspc) + abtoler1)
            dely(kloop) = dely(kloop) + (errymax * errymax)
          end do
        end do

        do kloop = 1, ktloop
          errmx2(jlooplo+kloop) = dely(kloop)
        end do
   print*, "leaving smvgear"
        return



!     ===============
      end if IREORDIF
!     ===============

! routine end restartTimeInterval




! routine start calcInitialTimeStepSize
!     --------------------------------------
!     Calculate initial time step size (s).
!
!     Sqrt (dely / [initialError * order]) =
!       rmsnorm of error scaled to initialError *
!       cnew + abtol / reltol
!     --------------------------------------
! MRD: leave in solver
! MRD: adaptivity manager (time step and order)
! MRD: this is a guess, later it will adapt
! MRD: client will call this

   print*, "doing time step stuff"
      rmstop = 0.0d0

      do kloop = 1, ktloop
        if (dely(kloop) > rmstop) then
          rmstop = dely(kloop)
        end if
      end do

      delt1 = Sqrt (initialError / (abst2(ncs) + (rmstop * order_inv)))
      delt  = Max  (Min (delt1, timremain, maxTimeStep), HMIN)


!     -----------------------
!     Set initial order to 1.
!     -----------------------

      nqqold = 0
      nqq    = 1
      jeval  = 1
      rdelt  = 1.0d0


!     --------------------------------------------------------------
!     Store initial concentration and first derivatives x time step.
!     --------------------------------------------------------------

      do jspc = 1, num1stOEqnsSolve

        j = jspc + num1stOEqnsSolve

        do kloop = 1, ktloop
          conc(kloop,jspc) = cnew(kloop,jspc)
          conc(kloop,j)    = delt * gloss(kloop,jspc)
        end do

      end do
! routine start calcInitialTimeStepSize


! routine start updateLimitTighten
!     ========
 200  continue
!     ========

print*, "in continue 200"


!     -------------------------------------------------------------------
!     Update coefficients of (for?) the order; note that pertst2 is the original
!     pertst^2.
!     -------------------------------------------------------------------
! MRD: Gear can be 1st, 2nd, etc. order
! MRD: this could be where it changes it order
! MRD: track these down, figure out what they are
      if (nqq /= nqqold) then

        nqqold = nqq
        kstep  = nqq + 1
        hratio = hratio * aset(nqq,1) / asn1
        asn1   = aset(nqq,1)
        enqq   = pertst2(nqq,1) * order
        eup    = pertst2(nqq,2) * order
        edwn   = pertst2(nqq,3) * order
        conp3  = 1.4d0 /  (eup**enqq3(nqq))
        conp2  = 1.2d0 / (enqq**enqq2(nqq))
        conp1  = 1.3d0 / (edwn**enqq1(nqq))
        nqqisc = nqq * num1stOEqnsSolve

      end if

! MRD: Below should go into adaptivity manager
!     ----------------------------------------------------------------
!     Limit size of rdelt, then recalculate new time step and update
!     hratio.  Use hratio to determine whether or not the predictor
!     should be updated. (edit by MRD on 2/27/2013)
!     ----------------------------------------------------------------

      hmtim  = Min (maxTimeStep, timremain)
      rdelt  = Min (rdelt, rdelmax, hmtim/delt)
      delt   = delt   * rdelt
      hratio = hratio * rdelt
      xelaps = xelaps + delt

      if ((Abs (hratio-1.0d0) > MAX_REL_CHANGE) .or. (numSuccessTdt >= nslp)) then
        jeval = 1 ! MRD: could be a boolean; this is signifying to whether or not to update Jacobian
      end if


!     ----------------------------------------------------------
!     If time step < HMIN, tighten absoloute error tolerance and
!     restart integration at beginning of time interval.
!     ----------------------------------------------------------
   print*, "tightning absolute error tolerance"
      if (delt < HMIN) then

        if (pr_smv2) then
          Write (lunsmv,950) delt, timremain, failureFraction, errmax(ncs)
        end if

 950    format ('Smvgear:  delt      = ', 1pe9.3, /,  &
     &          '          timremain = ', 1pe9.3, /,  &
     &          '          failureFraction      = ', 1pe9.3, /,  &
     &          '          errmax    = ', 1pe9.3)

        numErrTolDecreases = numErrTolDecreases + 1
        failureFraction     = failureFraction * 0.01d0

        if (numErrTolDecreases == 10) then

          if (pr_smv2) then
            Write (lunsmv,960)
          end if

 960      format ('Smvgear:  too many decreases of failureFraction.')

        call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

        end if


!       =========
        go to 100 ! routine start startTimeInterval
!       =========

      end if

    ! routine end updateLimitTighten


! routine start scalingDerivatives
!     -------------------------------------------------------------------
!     If the delt is different than during the last step (if rdelt /= 1),
!     then scale the derivatives.
!     -------------------------------------------------------------------
   print*, "scaling derivatives"
      if (rdelt /= 1.0d0) then

        rdelta = 1.0d0
        i1     = 1

        do j = 2, kstep

          rdelta = rdelta * rdelt
          i1     = i1 + num1stOEqnsSolve

          do i = i1, i1 + (num1stOEqnsSolve-1)
            do kloop = 1, ktloop
              conc(kloop,i) = conc(kloop,i) * rdelta
            end do
          end do

        end do

      end if

! routine end scalingDerivatives



! routine start resetDelMaxUpdateWithCnew

!     --------------------------------------------------------------
!     If the last step was successful, reset rdelmax = 10 and update
!     the chold array with current values of cnew.
!     --------------------------------------------------------------
   print*, "last time step was successful"
!     ================================
      IFSUCCESSIF: if (ifsuccess == 1) then
!     ================================

        rdelmax = 10.0d0

!       ---------------------------------------
!       Determine new absolute error tolerance.
!       ---------------------------------------

        if (Mod (numSuccessTdt, 3) == 2) then

          do k = 1, 5
            do kloop = 1, ktloop
              kgrp(kloop,k) = 0
            end do
          end do

          do jspc = 1, num1stOEqnsSolve
            do kloop = 1, ktloop

              cnw = cnew(kloop,jspc)

              if (cnw > abtol(1,ncs)) then
                kgrp(kloop,1) = kgrp(kloop,1) + 1
              else if (cnw > abtol(2,ncs)) then
                kgrp(kloop,2) = kgrp(kloop,2) + 1
              else if (cnw > abtol(3,ncs)) then
                kgrp(kloop,3) = kgrp(kloop,3) + 1
              else if (cnw > abtol(4,ncs)) then
                kgrp(kloop,4) = kgrp(kloop,4) + 1
              else if (cnw > abtol(5,ncs)) then
                kgrp(kloop,5) = kgrp(kloop,5) + 1
              end if

            end do
          end do

          do kloop = 1, ktloop

            k1 = kgrp(kloop,1)
            k2 = kgrp(kloop,2) + k1
            k3 = kgrp(kloop,3) + k2
            k4 = kgrp(kloop,4) + k3
            k5 = kgrp(kloop,5) + k4

            if (k1 > iabove) then
              yabst(kloop) = abtol(1,ncs)
            else if (k2 > iabove) then
              yabst(kloop) = abtol(2,ncs)
            else if (k3 > iabove) then
              yabst(kloop) = abtol(3,ncs)
            else if (k4 > iabove) then
              yabst(kloop) = abtol(4,ncs)
            else if (k5 > iabove) then
              yabst(kloop) = abtol(5,ncs)
            else
              yabst(kloop) = abtol(6,ncs)
            end if

          end do

        end if

!c
        do kloop = 1, ktloop
          do jspc = 1, num1stOEqnsSolve

            chold(kloop,jspc) =  &
     &        reltol3 /  &
     &        (Max (cnew(kloop,jspc), 0.0d0) +  &
     &         (yabst(kloop) * reltol2))

          end do
        end do

!     ==================
      end if IFSUCCESSIF
!     ==================

! routine end resetDelMaxUpdateWithCnew


   ! routine start computePredictConcPascal
!     ------------------------------------------------------------------
!     Compute the predicted concentration and derivatives by multiplying
!     previous values by the pascal triangle matrix.
!     ------------------------------------------------------------------
   print*, "computing predicted conc and derivatives using pascal triangle matrix"
      i1 = nqqisc + 1

      do jb = 1, nqq - 1

        i1 = i1 - num1stOEqnsSolve

        do i = i1,  nqqisc

          j = i + num1stOEqnsSolve

          do kloop = 1, ktloop
            conc(kloop,i)  = conc(kloop,i) + conc(kloop,j)
          end do

        end do

      end do

      do jspc = 1,  num1stOEqnsSolve

        j = jspc + num1stOEqnsSolve

        do kloop = 1, ktloop
          conc  (kloop,jspc) = conc(kloop,jspc) + conc(kloop,j)
          explic(kloop,jspc) = conc(kloop,j)
        end do

      end do

      do i = num1stOEqnsSolve + 1, nqqisc

        j = i + num1stOEqnsSolve

        do kloop = 1, ktloop
          conc(kloop,i) = conc(kloop,i) + conc(kloop,j)
        end do

      end do

   ! routine end computePredictConcPascal
!     ---------------------------------------



   ! routine start correctionLoop
!     ---------------------------------------

!     -------------------------------------------------------------------
!     Correction loop.
!
!     Take up to 3 corrector iterations.  Test convergence by requiring
!     that changes be less than the rms norm weighted by chold.
!     Accumulate the correction in the array dtlos.  It equals the
!     jth derivative of concentration multiplied by delt^kstep /
!     (factorial(kstep-1) * aset(kstep)); thus, it is proportional to the
!     actual errors to the lowest power of delt present (delt^kstep).
!     -------------------------------------------------------------------


!     ========
 250  continue ! correctionLoop
!     ========
   print*, "in 250 or correction loops"

      l3 = 0

      do jspc = 1, num1stOEqnsSolve
        do kloop = 1, ktloop
          cnew (kloop,jspc) = conc(kloop,jspc)
          dtlos(kloop,jspc) = 0.0d0
        end do
      end do
   ! routine stop correctionLoop


   ! routine start reEvalPredictor
!     ---------------------------------------

!     ------------------------------------------------------------------
!     If jeval = 1, re-evaluate predictor matrix P = I - H * aset(1) * J
!     before starting the corrector iteration.  After calling Pderiv,
!     set jeval = -1 to prevent recalling Pderiv unless necessary later.
!     Call Decomp to decompose the matrix.
!     ------------------------------------------------------------------

      if (jeval == 1) then
         print*, "re-evalulate predictor matrix"
        r1delt = -asn1 * delt

         ! MRD: The below functionality has replaced the subroutine Pderiv
         !       ===========
         !        call Pderiv  &
         !       ===========
         !         (mechanismObject, num1stOEqnsSolve, ncsp, nfdh2, nfdh3, nfdl1, nfdl2, irma, irmb,  &
         !     &   irmc, r1delt, cnew, rrate, npderiv, cc2, urate, nfdh1)


         nondiag  = iarray(ncsp) - num1stOEqnsSolve ! iarray is in common block
         mechanismObject%numRxns1 = nfdh2 + ioner(ncsp)

         call calculateTermOfJacobian (mechanismObject, cnew, urate)

         numCallsPredict  = numCallsPredict + 1
         ! MRD: iarray, npdhi, and npdlo are in common blocks
         call calculatePredictor (nondiag, iarray(ncsp), mechanismObject%numGridCellsInBlock, &
            &  npdhi(ncsp), npdlo(ncsp), r1delt, urate, cc2)

      ! MRD: End block of code that was in Pderiv

!DIR$   INLINE
!       ===========
        call Decomp  &
!       ===========
     &    (num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag)
!DIR$   NOINLINE

        jeval  = -1
        hratio = 1.0d0
        nslp   = numSuccessTdt + MBETWEEN
        drate  = 0.7d0

      end if

   ! routine end reEvalPredictor
!     ---------------------------------------



   ! routine start evalFirstDerivative
!     ---------------------------------------

!     -------------------------------------------------------------
!     Evaluate the first derivative using corrected values of cnew.
!     -------------------------------------------------------------

!     ========
 300  continue
!     ========

   print*, "in 300, evaluating first derivative"

     call velocity (mechanismObject, num1stOEqnsSolve, ncsp, cnew, gloss, trate, nfdh1)
     numCallsVelocity = numCallsVelocity + 1

!     ----------------------------------------------------------------
!     Zero first derviatives in surface zones for species with fixed
!     concentration boundary conditions.  Include surface emissions in
!     first derviatives.
!
!     Species with non-zero fluxes (LLNL addition, PSC, 7/24/96).
!     ----------------------------------------------------------------

      do kloop = 1, ktloop

        if (jreorder(jlooplo+kloop) <= (ilat*ilong)) then

          do jspc = 1, ntspec(ncs)

            jgas = inewold(jspc,1)

            if (ibcb(jgas) == 1) then

              gloss(kloop,jspc) = 0.0d0

            else

              if (do_semiss_inchem) then
                gloss(kloop,jspc) = gloss(kloop,jspc) +  &
     &                              yemis(jreorder(jlooplo+kloop),jgas)
              end if

            end if

          end do

        end if

      end do


!     ---------------------------------------------------------------
!     In the case of the chord method, compute error (gloss) from the
!     corrected calculation of the first derivative.
!     ---------------------------------------------------------------

      do jspc = 1, num1stOEqnsSolve

        j = jspc + num1stOEqnsSolve

        do kloop = 1, ktloop
          gloss(kloop,jspc) = (delt * gloss(kloop,jspc)) -  &
     &                        (conc(kloop,j) + dtlos(kloop,jspc))
        end do

      end do


   ! routine end evalFirstDerivative
!     ---------------------------------------

!     --------------------------------------------------------------
!     Solve the linear system of equations with the corrector error;
!     Backsub solves backsubstitution over matrix of partial derivs.
!     --------------------------------------------------------------
   print*, "solve system of equations"
!     ============
      call Backsub  &
!     ============
     &  (num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag, gloss)



   ! routine start sumAccumError
!     ---------------------------------------
!     ----------------------------------------------------------------
!     Sum up the accumulated error, correct the concentration with the
!     error, and begin to calculate the rmsnorm of the error relative
!     to chold.
!     ----------------------------------------------------------------
      print*, "summing up aculated error"
      do kloop = 1, ktloop
        dely(kloop) = 0.0d0
      end do

      if (asn1 == 1.0d0) then

        do i = 1, num1stOEqnsSolve
          do kloop = 1, ktloop

            dtlos(kloop,i) = dtlos(kloop,i) + gloss(kloop,i)
            cnew(kloop,i)  = conc(kloop,i)  + dtlos(kloop,i)
            errymax        = gloss(kloop,i) * chold(kloop,i)
            dely(kloop)    = dely(kloop)    + (errymax * errymax)

          end do
        end do

      else

        do i = 1, num1stOEqnsSolve
          do kloop = 1, ktloop

            dtlos(kloop,i) = dtlos(kloop,i) + gloss(kloop,i)
            cnew(kloop,i)  = conc(kloop,i)  + (asn1 * dtlos(kloop,i))
            errymax        = gloss(kloop,i) * chold(kloop,i)
            dely(kloop)    = dely(kloop)    + (errymax * errymax)

          end do
        end do

      end if


   ! routine end sumAccumError
!     ---------------------------------------


   ! routine start calculateRmsError
!     ---------------------------------------

!     ------------------------------------------------------------------
!     Set the previous rms error and calculate the new rms error.
!     If dcon < 1, then sufficient convergence has occurred.  Otherwise,
!     if the ratio of the current to previous rmsError is decreasing,
!     iterate more.  If it is not, then the convergence test failed.
!     ------------------------------------------------------------------

      rmsErrorPrevious = rmsError
      der2max = 0.0d0

      print*, "calculating new rms error"
      do kloop = 1, ktloop

        if (dely(kloop) > der2max) then
          der2max = dely(kloop)
        end if

      end do


      rmsError = Sqrt (der2max * order_inv)


      l3 = l3 + 1

      if (l3 > 1) then
        rmsrat = rmsError / rmsErrorPrevious
        drate  = Max (0.2d0*drate, rmsrat)
      else
        rmsrat = 1.0d0
      end if

      dcon = rmsError * Min (conpst(nqq), conp15(nqq)*drate)


   ! routine end calculateRmsError
!     ---------------------------------------


   ! routine start checkAccumulatedError
!     ---------------------------------------


!     --------------------------------------------------------
!     If convergence occurs, go on to check accumulated error.
!     --------------------------------------------------------

      if (dcon > 1.0d0) then
      print*, "convergence"
!       -------------------------------------------------------------------
!       If nonconvergence after one step, re-evaluate first derivative with
!       new values of cnew.
!       -------------------------------------------------------------------

        if (l3 == 1) then
        print*, "re evaulate first deriviatives"

!         =========
          go to 300 !evalFirstDerivative
!         =========

!         ----------------------------------------------------------------
!         The corrector iteration failed to converge.
!
!         If the Jacobian matrix is more than one step old, update the
!         Jacobian and try convergence again.  If the Jacobian is current,
!         then reduce the time step, reset the accumulated derivatives to
!         their values before the failed step, and retry with the smaller
!         step.
!         ----------------------------------------------------------------

        else if (jeval == 0) then
         print*, "conversion failure"
          numFailOldJacobian = numFailOldJacobian + 1
          jeval = 1

!         =========
          go to  250 ! correctionLoop
!         =========

        end if

           ! routine stop checkAccumulatedError
!     ---------------------------------------

        numFailAfterPredict     = numFailAfterPredict + 1
        rdelmax   = 2.0d0
        jeval     = 1
        ifsuccess = 0
        xelaps    = told
        rdelt     = fracdec

        i1 = nqqisc + 1

        do jb = 1, nqq

          i1 = i1 - num1stOEqnsSolve

          do i = i1, nqqisc

            j = i + num1stOEqnsSolve

            do kloop = 1, ktloop
              conc(kloop,i) = conc(kloop,i) - conc(kloop,j)
            end do

          end do

        end do

!       =========
        go to 200 ! updateLimitTighten
!       =========

      end if


   ! routine end checkAccumulatedError
!     ---------------------------------------


          ! routine start correctorIterationConveraged
!     -------------------------------------------------------------------
!     The corrector iteration converged.
!
!     Set jeval = 0, so that it does not need to be called the next step.
!     If all else goes well.  Next, test the accumulated error from the
!     convergence process above.
!     -------------------------------------------------------------------
   print*, "corrector iteration converaged"
      jeval = 0

      if (l3 > 1) then

        do kloop = 1, ktloop
          dely(kloop) = 0.0d0
        end do

!c
        do kloop = 1, ktloop
          do jspc = 1, num1stOEqnsSolve
            errymax     = dtlos(kloop,jspc) * chold(kloop,jspc)
            dely(kloop) = dely(kloop) + errymax * errymax
          end do
        end do

        der2max = 0.0d0

        do kloop = 1, ktloop

          if (dely(kloop) > der2max) then
            der2max = dely(kloop)
          end if

        end do

      end if

          ! routine end correctorIterationConveraged
!     -------------------------------------------------------------------


          ! routine start accumulatedErrorTestFailed
!     -------------------------------------------------------------------

!     ----------------------------------------------------------------
!     The accumulated error test failed.
!
!     In all cases, reset the derivatives to their values before the
!     last time step.  Next:
!       (a) re-estimate a time step at the same or one lower order and
!           retry the step;
!       (b) if the first attempts fail, retry the step at fracdec the
!           the prior step;
!       (c) iF this fails, reset the order to 1 and go back to the
!           beginning, at order = 1, because errors of the wrong order
!           have accumulated.
!     ----------------------------------------------------------------

!     ==============================
      DER2MAXIF: if (der2max > enqq) then
!     ==============================
      print*, "der2max > enqq"
        xelaps = told
        numFailErrorTest  = numFailErrorTest + 1
        jfail  = jfail  + 1
        i1     = nqqisc + 1

        do jb = 1, nqq

          i1 = i1 - num1stOEqnsSolve

          do i = i1, nqqisc

            j = i + num1stOEqnsSolve

            do kloop = 1, ktloop
              conc(kloop,i) = conc(kloop,i) - conc(kloop,j)
            end do

          end do

        end do

        rdelmax = 2.0d0

        if (jfail <= 6) then

          ifsuccess = 0
          rdeltup   = 0.0d0
      print*, "jfail <= 6"
!         =========
          go to 400
!         =========

        else if (jfail <= 20) then
      print*, "jfail <= 20"
          ifsuccess = 0
          rdelt     = fracdec

!         =========
          go to 200
!         =========

        else

          delt    = delt * 0.1d0
          rdelt   = 1.0d0
          jfail   = 0
          jrestar = jrestar + 1
          idoub   = 5

          do jspc = 1, num1stOEqnsSolve
            do kloop = 1, ktloop
              cnew(kloop,jspc) = conc(kloop,jspc)
            end do
          end do

          if (pr_smv2) then
            Write (lunsmv,970) delt, xelaps
          end if

 970      format ('delt dec to ', e13.5, ' at time ', e13.5,  &
     &            ' because of excessive errors.')

          if (jrestar == 100) then

            if (pr_smv2) then
              Write (lunsmv,980)
            end if

 980        format ('Smvgear:  Stopping because of excessive errors.')

            call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)

          end if

!         =========
          go to 150
!         =========

        end if

                  ! routine end accumulatedErrorTestFailed
!     -------------------------------------------------------------------

                          ! routine start successfulStep
!     -------------------------------------------------------------------


!     ====
      else
!     ====

!       -------------------------------------------------------------
!       All successful steps come through here.
!
!       After a successful step, update the concentration and all
!       derivatives, reset told, set ifsuccess = 1, increment numSuccessTdt,
!       and reset jfail = 0.
!       -------------------------------------------------------------
   print*, "successful stepping"

        if (pr_qqjk .and. do_qqjk_inchem) then
          xtimestep = xelaps - told

!         =================
          call Do_Smv2_Diag  &
!         =================
     &      (jlooplo, ktloop, pr_nc_period, tdt, told, do_cell_chem,  &
     &       jreorder, inewold, denair, cnew, xtimestep, &
     &       yda, qqkda, qqjda, qkgmi, qjgmi, &
     &       ilong, ilat, ivert, itloop, &
     &       CTMi1, CTMi2, CTMju1, CTMj2, CTMk1, CTMk2, &
     &       num_qjo, num_qks, num_qjs, num_active)
        end if


        jfail     = 0
        ifsuccess = 1
        numSuccessTdt    = numSuccessTdt + 1
        told      = xelaps

        i1 = 1

        do j = 2, kstep

          i1 = i1 + num1stOEqnsSolve

          asnqqj = aset(nqq,j)

          do jspc = 1, num1stOEqnsSolve

            i = jspc + i1 - 1

            do kloop = 1, ktloop
              conc(kloop,i) =  &
     &          conc(kloop,i) + (asnqqj * dtlos(kloop,jspc))
            end do

          end do

        end do

                                  ! routine end successfulStep
!     -------------------------------------------------------------------


                          ! routine start updateChemistryMassBalance
!     -------------------------------------------------------------------
!       ------------------------------
!       Update chemistry mass balance.
!       ------------------------------
print*, "update chemistry mass balance"

        if (asn1 == 1.0d0) then

          do i = 1, num1stOEqnsSolve
            do kloop = 1, ktloop

              smvdm(kloop,i) =  &
     &          smvdm(kloop,i) + dtlos(kloop,i) + explic(kloop,i)

              conc(kloop,i) = conc(kloop,i) + dtlos(kloop,i)

            end do
          end do

        else

          do i = 1, num1stOEqnsSolve
            do kloop = 1, ktloop

              dtasn1         = asn1 * dtlos(kloop,i)
              smvdm(kloop,i) = smvdm(kloop,i) + dtasn1 + explic(kloop,i)
              conc (kloop,i) = conc (kloop,i) + dtasn1

            end do
          end do

        end if

! routine end updateChemistryMassBalance
!     -------------------------------------------------------------------


!       ---------------------------------------------------
!       Exit smvgear if a time interval has been completed.
!       ---------------------------------------------------

        timremain = chemTimeInterval - xelaps
   print*, "checking time interval"
        if (timremain <= 1.0d-06) return

!       -------------------------------------------------------------------
!       idoub counts the number of successful steps before re-testing the
!       step-size and order:
!         if idoub > 1, decrease idoub and go on to the next time step with
!                       the current step-size and order;
!         if idoub = 1, store the value of the error (dtlos) for the time
!                       step prediction, which will occur when idoub = 0,
!                       but go on to the next step with the current step
!                       size and order;
!         if idoub = 0, test the time step and order for a change.
!       -------------------------------------------------------------------
! routine start goToNextTimeStep
!     -------------------------------------------------------------------
!       ------------------------------
        if (idoub > 1) then

          idoub = idoub - 1

          if (idoub == 1) then

            do jspc = 1, num1stOEqnsSolve, 2

              jg1 = jspc + 1

              do kloop = 1, ktloop
                cest(kloop,jspc) = dtlos(kloop,jspc)
                cest(kloop,jg1)  = dtlos(kloop,jg1)
              end do

            end do

          end if

          rdelt = 1.0d0
         print*, "Going to 200"
!         =========
          go to 200
!         =========

        end if

        ! routine end goToNextTimeStep
!     -------------------------------------------------------------------
!       ------------------------------

!     ================
      end if DER2MAXIF
!     ================


!     ------------------------------------------------------------------
!     Test whether to change the step-size and order.
!
!     Determine the time step at (a) one order lower than, (b) the same
!     order as, and (c) one order higher than the current order.  In the
!     case of multiple grid-cells in a grid-block, find the minimum
!     step size among all the cells for each of the orders.  Then, in
!     all cases, choose the longest time step among the three steps
!     paired with orders, and choose the order allowing this longest
!     step.
!     ------------------------------------------------------------------

!     ---------------------------------------------------------------
!     Estimate the time step ratio (rdeltup) at one order higher than
!     the current order.  If nqq >= MAXORD, then we do not allow the
!     order to increase.
!     ---------------------------------------------------------------
        ! routine start changeStepSizeAndOrderTest
!     -------------------------------------------------------------------
   print*, "testing whether or not to change time order"
      if (nqq < MAXORD) then

        do kloop = 1, ktloop
          dely(kloop) = 0.0d0
        end do

        do jspc = 1, num1stOEqnsSolve
          do kloop = 1, ktloop
            errymax     = (dtlos(kloop,jspc) - cest(kloop,jspc)) *  &
     &                    chold(kloop,jspc)
            dely(kloop) = dely(kloop) + (errymax * errymax)
          end do
        end do

        der3max = 0.0d0

        do kloop = 1, ktloop

          if (dely(kloop) > der3max) then
            der3max = dely(kloop)
          end if

        end do

        rdeltup = 1.0d0 / ((conp3 * der3max**enqq3(nqq)) + 1.4d-6)

      else

        rdeltup = 0.0d0

      end if

        ! routine start changeStepSizeAndOrderTest
!     -------------------------------------------------------------------

!     ========
 400  continue
!     ========
   print*, "in 400"

!     ------------------------------------------------------------
!     Estimate the time step ratio (rdeltsm) at the current order.
!     der2max was calculated during the error tests earlier.
!     ------------------------------------------------------------

      rdeltsm = 1.0d0 / ((conp2 * der2max**enqq2(nqq)) + 1.2d-6)


!     ------------------------------------------------------------------
!     Estimate the time step ratio (rdeltdn) at one order lower than
!     the current order.  if nqq = 1, then we cannot test a lower order.
!     ------------------------------------------------------------------
        ! routine start estimateTimeStepRatioLowerOrder
!     -------------------------------------------------------------------

      if (nqq > 1) then

        do kloop = 1, ktloop
          dely(kloop) = 0.0d0
        end do

        kstepisc = (kstep - 1) * num1stOEqnsSolve

!c
        do kloop = 1, ktloop
          do jspc = 1, num1stOEqnsSolve

            i = jspc + kstepisc

            errymax     = conc(kloop,i) * chold(kloop,jspc)
            dely(kloop) = dely(kloop) + (errymax * errymax)

          end do

        end do

        der1max = 0.0d0

        do kloop = 1, ktloop
          if (dely(kloop) > der1max) then
            der1max = dely(kloop)
          end if
        end do

        rdeltdn = 1.0d0 / ((conp1 * der1max**enqq1(nqq)) + 1.3d-6)

      else

        rdeltdn = 0.0d0

      end if


!     -----------------------------------------------------------------
!     Find the largest of the predicted time step ratios of each order.
!     -----------------------------------------------------------------

      rdelt = Max (rdeltup, rdeltsm, rdeltdn)

      ! routine start estimateTimeStepRatioLowerOrder
!     -------------------------------------------------------------------



!     ---------------------------------------------------------------
!     If the last step was successful and rdelt is small, keep the
!     current step and order, and allow three successful steps before
!     re-checking the time step and order.
!     ---------------------------------------------------------------

      if ((rdelt < 1.1d0) .and. (ifsuccess == 1)) then

        idoub = 3

!       =========
        go to 200
!       =========

!       --------------------------------------------------------------
!       If the maximum time step ratio is that of one order lower than
!       the current order, decrease the order.  Do not minimize rdelt
!       to <= 1, when ifsuccess = 0 since this is less efficient.
!       --------------------------------------------------------------

      else if (rdelt == rdeltdn) then

        nqq = nqq - 1

!       ---------------------------------------------------------------
!       If the maximum time step ratio is that of one order higher than
!       the current order, increase the order and add a derivative term
!       for the higher order.
!       ---------------------------------------------------------------

      else if (rdelt == rdeltup) then

! routine start increareOrderAddDerTerm
!     -------------------------------------------------------------------

        real_kstep = kstep
        consmult   = aset(nqq,kstep) / real_kstep
        nqq        = kstep
        nqisc      = nqq * num1stOEqnsSolve

        do jspc = 1, num1stOEqnsSolve, 2

          jg1 = jspc + 1
          i1  = jspc + nqisc
          i2  = jg1  + nqisc

          do kloop = 1, ktloop
            conc(kloop,i1) = dtlos(kloop,jspc) * consmult
            conc(kloop,i2) = dtlos(kloop,jg1)  * consmult
          end do

        end do

! routine end increareOrderAddDerTerm
!     -------------------------------------------------------------------

      end if

!     ----------------------------------------------------------------
!     If the last two steps have failed, re-set idoub to the current
!     order + 1.  Do not minimize rdelt if jfail >= 2 since tests show
!     that this merely leads to additional computations.
!     ----------------------------------------------------------------

      idoub = nqq + 1

!     =========
      go to 200
!     =========

      return

      end subroutine Smvgear


      subroutine initializeGear ()

      end subroutine initializeGear


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

      implicit none

#     include "smv2chem_par.h"
#     include "smv2chem2.h"


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
!   cc2    : array of iarray units holding values of each matrix
!            position actually used; originally,
!            cc2 = P = I - delt * aset(nqq,1) * partial_derivatives;
!            however, cc2 is decomposed here
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!
!-----------------------------------------------------------------------------
! MRD: LU Decomp  - should go into sparseMatrix module
      subroutine Decomp  &
     &  (num1stOEqnsSolve, ktloop, ncsp, cc2, vdiag)

      implicit none

#     include "smv2chem_par.h"
#     include "smv2chem2.h"


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

        iar = jarrdiag(j,ncsp)

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


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Update
!
! DESCRIPTION
!   This routine updates photodissociation rates.
!
!   Photorates are included in first and partial derivative equations.
!
! ARGUMENTS
!   ktloop   : # of grid-cells in a grid-block
!   nallr    : # of active rxns
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!   ncsp     : ncs       => for daytime   gas chemistry
!              ncs + ICS => for nighttime gas chemistry
!   jphotrat : tbd
!   ptratk1  : tbd
!   rrate    : rate constants
!   trate    : rxn rate (moles l^-1-h2o s^-1 or # cm^-3 s^-1 (?))
!
!-----------------------------------------------------------------------------
! MRD: Update is probably setPhotolysisCoeffs and computeRateCoeffs
! MRD: Is going into the mechanism
! MRD: Solver will not call it

      subroutine Update  &
     &  (ktloop, nallr, ncs, ncsp, jphotrat, pratk1, rrate, trate)

      implicit none

#     include "smv2chem_par.h"
#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: ktloop
      integer, intent(in)  :: nallr
      integer, intent(in)  :: ncs
      integer, intent(in)  :: ncsp
      integer, intent(in)  :: jphotrat(ICS)
      real*8,  intent(in)  :: pratk1  (KBLOOP, IPHOT)

      real*8,  intent(inout) :: rrate(KBLOOP, NMTRATE)
      real*8,  intent(inout) :: trate(KBLOOP, NMTRATE*2)


!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: i, j
      integer :: kloop
      integer :: nh, nk, nkn


!     ----------------
!     Begin execution.
!     ----------------

!c    Write (6,*) 'Update called.'


!     -------------------------------
!     Load photolysis rate constants.
!     -------------------------------

      do j = 1, jphotrat(ncs)

        nkn = nknphotrt(j,ncs)

        do kloop = 1, ktloop
          rrate(kloop,nkn) = pratk1(kloop,j)
        end do

      end do


!     ------------------------------------------------------
!     Set rates where photoreaction has no active loss term.
!     ------------------------------------------------------

      do i = 1, nolosp(ncsp)

        nk  = nknlosp(i,ncs)
        nkn = newfold(nk,ncs)
        nh  = nkn + nallr

        do kloop = 1, ktloop
          trate(kloop,nkn) =  rrate(kloop,nkn)
          trate(kloop,nh)  = -trate(kloop,nkn)
        end do

      end do


      return

      end subroutine Update


