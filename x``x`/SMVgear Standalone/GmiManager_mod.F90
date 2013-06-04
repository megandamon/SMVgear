! MRD: Notes from Tom
! two levels of management
! top level : GMI specific reordering
! below: list of boxes - focused on managing gear updates its order and timestep


module GmiManager_mod

   use Smv2Chem2_mod

   implicit none
   private

#     include "smv2chem_par.h"

   public :: Manager_type
   public :: resetGear
   public :: startTimeInterval
   public :: calculateErrorTolerances
   public :: calcInitialTimeStepSize
   public :: updateCoefficients
   public :: calculateTimeStep
   public :: tightenErrorTolerance
   public :: calculateNewRmsError
   public :: testAccumulatedError
   public :: estimateTimeStepRatio
   public :: resetBeforeUpdate
   public :: calcNewAbsoluteErrorTolerance
   public :: scaleDerivatives
   public :: updateChold
   public :: predictConcAndDerivatives
   public :: resetCnewDerivatives
   public :: updateDerivatives
   public :: determineInitialAbTol
   public :: setInitialOrder
   public :: initCorrector
   public :: setConvergenceTerms
   public :: sumAccumulatedError
   public :: updateAfterNonConvTightenLimits
   public :: updateAfterAccumErrorTestFails
   public :: resetTermsBeforeStartingOver
   public :: updateAndResetAfterSucessfulStep
   public :: storeAccumErrorAndSetTimeStepRatio
   public :: estimateTimeStepRatioOneOrderHigher
   public :: increaseOrderAndAddDerivativeTerm

   public :: LIMIT_EXCESSIVE_FAILURES
   public :: REORDER_GRID_CELLS, SOLVE_CHEMISTRY
   public :: EVAL_PREDICTOR, DO_NOT_EVAL_PREDICTOR, PREDICTOR_JUST_CALLED
   public :: STEP_SUCCESS, STEP_FAILURE

   integer, parameter :: REORDER_GRID_CELLS = 1
   integer, parameter :: SOLVE_CHEMISTRY = 2
   integer, parameter :: STEP_SUCCESS = 1
   integer, parameter :: STEP_FAILURE = 0
   integer, parameter :: LIMIT_EXCESSIVE_FAILURES = 100

   ! MRD: do these belong here?
   integer, parameter :: EVAL_PREDICTOR = 1
   integer, parameter :: DO_NOT_EVAL_PREDICTOR = 0
   integer, parameter :: PREDICTOR_JUST_CALLED = -1



! MRD: add type bound procedures here
! can remove "_type"
! need a constructor
   type Manager_type

       ! private
       integer :: numErrTolDecreases
       integer :: numFailOldJacobian ! of times corrector failed to converge while the Jacobian was old
       integer :: numFailuresAfterVelocity ! of times correcter failed to converge after old "Pderiv" was called
       integer :: numFailErrorTest
       integer :: numFailAfterPredict
       integer :: numCallsPredict ! total # of times predictor is called
       integer :: numSuccessTdt ! numSuccessTdt    : total # of successful time steps taken
       integer :: numCallsVelocity ! total # of times velocity is called
       integer :: num1stOEqnsSolve !# of first-order eqns to solve, = # of spc = order of
               ! original matrix; num1stOEqnsSolve has a different value for day and
               ! night and for gas- and aqueous-phase chemistry;
               !  # spc with prod or loss terms in Smvgear (?)
       real*8  :: abtoler1
       real*8  :: order, order_inv
       real*8  :: chemTimeInterval !total chem time interval; same as chemintv (s)
       real*8  :: maxTimeStep ! max time step at a given time (s)
       real*8  :: failureFraction ! = 1 originially, but is decreased if excessive failures
               ! occur in order to reduce absolute error tolerance
       real*8  :: timeRemainingInChemInterval ! remaining time in an chem interval (s)
       real*8  :: iabove
       real*8  :: initialError, initialError_inv
       real*8  :: errmax_ncs_inv
       real*8  :: elapsedTimeInChemInterval ! elapsed time in chem interval (s)
       real*8  :: told !stores last value of elapsedTimeInChemInterval in case current step fails
       real*8  :: reltol1, reltol2, reltol3
       real*8  :: rmsError
       integer :: idoub ! records # of steps since the last change in step size or
               ! order; it must be at least oneHigherOrderIntegration = orderOfIntegrationMethod+1 before doubling is
               ! allowed
       integer :: nslp ! last time step # during which "Pderiv" was called
       integer :: numExcessiveFailures ! counts # of times Smvgear starts over at order 1 because of
               ! excessive failures

       integer :: oldOrderOfIntegrationMethod ! value of orderOfIntegrationMethod during last time step
       integer :: orderOfIntegrationMethod ! order of integration method; varies between 1 and MAXORD
       integer :: oneHigherOrderIntegration ! orderOfIntegrationMethod + 1
       real*8  :: hratio ! relative change in currentTimeStep*coeffsForIntegrationOrder(1) each change in step or order
               ! when Abs(hratio-1) > MAX_REL_CHANGE, reset jeval = 1 to call Pderiv
       real*8 :: integrationOrderCoeff1 ! value of coeffsForIntegrationOrder(orderOfIntegrationMethod,1)
       real*8 :: enqq ! pertst^2*order for current order
       real*8  :: conp1, conp2, conp3
       integer :: nqqisc ! orderOfIntegrationMethod * num1stOEqnsSolve
       real*8 :: timeStepRatio ! factor by which currentTimeStep is increased or decreased
       real*8 :: rdelmax ! max factor by which currentTimeStep can be increased in a single step;
               !                 as in Lsodes, set it to 1d4 initially to compensate for the
               !                 small initial currentTimeStep, but then set it to 10 after successful
               !                 steps and to 2 after unsuccessful steps
       real*8 :: der2max
       real*8  :: drate ! parameter which is used to determine whether convergence
               !                 has occurred
       real*8  :: dcon
       real*8  :: accumulatedError (KBLOOP, MXGSAER)! an array of length num1stOEqnsSolve, used for the accumulated corrections;
               ! on a successful return; accumulatedError(kloop,i) contains the estimated one step local error in cnew
      real*8  :: chold (KBLOOP, MXGSAER) ! 1 / (reltol * cnew + abtol); multiply chold by local errors in different error tests
      real*8  :: timeStepRatioLowerOrder ! time step ratio at one order lower  than current order
      real*8  :: timeStepRatioHigherOrder ! time step ratio at one order higher than current order
      integer :: ifsuccess ! identifies whether step is successful (=1) or not (=0)
      real*8 :: tolerance (KBLOOP)
      integer :: correctorIterations
      real*8  :: currentTimeStep
    end type Manager_type

contains



!-----------------------------------------------------------------------------
!
! ROUTINE
!   increaseOrderAndAddDerivativeTerm
! DESCRIPTION
! Increase the order and add a derivative term for the higher order.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
      subroutine increaseOrderAndAddDerivativeTerm (this, cnewDerivatives, &
                                    & ktloop)
         implicit none

         type(Manager_type) :: this
         real*8, intent(out) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer, intent(in) :: ktloop

         integer :: kloop, jspc
         integer :: jg1, i1, i2
         integer :: nqisc
         real*8  :: consmult
         real*8 :: oneHigherOrderIntegrationSave

         oneHigherOrderIntegrationSave = this%oneHigherOrderIntegration
         consmult   = coeffsForIntegrationOrder(this%orderOfIntegrationMethod,this%oneHigherOrderIntegration) / oneHigherOrderIntegrationSave
         this%orderOfIntegrationMethod = this%oneHigherOrderIntegration
         nqisc      = this%orderOfIntegrationMethod * this%num1stOEqnsSolve

         do jspc = 1, this%num1stOEqnsSolve, 2

            jg1 = jspc + 1
            i1  = jspc + nqisc
            i2  = jg1  + nqisc

            do kloop = 1, ktloop
              cnewDerivatives(kloop,i1) = this%accumulatedError(kloop,jspc) * consmult
              cnewDerivatives(kloop,i2) = this%accumulatedError(kloop,jg1)  * consmult
            end do

          end do

      end subroutine increaseOrderAndAddDerivativeTerm

!-----------------------------------------------------------------------------
!
! ROUTINE
!   estimateTimeStepRatioOneOrderHigher
! DESCRIPTION
! Estimate the time step ratio (timeStepRatioHigherOrder) at one order higher than
! the current order.  If orderOfIntegrationMethod >= MAXORD, then we do not allow the
! order to increase.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------

      subroutine estimateTimeStepRatioOneOrderHigher (this, maxOrder, ktloop, &
                  dely, cest)
         implicit none

         type(Manager_type) :: this
         integer, intent(in) :: maxOrder
         integer, intent(in) :: ktloop

         real*8, intent(inout) :: dely(KBLOOP)
         real*8, intent(in) :: cest  (KBLOOP, MXGSAER)

         integer :: kloop, jspc
         real*8 :: delyMax
         real*8 :: errymax

         if (this%orderOfIntegrationMethod < maxOrder) then

           do kloop = 1, ktloop
             dely(kloop) = 0.0d0
           end do

           do jspc = 1, this%num1stOEqnsSolve
             do kloop = 1, ktloop
               errymax     = (this%accumulatedError(kloop,jspc) - cest(kloop,jspc)) *  &
        &                    this%chold(kloop,jspc)
               dely(kloop) = dely(kloop) + (errymax * errymax)
             end do
           end do

           delyMax = 0.0d0
           do kloop = 1, ktloop
             if (dely(kloop) > delyMax) then
               delyMax = dely(kloop)
             end if
           end do

           this%timeStepRatioHigherOrder = 1.0d0 / &
            ((this%conp3 * delyMax**enqq3(this%orderOfIntegrationMethod)) + 1.4d-6)

         else
           this%timeStepRatioHigherOrder = 0.0d0
         end if

      end subroutine estimateTimeStepRatioOneOrderHigher

!-----------------------------------------------------------------------------
!
! ROUTINE
!   storeAccumErrorAndStepTimeStepRatio
! DESCRIPTION
! if idoub = 1, store the value of the error (accumulatedError) for the time
! step prediction, which will occur when idoub = 0,
! but go on to the next step with the current step size and order
! if idoub = 0, test the time step and order for a change.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------

      subroutine storeAccumErrorAndSetTimeStepRatio (this, accumErrorStorage, ktloop)
         implicit none

         type(Manager_type) :: this
         integer, intent(in) :: ktloop

         real*8, intent(out) :: accumErrorStorage(KBLOOP, MXGSAER)
         integer :: jg1, jspc, kloop

         this%idoub = this%idoub - 1
         if (this%idoub == 1) then
           do jspc = 1, this%num1stOEqnsSolve, 2
             jg1 = jspc + 1
             do kloop = 1, ktloop
               accumErrorStorage(kloop,jspc) = this%accumulatedError(kloop,jspc)
               accumErrorStorage(kloop,jg1)  = this%accumulatedError(kloop,jg1)
             end do
           end do
         end if

         this%timeStepRatio = 1.0d0

      end subroutine storeAccumErrorAndSetTimeStepRatio

!-----------------------------------------------------------------------------
!
! ROUTINE
!   resetTermsBeforeStartingOver
! DESCRIPTION
! After a successful step, update the concentration and all
! derivatives, reset told, set ifsuccess = 1, increment numSuccessTdt,
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
      subroutine updateAndResetAfterSucessfulStep (this, cnewDerivatives, ktloop)

         implicit none

         type(manager_type) :: this
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer, intent(in) :: ktloop

         this%numFailuresAfterVelocity     = 0
         this%ifsuccess = STEP_SUCCESS
         this%numSuccessTdt    = this%numSuccessTdt + 1
         this%told      = this%elapsedTimeInChemInterval

         call updateDerivatives(this, cnewDerivatives, ktloop)

      end subroutine updateAndResetAfterSucessfulStep

!-----------------------------------------------------------------------------
!
! ROUTINE
!   resetTermsBeforeStartingOver
! DESCRIPTION
! iF this fails, reset the order to 1 and go back to the
! beginning, at order = 1, because errors of the wrong order have accumulated.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
      subroutine resetTermsBeforeStartingOver (this, cnew, cnewDerivatives, &
                                          & ktloop, lunsmv, pr_smv2)

         implicit none
         type(manager_type) :: this
         real*8, intent(out) :: cnew(KBLOOP, MXGSAER)
         real*8, intent(in) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer, intent(in) :: ktloop
         integer, intent(in) :: lunsmv
         logical, intent(in) :: pr_smv2

         integer :: jspc, kloop

         this%currentTimeStep    = this%currentTimeStep * 0.1d0
         this%timeStepRatio   = 1.0d0
         this%numFailuresAfterVelocity   = 0
         this%numExcessiveFailures = this%numExcessiveFailures + 1
         this%idoub   = 5

         do jspc = 1, this%num1stOEqnsSolve
           do kloop = 1, ktloop
             cnew(kloop,jspc) = cnewDerivatives(kloop,jspc)
           end do
         end do

         if (pr_smv2) then
           Write (lunsmv,970) this%currentTimeStep, this%elapsedTimeInChemInterval
         end if

970      format ('currentTimeStep dec to ', e13.5, ' at time ', e13.5,  &
         ' because of excessive errors.')

      end subroutine resetTermsBeforeStartingOver

!-----------------------------------------------------------------------------
!
! ROUTINE
!   updateAfterNonConvTightenLimits
! DESCRIPTION
! if the Jacobian is current, then reduce the time step,
! reset the accumulated derivatives to their values before the failed step,
! and retry with the smaller step.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine updateAfterNonConvTightenLimits (this, maxFactorTimeStepIncrease, &
                              & elapsedTime, timeStepRatio)

      implicit none

      type (Manager_type) :: this
      real*8, intent(in) :: maxFactorTimeStepIncrease
      real*8, intent(in) :: elapsedTime
      real*8, intent(in) :: timeStepRatio

      this%numFailAfterPredict     = this%numFailAfterPredict + 1
      this%rdelmax   = maxFactorTimeStepIncrease
      this%ifsuccess = STEP_SUCCESS
      this%elapsedTimeInChemInterval    = elapsedTime
      this%timeStepRatio     = timeStepRatio

   end subroutine updateAfterNonConvTightenLimits

!-----------------------------------------------------------------------------
!
! ROUTINE
! The accumulated error test failed.
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine updateAfterAccumErrorTestFails (this)

      implicit none

      type (Manager_type) :: this
      this%elapsedTimeInChemInterval = this%told
      this%numFailErrorTest  = this%numFailErrorTest + 1
      this%numFailuresAfterVelocity  = this%numFailuresAfterVelocity  + 1
      this%rdelmax = 2.0d0

   end subroutine updateAfterAccumErrorTestFails

!-----------------------------------------------------------------------------
!
! ROUTINE
!   sumAccumulatedError
! DESCRIPTION
! Sum up the accumulated error, correct the concentration with the
! error, and begin to calculate the rmsnorm of the error relative
! to chold
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine sumAccumulatedError(this, cnew, cnewDerivatives, dely, gloss, &
                                 ktloop)
      implicit none

      type (Manager_type) :: this
      real*8, intent(out) :: cnew(KBLOOP, MXGSAER)
      real*8, intent(in) :: cnewDerivatives(KBLOOP, MXGSAER*7)
      real*8, intent(inout) :: dely(KBLOOP)
      real*8, intent(in) :: gloss(KBLOOP, MXGSAER)
      integer, intent(in) :: ktloop

      integer :: i, kloop
      real*8 :: errymax

      do kloop = 1, ktloop
        dely(kloop) = 0.0d0
      end do

      do i = 1, this%num1stOEqnsSolve
         do kloop = 1, ktloop
            this%accumulatedError(kloop,i) = this%accumulatedError(kloop,i) + gloss(kloop,i) !*
            cnew(kloop,i)  = cnewDerivatives(kloop,i)  + (this%integrationOrderCoeff1 * this%accumulatedError(kloop,i))
            errymax        = gloss(kloop,i) * this%chold(kloop,i) !*
            dely(kloop)    = dely(kloop)    + (errymax * errymax) !*
         end do
      end do

   end subroutine sumAccumulatedError


!-----------------------------------------------------------------------------
!
! ROUTINE
!   setConvergenceTerms
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine setConvergenceTerms(this, maxAllowableSteps)
      implicit none

      type (Manager_type) :: this
      integer, intent(in) :: maxAllowableSteps

       this%hratio = 1.0d0
       this%nslp   = this%numSuccessTdt + maxAllowableSteps
       this%drate  = 0.7d0

   end subroutine setConvergenceTerms

!-----------------------------------------------------------------------------
!
! ROUTINE
!   initCorrector
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine initCorrector(this, ktloop, concentrationsNew, cnewDerivatives)
      implicit none

      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8, intent(out) :: concentrationsNew(KBLOOP, MXGSAER)
      real*8, intent(in) :: cnewDerivatives(KBLOOP, MXGSAER*7)

      integer :: jspc, kloop

      this%correctorIterations = 0
      do jspc = 1, this%num1stOEqnsSolve
        do kloop = 1, ktloop
          concentrationsNew (kloop,jspc) = cnewDerivatives(kloop,jspc)
          this%accumulatedError(kloop,jspc) = 0.0d0
        end do
      end do
   end subroutine initCorrector

!-----------------------------------------------------------------------------
!
! ROUTINE
!   setInitialOrder
! DESCRIPTION
! Set initial order to 1
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine setInitialOrder(this, evaluatePredictor)
      implicit none

      type (Manager_type) :: this
      integer, intent(out) :: evaluatePredictor

      this%oldOrderOfIntegrationMethod = 0
      this%orderOfIntegrationMethod    = 1
      this%timeStepRatio  = 1.0d0

      evaluatePredictor  = EVAL_PREDICTOR
   end subroutine setInitialOrder


!-----------------------------------------------------------------------------
!
! ROUTINE
!   determineInitialAbTol
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
      subroutine determineInitialAbTol(this, concentrationsNew, concAboveAbtolCount, &
                                       & ireord, ktloop, ncs, yabst)
         implicit none

         type (Manager_type) :: this
         real*8, intent(in) :: concentrationsNew(KBLOOP, MXGSAER)
         integer, intent(inout) :: concAboveAbtolCount(KBLOOP, 5)
         integer, intent(in) :: ireord
         integer, intent(in) :: ktloop
         integer, intent(in) :: ncs
         real*8, intent(inout) :: yabst(KBLOOP)

         integer :: kloop

         if (ireord == SOLVE_CHEMISTRY) then
            call calcNewAbsoluteErrorTolerance (this, concentrationsNew, concAboveAbtolCount, ktloop, yabst, ncs)
            do kloop = 1, ktloop
               this%tolerance(kloop) = yabst(kloop) * this%reltol1
            end do
         else
            do kloop = 1, ktloop
               this%tolerance(kloop) = this%abtoler1
            end do
         end if

      end subroutine determineInitialAbTol

!-----------------------------------------------------------------------------
!
! ROUTINE
!   updateDerivatives
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------

      subroutine updateDerivatives(this, cnewDerivatives, ktloop)
         implicit none
         type (Manager_type) :: this
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer, intent(in) :: ktloop

         integer :: i, i1, j, jspc, kloop
         real*8  :: asnqqj

         i1 = 1
         do j = 2, this%oneHigherOrderIntegration
           i1 = i1 + this%num1stOEqnsSolve
           asnqqj = coeffsForIntegrationOrder(this%orderOfIntegrationMethod,j)
           do jspc = 1, this%num1stOEqnsSolve
             i = jspc + i1 - 1
             do kloop = 1, ktloop
               cnewDerivatives(kloop,i) =  cnewDerivatives(kloop,i) + (asnqqj * this%accumulatedError(kloop,jspc))
             end do
           end do
         end do
      end subroutine updateDerivatives

      subroutine resetCnewDerivatives(this, cnewDerivatives, ktloop)
         implicit none

         type (Manager_type) :: this
         real*8, intent(inout) :: cnewDerivatives(KBLOOP, MXGSAER*7)
         integer, intent(in) :: ktloop
         integer :: i,i1,j,jb,kloop

         i1 = this%nqqisc + 1
         j = 0
         do jb = 1, this%orderOfIntegrationMethod
            i1 = i1 - this%num1stOEqnsSolve
            do i = i1, this%nqqisc
               j = i + this%num1stOEqnsSolve
               do kloop = 1, ktloop
                  cnewDerivatives(kloop,i) = cnewDerivatives(kloop,i) - cnewDerivatives(kloop,j)
             end do
           end do
         end do

      end subroutine resetCnewDerivatives


!-----------------------------------------------------------------------------
!
! ROUTINE
!   predictConcAndDerivatives
! DESCRIPTION
!     Compute the predicted concentration and derivatives by multiplying
!     previous values by the pascal triangle matrix.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
      subroutine predictConcAndDerivatives(this, conc, explic, ktloop, prDiag)

         implicit none

         type (Manager_type) :: this
         real*8, intent(inout) :: conc(KBLOOP, MXGSAER*7)
         real*8, intent(out) :: explic(KBLOOP, MXGSAER)
         integer, intent(in) :: ktloop
         logical, intent(in) :: prDiag

         integer :: i,i1,j,jb,jspc,kloop

         if (prDiag) Write(*,*) "Computing predicted conc and derivatives using pascal triangle matrix"
         i1 = this%nqqisc + 1

         do jb = 1, this%orderOfIntegrationMethod - 1
           i1 = i1 - this%num1stOEqnsSolve
           do i = i1,  this%nqqisc
             j = i + this%num1stOEqnsSolve
             do kloop = 1, ktloop
               conc(kloop,i)  = conc(kloop,i) + conc(kloop,j)
             end do
           end do
         end do

         do jspc = 1,  this%num1stOEqnsSolve
           j = jspc + this%num1stOEqnsSolve
           do kloop = 1, ktloop
             conc  (kloop,jspc) = conc(kloop,jspc) + conc(kloop,j)
             explic(kloop,jspc) = conc(kloop,j)
           end do

         end do

         do i = this%num1stOEqnsSolve + 1, this%nqqisc
           j = i + this%num1stOEqnsSolve
           do kloop = 1, ktloop
             conc(kloop,i) = conc(kloop,i) + conc(kloop,j)
           end do
         end do

      end subroutine predictConcAndDerivatives


!-----------------------------------------------------------------------------
!
! ROUTINE
!   updateChold
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine updateChold (this, ktloop, cnew, yabst)
      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8,  intent(in) :: cnew  (KBLOOP, MXGSAER)
      real*8, intent(in)  :: yabst (KBLOOP)

      integer :: kloop, jspc

      do kloop = 1, ktloop
       do jspc = 1, this%num1stOEqnsSolve

         this%chold(kloop,jspc) =  &
  &        this%reltol3 /  &
  &        (Max (cnew(kloop,jspc), 0.0d0) +  &
  &         (yabst(kloop) * this%reltol2))

       end do
     end do

   end subroutine updateChold

!-----------------------------------------------------------------------------
!
! ROUTINE
!   scaleDerivatives
! DESCRIPTION
! Created by: Megan Rose Damon
! cnewDerivatives   : an array of length num1stOEqnsSolve*(MAXORD+1) that carries the
!   derivatives of cnew, scaled by currentTimeStep^j/factorial(j), where j is
!   the jth derivative; j varies from 1 to orderOfIntegrationMethod; e.g., conc(jspc,2)
!   stores currentTimeStep*y' (estimated)
!-----------------------------------------------------------------------------
   subroutine scaleDerivatives (this, ktloop, cnewDerivatives)
      type (Manager_type) :: this
      integer, intent(in)  :: ktloop
      real*8, intent(inout)  :: cnewDerivatives  (KBLOOP, MXGSAER*7)

      real*8  :: rdelta
      integer :: i1, j, i, kloop

      rdelta = 1.0d0
      i1     = 1

      do j = 2, this%oneHigherOrderIntegration
         rdelta = rdelta * this%timeStepRatio
         i1 = i1 + this%num1stOEqnsSolve
          do i = i1, i1 + (this%num1stOEqnsSolve-1)
            do kloop = 1, ktloop
              cnewDerivatives(kloop,i) = cnewDerivatives(kloop,i) * rdelta
            end do
          end do
      end do

   end subroutine scaleDerivatives

!-----------------------------------------------------------------------------
!
! ROUTINE
!   calcNewAbsoluteErrorTolerance
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine  calcNewAbsoluteErrorTolerance (this, cnew, concAboveAbtolCount, ktloop, yabst, ncs)
      type (Manager_type) :: this
      real*8,  intent(in) :: cnew  (KBLOOP, MXGSAER)
      integer, intent(inout) :: concAboveAbtolCount(KBLOOP, 5)
      integer, intent(in)  :: ktloop
      real*8, intent(out)  :: yabst (KBLOOP)
      integer, intent(in)  :: ncs

      integer :: jspc, kloop
      integer :: k1, k2, k3, k4, k5, k
      real*8  :: cnw

      do k = 1, 5
          do kloop = 1, ktloop
            concAboveAbtolCount(kloop,k) = 0
          end do
      end do

      do jspc = 1, this%num1stOEqnsSolve
          do kloop = 1, ktloop
            cnw = cnew(kloop,jspc)
            do k = 1, 5
               if (cnw > absoluteErrorTolerance(k,ncs)) then
                 concAboveAbtolCount(kloop,k) = concAboveAbtolCount(kloop,k) + 1
                 exit
               end if
            end do
          end do
        end do

        do kloop = 1, ktloop

          k1 = concAboveAbtolCount(kloop,1)
          k2 = concAboveAbtolCount(kloop,2) + k1
          k3 = concAboveAbtolCount(kloop,3) + k2
          k4 = concAboveAbtolCount(kloop,4) + k3
          k5 = concAboveAbtolCount(kloop,5) + k4

          if (k1 > this%iabove) then
            yabst(kloop) = absoluteErrorTolerance(1,ncs) ! MRD: these yabst should be passed in
          else if (k2 > this%iabove) then    ! does the driver pass them in?
            yabst(kloop) = absoluteErrorTolerance(2,ncs) ! or does the mechanism specify them
          else if (k3 > this%iabove) then    ! tabled for now
            yabst(kloop) = absoluteErrorTolerance(3,ncs)
          else if (k4 > this%iabove) then
            yabst(kloop) = absoluteErrorTolerance(4,ncs)
          else if (k5 > this%iabove) then
            yabst(kloop) = absoluteErrorTolerance(5,ncs)
          else
            yabst(kloop) = absoluteErrorTolerance(6,ncs)
          end if

        end do

        end subroutine calcNewAbsoluteErrorTolerance
!-----------------------------------------------------------------------------
!
! ROUTINE
!   resetBeforeUpdate
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine resetBeforeUpdate (this)
      type (Manager_type) :: this

      this%hratio    = 0.0d0
      this%integrationOrderCoeff1      = 1.0d0
      this%ifsuccess = 1
      this%rdelmax   = 1.0d4

   end subroutine resetBeforeUpdate

!-----------------------------------------------------------------------------
!
! ROUTINE
!   estimateTimeStepRatio
! DESCRIPTION
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine estimateTimeStepRatio (this, ktloop, dely, conc)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8, intent(inout)  :: dely  (KBLOOP)
      real*8, intent(in)  :: conc  (KBLOOP, MXGSAER*7)

      integer :: kloop, kstepisc, jspc, i
      real*8  :: errymax, der1max
      real*8  :: rdeltsm   ! time step ratio at current order

      !     Estimate the time step ratio (rdeltsm) at the current order.
      !     der2max was calculated during the error tests earlier.
      rdeltsm = 1.0d0 / ((this%conp2 * this%der2max**enqq2(this%orderOfIntegrationMethod)) + 1.2d-6)

      !     Estimate the time step ratio (timeStepRatioLowerOrder) at one order lower than
      !     the current order.  if orderOfIntegrationMethod = 1, then we cannot test a lower order
      if (this%orderOfIntegrationMethod > 1) then

         do kloop = 1, ktloop
            dely(kloop) = 0.0d0
         end do

         kstepisc = (this%oneHigherOrderIntegration - 1) * this%num1stOEqnsSolve

         do kloop = 1, ktloop
            do jspc = 1, this%num1stOEqnsSolve
               i = jspc + kstepisc
               errymax     = conc(kloop,i) * this%chold(kloop,jspc)
               dely(kloop) = dely(kloop) + (errymax * errymax)
            end do
         end do

        der1max = 0.0d0

        do kloop = 1, ktloop
          if (dely(kloop) > der1max) then
            der1max = dely(kloop)
          end if
        end do

        this%timeStepRatioLowerOrder = 1.0d0 / ((this%conp1 * der1max**enqq1(this%orderOfIntegrationMethod)) + 1.3d-6)

      else
        this%timeStepRatioLowerOrder = 0.0d0
      end if

      !     Find the largest of the predicted time step ratios of each order.
      this%timeStepRatio = Max (this%timeStepRatioHigherOrder, rdeltsm, this%timeStepRatioLowerOrder)


   end subroutine estimateTimeStepRatio

!-----------------------------------------------------------------------------
!
! ROUTINE
!   testAccumulatedError
! DESCRIPTION
!     Test the accumulated error from the convergence process.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine testAccumulatedError (this, ktloop, dely)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8, intent(inout)  :: dely  (KBLOOP)

      integer:: jspc, kloop
      real*8  :: errymax

      do kloop = 1, ktloop
         dely(kloop) = 0.0d0
      end do

      do kloop = 1, ktloop
         do jspc = 1, this%num1stOEqnsSolve
            errymax     = this%accumulatedError(kloop,jspc) * this%chold(kloop,jspc)
            dely(kloop) = dely(kloop) + errymax * errymax
          end do
      end do

      this%der2max = 0.0d0

      do kloop = 1, ktloop
         if (dely(kloop) > this%der2max) then
            this%der2max = dely(kloop)
         end if
      end do

   end subroutine testAccumulatedError

!-----------------------------------------------------------------------------
! ROUTINE
!   calculateNewRmsError
! DESCRIPTION
!     Set the previous rms error and calculate the new rms error.
!     If dcon < 1, then sufficient convergence has occurred.  Otherwise,
!     if the ratio of the current to previous rmsError is decreasing,
!     iterate more.  If it is not, then the convergence test failed.
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine calculateNewRmsError (this, ktloop, dely, l3)
      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in) :: ktloop
      real*8, intent(in)  :: dely  (KBLOOP)
      integer, intent(inout) :: l3

      real*8  :: rmsErrorPrevious, rmsrat
      integer :: kloop

      rmsErrorPrevious = this%rmsError
      this%der2max = 0.0d0

      ! make this a one line using maxval
      do kloop = 1, ktloop
        if (dely(kloop) > this%der2max) then
          this%der2max = dely(kloop)
        end if
      end do

      this%rmsError = Sqrt (this%der2max * this%order_inv)
      l3 = l3 + 1

      if (l3 > 1) then
        rmsrat = this%rmsError / rmsErrorPrevious
        this%drate  = Max (0.2d0*this%drate, rmsrat)
      else
        rmsrat = 1.0d0
      end if

      if (this%orderOfIntegrationMethod /= 0) then
         this%dcon = this%rmsError * Min (conpst(this%orderOfIntegrationMethod), conp15(this%orderOfIntegrationMethod)*this%drate)
      else
         this%dcon = 0
      endif

   end subroutine calculateNewRmsError

!-----------------------------------------------------------------------------
! ROUTINE
!   tightenErrorTolerance
! DESCRIPTION
!     tighten absoloute error tolerance and restart integration
!     at beginning of time interval
! Created by: Megan Rose Damon
!   pr_smv2  : should the SmvgearII     output file be written
!   lunsmv   : logical unit number to write to when pr_smv2 is true
!   ncs      : identifies gas chemistry type (1..NCSGAS)
! WARNING: this routine may not have been adequately tested
!-----------------------------------------------------------------------------
   subroutine tightenErrorTolerance (this, pr_smv2, lunsmv, ncs)

      use GmiPrintError_mod, only : GmiPrintError

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this

      logical, intent(in)  :: pr_smv2
      integer, intent(in)  :: lunsmv
      integer, intent(in)  :: ncs ! ncs is argument to Smvgear

      if (pr_smv2) then
         Write (lunsmv,950) this%currentTimeStep, this%timeRemainingInChemInterval, this%failureFraction, relativeErrorTolerance(ncs)
      end if

      950    format ('Smvgear:  currentTimeStep      = ', 1pe9.3, /,  '          timremain = ', 1pe9.3, /,  &
         &          '          failureFraction      = ', 1pe9.3, /,  '          errmax    = ', 1pe9.3)

      this%numErrTolDecreases = this%numErrTolDecreases + 1
      this%failureFraction     = this%failureFraction * 0.01d0

      ! handle magic number
      if (this%numErrTolDecreases == 10) then
         if (pr_smv2) then
            Write (lunsmv,960)
         end if

         960      format ('Smvgear:  too many decreases of failureFraction.')

         call GmiPrintError ('Problem in Smvgear', .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
      end if

   end subroutine tightenErrorTolerance

!-----------------------------------------------------------------------------
! ROUTINE
!   calculateTimeStep
! DESCRIPTION
!     Limit size of timeStepRatio, then recalculate new time step and update
!     hratio.  Use hratio to determine whether or not the predictor
!     should be updated. (edit by MRD on 2/27/2013)
! Created by: Megan Rose Damon
! Unit testing ideas: can't take a step bigger than what we have remaining
! Make this a method on a class?
! two routines: update time step and determine Jacobian (something like this)
!-----------------------------------------------------------------------------
   subroutine calculateTimeStep (this, jeval, maxRelChange)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(out) :: jeval
      real*8, intent(in) :: maxRelChange

      real*8  :: hmtim

      hmtim  = Min (this%maxTimeStep, this%timeRemainingInChemInterval)
      this%timeStepRatio  = Min (this%timeStepRatio, this%rdelmax, hmtim/this%currentTimeStep)
      this%currentTimeStep   = this%currentTimeStep   * this%timeStepRatio

      this%hratio = this%hratio * this%timeStepRatio
      this%elapsedTimeInChemInterval = this%elapsedTimeInChemInterval + this%currentTimeStep
      ! rename nslp
      if ((Abs (this%hratio-1.0d0) > maxRelChange) .or. (this%numSuccessTdt >= this%nslp)) then
        jeval = 1 ! MRD: could be a boolean; this is signifying to whether or not to update Jacobian
      end if

   end subroutine calculateTimeStep


!-----------------------------------------------------------------------------
!
! ROUTINE
!   updateCoefficients
!
! DESCRIPTION
!     Update coefficients of (for?) the order; note that pertst2 is the original
!     pertst^2.
! MRD: Gear can be 1st, 2nd, etc. order
! MRD: this could be where it changes it order
! MRD: track these down, figure out what they are
! Created by: Megan Rose Damon
!-----------------------------------------------------------------------------
   subroutine updateCoefficients (this)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this

      real*8 :: eup ! pertst^2*order for one order higher than current order
      real*8  :: edwn ! pertst^2*order for one order lower  than current order

      this%oldOrderOfIntegrationMethod = this%orderOfIntegrationMethod
      this%oneHigherOrderIntegration  = this%orderOfIntegrationMethod + 1
      this%hratio = this%hratio * coeffsForIntegrationOrder(this%orderOfIntegrationMethod,1) / this%integrationOrderCoeff1
      this%integrationOrderCoeff1   = coeffsForIntegrationOrder(this%orderOfIntegrationMethod,1)
      this%enqq   = coeffsForSelectingStepAndOrder(this%orderOfIntegrationMethod,1) * this%order
      eup    = coeffsForSelectingStepAndOrder(this%orderOfIntegrationMethod,2) * this%order
      edwn   = coeffsForSelectingStepAndOrder(this%orderOfIntegrationMethod,3) * this%order
      this%conp3  = 1.4d0 /  (eup**enqq3(this%orderOfIntegrationMethod)) !eup is zero
      this%conp2  = 1.2d0 / (this%enqq**enqq2(this%orderOfIntegrationMethod)) !enqq is zero
      this%conp1  = 1.3d0 / (edwn**enqq1(this%orderOfIntegrationMethod)) !edwn is zero
      this%nqqisc = this%orderOfIntegrationMethod * this%num1stOEqnsSolve

   end subroutine updateCoefficients

!-----------------------------------------------------------------------------
!
! ROUTINE
!   calcInitialTimeStepSize
!
! DESCRIPTION
!     Calculate initial time step size (s).
!     Sqrt (dely / [initialError * order]) =
!       rmsnorm of error scaled to initialError *
!       cnew + abtol / reltol
! This is a guess, later it will adapt
! Created by: Megan Rose Damon
!   ktloop   : # of grid-cells in a grid-block
!   dely     : TBD
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!-----------------------------------------------------------------------------
   subroutine calcInitialTimeStepSize (this, ktloop, dely, ncs)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in)  :: ktloop
      real*8, intent(in)  :: dely  (KBLOOP)
      integer, intent(in)  :: ncs ! ncs is argument to Smvgear

      integer :: kloop
      real*8  :: rmstop
      real*8  :: delt1

      rmstop = 0.0d0
      do kloop = 1, ktloop
         if (dely(kloop) > rmstop) rmstop = dely(kloop)
      end do

      delt1 = Sqrt (this%initialError / (abst2(ncs) + (rmstop * this%order_inv)))
      this%currentTimeStep  = Max  (Min (delt1, this%timeRemainingInChemInterval, this%maxTimeStep), HMIN)

   end subroutine calcInitialTimeStepSize

!-----------------------------------------------------------------------------
!
! ROUTINE
!   calculateErrorTolerances
!
! DESCRIPTION
! Use lowest absolute error tolerance when reordering.
! It is conceviable that there could be different norms or error criteria
! Created by: Megan Rose Damon
!   ktloop   : # of grid-cells in a grid-block
!   cnew     : stores conc (y (estimated)) (molec/cm^3)
!   gloss    : value of first derivatives on output from velocity; right-side
!              of eqn on input to Backsub; error term (solution from Backsub)
!              on output from Backsub
!   dely     : TBD
!-----------------------------------------------------------------------------
   subroutine calculateErrorTolerances (this, ktloop, cnew, gloss, dely)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in)  :: ktloop
      real*8, intent(in) :: cnew  (KBLOOP, MXGSAER)
      real*8, intent(in) :: gloss (KBLOOP, MXGSAER)
      real*8, intent(inout)  :: dely  (KBLOOP)

      integer :: kloop
      integer :: jspc
      real*8  :: errymax

      do kloop = 1, ktloop
         do jspc = 1, this%num1stOEqnsSolve
            errymax     = gloss(kloop,jspc) / (cnew(kloop,jspc) + this%tolerance(kloop))
            dely(kloop) = dely(kloop) + (errymax * errymax)
         end do
       end do

   end subroutine calculateErrorTolerances

!-----------------------------------------------------------------------------
!
! ROUTINE
!   startTimeInterval
!
! DESCRIPTION
!   This routine starts time interval or is called after total failure.
! Created by: Megan Rose Damon
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!-----------------------------------------------------------------------------
   subroutine startTimeInterval (this, ncs)

      ! ----------------------
      ! Argument declarations.
      ! ----------------------
      type (Manager_type) :: this
      integer, intent(in)  :: ncs ! ncs is argument to Smvgear

      this%idoub     = 2
      this%nslp      = MBETWEEN
      this%numExcessiveFailures   = 0
      this%elapsedTimeInChemInterval    = 0.0d0
      this%told      = 0.0d0
      this%timeRemainingInChemInterval = this%chemTimeInterval

      this%reltol1   = this%failureFraction * this%initialError_inv

      this%reltol2   = this%failureFraction * this%errmax_ncs_inv

      this%reltol3   = this%errmax_ncs_inv

      ! MRD: abtol in common block?
      ! put abtol in another structure
      this%abtoler1  = absoluteErrorTolerance(6,ncs) * this%reltol1


   end subroutine startTimeInterval
!-----------------------------------------------------------------------------
!
! ROUTINE
!   resetGear
!
! DESCRIPTION
!   This routine is called once at the beginning of Smvgear
!
! Created by: Megan Rose Damon
!   ncs      : identifies gas chemistry type (1..NCSGAS)
!-----------------------------------------------------------------------------
      subroutine resetGear (this, ncsp, ncs, ifsun, hmaxnit)

         use Smv2Chem2_mod

         ! ----------------------
         ! Argument declarations.
         ! ----------------------
         type (Manager_type) :: this
         integer, intent(out) :: ncsp
         integer, intent(in)  :: ncs ! ncs is argument to Smvgear
         integer, intent(in)  :: ifsun ! ifsun is an argument to Smvgear
         real*8,  intent(in)  :: hmaxnit

         this%numFailOldJacobian     = 0
         this%numFailuresAfterVelocity     = 0
         this%numFailErrorTest     = 0
         this%numFailAfterPredict     = 0
         this%numCallsPredict   = 0
         this%numSuccessTdt    = 0
         this%numCallsVelocity   = 0
         this%numErrTolDecreases  = 0
         this%dcon = 0.0d0

         ! MAX_REL_CHANGE moved to parameter
         this%rmsError    = 1.0d0

         ! MRD: derived from common block
         this%num1stOEqnsSolve    = numOrigSpcGtrOrEql1PdTerm(ncs)

         this%order     = this%num1stOEqnsSolve
         this%order_inv = 1.0d0 / this%num1stOEqnsSolve

         ! MRD: derived from common bloc
         this%chemTimeInterval = timeintv(ncs)

         ! MRD: ICS is from common block
         ncsp      = (ifsun - 1) * ICS + ncs

         this%maxTimeStep = hmaxnit
         if (ifsun == 1) this%maxTimeStep = hmaxday(ncs)

         this%failureFraction   = 1.0d0
         this%iabove = this%order * 0.4d0
         this%initialError     = Min (relativeErrorTolerance(ncs), 1.0d-03)
         this%initialError_inv = 1.0d0 / this%initialError
         this%errmax_ncs_inv = 1.0d0 / relativeErrorTolerance(ncs)

      end subroutine resetGear


end module GmiManager_mod
