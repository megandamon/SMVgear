
!=============================================================================
!
! $Id: smv2chem_diag.F90,v 1.1.1.1 2008-02-12 16:06:36 trayanov Exp $
!
! CODE DEVELOPER
!   Original code from Peter Connell, LLNL
!   Gmimod modifications:  John Tannahill
!                          jrt@llnl.gov
!
! FILE
!   smv2chem_diag.F
!
! ROUTINES
!   Do_Smv2_Diag
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Smv2_Diag
!
! DESCRIPTION
!   This routine collects the Smvgear II chemical diagnostics.  It accumulates
!   averages for species and reaction rates.
!
! ARGUMENTS
!   jlooplo      : low ntloop grid-cell - 1 in a grid-block
!   ktloop       : number of grid-cells     in a grid-block
!   pr_nc_period : NetCDF output period
!   tdt          : model time step (s)
!   told         : stores last value of xelaps in case current step fails
!   do_cell_chem : do chemistry for a particular cell?
!   jreorder     : gives original grid-cell from re-ordered grid-cell
!   inewold      : original spc # of each new jnew spc
!   denair       : density of air (molec/cm^3)
!   cnew         : init (and final) spc conc
!                  (# cm^-3-air or moles l^-1-h2o (?))
!   xtimestep    : xelaps - told
!
!-----------------------------------------------------------------------------

      subroutine Do_Smv2_Diag  &
     &  (jlooplo, ktloop, pr_nc_period, tdt, told, do_cell_chem,  &
     &   jreorder, inewold, denair, cnew, xtimestep, &
     &   yda, qqkda, qqjda, qkgmi, qjgmi, &
     &   ilong, ilat, ivert, itloop, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   num_qjo, num_qks, num_qjs, num_active)

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilong, ilat, ivert, itloop
      integer, intent(in) :: num_qjo, num_qks, num_qjs, num_active
      real*8 , intent(in   ) :: qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)
      real*8 , intent(in   ) :: qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)
      real*8 , intent(inout) :: qqkda(i1:i2, ju1:j2, k1:k2, num_qks)
      real*8 , intent(inout) :: yda  (i1:i2, ju1:j2, k1:k2, num_active)

      integer, intent(in)  :: jlooplo
      integer, intent(in)  :: ktloop
      real*8,  intent(in)  :: pr_nc_period
      real*8,  intent(in)  :: tdt
      real*8,  intent(in)  :: told
      logical, intent(in)  :: do_cell_chem(ilong, ilat, ivert)
      integer, intent(in)  :: jreorder(itloop)
      integer, intent(in)  :: inewold (MXGSAER, ICS)
      real*8,  intent(in)  :: denair  (ktloop)
      real*8,  intent(in)  :: cnew    (KBLOOP, MXGSAER)

      real*8,  intent(inout) :: xtimestep


!     ----------------------
!     Variable declarations.
!     ----------------------

      logical :: end_chem_tstep

      logical, save :: end_period = .true.
      logical, save :: first      = .true.

      integer :: ic, ix
      integer :: il, ij, ik
      integer :: ind
      integer :: jloop, kloop
      integer :: smv2di

      integer, save :: nsteps            = 0
      integer, save :: nsteps_per_period = 0

      integer, allocatable, save :: lonloop(:)
      integer, allocatable, save :: latloop(:)
      integer, allocatable, save :: altloop(:)

      integer, allocatable, save :: isteps (:,:,:)

      real*8  :: rsteps_per_period
      real*8  :: tacc_tdt
      real*8  :: tfrac

      real*8  :: qjblock (ktloop, num_qjs)
      real*8  :: qqjblock(ktloop, num_qjs)

      real*8  :: qkblock (ktloop, num_qks)
      real*8  :: qqkblock(ktloop, num_qks)

      real*8  :: yblock  (ktloop, num_active)

      real*8, allocatable, save :: taccum(:,:,:)

      real*8, allocatable, save :: qqjts (:,:,:,:)
      real*8, allocatable, save :: qqkts (:,:,:,:)
      real*8, allocatable, save :: yts   (:,:,:,:)


!     ----------------
!     Begin execution.
!     ----------------

!     ==========
      if (first) then
!     ==========

        first = .false.

        Allocate (lonloop(itloop))
        Allocate (latloop(itloop))
        Allocate (altloop(itloop))
        lonloop = 0; latloop = 0; altloop = 0

        Allocate (isteps(ilong, ilat, ivert))
        isteps = 0

        Allocate (taccum(ilong, ilat, ivert))
        taccum = 0.0d0

        Allocate (qqjts(ilong, ilat, ivert, num_qjs))
        Allocate (qqkts(ilong, ilat, ivert, num_qks))
        Allocate (yts  (ilong, ilat, ivert, num_active))
        qqjts = 0.0d0; qqkts = 0.0d0; yts = 0.0d0


        ind = 0

        do ik = k1, k2
          do ij = ju1, j2
            do il = i1, i2

              ind = ind + 1

              lonloop(ind) = il
              latloop(ind) = ij
              altloop(ind) = ik

            end do
          end do
        end do

        nsteps_per_period = Nint (pr_nc_period / tdt)

      end if


!     ===============
      if (end_period) then
!     ===============

        end_period = .false.

!       -------------------------------------
!       Initialize arrays at start of period.
!       -------------------------------------

        isteps (:,:,:) = 0

        taccum (:,:,:) = 0.0d0

        qqjda(:,:,:,:) = 0.0d0
        qqkda(:,:,:,:) = 0.0d0
        yda  (:,:,:,:) = 0.0d0

        qqjts(:,:,:,:) = 0.0d0
        qqkts(:,:,:,:) = 0.0d0
        yts  (:,:,:,:) = 0.0d0

      end if


      do kloop = 1, ktloop

        jloop  = jlooplo + kloop
        smv2di = jreorder(jloop)

        il = lonloop(smv2di)
        ij = latloop(smv2di)
        ik = altloop(smv2di)

        qjblock(kloop,:) = qjgmi(il,ij,ik,:)
        qkblock(kloop,:) = qkgmi(il,ij,ik,:)

        do ic = 1, num_active
          ix = inewold(ic,1)
          yblock(kloop,ix) = cnew(kloop,ic)
        end do

      end do


!     =====================
      call Calc_Rate_Setkin  &
!     =====================
     &  (ktloop, num_qjs, num_qks, num_active, qkblock, qjblock,  &
     &   yblock, qqkblock, qqjblock)


!     -------------------------------------------------
!     Test for complete restart or time step overshoot.
!     -------------------------------------------------

      jloop  = jlooplo + 1
      smv2di = jreorder(jloop)

      il = lonloop(smv2di) - i1  + 1
      ij = latloop(smv2di) - ju1 + 1
      ik = altloop(smv2di)

      tacc_tdt = taccum(il,ij,ik)


      if ((told == 0.0d0) .and. (tacc_tdt /= 0.0d0)) then

!       ------------------------------------------------
!       SmvgearII has restarted timestep for this block.
!       ------------------------------------------------

        tacc_tdt = 0.0d0

        do kloop = 1, ktloop

          jloop  = jlooplo + kloop
          smv2di = jreorder(jloop)

          il = lonloop(smv2di) - i1  + 1
          ij = latloop(smv2di) - ju1 + 1
          ik = altloop(smv2di)

          taccum(il,ij,ik)   = 0.0d0

          qqjts (il,ij,ik,:) = 0.0d0
          qqkts (il,ij,ik,:) = 0.0d0
          yts   (il,ij,ik,:) = 0.0d0

        end do

      end if


      if ((Abs ((tacc_tdt + xtimestep) - tdt) <= 1.0d-06) .or.  &
     &    ((tacc_tdt + xtimestep) > tdt )) then

!       -----------------------------------------------
!       SmvgearII has overshot timestep for this block.
!       -----------------------------------------------

        xtimestep = tdt - tacc_tdt

      end if

      tfrac = xtimestep / pr_nc_period


!     ------------------------------------
!     Accumulate rates and concentrations.
!     ------------------------------------

      do kloop = 1, ktloop

        jloop  = jlooplo + kloop
        smv2di = jreorder(jloop)

        il = lonloop(smv2di) - i1  + 1
        ij = latloop(smv2di) - ju1 + 1
        ik = altloop(smv2di)

        where (qqjblock(kloop,:) > 0.0d0) qqjts(il,ij,ik,:)  =  &
     &    qqjts(il,ij,ik,:) +  &
     &    ((qqjblock(kloop,:)  / denair(kloop)) * tfrac)

        where (qqkblock(kloop,:) > 0.0d0) qqkts(il,ij,ik,:)  =  &
     &    qqkts(il,ij,ik,:) +  &
     &    ((qqkblock(kloop,:)  / denair(kloop)) * tfrac)

        yts(il,ij,ik,:) =  &
     &    yts(il,ij,ik,:) +  &
     &    (yblock(kloop,:) * tfrac)

        taccum(il,ij,ik) = taccum(il,ij,ik) + xtimestep

      end do


      if (Any (do_cell_chem(:,:,:) .and.  &
     &         (tdt - taccum(:,:,:) > 1.0d-06))) then
        end_chem_tstep = .false.
      else
        end_chem_tstep = .true.
      end if


!     ===================
      if (end_chem_tstep) then
!     ===================

        nsteps = nsteps + 1

        taccum(:,:,:) = 0.0d0

!       ----------------------------------------------------------------
!       Keep count of number of time steps cells have calculated values.
!       ----------------------------------------------------------------

        where (do_cell_chem(:,:,:)) isteps(:,:,:) = isteps(:,:,:) + 1


        qqjda(i1:i2,ju1:j2,k1:k2,:) =  &
     &    qqjda(i1:i2,ju1:j2,k1:k2,:) + qqjts(1:ilong,1:ilat,1:ivert,:)

        qqkda(i1:i2,ju1:j2,k1:k2,:) =  &
     &    qqkda(i1:i2,ju1:j2,k1:k2,:) + qqkts(1:ilong,1:ilat,1:ivert,:)

        yda  (i1:i2,ju1:j2,k1:k2,:) =  &
     &    yda  (i1:i2,ju1:j2,k1:k2,:) + yts  (1:ilong,1:ilat,1:ivert,:)


        qqjts(:,:,:,:) = 0.0d0
        qqkts(:,:,:,:) = 0.0d0
        yts  (:,:,:,:) = 0.0d0

      end if


!     ================================
      if (nsteps == nsteps_per_period) then
!     ================================

        end_period = .true.

        nsteps = 0

        rsteps_per_period = nsteps_per_period


        do ic = 1, num_qjs
          where (do_cell_chem(1:ilong,1:ilat,1:ivert))

            qqjda(i1:i2,ju1:j2,k1:k2,ic) =  &
     &        qqjda(i1:i2,ju1:j2,k1:k2,ic) *  &
     &        (rsteps_per_period / isteps(1:ilong,1:ilat,1:ivert))

          end where
        end do

        do ic = 1, num_qks
          where (do_cell_chem(1:ilong,1:ilat,1:ivert))

            qqkda(i1:i2,ju1:j2,k1:k2,ic) =  &
     &        qqkda(i1:i2,ju1:j2,k1:k2,ic) *  &
     &        (rsteps_per_period / isteps(1:ilong,1:ilat,1:ivert))

          end where
        end do

        do ic = 1, num_active
          where (do_cell_chem(1:ilong,1:ilat,1:ivert))

            yda(i1:i2,ju1:j2,k1:k2,ic) =  &
     &        yda(i1:i2,ju1:j2,k1:k2,ic) *  &
     &        (rsteps_per_period / isteps(1:ilong,1:ilat,1:ivert))

!c?       elsewhere

!c          yda(:,:,:,ic) = yinit(:,:,:,ic)

          end where
        end do

      end if


      return

      end subroutine Do_Smv2_Diag

