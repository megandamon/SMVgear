!-----------------------------------------------------------------------------
!
! ROUTINE
!   doSmv2Solver
!
! DESCRIPTION
!   This is the main control routine for the ordinary differential equation
!   solver, "Smvgear II" (Sparse Matrix Vectorized Gear-type code).
!
! ARGUMENTS
!   do_qqjk_inchem   : if pr_qqjk is on, should qqj's & qqk's be determined
!                      inside the chemistry solver, or outside?
!   do_semiss_inchem : do surface emissions inside the chemistry solver, or
!                      outside?
!   pr_diag          : print some diagnostic output to screen?
!   pr_qqjk          : should the periodic qqjk output file be written?
!   pr_smv2          : should the SmvgearII     output file be written
!                      (non-parallel mode only)?
!   loc_proc         : local processor #
!   ilat             : # of latitudes
!   ilong            : # of longitudes
!   ivert            : # of vertical layers
!   itloop           : # of zones (ilong * ilat * ivert)
!   pr_nc_period     : NetCDF output period
!   tdt              : model time step (s)
!   do_cell_chem     : do chemistry for a particular cell?
!   arate            : thermal    rate constants (units vary)
!   prate            : photolysis rate constants (s^-1)
!   yemis            : surface emissions (molec/cm^3/s)
!   cx               : spc conc (molec/cm^3)
!
!-----------------------------------------------------------------------------

      program doSmv2Solver  !&
!     &  (do_qqjk_inchem, do_semiss_inchem, pr_diag, pr_qqjk, pr_smv2,  &
!     &   loc_proc, ilat, ilong, ivert, itloop, pr_nc_period, tdt,  &
!     &   do_cell_chem, arate, prate, yemis, cx, &
!     &   yda, qqkda, qqjda, qkgmi, qjgmi, &
!     &   i1, i2, ju1, j2, k1, k2, &
!     &   num_qjo, num_qks, num_qjs, num_active)

      use timing_mod
      implicit none

#     include "smv2chem_par.h"
#     include "smv2chem1.h"
#     include "smv2chem2.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: i1, i2, ju1, j2, k1, k2
      integer :: numQjo, numQks, numQjs, numActive
      real*8, allocatable :: qjGmi(:, :, :, :)
      real*8, allocatable :: qkGmi(:, :, :, :)
      real*8, allocatable :: qqjda(:, :, :, :)
      real*8, allocatable :: qqkda(:, :, :, :)
      real*8, allocatable :: yda  (:, :, :, :)
      logical  :: doQqjkInchem
      logical  :: doSurfEmissInChem
      logical  :: prDiag
      logical  :: prQqjk
      logical  :: prSmv2
      integer  :: localProc
      integer  :: numLat, numLong, numVert
      integer  :: numZones
      real*8  :: prNcPeriod
      real*8  :: timeStep
      logical, allocatable  :: doCellChem(:)
      real*8, allocatable :: thermalRateConstants(:, :)
      real*8, allocatable :: photolysisRateConstants(:, :)
      real*8, allocatable :: surfaceEmissions(:, :)

      real*8, allocatable :: speciesConst(:, :)


!     ----------------------
!     Variable declarations.
!     ----------------------

      logical, save :: first = .true.

      integer, allocatable :: jReOrder(:)
      integer, allocatable :: lReOrder(:)

      real*8, allocatable  :: errorMx2  (:)

      real*8, allocatable, save :: cSumA(:)
      real*8, allocatable, save :: cSumB(:)
 
      integer :: rank,errorInt
      character(len=100) :: smv2Chem1Entry
      character(len=100) :: smv2Chem1Exit
      character(len=100) :: smv2Chem2Entry
      character(len=100) :: smv2Chem2Exit
      character(len=100) :: physProcEntry
      character(len=100) :: physProcExit

      character(len=128) :: tempText

!     ----------------
!     Begin execution.
!     ----------------

!      call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
      call timingInit

      rank = 17
      if (prDiag) then
        Write (6,*) 'doSmv2Solver called by ', localProc
      end if


      write(smv2Chem1Entry,1001) rank
 1001 format('smv2chem1_entry.proc',i4.4)
      write(smv2Chem1Exit,1002) rank
 1002 format('smv2chem1_exit.proc',i4.4)
      write(smv2Chem2Entry,1003) rank
 1003 format('smv2chem2_entry.proc',i4.4)
      write(smv2Chem2Exit, 1004) rank
 1004 format('smv2chem2_exit.proc',i4.4)
      write(physProcEntry, 1005) rank
 1005 format('physproc_entry.proc',i4.4)
      write(physProcExit, 1006) rank
 1006 format('physproc_exit.proc',i4.4)

!read smv2chem1 on entry to physproc to be used for standAlone code
      open(file=trim(smv2Chem1Entry),unit=23,form="formatted")
      read(23,*) 
      read(23,*) ifreord
      read(23,*) 
      read(23,*) ih2o
      read(23,*) 
      read(23,*) imgas
      read(23,*) 
      read(23,*) initrogen
      read(23,*) 
      read(23,*) ioxygen
      read(23,*) 
      read(23,*) kuloop
      read(23,*) 
      read(23,*) lunsmv
      read(23,*) 
      read(23,*) ncs
      read(23,*) 
      read(23,*) jphotrat
      read(23,*) 
      read(23,*) nrates
      read(23,*) 
      read(23,*) ntloopncs
      read(23,*) 
      read(23,*) ntspec
      read(23,*) 
      read(23,*) inewold
      read(23,*) 
      read(23,*) npphotrat
      read(23,*) 
      read(23,*) fracdec
      read(23,*) 
      read(23,*) hmaxnit

!read smv2chem2 on entry to physproc to be used for standAlone code
      open(file=trim(smv2Chem2Entry),unit=25,form="formatted")
      read(25,*) 
      read(25,*) ioner
      read(25,*) 
      read(25,*) nallrat
      read(25,*) 
      read(25,*) inorep
      read(25,*) 
      read(25,*) ithrr
      read(25,*) 
      read(25,*) itwor
      read(25,*) 
      read(25,*) nm3bod
      read(25,*) 
      read(25,*) nmair
      read(25,*) 
      read(25,*) nmn2
      read(25,*) 
      read(25,*) nmo2
      read(25,*) 
      read(25,*) nmoth
      read(25,*) 
      read(25,*) ntrates
      read(25,*) 
      read(25,*) mappl
      read(25,*) 
      read(25,*) lgasbino
      read(25,*) 
      read(25,*) nreacoth
      read(25,*) 
      read(25,*) lgas3bod
      read(25,*) 
      read(25,*) losinacp
      read(25,*) 
      read(25,*) nreac3b
      read(25,*) 
      read(25,*) nreacair
      read(25,*) 
      read(25,*) nreacn2
      read(25,*) 
      read(25,*) nreaco2
      read(25,*) 
      read(25,*) jphotnk
      read(25,*) 
      read(25,*) noldfnew
      read(25,*) 
      read(25,*) irm2
      read(25,*) 
      read(25,*) ischang
      read(25,*) 
      read(25,*) kzthi
      read(25,*) 
      read(25,*) kztlo
      read(25,*) 
      read(25,*) ikztot
      read(25,*) 
      read(25,*) kbh1
      read(25,*) 
      read(25,*) kbh2
      read(25,*) 
      read(25,*) kbh3
      read(25,*) 
      read(25,*) kbh4
      read(25,*) 
      read(25,*) kbh5
      read(25,*) 
      read(25,*) kbl1
      read(25,*) 
      read(25,*) kbl2
      read(25,*) 
      read(25,*) kbl3
      read(25,*) 
      read(25,*) kbl4
      read(25,*) 
      read(25,*) kbl5
      read(25,*) 
      read(25,*) mbh1
      read(25,*) 
      read(25,*) mbh2
      read(25,*) 
      read(25,*) mbh3
      read(25,*) 
      read(25,*) mbh4
      read(25,*) 
      read(25,*) mbh5
      read(25,*) 
      read(25,*) mbl1
      read(25,*) 
      read(25,*) mbl2
      read(25,*) 
      read(25,*) mbl3
      read(25,*) 
      read(25,*) mbl4
      read(25,*) 
      read(25,*) mbl5
      read(25,*) 
      read(25,*) kzeroa
      read(25,*) 
      read(25,*) kzerob
      read(25,*) 
      read(25,*) kzeroc
      read(25,*) 
      read(25,*) kzerod
      read(25,*) 
      read(25,*) kzeroe
      read(25,*) 
      read(25,*) mzeroa
      read(25,*) 
      read(25,*) mzerob
      read(25,*) 
      read(25,*) mzeroc
      read(25,*) 
      read(25,*) mzerod
      read(25,*) 
      read(25,*) mzeroe
      read(25,*) 
      read(25,*) imztot
      read(25,*) 
      read(25,*) ijval
      read(25,*) 
      read(25,*) jzeroa
      read(25,*) 
      read(25,*) idh1
      read(25,*) 
      read(25,*) idh2
      read(25,*) 
      read(25,*) idh3
      read(25,*) 
      read(25,*) idh4
      read(25,*) 
      read(25,*) idh5
      read(25,*) 
      read(25,*) idl1
      read(25,*) 
      read(25,*) idl2
      read(25,*) 
      read(25,*) idl3
      read(25,*) 
      read(25,*) idl4
      read(25,*) 
      read(25,*) idl5
      read(25,*) 
      read(25,*) ikdeca
      read(25,*) 
      read(25,*) ikdecb
      read(25,*) 
      read(25,*) ikdecc
      read(25,*) 
      read(25,*) ikdecd
      read(25,*) 
      read(25,*) ikdece
      read(25,*) 
      read(25,*) kjdeca
      read(25,*) 
      read(25,*) kjdecb
      read(25,*) 
      read(25,*) kjdecc
      read(25,*) 
      read(25,*) kjdecd
      read(25,*) 
      read(25,*) kjdece
      read(25,*) 
      read(25,*) ijthi
      read(25,*) 
      read(25,*) ijtlo
      read(25,*) 
      read(25,*) jarrdiag
      read(25,*) 
      read(25,*) jhiz1
      read(25,*) 
      read(25,*) jloz1
      read(25,*) 
      read(25,*) iarray
      read(25,*) 
      read(25,*) npdhi
      read(25,*) 
      read(25,*) npdlo
      read(25,*) 
      read(25,*) iialpd
      read(25,*) 
      read(25,*) ipospd
      read(25,*) 
      read(25,*) nkpdterm
      read(25,*) 
      read(25,*) nfrhi
      read(25,*) 
      read(25,*) nfrlo
      read(25,*) 
      read(25,*) nplhi
      read(25,*) 
      read(25,*) npllo
      read(25,*) 
      read(25,*) jspcnfr
      read(25,*) 
      read(25,*) jspnpl
      read(25,*) 
      read(25,*) nknfr
      read(25,*) 
      read(25,*) lossra
      read(25,*) 
      read(25,*) lossrb
      read(25,*) 
      read(25,*) lossrc
      read(25,*) 
      read(25,*) lossrd
      read(25,*) 
      read(25,*) lossre
      read(25,*) 
      read(25,*) nph1
      read(25,*) 
      read(25,*) nph2
      read(25,*) 
      read(25,*) nph3
      read(25,*) 
      read(25,*) nph4
      read(25,*) 
      read(25,*) nph5
      read(25,*) 
      read(25,*) npl1
      read(25,*) 
      read(25,*) npl2
      read(25,*) 
      read(25,*) npl3
      read(25,*) 
      read(25,*) npl4
      read(25,*) 
      read(25,*) npl5
      read(25,*) 
      read(25,*) nolosp
      read(25,*) 
      read(25,*) newfold
      read(25,*) 
      read(25,*) nknlosp
      read(25,*) 
      read(25,*) nknphotrt
      read(25,*) 
      read(25,*) abst2
      read(25,*) 
      read(25,*) errmax
      read(25,*) 
      read(25,*) hmaxday
      read(25,*) 
      read(25,*) timeintv
      read(25,*) 
      read(25,*) abtol
      read(25,*) 
      read(25,*) enqq1
      read(25,*) 
      read(25,*) enqq2
      read(25,*) 
      read(25,*) enqq3
      read(25,*) 
      read(25,*) conp15
      read(25,*) 
      read(25,*) conpst
      read(25,*) 
      read(25,*) pertst2
      read(25,*) 
      read(25,*) aset
      read(25,*) 
      read(25,*) fracpl
      read(25,*) 
      read(25,*) fracnfr

      open(file=trim(physProcEntry),unit=27,form="formatted")
      read(27,*)
      read(27,*) i1, i2, ju1, j2, k1, k2
      read(27,*)
      read(27,*) numQjo, numQks, numQjs, numActive
      allocate(qjGmi(i1:i2, ju1:j2, k1:k2, numQjo))
      allocate(qkGmi(i1:i2, ju1:j2, k1:k2, numQks))
      allocate(qqjda(i1:i2, ju1:j2, k1:k2, numQjs))
      allocate(qqkda(i1:i2, ju1:j2, k1:k2, numQks))
      allocate(yda  (i1:i2, ju1:j2, k1:k2, numActive))
      read(27,*)
      read(27,*) qjGmi(i1:i2, ju1:j2, k1:k2, 1:numQjo)
      read(27,*)
      read(27,*) qkGmi(i1:i2, ju1:j2, k1:k2, 1:numQks)
      read(27,*)
      read(27,*) qqjda(i1:i2, ju1:j2, k1:k2, 1:numQjs)
      read(27,*)
      read(27,*) qqkda(i1:i2, ju1:j2, k1:k2, 1:numQks)
      read(27,*)
      read(27,*) yda  (i1:i2, ju1:j2, k1:k2, 1:numActive)
      read(27,*)
      read(27,*) doQqjkInchem
      read(27,*)
      read(27,*) doSurfEmissInChem
      read(27,*)
      read(27,*) prDiag
      read(27,*)
      read(27,*) prQqjk
      read(27,*)
      read(27,*) prSmv2
      read(27,*)
      read(27,*) localProc
      read(27,*)
      read(27,*) numLat, numLong, numVert
      read(27,*)
      read(27,*) numZones
      read(27,*)
      read(27,*) prNcPeriod
      read(27,*)
      read(27,*) timeStep
      allocate(doCellChem(numZones))
      allocate(thermalRateConstants(numZones, ITHERM))
      allocate(photolysisRateConstants(numZones, IPHOT))
      allocate(surfaceEmissions(numLat*numLong, IGAS))
      allocate(speciesConst(numZones, IGAS))
      read(27,*)
      read(27,*) doCellChem(1:numZones)
      read(27,*)
      read(27,*) thermalRateConstants(1:numZones, 1:ITHERM)
      read(27,*)
      read(27,*) photolysisRateConstants(1:numZones, 1:IPHOT)
      read(27,*)
      read(27,*) surfaceEmissions(1:numLat*numLong, 1:IGAS)
      read(27,*)
      read(27,*) speciesConst(1:numZones, 1:IGAS)

!     ==========
      if (first) then
!     ==========

        first = .false.

        Allocate (cSumA(numZones))
        Allocate (cSumB(numZones))
        cSumA = 0.0d0; cSumB = 0.0d0

      end if

      allocate (jReOrder(numZones))
      allocate (lReOrder(numZones))

      allocate (errorMx2(numZones))

      jReOrder(:) = 0; lReOrder(:) = 0

      errorMx2  (:) = 0.0d0

      call timingOn("Physproc")
!     =============
      call physProc  &
!     =============
     &  (doQqjkInchem, doSurfEmissInChem, prQqjk, prSmv2, numLat,  &
     &   numLong, numVert, ifreord, imgas, initrogen, ioxygen, numZones,  &
     &   kuloop, lunsmv, ncs, fracdec, hmaxnit, prNcPeriod, timeStep,  &
     &   doCellChem, jphotrat, nrates, ntloopncs, ntspec, inewold,  &
     &   npphotrat, thermalRateConstants, photolysisRateConstants, surfaceEmissions, jReOrder, lReOrder, cSumA,  &
     &   cSumB, errorMx2, speciesConst, &
     &   yda, qqkda, qqjda, qkGmi, qjGmi, &
     &   i1, i2, ju1, j2, k1, k2, &
     &   numQjo, numQks, numQjs, numActive)

      call timingOff("Physproc")

!dump smv2chem1 on exit from physproc to be used for standAlone code
      open(file=trim(smv2Chem1Exit),unit=24,form="formatted")
      write(24,*) "ifreord   "
      write(24,*) ifreord
      write(24,*) "ih2o      "
      write(24,*) ih2o
      write(24,*) "imgas     "
      write(24,*) imgas
      write(24,*) "initrogen "
      write(24,*) initrogen
      write(24,*) "ioxygen   "
      write(24,*) ioxygen
      write(24,*) "kuloop    "
      write(24,*) kuloop
      write(24,*) "lunsmv    "
      write(24,*) lunsmv
      write(24,*) "ncs       "
      write(24,*) ncs
      write(24,*) "jphotrat (ICS)            "
      write(24,*) jphotrat
      write(24,*) "nrates   (ICS)            "
      write(24,*) nrates
      write(24,*) "ntloopncs(ICS)            "
      write(24,*) ntloopncs
      write(24,*) "ntspec   (ICS)            "
      write(24,*) ntspec
      write(24,*) "inewold  (MXGSAER, ICS)   "
      write(24,*) inewold
      write(24,*) "npphotrat(IPHOT,   ICS)   "
      write(24,*) npphotrat
      write(24,*) "fracdec   "
      write(24,*) fracdec
      write(24,*) "hmaxnit   "
      write(24,*) hmaxnit

!dump smv2chem2 on exit from physproc to be used for standAlone code
      open(file=trim(smv2Chem2Exit),unit=26,form="formatted")
      write(26,*) "ioner  (ICP)"
      write(26,*) ioner
      write(26,*) "nallrat(ICP)"
      write(26,*) nallrat
      write(26,*) "inorep (ICS)"
      write(26,*) inorep
      write(26,*) "ithrr  (ICS)"
      write(26,*) ithrr
      write(26,*) "itwor  (ICS)"
      write(26,*) itwor
      write(26,*) "nm3bod (ICS)"
      write(26,*) nm3bod
      write(26,*) "nmair  (ICS)"
      write(26,*) nmair
      write(26,*) "nmn2   (ICS)"
      write(26,*) nmn2
      write(26,*) "nmo2 (ICS)"
      write(26,*) nmo2
      write(26,*) "nmoth  (ICS)"
      write(26,*) nmoth
      write(26,*) "ntrates(ICS)"
      write(26,*) ntrates
      write(26,*) "mappl   (MXGSAER, ICS)"
      write(26,*) mappl
      write(26,*) "lgasbino(MAXGL2,  ICS)"
      write(26,*) lgasbino
      write(26,*) "nreacoth(MAXGL2,  ICS)"
      write(26,*) nreacoth
      write(26,*) "lgas3bod(MAXGL3,  ICS)"
      write(26,*) lgas3bod
      write(26,*) "losinacp(MAXGL3,  ICS)"
      write(26,*) losinacp
      write(26,*) "nreac3b (MAXGL3,  ICS)"
      write(26,*) nreac3b
      write(26,*) "nreacair(MAXGL3,  ICS)"
      write(26,*) nreacair
      write(26,*) "nreacn2 (MAXGL3,  ICS)"
      write(26,*) nreacn2
      write(26,*) "nreaco2 (MAXGL3,  ICS)"
      write(26,*) nreaco2
      write(26,*) "jphotnk (NMTRATE, ICS)"
      write(26,*) jphotnk
      write(26,*) "noldfnew(NMTRATE, ICS)"
      write(26,*) noldfnew
      write(26,*) "irm2(NMRPROD, NMTRATE, ICS)"
      write(26,*) irm2
      write(26,*) "ischang(ICS)"
      write(26,*) ischang
      write(26,*) "kzthi(ICP)"
      write(26,*) kzthi
      write(26,*) "kztlo(ICP)"
      write(26,*) kztlo
      write(26,*) "ikztot(MXCOUNT4)"
      write(26,*) ikztot
      write(26,*) "kbh1(MXCOUNT4)"
      write(26,*) kbh1
      write(26,*) "kbh2(MXCOUNT4)"
      write(26,*) kbh2
      write(26,*) "kbh3(MXCOUNT4)"
      write(26,*) kbh3
      write(26,*) "kbh4(MXCOUNT4)"
      write(26,*) kbh4
      write(26,*) "kbh5(MXCOUNT4)"
      write(26,*) kbh5
      write(26,*) "kbl1(MXCOUNT4)"
      write(26,*) kbl1
      write(26,*) "kbl2(MXCOUNT4)"
      write(26,*) kbl2
      write(26,*) "kbl3(MXCOUNT4)"
      write(26,*) kbl3
      write(26,*) "kbl4(MXCOUNT4)"
      write(26,*) kbl4
      write(26,*) "kbl5(MXCOUNT4)"
      write(26,*) kbl5
      write(26,*) "mbh1(MXCOUNT4)"
      write(26,*) mbh1
      write(26,*) "mbh2(MXCOUNT4)"
      write(26,*) mbh2
      write(26,*) "mbh3(MXCOUNT4)"
      write(26,*) mbh3
      write(26,*) "mbh4(MXCOUNT4)"
      write(26,*) mbh4
      write(26,*) "mbh5(MXCOUNT4)"
      write(26,*) mbh5
      write(26,*) "mbl1(MXCOUNT4)"
      write(26,*) mbl1
      write(26,*) "mbl2(MXCOUNT4)"
      write(26,*) mbl2
      write(26,*) "mbl3(MXCOUNT4)"
      write(26,*) mbl3
      write(26,*) "mbl4(MXCOUNT4)"
      write(26,*) mbl4
      write(26,*) "mbl5(MXCOUNT4)"
      write(26,*) mbl5
      write(26,*) "kzeroa(MXCOUNT4)"
      write(26,*) kzeroa
      write(26,*) "kzerob(MXCOUNT4)"
      write(26,*) kzerob
      write(26,*) "kzeroc(MXCOUNT4)"
      write(26,*) kzeroc
      write(26,*) "kzerod(MXCOUNT4)"
      write(26,*) kzerod
      write(26,*) "kzeroe(MXCOUNT4)"
      write(26,*) kzeroe
      write(26,*) "mzeroa(MXCOUNT4)"
      write(26,*) mzeroa
      write(26,*) "mzerob(MXCOUNT4)"
      write(26,*) mzerob
      write(26,*) "mzeroc(MXCOUNT4)"
      write(26,*) mzeroc
      write(26,*) "mzerod(MXCOUNT4)"
      write(26,*) mzerod
      write(26,*) "mzeroe(MXCOUNT4)"
      write(26,*) mzeroe
      write(26,*) "imztot(MXGSAER, ICP)"
      write(26,*) imztot
      write(26,*) "ijval (MXCOUNT3)"
      write(26,*) ijval
      write(26,*) "jzeroa(MXCOUNT3)"
      write(26,*) jzeroa
      write(26,*) "idh1  (MXCOUNT3)"
      write(26,*) idh1
      write(26,*) "idh2  (MXCOUNT3)"
      write(26,*) idh2
      write(26,*) "idh3  (MXCOUNT3)"
      write(26,*) idh3
      write(26,*) "idh4  (MXCOUNT3)"
      write(26,*) idh4
      write(26,*) "idh5  (MXCOUNT3)"
      write(26,*) idh5
      write(26,*) "idl1  (MXCOUNT3)"
      write(26,*) idl1
      write(26,*) "idl2  (MXCOUNT3)"
      write(26,*) idl2
      write(26,*) "idl3  (MXCOUNT3)"
      write(26,*) idl3
      write(26,*) "idl4  (MXCOUNT3)"
      write(26,*) idl4
      write(26,*) "idl5  (MXCOUNT3)"
      write(26,*) idl5
      write(26,*) "ikdeca(MXCOUNT3)"
      write(26,*) ikdeca
      write(26,*) "ikdecb(MXCOUNT3)"
      write(26,*) ikdecb
      write(26,*) "ikdecc(MXCOUNT3)"
      write(26,*) ikdecc
      write(26,*) "ikdecd(MXCOUNT3)"
      write(26,*) ikdecd
      write(26,*) "ikdece(MXCOUNT3)"
      write(26,*) ikdece
      write(26,*) "kjdeca(MXCOUNT3)"
      write(26,*) kjdeca
      write(26,*) "kjdecb(MXCOUNT3)"
      write(26,*) kjdecb
      write(26,*) "kjdecc(MXCOUNT3)"
      write(26,*) kjdecc
      write(26,*) "kjdecd(MXCOUNT3)"
      write(26,*) kjdecd
      write(26,*) "kjdece(MXCOUNT3)"
      write(26,*) kjdece
      write(26,*) "ijthi   (MXGSAER, ICP)"
      write(26,*) ijthi
      write(26,*) "ijtlo(MXGSAER, ICP)"
      write(26,*) ijtlo
      write(26,*) "jarrdiag(MXGSAER, ICP)"
      write(26,*) jarrdiag
      write(26,*) "jhiz1   (MXGSAER, ICP)"
      write(26,*) jhiz1
      write(26,*) "jloz1(MXGSAER, ICP)"
      write(26,*) jloz1
      write(26,*) "iarray(ICP)"
      write(26,*) iarray
      write(26,*) "npdhi (ICP)"
      write(26,*) npdhi
      write(26,*) "npdlo(ICP)"
      write(26,*) npdlo
      write(26,*) "iialpd  (MXCOUNT2)"
      write(26,*) iialpd
      write(26,*) "ipospd  (MXCOUNT2)"
      write(26,*) ipospd
      write(26,*) "nkpdterm(MXCOUNT2)"
      write(26,*) nkpdterm
      write(26,*) "nfrhi(ICP)"
      write(26,*) nfrhi
      write(26,*) "nfrlo(ICP)"
      write(26,*) nfrlo
      write(26,*) "nplhi(ICP)"
      write(26,*) nplhi
      write(26,*) "npllo(ICP)"
      write(26,*) npllo
      write(26,*) "jspcnfr(MXCOUNT4)"
      write(26,*) jspcnfr
      write(26,*) "jspnpl(MXCOUNT4)"
      write(26,*) jspnpl
      write(26,*) "nknfr  (MXCOUNT4)"
      write(26,*) nknfr
      write(26,*) "lossra (MXCOUNT4)"
      write(26,*) lossra
      write(26,*) "lossrb (MXCOUNT4)"
      write(26,*) lossrb
      write(26,*) "lossrc(MXCOUNT4)"
      write(26,*) lossrc
      write(26,*) "lossrd (MXCOUNT4)"
      write(26,*) lossrd
      write(26,*) "lossre(MXCOUNT4)"
      write(26,*) lossre
      write(26,*) "nph1(MXCOUNT4)"
      write(26,*) nph1
      write(26,*) "nph2(MXCOUNT4)"
      write(26,*) nph2
      write(26,*) "nph3(MXCOUNT4)"
      write(26,*) nph3
      write(26,*) "nph4(MXCOUNT4)"
      write(26,*) nph4
      write(26,*) "nph5(MXCOUNT4)"
      write(26,*) nph5
      write(26,*) "npl1(MXCOUNT4)"
      write(26,*) npl1
      write(26,*) "npl2(MXCOUNT4)"
      write(26,*) npl2
      write(26,*) "npl3(MXCOUNT4)"
      write(26,*) npl3
      write(26,*) "npl4(MXCOUNT4)"
      write(26,*) npl4
      write(26,*) "npl5(MXCOUNT4)"
      write(26,*) npl5
      write(26,*) "nolosp(ICP)"
      write(26,*) nolosp
      write(26,*) "newfold(NMTRATE*2, ICS)"
      write(26,*) newfold
      write(26,*) "nknlosp(MAXGL3, ICS)"
      write(26,*) nknlosp
      write(26,*) "nknphotrt(IPHOT,ICS)"
      write(26,*) nknphotrt
      write(26,*) "abst2   (ICS)"
      write(26,*) abst2
      write(26,*) "errmax  (ICS)"
      write(26,*) errmax
      write(26,*) "hmaxday (ICS)"
      write(26,*) hmaxday
      write(26,*) "timeintv(ICS)"
      write(26,*) timeintv
      write(26,*) "abtol(6, ICS)"
      write(26,*) abtol
      write(26,*) "enqq1 (MORDER)"
      write(26,*) enqq1
      write(26,*) "enqq2 (MORDER)"
      write(26,*) enqq2
      write(26,*) "enqq3 (MORDER)"
      write(26,*) enqq3
      write(26,*) "conp15(MORDER)"
      write(26,*) conp15
      write(26,*) "conpst(MORDER)"
      write(26,*) conpst
      write(26,*) "pertst2(MORDER, 3)"
      write(26,*) pertst2
      write(26,*) "aset(10, 8)"
      write(26,*) aset
      write(26,*) "fracpl (MXCOUNT2)"
      write(26,*) fracpl
      write(26,*) "fracnfr(MXCOUNT4)"
      write(26,*) fracnfr

      open(file=trim(physProcExit),unit=28,form="formatted")
      write(28,*) "i1, i2, ju1, j2, k1, k2"
      write(28,*) i1, i2, ju1, j2, k1, k2
      write(28,*) "num_qjo, num_qks, num_qjs, num_active"
      write(28,*) numQjo, numQks, numQjs, numActive
      write(28,*) "qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)"
      write(28,*) qjGmi(i1:i2, ju1:j2, k1:k2, 1:numQjo)
      write(28,*) "qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)"
      write(28,*) qkGmi(i1:i2, ju1:j2, k1:k2, 1:numQks)
      write(28,*) "qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)"
      write(28,*) qqjda(i1:i2, ju1:j2, k1:k2, 1:numQjs)
      write(28,*) "qqkda(i1:i2, ju1:j2, k1:k2, num_qks)"
      write(28,*) qqkda(i1:i2, ju1:j2, k1:k2, 1:numQks)
      write(28,*) "yda  (i1:i2, ju1:j2, k1:k2, num_active)"
      write(28,*) yda  (i1:i2, ju1:j2, k1:k2, 1:numActive)
      write(28,*) "do_qqjk_inchem"
      write(28,*) doQqjkInchem
      write(28,*) "do_semiss_inchem"
      write(28,*) doSurfEmissInChem
      write(28,*) "pr_diag"
      write(28,*) prDiag
      write(28,*) "pr_qqjk"
      write(28,*) prQqjk
      write(28,*) "pr_smv2"
      write(28,*) prSmv2
      write(28,*) "loc_proc"
      write(28,*) localProc
      write(28,*) "ilat, ilong, ivert"
      write(28,*) numLat, numLong, numVert
      write(28,*) "itloop"
      write(28,*) numZones
      write(28,*) "pr_nc_period"
      write(28,*) prNcPeriod
      write(28,*) "tdt"
      write(28,*) timeStep
      write(28,*) "do_cell_chem(itloop)"
      write(28,*) doCellChem(1:numZones)
      write(28,*) "arate(itloop, ITHERM)"
      write(28,*) thermalRateConstants(1:numZones, 1:ITHERM)
      write(28,*) "prate(itloop, IPHOT)"
      write(28,*) photolysisRateConstants(1:numZones, 1:IPHOT)
      write(28,*) "yemis(ilat*ilong, IGAS)"
      write(28,*) surfaceEmissions(1:numLat*numLong, 1:IGAS)
      write(28,*) "cx(itloop, IGAS)"
      write(28,*) speciesConst(1:numZones, 1:IGAS)

      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
      close(28)

      deallocate (jReOrder)
      deallocate (lReOrder)
      deallocate (errorMx2)
      deallocate(qjGmi)
      deallocate(qkGmi)
      deallocate(qqjda)
      deallocate(qqkda)
      deallocate(yda)
      deallocate(doCellChem)
      deallocate(thermalRateConstants)
      deallocate(photolysisRateConstants)
      deallocate(surfaceEmissions)
      deallocate(speciesConst)

      call timingPrint

      print*, "Exiting doSmv2Solver"

      end program doSmv2Solver

