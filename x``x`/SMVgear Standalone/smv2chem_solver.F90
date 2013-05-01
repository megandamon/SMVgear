!-----------------------------------------------------------------------------
!
! ROUTINE
!   doSmv2Solver
!
! DESCRIPTION
!   This is the main control routine for the ordinary differential equation
!   solver, "Smvgear II" (Sparse Matrix Vectorized Gear-type code).
!-----------------------------------------------------------------------------

   program doSmv2Solver
      use pFUnit
      use timing_mod
      use SmvChem_mod
      implicit none
      external testFailOnPurpose
      external testSpeciesConst

#     include "smv2chem_par.h"
#     include "smv2chem2.h"
#     include "smv2chem1.h"

      logical, save :: first = .true.
      integer :: rank,errorInt
      integer, allocatable :: jReOrder(:)
      integer, allocatable :: lReOrder(:)
      real*8, allocatable  :: errorMx2  (:)
      real*8, allocatable, save :: cSumA(:)
      real*8, allocatable, save :: cSumB(:)
      character(len=100) :: smv2Chem1Entry
      character(len=100) :: smv2Chem1Exit
      character(len=100) :: smv2Chem2Entry
      character(len=100) :: smv2Chem2Exit
      character(len=100) :: physProcEntry
      character(len=100) :: physProcExit
      type (SmvChem_type) :: chemObject
      type (TestSuite_type) :: suite
      type (TestResult_type) :: result
      character(len=100) :: summary_statement

      call pFUnit_init()

!      call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
      call timingInit

      rank = 17
      if (chemObject%prDiag) then
        Write (6,*) 'doSmv2Solver called by ', chemObject%localProc
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

      call readSmv2Chem1Entry (smv2Chem1Entry)
      call readSmv2Chem2Entry (smv2Chem2Entry)
      call readPhysProc (chemObject, physProcEntry)

      print*, "read from: ", trim(physProcEntry), " prDiag = ", chemObject%prDiag
      chemObject%prDiag = 1

      if (first) then
         first = .false.
        Allocate (cSumA(chemObject%numZones))
        Allocate (cSumB(chemObject%numZones))
        cSumA = 0.0d0; cSumB = 0.0d0
      end if

      allocate (jReOrder(chemObject%numZones))
      allocate (lReOrder(chemObject%numZones))
      allocate (errorMx2(chemObject%numZones))
      jReOrder(:) = 0; lReOrder(:) = 0
      errorMx2  (:) = 0.0d0

      suite = TestSuite('smvgear tests')
      call add(suite, TestCase1Step('testFailOnPurpose', testFailOnPurpose))


      call timingOn("Physproc")
      call physProc  &
     &  (chemObject%doQqjkInchem, chemObject%doSurfEmissInChem, chemObject%prQqjk, chemObject%prSmv2, chemObject%numLat,  &
     &   chemObject%numLong, chemObject%numVert, ifreord, imgas, initrogen, ioxygen, chemObject%numZones,  &
     &   kuloop, lunsmv, ncs, fracdec, hmaxnit, chemObject%prNcPeriod, chemObject%timeStep,  &
     &   chemObject%doCellChem, jphotrat, nrates, ntloopncs, ntspec, inewold,  &
     &   npphotrat, chemObject%thermalRateConstants, chemObject%photolysisRateConstants, chemObject%surfaceEmissions, jReOrder, lReOrder, cSumA,  &
     &   cSumB, errorMx2, chemObject%speciesConst, &
     &   chemObject%yda, chemObject%qqkda, chemObject%qqjda, chemObject%qkGmi, chemObject%qjGmi, &
     &   chemObject%i1, chemObject%i2, chemObject%ju1, chemObject%j2, chemObject%k1, chemObject%k2, &
     &   chemObject%numQjo, chemObject%numQks, chemObject%numQjs, chemObject%numActive)
      call timingOff("Physproc")

      !call add(suite, TestCase1Step('testSpeciesConst'), testSpeciesConst, chemObject%speciesConst(1,1), 1077508322.73440)
      !speciesConst(1,1):   1077508322.73440
      !speciesConst(numZones,IGAS):   540086297554450.

      call writeSmv2Chem1Exit (smv2Chem1Exit)
      call writeSmv2Chem2Exit (smv2Chem2Exit)
      call writePhysProc (chemObject, physProcEntry)

      deallocate (jReOrder)
      deallocate (lReOrder)
      deallocate (errorMx2)
      call deallocateVariables(chemObject)

      ! Run the tests and accumulate the results in "result"
      result = newTestResult(mode=MODE_USE_STDOUT)
      call Run(suite, result)
      summary_statement=Summary(result)
      print*,trim(summary_statement)

      call clean(result)
      call clean(suite)
      call pFUnit_finalize()

      ! call timingPrint
      print*, "Exiting doSmv2Solver"

   end program doSmv2Solver






      subroutine readSmv2Chem1Entry (fileName)
         implicit none

#     include "smv2chem_par.h"
#     include "smv2chem1.h"

         ! Arguments
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Reading from: ", fileName
         open(newunit=fileNumber, file=trim(fileName),form="formatted")

         read(fileNumber,'(a)')
         read(fileNumber,*) ifreord
         read(fileNumber,'(a)')
         read(fileNumber,*) ih2o
         read(fileNumber,'(a)')
         read(fileNumber,*) imgas
         read(fileNumber,'(a)')
         read(fileNumber,*) initrogen
         read(fileNumber,'(a)')
         read(fileNumber,*) ioxygen
         read(fileNumber,'(a)')
         read(fileNumber,*) kuloop
         read(fileNumber,'(a)')
         read(fileNumber,*) lunsmv
         read(fileNumber,'(a)')
         read(fileNumber,*) ncs
         read(fileNumber,'(a)')
         read(fileNumber,*) jphotrat
         read(fileNumber,'(a)')
         read(fileNumber,*) nrates
         read(fileNumber,'(a)')
         read(fileNumber,*) ntloopncs
         read(fileNumber,'(a)')
         read(fileNumber,*) ntspec
         read(fileNumber,'(a)')
         read(fileNumber,*) inewold
         read(fileNumber,'(a)')
         read(fileNumber,*) npphotrat
         read(fileNumber,'(a)')
         read(fileNumber,*) fracdec
         read(fileNumber,'(a)')
         read(fileNumber,*) hmaxnit

         close(fileNumber)

      end subroutine readSmv2Chem1Entry


      subroutine readSmv2Chem2Entry (fileName)
         implicit none

#     include "smv2chem_par.h"
#     include "smv2chem2.h"

         ! Arguments
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Reading from: ", fileName
         open(newunit=fileNumber, file=trim(fileName),form="formatted")

         read(fileNumber,'(a)')
         read(fileNumber,*) ioner
         read(fileNumber,'(a)')
         read(fileNumber,*) nallrat
         read(fileNumber,'(a)')
         read(fileNumber,*) inorep
         read(fileNumber,'(a)')
         read(fileNumber,*) ithrr
         read(fileNumber,'(a)')
         read(fileNumber,*) itwor
         read(fileNumber,'(a)')
         read(fileNumber,*) nm3bod
         read(fileNumber,'(a)')
         read(fileNumber,*) nmair
         read(fileNumber,'(a)')
         read(fileNumber,*) nmn2
         read(fileNumber,'(a)')
         read(fileNumber,*) nmo2
         read(fileNumber,'(a)')
         read(fileNumber,*) nmoth
         read(fileNumber,'(a)')
         read(fileNumber,*) ntrates
         read(fileNumber,'(a)')
         read(fileNumber,*) mappl
         read(fileNumber,'(a)')
         read(fileNumber,*) lgasbino
         read(fileNumber,'(a)')
         read(fileNumber,*) nreacoth
         read(fileNumber,'(a)')
         read(fileNumber,*) lgas3bod
         read(fileNumber,'(a)')
         read(fileNumber,*) losinacp
         read(fileNumber,'(a)')
         read(fileNumber,*) nreac3b
         read(fileNumber,'(a)')
         read(fileNumber,*) nreacair
         read(fileNumber,'(a)')
         read(fileNumber,*) nreacn2
         read(fileNumber,'(a)')
         read(fileNumber,*) nreaco2
         read(fileNumber,'(a)')
         read(fileNumber,*) jphotnk
         read(fileNumber,'(a)')
         read(fileNumber,*) noldfnew
         read(fileNumber,'(a)')
         read(fileNumber,*) irm2
         read(fileNumber,'(a)')
         read(fileNumber,*) ischang
         read(fileNumber,'(a)')
         read(fileNumber,*) kzthi
         read(fileNumber,'(a)')
         read(fileNumber,*) kztlo
         read(fileNumber,'(a)')
         read(fileNumber,*) ikztot
         read(fileNumber,'(a)')
         read(fileNumber,*) kbh1
         read(fileNumber,'(a)')
         read(fileNumber,*) kbh2
         read(fileNumber,'(a)')
         read(fileNumber,*) kbh3
         read(fileNumber,'(a)')
         read(fileNumber,*) kbh4
         read(fileNumber,'(a)')
         read(fileNumber,*) kbh5
         read(fileNumber,'(a)')
         read(fileNumber,*) kbl1
         read(fileNumber,'(a)')
         read(fileNumber,*) kbl2
         read(fileNumber,'(a)')
         read(fileNumber,*) kbl3
         read(fileNumber,'(a)')
         read(fileNumber,*) kbl4
         read(fileNumber,'(a)')
         read(fileNumber,*) kbl5
         read(fileNumber,'(a)')
         read(fileNumber,*) mbh1
         read(fileNumber,'(a)')
         read(fileNumber,*) mbh2
         read(fileNumber,'(a)')
         read(fileNumber,*) mbh3
         read(fileNumber,'(a)')
         read(fileNumber,*) mbh4
         read(fileNumber,'(a)')
         read(fileNumber,*) mbh5
         read(fileNumber,'(a)')
         read(fileNumber,*) mbl1
         read(fileNumber,'(a)')
         read(fileNumber,*) mbl2
         read(fileNumber,'(a)')
         read(fileNumber,*) mbl3
         read(fileNumber,'(a)')
         read(fileNumber,*) mbl4
         read(fileNumber,'(a)')
         read(fileNumber,*) mbl5
         read(fileNumber,'(a)')
         read(fileNumber,*) kzeroa
         read(fileNumber,'(a)')
         read(fileNumber,*) kzerob
         read(fileNumber,'(a)')
         read(fileNumber,*) kzeroc
         read(fileNumber,'(a)')
         read(fileNumber,*) kzerod
         read(fileNumber,'(a)')
         read(fileNumber,*) kzeroe
         read(fileNumber,'(a)')
         read(fileNumber,*) mzeroa
         read(fileNumber,'(a)')
         read(fileNumber,*) mzerob
         read(fileNumber,'(a)')
         read(fileNumber,*) mzeroc
         read(fileNumber,'(a)')
         read(fileNumber,*) mzerod
         read(fileNumber,'(a)')
         read(fileNumber,*) mzeroe
         read(fileNumber,'(a)')
         read(fileNumber,*) imztot
         read(fileNumber,'(a)')
         read(fileNumber,*) ijval
         read(fileNumber,'(a)')
         read(fileNumber,*) jzeroa
         read(fileNumber,'(a)')
         read(fileNumber,*) idh1
         read(fileNumber,'(a)')
         read(fileNumber,*) idh2
         read(fileNumber,'(a)')
         read(fileNumber,*) idh3
         read(fileNumber,'(a)')
         read(fileNumber,*) idh4
         read(fileNumber,'(a)')
         read(fileNumber,*) idh5
         read(fileNumber,'(a)')
         read(fileNumber,*) idl1
         read(fileNumber,'(a)')
         read(fileNumber,*) idl2
         read(fileNumber,'(a)')
         read(fileNumber,*) idl3
         read(fileNumber,'(a)')
         read(fileNumber,*) idl4
         read(fileNumber,'(a)')
         read(fileNumber,*) idl5
         read(fileNumber,'(a)')
         read(fileNumber,*) ikdeca
         read(fileNumber,'(a)')
         read(fileNumber,*) ikdecb
         read(fileNumber,'(a)')
         read(fileNumber,*) ikdecc
         read(fileNumber,'(a)')
         read(fileNumber,*) ikdecd
         read(fileNumber,'(a)')
         read(fileNumber,*) ikdece
         read(fileNumber,'(a)')
         read(fileNumber,*) kjdeca
         read(fileNumber,'(a)')
         read(fileNumber,*) kjdecb
         read(fileNumber,'(a)')
         read(fileNumber,*) kjdecc
         read(fileNumber,'(a)')
         read(fileNumber,*) kjdecd
         read(fileNumber,'(a)')
         read(fileNumber,*) kjdece
         read(fileNumber,'(a)')
         read(fileNumber,*) ijthi
         read(fileNumber,'(a)')
         read(fileNumber,*) ijtlo
         read(fileNumber,'(a)')
         read(fileNumber,*) jarrdiag
         read(fileNumber,'(a)')
         read(fileNumber,*) jhiz1
         read(fileNumber,'(a)')
         read(fileNumber,*) jloz1
         read(fileNumber,'(a)')
         read(fileNumber,*) iarray
         read(fileNumber,'(a)')
         read(fileNumber,*) npdhi
         read(fileNumber,'(a)')
         read(fileNumber,*) npdlo
         read(fileNumber,'(a)')
         read(fileNumber,*) iialpd
         read(fileNumber,'(a)')
         read(fileNumber,*) ipospd
         read(fileNumber,'(a)')
         read(fileNumber,*) nkpdterm
         read(fileNumber,'(a)')
         read(fileNumber,*) nfrhi
         read(fileNumber,'(a)')
         read(fileNumber,*) nfrlo
         read(fileNumber,'(a)')
         read(fileNumber,*) nplhi
         read(fileNumber,'(a)')
         read(fileNumber,*) npllo
         read(fileNumber,'(a)')
         read(fileNumber,*) jspcnfr
         read(fileNumber,'(a)')
         read(fileNumber,*) jspnpl
         read(fileNumber,'(a)')
         read(fileNumber,*) nknfr
         read(fileNumber,'(a)')
         read(fileNumber,*) lossra
         read(fileNumber,'(a)')
         read(fileNumber,*) lossrb
         read(fileNumber,'(a)')
         read(fileNumber,*) lossrc
         read(fileNumber,'(a)')
         read(fileNumber,*) lossrd
         read(fileNumber,'(a)')
         read(fileNumber,*) lossre
         read(fileNumber,'(a)')
         read(fileNumber,*) nph1
         read(fileNumber,'(a)')
         read(fileNumber,*) nph2
         read(fileNumber,'(a)')
         read(fileNumber,*) nph3
         read(fileNumber,'(a)')
         read(fileNumber,*) nph4
         read(fileNumber,'(a)')
         read(fileNumber,*) nph5
         read(fileNumber,'(a)')
         read(fileNumber,*) npl1
         read(fileNumber,'(a)')
         read(fileNumber,*) npl2
         read(fileNumber,'(a)')
         read(fileNumber,*) npl3
         read(fileNumber,'(a)')
         read(fileNumber,*) npl4
         read(fileNumber,'(a)')
         read(fileNumber,*) npl5
         read(fileNumber,'(a)')
         read(fileNumber,*) nolosp
         read(fileNumber,'(a)')
         read(fileNumber,*) newfold
         read(fileNumber,'(a)')
         read(fileNumber,*) nknlosp
         read(fileNumber,'(a)')
         read(fileNumber,*) nknphotrt
         read(fileNumber,'(a)')
         read(fileNumber,*) abst2
         read(fileNumber,'(a)')
         read(fileNumber,*) errmax
         read(fileNumber,'(a)')
         read(fileNumber,*) hmaxday
         read(fileNumber,'(a)')
         read(fileNumber,*) timeintv
         read(fileNumber,'(a)')
         read(fileNumber,*) abtol
         read(fileNumber,'(a)')
         read(fileNumber,*) enqq1
         read(fileNumber,'(a)')
         read(fileNumber,*) enqq2
         read(fileNumber,'(a)')
         read(fileNumber,*) enqq3
         read(fileNumber,'(a)')
         read(fileNumber,*) conp15
         read(fileNumber,'(a)')
         read(fileNumber,*) conpst
         read(fileNumber,'(a)')
         read(fileNumber,*) pertst2
         read(fileNumber,'(a)')
         read(fileNumber,*) aset
         read(fileNumber,'(a)')
         read(fileNumber,*) fracpl
         read(fileNumber,'(a)')
         read(fileNumber,*) fracnfr

         close(fileNumber)

      end subroutine readSmv2Chem2Entry




      subroutine writeSmv2Chem1Exit (fileName)
         implicit none

#     include "smv2chem_par.h"
#     include "smv2chem1.h"

         ! Arguments
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Writing to: ", fileName
         open(newunit=fileNumber, file=trim(fileName),form="formatted")

         write(fileNumber,*) "ifreord   "
         write(fileNumber,*) ifreord
         write(fileNumber,*) "ih2o      "
         write(fileNumber,*) ih2o
         write(fileNumber,*) "imgas     "
         write(fileNumber,*) imgas
         write(fileNumber,*) "initrogen "
         write(fileNumber,*) initrogen
         write(fileNumber,*) "ioxygen   "
         write(fileNumber,*) ioxygen
         write(fileNumber,*) "kuloop    "
         write(fileNumber,*) kuloop
         write(fileNumber,*) "lunsmv    "
         write(fileNumber,*) lunsmv
         write(fileNumber,*) "ncs       "
         write(fileNumber,*) ncs
         write(fileNumber,*) "jphotrat (ICS)            "
         write(fileNumber,*) jphotrat
         write(fileNumber,*) "nrates   (ICS)            "
         write(fileNumber,*) nrates
         write(fileNumber,*) "ntloopncs(ICS)            "
         write(fileNumber,*) ntloopncs
         write(fileNumber,*) "ntspec   (ICS)            "
         write(fileNumber,*) ntspec
         write(fileNumber,*) "inewold  (MXGSAER, ICS)   "
         write(fileNumber,*) inewold
         write(fileNumber,*) "npphotrat(IPHOT,   ICS)   "
         write(fileNumber,*) npphotrat
         write(fileNumber,*) "fracdec   "
         write(fileNumber,*) fracdec
         write(fileNumber,*) "hmaxnit   "
         write(fileNumber,*) hmaxnit

         close(fileNumber)

      end subroutine writeSmv2Chem1Exit


      subroutine writeSmv2Chem2Exit (fileName)
         implicit none

#     include "smv2chem_par.h"
#     include "smv2chem2.h"

         ! Arguments
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Writing to: ", fileName
         open(newunit=fileNumber, file=trim(fileName),form="formatted")

         write(fileNumber,*) "ioner  (ICP)"
         write(fileNumber,*) ioner
         write(fileNumber,*) "nallrat(ICP)"
         write(fileNumber,*) nallrat
         write(fileNumber,*) "inorep (ICS)"
         write(fileNumber,*) inorep
         write(fileNumber,*) "ithrr  (ICS)"
         write(fileNumber,*) ithrr
         write(fileNumber,*) "itwor  (ICS)"
         write(fileNumber,*) itwor
         write(fileNumber,*) "nm3bod (ICS)"
         write(fileNumber,*) nm3bod
         write(fileNumber,*) "nmair  (ICS)"
         write(fileNumber,*) nmair
         write(fileNumber,*) "nmn2   (ICS)"
         write(fileNumber,*) nmn2
         write(fileNumber,*) "nmo2 (ICS)"
         write(fileNumber,*) nmo2
         write(fileNumber,*) "nmoth  (ICS)"
         write(fileNumber,*) nmoth
         write(fileNumber,*) "ntrates(ICS)"
         write(fileNumber,*) ntrates
         write(fileNumber,*) "mappl   (MXGSAER, ICS)"
         write(fileNumber,*) mappl
         write(fileNumber,*) "lgasbino(MAXGL2,  ICS)"
         write(fileNumber,*) lgasbino
         write(fileNumber,*) "nreacoth(MAXGL2,  ICS)"
         write(fileNumber,*) nreacoth
         write(fileNumber,*) "lgas3bod(MAXGL3,  ICS)"
         write(fileNumber,*) lgas3bod
         write(fileNumber,*) "losinacp(MAXGL3,  ICS)"
         write(fileNumber,*) losinacp
         write(fileNumber,*) "nreac3b (MAXGL3,  ICS)"
         write(fileNumber,*) nreac3b
         write(fileNumber,*) "nreacair(MAXGL3,  ICS)"
         write(fileNumber,*) nreacair
         write(fileNumber,*) "nreacn2 (MAXGL3,  ICS)"
         write(fileNumber,*) nreacn2
         write(fileNumber,*) "nreaco2 (MAXGL3,  ICS)"
         write(fileNumber,*) nreaco2
         write(fileNumber,*) "jphotnk (NMTRATE, ICS)"
         write(fileNumber,*) jphotnk
         write(fileNumber,*) "noldfnew(NMTRATE, ICS)"
         write(fileNumber,*) noldfnew
         write(fileNumber,*) "irm2(NMRPROD, NMTRATE, ICS)"
         write(fileNumber,*) irm2
         write(fileNumber,*) "ischang(ICS)"
         write(fileNumber,*) ischang
         write(fileNumber,*) "kzthi(ICP)"
         write(fileNumber,*) kzthi
         write(fileNumber,*) "kztlo(ICP)"
         write(fileNumber,*) kztlo
         write(fileNumber,*) "ikztot(MXCOUNT4)"
         write(fileNumber,*) ikztot
         write(fileNumber,*) "kbh1(MXCOUNT4)"
         write(fileNumber,*) kbh1
         write(fileNumber,*) "kbh2(MXCOUNT4)"
         write(fileNumber,*) kbh2
         write(fileNumber,*) "kbh3(MXCOUNT4)"
         write(fileNumber,*) kbh3
         write(fileNumber,*) "kbh4(MXCOUNT4)"
         write(fileNumber,*) kbh4
         write(fileNumber,*) "kbh5(MXCOUNT4)"
         write(fileNumber,*) kbh5
         write(fileNumber,*) "kbl1(MXCOUNT4)"
         write(fileNumber,*) kbl1
         write(fileNumber,*) "kbl2(MXCOUNT4)"
         write(fileNumber,*) kbl2
         write(fileNumber,*) "kbl3(MXCOUNT4)"
         write(fileNumber,*) kbl3
         write(fileNumber,*) "kbl4(MXCOUNT4)"
         write(fileNumber,*) kbl4
         write(fileNumber,*) "kbl5(MXCOUNT4)"
         write(fileNumber,*) kbl5
         write(fileNumber,*) "mbh1(MXCOUNT4)"
         write(fileNumber,*) mbh1
         write(fileNumber,*) "mbh2(MXCOUNT4)"
         write(fileNumber,*) mbh2
         write(fileNumber,*) "mbh3(MXCOUNT4)"
         write(fileNumber,*) mbh3
         write(fileNumber,*) "mbh4(MXCOUNT4)"
         write(fileNumber,*) mbh4
         write(fileNumber,*) "mbh5(MXCOUNT4)"
         write(fileNumber,*) mbh5
         write(fileNumber,*) "mbl1(MXCOUNT4)"
         write(fileNumber,*) mbl1
         write(fileNumber,*) "mbl2(MXCOUNT4)"
         write(fileNumber,*) mbl2
         write(fileNumber,*) "mbl3(MXCOUNT4)"
         write(fileNumber,*) mbl3
         write(fileNumber,*) "mbl4(MXCOUNT4)"
         write(fileNumber,*) mbl4
         write(fileNumber,*) "mbl5(MXCOUNT4)"
         write(fileNumber,*) mbl5
         write(fileNumber,*) "kzeroa(MXCOUNT4)"
         write(fileNumber,*) kzeroa
         write(fileNumber,*) "kzerob(MXCOUNT4)"
         write(fileNumber,*) kzerob
         write(fileNumber,*) "kzeroc(MXCOUNT4)"
         write(fileNumber,*) kzeroc
         write(fileNumber,*) "kzerod(MXCOUNT4)"
         write(fileNumber,*) kzerod
         write(fileNumber,*) "kzeroe(MXCOUNT4)"
         write(fileNumber,*) kzeroe
         write(fileNumber,*) "mzeroa(MXCOUNT4)"
         write(fileNumber,*) mzeroa
         write(fileNumber,*) "mzerob(MXCOUNT4)"
         write(fileNumber,*) mzerob
         write(fileNumber,*) "mzeroc(MXCOUNT4)"
         write(fileNumber,*) mzeroc
         write(fileNumber,*) "mzerod(MXCOUNT4)"
         write(fileNumber,*) mzerod
         write(fileNumber,*) "mzeroe(MXCOUNT4)"
         write(fileNumber,*) mzeroe
         write(fileNumber,*) "imztot(MXGSAER, ICP)"
         write(fileNumber,*) imztot
         write(fileNumber,*) "ijval (MXCOUNT3)"
         write(fileNumber,*) ijval
         write(fileNumber,*) "jzeroa(MXCOUNT3)"
         write(fileNumber,*) jzeroa
         write(fileNumber,*) "idh1  (MXCOUNT3)"
         write(fileNumber,*) idh1
         write(fileNumber,*) "idh2  (MXCOUNT3)"
         write(fileNumber,*) idh2
         write(fileNumber,*) "idh3  (MXCOUNT3)"
         write(fileNumber,*) idh3
         write(fileNumber,*) "idh4  (MXCOUNT3)"
         write(fileNumber,*) idh4
         write(fileNumber,*) "idh5  (MXCOUNT3)"
         write(fileNumber,*) idh5
         write(fileNumber,*) "idl1  (MXCOUNT3)"
         write(fileNumber,*) idl1
         write(fileNumber,*) "idl2  (MXCOUNT3)"
         write(fileNumber,*) idl2
         write(fileNumber,*) "idl3  (MXCOUNT3)"
         write(fileNumber,*) idl3
         write(fileNumber,*) "idl4  (MXCOUNT3)"
         write(fileNumber,*) idl4
         write(fileNumber,*) "idl5  (MXCOUNT3)"
         write(fileNumber,*) idl5
         write(fileNumber,*) "ikdeca(MXCOUNT3)"
         write(fileNumber,*) ikdeca
         write(fileNumber,*) "ikdecb(MXCOUNT3)"
         write(fileNumber,*) ikdecb
         write(fileNumber,*) "ikdecc(MXCOUNT3)"
         write(fileNumber,*) ikdecc
         write(fileNumber,*) "ikdecd(MXCOUNT3)"
         write(fileNumber,*) ikdecd
         write(fileNumber,*) "ikdece(MXCOUNT3)"
         write(fileNumber,*) ikdece
         write(fileNumber,*) "kjdeca(MXCOUNT3)"
         write(fileNumber,*) kjdeca
         write(fileNumber,*) "kjdecb(MXCOUNT3)"
         write(fileNumber,*) kjdecb
         write(fileNumber,*) "kjdecc(MXCOUNT3)"
         write(fileNumber,*) kjdecc
         write(fileNumber,*) "kjdecd(MXCOUNT3)"
         write(fileNumber,*) kjdecd
         write(fileNumber,*) "kjdece(MXCOUNT3)"
         write(fileNumber,*) kjdece
         write(fileNumber,*) "ijthi   (MXGSAER, ICP)"
         write(fileNumber,*) ijthi
         write(fileNumber,*) "ijtlo(MXGSAER, ICP)"
         write(fileNumber,*) ijtlo
         write(fileNumber,*) "jarrdiag(MXGSAER, ICP)"
         write(fileNumber,*) jarrdiag
         write(fileNumber,*) "jhiz1   (MXGSAER, ICP)"
         write(fileNumber,*) jhiz1
         write(fileNumber,*) "jloz1(MXGSAER, ICP)"
         write(fileNumber,*) jloz1
         write(fileNumber,*) "iarray(ICP)"
         write(fileNumber,*) iarray
         write(fileNumber,*) "npdhi (ICP)"
         write(fileNumber,*) npdhi
         write(fileNumber,*) "npdlo(ICP)"
         write(fileNumber,*) npdlo
         write(fileNumber,*) "iialpd  (MXCOUNT2)"
         write(fileNumber,*) iialpd
         write(fileNumber,*) "ipospd  (MXCOUNT2)"
         write(fileNumber,*) ipospd
         write(fileNumber,*) "nkpdterm(MXCOUNT2)"
         write(fileNumber,*) nkpdterm
         write(fileNumber,*) "nfrhi(ICP)"
         write(fileNumber,*) nfrhi
         write(fileNumber,*) "nfrlo(ICP)"
         write(fileNumber,*) nfrlo
         write(fileNumber,*) "nplhi(ICP)"
         write(fileNumber,*) nplhi
         write(fileNumber,*) "npllo(ICP)"
         write(fileNumber,*) npllo
         write(fileNumber,*) "jspcnfr(MXCOUNT4)"
         write(fileNumber,*) jspcnfr
         write(fileNumber,*) "jspnpl(MXCOUNT4)"
         write(fileNumber,*) jspnpl
         write(fileNumber,*) "nknfr  (MXCOUNT4)"
         write(fileNumber,*) nknfr
         write(fileNumber,*) "lossra (MXCOUNT4)"
         write(fileNumber,*) lossra
         write(fileNumber,*) "lossrb (MXCOUNT4)"
         write(fileNumber,*) lossrb
         write(fileNumber,*) "lossrc(MXCOUNT4)"
         write(fileNumber,*) lossrc
         write(fileNumber,*) "lossrd (MXCOUNT4)"
         write(fileNumber,*) lossrd
         write(fileNumber,*) "lossre(MXCOUNT4)"
         write(fileNumber,*) lossre
         write(fileNumber,*) "nph1(MXCOUNT4)"
         write(fileNumber,*) nph1
         write(fileNumber,*) "nph2(MXCOUNT4)"
         write(fileNumber,*) nph2
         write(fileNumber,*) "nph3(MXCOUNT4)"
         write(fileNumber,*) nph3
         write(fileNumber,*) "nph4(MXCOUNT4)"
         write(fileNumber,*) nph4
         write(fileNumber,*) "nph5(MXCOUNT4)"
         write(fileNumber,*) nph5
         write(fileNumber,*) "npl1(MXCOUNT4)"
         write(fileNumber,*) npl1
         write(fileNumber,*) "npl2(MXCOUNT4)"
         write(fileNumber,*) npl2
         write(fileNumber,*) "npl3(MXCOUNT4)"
         write(fileNumber,*) npl3
         write(fileNumber,*) "npl4(MXCOUNT4)"
         write(fileNumber,*) npl4
         write(fileNumber,*) "npl5(MXCOUNT4)"
         write(fileNumber,*) npl5
         write(fileNumber,*) "nolosp(ICP)"
         write(fileNumber,*) nolosp
         write(fileNumber,*) "newfold(NMTRATE*2, ICS)"
         write(fileNumber,*) newfold
         write(fileNumber,*) "nknlosp(MAXGL3, ICS)"
         write(fileNumber,*) nknlosp
         write(fileNumber,*) "nknphotrt(IPHOT,ICS)"
         write(fileNumber,*) nknphotrt
         write(fileNumber,*) "abst2   (ICS)"
         write(fileNumber,*) abst2
         write(fileNumber,*) "errmax  (ICS)"
         write(fileNumber,*) errmax
         write(fileNumber,*) "hmaxday (ICS)"
         write(fileNumber,*) hmaxday
         write(fileNumber,*) "timeintv(ICS)"
         write(fileNumber,*) timeintv
         write(fileNumber,*) "abtol(6, ICS)"
         write(fileNumber,*) abtol
         write(fileNumber,*) "enqq1 (MORDER)"
         write(fileNumber,*) enqq1
         write(fileNumber,*) "enqq2 (MORDER)"
         write(fileNumber,*) enqq2
         write(fileNumber,*) "enqq3 (MORDER)"
         write(fileNumber,*) enqq3
         write(fileNumber,*) "conp15(MORDER)"
         write(fileNumber,*) conp15
         write(fileNumber,*) "conpst(MORDER)"
         write(fileNumber,*) conpst
         write(fileNumber,*) "pertst2(MORDER, 3)"
         write(fileNumber,*) pertst2
         write(fileNumber,*) "aset(10, 8)"
         write(fileNumber,*) aset
         write(fileNumber,*) "fracpl (MXCOUNT2)"
         write(fileNumber,*) fracpl
         write(fileNumber,*) "fracnfr(MXCOUNT4)"
         write(fileNumber,*) fracnfr

         close(fileNumber)

      end subroutine writeSmv2Chem2Exit
