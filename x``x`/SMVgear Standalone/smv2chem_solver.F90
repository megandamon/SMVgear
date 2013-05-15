!-----------------------------------------------------------------------------
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
      use Smv2Chem1_mod
      use Smv2Chem2_mod
		use ChemTable_mod

      implicit none
      external testFailOnPurpose
      external testSpeciesConst

#     include "smv2chem_par.h"

      logical, save :: first = .true.
      integer :: Tinit, Tfin, Tclockrate
      integer :: rank,errorInt, iCell
      integer, allocatable :: jReOrder(:)
      integer, allocatable :: lReOrder(:)
      real*8, allocatable  :: errorMx2  (:)
      real*8, allocatable, save :: cSumA(:)
      real*8, allocatable, save :: cSumB(:)
      real*8 :: MassInit,MassFin
      character(len=100) :: smv2Chem1Entry
      character(len=100) :: smv2Chem1Exit
      character(len=100) :: smv2Chem2Entry
      character(len=100) :: smv2Chem2Exit
      character(len=100) :: physProcEntry
      character(len=100) :: physProcExit
      character(len=100) :: summary_statement
      character(len=128) :: tempText
      type (SmvChem_type) :: chemObject
      type (Smv2Chem1_type) :: smv2chem1Object
      type (TestSuite_type) :: suite
      type (TestResult_type) :: result

      !call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
      !call timingInit
      call pFUnit_init()

      chemObject%prDiag = .true.

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

      call readSmv2Chem1Entry (smv2Chem1Entry, smv2chem1Object)
      call readSmv2Chem2Entry (smv2Chem2Entry)
      call readPhysProc (chemObject, physProcEntry)

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

      MassInit = sum(chemObject%speciesConst)
      call CalcTabl(GenChem,1,1,1)

      smv2chem1Object%intendedNumGridCellsInBlock = BLOCKSIZE
      smv2chem1Object%reorderGridCellsStiffness = DOREORD

      Call system_clock(Tinit)
      !call timingOn("Physproc")

      call physProc  &
     &  (chemObject%doQqjkInchem, chemObject%doSurfEmissInChem, chemObject%prQqjk, &
     &   chemObject%prSmv2, chemObject%numLat,  &
     &   chemObject%numLong, chemObject%numVert, smv2chem1Object%reorderGridCellsStiffness, &
     &   smv2chem1Object%imgas, smv2chem1Object%speciesNumOfNitrogen, &
     &   smv2chem1Object%speciesNumOfOxygen, chemObject%numZones,  &
     &   smv2chem1Object%intendedNumGridCellsInBlock, smv2chem1Object%unitNumberPrSmv2, &
     &   smv2chem1Object%gasChemistryType, smv2chem1Object%fractionDecreaseConvergence, &
     &   smv2chem1Object%maxTimeStepNight, chemObject%prNcPeriod, chemObject%timeStep,  &
     &   chemObject%doCellChem, smv2chem1Object%jphotrat, smv2chem1Object%numKineticRxns, &
     &   smv2chem1Object%ntloopncs, smv2chem1Object%numActiveAndInactiveGases, &
     &   smv2chem1Object%inewold,  &
     &   smv2chem1Object%npphotrat, chemObject%thermalRateConstants, &
     &   chemObject%photolysisRateConstants, chemObject%surfaceEmissions, &
     &   jReOrder, lReOrder, cSumA,  &
     &   cSumB, errorMx2, chemObject%speciesConst, &
     &   chemObject%yda, chemObject%qqkda, chemObject%qqjda, chemObject%qkGmi, chemObject%qjGmi, &
     &   chemObject%i1, chemObject%i2, chemObject%ju1, chemObject%j2, chemObject%k1, chemObject%k2, &
     &   chemObject%numQjo, chemObject%numQks, chemObject%numQjs, chemObject%numActive, chemObject%prDiag)

      !call timingOff("Physproc")
      Call system_clock(Tfin,Tclockrate)

      !call add(suite, TestCase1Step('testSpeciesConst'), testSpeciesConst, chemObject%speciesConst(1,1), 1077508322.73440)
      !speciesConst(1,1):   1077508322.73440
      !speciesConst(numZones,IGAS):   540086297554450.

      ! regression testing area
      call writeSmv2Chem1Exit (smv2Chem1Exit, smv2Chem1Object)
      call writeSmv2Chem2Exit (smv2Chem2Exit)
      call writePhysProc (chemObject, physProcExit)

		open(file=trim("cxfin"),unit=30,form="formatted")
		write(30,*) "cx(itloop,IGAS) = ", chemObject%numZones, " by ", IGAS
		do iCell=1,chemObject%numZones
			write(30,*) "Cell: ",iCell, "Max/Min = ",log10(maxval(chemObject%speciesConst(iCell,:))), &
			& "/",log10(minval(chemObject%speciesConst(iCell,:)))
			write(30,'(es22.15)') chemObject%speciesConst(iCell,:)
		end do
		close(30)

		MassFin = sum(chemObject%speciesConst)

		Write(*,*) 'Relative mass change = ', (MassFin-MassInit)/MassInit
      Write(*,*) 'Reordering = ', smv2chem1Object%reorderGridCellsStiffness
      Write(*,*) 'Blocksize = ', smv2chem1Object%intendedNumGridCellsInBlock
      Write(*,*) '(" Time taken = ", f9.6, " seconds.")', Real(Tfin-Tinit)/Real(Tclockrate)
		write(*,*) 'Species 61:'
		write(*,*) '82 - > ', smv2chem1Object%inewold(82,1)
		write(*,*) 'Inputs = ', count(GenChem%RxnIn==82)
		write(*,*) 'Outputs = ', count(GenChem%RxnOut==82)
      write(*,*) 'Frac-Out = ', count(GenChem%FracRxnOut==82)


      ! Run the tests and accumulate the results in "result"
      result = newTestResult(mode=MODE_USE_STDOUT)
      call Run(suite, result)
      summary_statement=Summary(result)
      print*,trim(summary_statement)

		deallocate (jReOrder)
      deallocate (lReOrder)
      deallocate (errorMx2)
      deallocate (smv2chem1Object%jphotrat)
      deallocate (smv2chem1Object%numKineticRxns)
      deallocate (smv2chem1Object%ntloopncs)
      deallocate (smv2chem1Object%numActiveAndInactiveGases)
      deallocate (smv2chem1Object%inewold)
      deallocate (smv2chem1Object%npphotrat)

      call deallocateVariables(chemObject)

	   call clean(result)
      call clean(suite)
      call pFUnit_finalize()

      print*, "Exiting doSmv2Solver"

   end program doSmv2Solver






      subroutine readSmv2Chem1Entry (fileName, smv2chem1Object)
         use Smv2Chem1_mod, only : Smv2Chem1_type
         implicit none

#     include "smv2chem_par.h"

         ! Arguments
         character(len=100), intent(in) :: fileName
         type(Smv2Chem1_type), intent(inout) :: smv2chem1Object

         integer :: fileNumber
         integer :: testVar
         character(len=128) :: testString

         print*, "Reading from: ", fileName
         open(newunit=fileNumber, file=trim(fileName),form="formatted")

         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%reorderGridCellsStiffness
         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%speciesNumOfWaterVapor
         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%imgas
         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%speciesNumOfNitrogen

         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%speciesNumOfOxygen
         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%intendedNumGridCellsInBlock
         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%unitNumberPrSmv2
         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%gasChemistryType
         read(fileNumber,'(a)')
         allocate (smv2chem1Object%jphotrat(ICS))
         read(fileNumber,*) smv2chem1Object%jphotrat
         read(fileNumber,'(a)')
         allocate (smv2chem1Object%numKineticRxns(ICS))
         read(fileNumber,*) smv2chem1Object%numKineticRxns
         read(fileNumber,'(a)')
         allocate (smv2chem1Object%ntloopncs(ICS))
         read(fileNumber,*) smv2chem1Object%ntloopncs
         read(fileNumber,'(a)')
         allocate (smv2chem1Object%numActiveAndInactiveGases(ICS))
         read(fileNumber,*) smv2chem1Object%numActiveAndInactiveGases
         read(fileNumber,'(a)')
         allocate (smv2chem1Object%inewold(MXGSAER,ICS))
         read(fileNumber,*) smv2chem1Object%inewold
         read(fileNumber,'(a)')
         allocate (smv2chem1Object%npphotrat(IPHOT,ICS))
         read(fileNumber,*) smv2chem1Object%npphotrat
         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%fractionDecreaseConvergence
         read(fileNumber,'(a)')
         read(fileNumber,*) smv2chem1Object%maxTimeStepNight

         close(fileNumber)

      end subroutine readSmv2Chem1Entry


      subroutine readSmv2Chem2Entry (fileName)
         use Smv2Chem2_mod
         implicit none

#     include "smv2chem_par.h"

         ! Arguments
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Reading from: ", fileName
         open(newunit=fileNumber, file=trim(fileName),form="formatted")

         allocate(numRxnsOneActiveReactant(ICP))
         allocate(nallrat(ICP))
         allocate(lastReorderedRxn(ICS))
         allocate(numRxnsThreeActiveReactant (ICS))
         allocate(numRxnsTwoActiveReactant (ICS))
         allocate(nm3bod(ICS))
         allocate(numRxns3rdSpcO2plusN2 (ICS))
         allocate(numRxns3rdSpcN2(ICS))
         allocate(numRxns3rdSpcO2(ICS))
         allocate(nmoth(ICS))
         allocate(numKinPhotoRxns(ICS))
         allocate(origSpcToReordSpcMap(MXGSAER, ICS))
         allocate(lgasbino(MAXGL2,ICS))
         allocate(nreacoth(MAXGL2,ICS))
         allocate(lgas3bod(MAXGL3,ICS))
         allocate(losinacp(MAXGL3,ICS))
         allocate(nreac3b (MAXGL3,ICS))
         allocate(nreacair(MAXGL3,ICS))
         allocate(nreacn2 (MAXGL3,ICS))
         allocate(nreaco2 (MAXGL3,ICS))
         allocate(jphotnk(NMTRATE, ICS))
         allocate(oldRxnRateNum(NMTRATE, ICS))
         allocate(newSpcNumActiveProduct(NMRPROD, NMTRATE, ICS))
         allocate(numOrigSpcGtrOrEql1PdTerm(ICS))
         allocate(kzthi(ICP))
         allocate(kztlo(ICP))

         allocate(ikztot(MXCOUNT4))

         allocate(kbh1(MXCOUNT4))
         allocate(kbh2(MXCOUNT4))
         allocate(kbh3(MXCOUNT4))
         allocate( kbh4(MXCOUNT4))
         allocate( kbh5(MXCOUNT4))
         allocate(kbl1(MXCOUNT4))
         allocate( kbl2(MXCOUNT4))
         allocate(kbl3(MXCOUNT4))
         allocate( kbl4(MXCOUNT4))
         allocate( kbl5(MXCOUNT4))

         allocate(mbh1(MXCOUNT4))
         allocate( mbh2(MXCOUNT4))
         allocate(mbh3(MXCOUNT4))
         allocate( mbh4(MXCOUNT4))
         allocate( mbh5(MXCOUNT4))
         allocate(mbl1(MXCOUNT4))
         allocate( mbl2(MXCOUNT4))
         allocate(mbl3(MXCOUNT4))
         allocate( mbl4(MXCOUNT4))
         allocate( mbl5(MXCOUNT4))

         allocate(kzeroa(MXCOUNT4))
         allocate(kzerob(MXCOUNT4))
         allocate( kzeroc(MXCOUNT4))
         allocate(kzerod(MXCOUNT4))
         allocate( kzeroe(MXCOUNT4))

         allocate(mzeroa(MXCOUNT4))
         allocate(mzerob(MXCOUNT4))
         allocate( mzeroc(MXCOUNT4))
         allocate(mzerod(MXCOUNT4))
         allocate( mzeroe(MXCOUNT4))
         allocate(imztot(MXGSAER, ICP))

         allocate(ijval (MXCOUNT3))
         allocate(jzeroa(MXCOUNT3))

         allocate(idh1  (MXCOUNT3))
         allocate(idh2  (MXCOUNT3))
         allocate(idh3  (MXCOUNT3))
         allocate( idh4  (MXCOUNT3))
         allocate( idh5(MXCOUNT3))
         allocate(idl1  (MXCOUNT3))
         allocate( idl2  (MXCOUNT3))
         allocate(idl3  (MXCOUNT3))
         allocate( idl4  (MXCOUNT3))
         allocate( idl5(MXCOUNT3))

         allocate(ikdeca(MXCOUNT3))
         allocate( ikdecb(MXCOUNT3))
         allocate(ikdecc(MXCOUNT3))
         allocate( ikdecd(MXCOUNT3))
         allocate(ikdece(MXCOUNT3))

         allocate(kjdeca(MXCOUNT3))
         allocate(kjdecb(MXCOUNT3))
         allocate( kjdecc(MXCOUNT3))
         allocate(kjdecd(MXCOUNT3))
         allocate( kjdece(MXCOUNT3))

         allocate(ijthi(MXGSAER, ICP))
         allocate(ijtlo(MXGSAER, ICP))
         allocate(diagonalTermDecomp(MXGSAER, ICP))
         allocate(jhiz1(MXGSAER, ICP))
         allocate(jloz1(MXGSAER, ICP))

         allocate(sparseMatrixDimension(ICP))
         allocate(npdhi(ICP))
         allocate(npdlo(ICP))
         allocate(nfrhi(ICP))
         allocate(nfrlo(ICP))
         allocate(nplhi(ICP))
         allocate(npllo(ICP))


         allocate(iialpd(MXCOUNT2))
         allocate(ipospd(MXCOUNT2))
         allocate(nkpdterm(MXCOUNT2))

         allocate(jspcnfr(MXCOUNT4))
         allocate(jspnpl(MXCOUNT4))
         allocate(nknfr(MXCOUNT4))

         allocate(reOrderedRxnRateEachLossA(MXCOUNT4))
         allocate(reOrderedRxnRateEachLossB(MXCOUNT4))
         allocate(reOrderedRxnRateEachLossC(MXCOUNT4))
         allocate(reOrderedRxnRateEachLossD(MXCOUNT4))
         allocate(reOrderedRxnRateEachLossE(MXCOUNT4))

         allocate(nph1(MXCOUNT4))
         allocate(nph2(MXCOUNT4))
         allocate(nph3(MXCOUNT4))
         allocate(nph4(MXCOUNT4))
         allocate(nph5(MXCOUNT4))
         allocate(npl1(MXCOUNT4))
         allocate(npl2(MXCOUNT4))
         allocate(npl3(MXCOUNT4))
         allocate(npl4(MXCOUNT4))
         allocate(npl5(MXCOUNT4))

         allocate(nolosp(ICP))

         allocate(newRxnRateNumber(NMTRATE*2, ICS))

         allocate(nknlosp(MAXGL3, ICS))

         allocate(nknphotrt(IPHOT,ICS))


         allocate(abst2(ICS))
         allocate(relativeErrorTolerance(ICS))
         allocate(hmaxday(ICS))
         allocate(timeintv(ICS))

         allocate(absoluteErrorTolerance(6, ICS))
         allocate(enqq1 (MORDER))
         allocate(enqq2 (MORDER))
         allocate( enqq3 (MORDER))
         allocate(conp15(MORDER))
         allocate( conpst(MORDER))

         allocate(coeffsForSelectingStepAndOrder(MORDER, 3))

         allocate(coeffsForIntegrationOrder(10, 8))
         allocate(fracpl (MXCOUNT2))
         allocate(fracnfr(MXCOUNT4))


         read(fileNumber,'(a)')
         read(fileNumber,*) numRxnsOneActiveReactant
         read(fileNumber,'(a)')
         read(fileNumber,*) nallrat
         read(fileNumber,'(a)')
         read(fileNumber,*) lastReorderedRxn
         read(fileNumber,'(a)')
         read(fileNumber,*) numRxnsThreeActiveReactant
         read(fileNumber,'(a)')
         read(fileNumber,*) numRxnsTwoActiveReactant
         read(fileNumber,'(a)')
         read(fileNumber,*) nm3bod
         read(fileNumber,'(a)')
         read(fileNumber,*) numRxns3rdSpcO2plusN2
         read(fileNumber,'(a)')
         read(fileNumber,*) numRxns3rdSpcN2
         read(fileNumber,'(a)')
         read(fileNumber,*) numRxns3rdSpcO2
         read(fileNumber,'(a)')
         read(fileNumber,*) nmoth
         read(fileNumber,'(a)')
         read(fileNumber,*) numKinPhotoRxns
         read(fileNumber,'(a)')
         read(fileNumber,*) origSpcToReordSpcMap
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
         read(fileNumber,*) oldRxnRateNum
         read(fileNumber,'(a)')
         read(fileNumber,*) newSpcNumActiveProduct
         read(fileNumber,'(a)')
         read(fileNumber,*) numOrigSpcGtrOrEql1PdTerm
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
         read(fileNumber,*) diagonalTermDecomp
         read(fileNumber,'(a)')
         read(fileNumber,*) jhiz1
         read(fileNumber,'(a)')
         read(fileNumber,*) jloz1
         read(fileNumber,'(a)')
         read(fileNumber,*) sparseMatrixDimension
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
         read(fileNumber,*) reOrderedRxnRateEachLossA
         read(fileNumber,'(a)')
         read(fileNumber,*) reOrderedRxnRateEachLossB
         read(fileNumber,'(a)')
         read(fileNumber,*) reOrderedRxnRateEachLossC
         read(fileNumber,'(a)')
         read(fileNumber,*) reOrderedRxnRateEachLossD
         read(fileNumber,'(a)')
         read(fileNumber,*) reOrderedRxnRateEachLossE
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
         read(fileNumber,*) newRxnRateNumber
         read(fileNumber,'(a)')
         read(fileNumber,*) nknlosp
         read(fileNumber,'(a)')
         read(fileNumber,*) nknphotrt
         read(fileNumber,'(a)')
         read(fileNumber,*) abst2
         read(fileNumber,'(a)')
         read(fileNumber,*) relativeErrorTolerance
         read(fileNumber,'(a)')
         read(fileNumber,*) hmaxday
         read(fileNumber,'(a)')
         read(fileNumber,*) timeintv
         read(fileNumber,'(a)')
         read(fileNumber,*) absoluteErrorTolerance
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
         read(fileNumber,*) coeffsForSelectingStepAndOrder
         read(fileNumber,'(a)')
         read(fileNumber,*) coeffsForIntegrationOrder
         read(fileNumber,'(a)')
         read(fileNumber,*) fracpl
         read(fileNumber,'(a)')
         read(fileNumber,*) fracnfr

         close(fileNumber)

      end subroutine readSmv2Chem2Entry




      subroutine writeSmv2Chem1Exit (fileName, smv2chem1Object)
         use Smv2Chem1_mod
         implicit none

#     include "smv2chem_par.h"

         ! Arguments
         character(len=100), intent(in) :: fileName
         type(Smv2Chem1_type) :: smv2chem1Object

         integer :: fileNumber

         print*, "Writing to: ", fileName
         open(newunit=fileNumber, file=trim(fileName),form="formatted")

         write(fileNumber,*) "ifreord   "
         write(fileNumber,*) smv2chem1Object%reorderGridCellsStiffness
         write(fileNumber,*) "ih2o      "
         write(fileNumber,*) smv2chem1Object%speciesNumOfWaterVapor
         write(fileNumber,*) "imgas     "
         write(fileNumber,*) smv2chem1Object%imgas
         write(fileNumber,*) "initrogen "
         write(fileNumber,*) smv2chem1Object%speciesNumOfNitrogen
         write(fileNumber,*) "ioxygen   "
         write(fileNumber,*) smv2chem1Object%speciesNumOfOxygen
         write(fileNumber,*) "kuloop    "
         write(fileNumber,*) smv2chem1Object%intendedNumGridCellsInBlock
         write(fileNumber,*) "lunsmv    "
         write(fileNumber,*) smv2chem1Object%unitNumberPrSmv2
         write(fileNumber,*) "ncs       "
         write(fileNumber,*) smv2chem1Object%gasChemistryType
         write(fileNumber,*) "jphotrat (ICS)            "
         write(fileNumber,*) smv2chem1Object%jphotrat
         write(fileNumber,*) "nrates   (ICS)            "
         write(fileNumber,*) smv2chem1Object%numKineticRxns
         write(fileNumber,*) "ntloopncs(ICS)            "
         write(fileNumber,*) smv2chem1Object%ntloopncs
         write(fileNumber,*) "ntspec   (ICS)            "
         write(fileNumber,*) smv2chem1Object%numActiveAndInactiveGases
         write(fileNumber,*) "inewold  (MXGSAER, ICS)   "
         write(fileNumber,*) smv2chem1Object%inewold
         write(fileNumber,*) "npphotrat(IPHOT,   ICS)   "
         write(fileNumber,*) smv2chem1Object%npphotrat
         write(fileNumber,*) "fracdec   "
         write(fileNumber,*) smv2chem1Object%fractionDecreaseConvergence
         write(fileNumber,*) "hmaxnit   "
         write(fileNumber,*) smv2chem1Object%maxTimeStepNight

         close(fileNumber)

      end subroutine writeSmv2Chem1Exit


      subroutine writeSmv2Chem2Exit (fileName)
         use Smv2Chem2_mod
         implicit none

#     include "smv2chem_par.h"

         ! Arguments
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Writing to: ", fileName
         open(newunit=fileNumber, file=trim(fileName),form="formatted")

         write(fileNumber,*) "ioner  (ICP)"
         write(fileNumber,*) numRxnsOneActiveReactant
         write(fileNumber,*) "nallrat(ICP)"
         write(fileNumber,*) nallrat
         write(fileNumber,*) "inorep (ICS)"
         write(fileNumber,*) lastReorderedRxn
         write(fileNumber,*) "ithrr  (ICS)"
         write(fileNumber,*) numRxnsThreeActiveReactant
         write(fileNumber,*) "itwor  (ICS)"
         write(fileNumber,*) numRxnsTwoActiveReactant
         write(fileNumber,*) "nm3bod (ICS)"
         write(fileNumber,*) nm3bod
         write(fileNumber,*) "nmair  (ICS)"
         write(fileNumber,*) numRxns3rdSpcO2plusN2
         write(fileNumber,*) "nmn2   (ICS)"
         write(fileNumber,*) numRxns3rdSpcN2
         write(fileNumber,*) "nmo2 (ICS)"
         write(fileNumber,*) numRxns3rdSpcO2
         write(fileNumber,*) "nmoth  (ICS)"
         write(fileNumber,*) nmoth
         write(fileNumber,*) "ntrates(ICS)"
         write(fileNumber,*) numKinPhotoRxns
         write(fileNumber,*) "mappl   (MXGSAER, ICS)"
         write(fileNumber,*) origSpcToReordSpcMap
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
         write(fileNumber,*) oldRxnRateNum
         write(fileNumber,*) "irm2(NMRPROD, NMTRATE, ICS)"
         write(fileNumber,*) newSpcNumActiveProduct
         write(fileNumber,*) "ischang(ICS)"
         write(fileNumber,*) numOrigSpcGtrOrEql1PdTerm
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
         write(fileNumber,*) diagonalTermDecomp
         write(fileNumber,*) "jhiz1   (MXGSAER, ICP)"
         write(fileNumber,*) jhiz1
         write(fileNumber,*) "jloz1(MXGSAER, ICP)"
         write(fileNumber,*) jloz1
         write(fileNumber,*) "iarray(ICP)"
         write(fileNumber,*) sparseMatrixDimension
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
         write(fileNumber,*) reOrderedRxnRateEachLossA
         write(fileNumber,*) "lossrb (MXCOUNT4)"
         write(fileNumber,*) reOrderedRxnRateEachLossB
         write(fileNumber,*) "lossrc(MXCOUNT4)"
         write(fileNumber,*) reOrderedRxnRateEachLossC
         write(fileNumber,*) "lossrd (MXCOUNT4)"
         write(fileNumber,*) reOrderedRxnRateEachLossD
         write(fileNumber,*) "lossre(MXCOUNT4)"
         write(fileNumber,*) reOrderedRxnRateEachLossE
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
         write(fileNumber,*) newRxnRateNumber
         write(fileNumber,*) "nknlosp(MAXGL3, ICS)"
         write(fileNumber,*) nknlosp
         write(fileNumber,*) "nknphotrt(IPHOT,ICS)"
         write(fileNumber,*) nknphotrt
         write(fileNumber,*) "abst2   (ICS)"
         write(fileNumber,*) abst2
         write(fileNumber,*) "errmax  (ICS)"
         write(fileNumber,*) relativeErrorTolerance
         write(fileNumber,*) "hmaxday (ICS)"
         write(fileNumber,*) hmaxday
         write(fileNumber,*) "timeintv(ICS)"
         write(fileNumber,*) timeintv
         write(fileNumber,*) "abtol(6, ICS)"
         write(fileNumber,*) absoluteErrorTolerance
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
         write(fileNumber,*) coeffsForSelectingStepAndOrder
         write(fileNumber,*) "aset(10, 8)"
         write(fileNumber,*) coeffsForIntegrationOrder
         write(fileNumber,*) "fracpl (MXCOUNT2)"
         write(fileNumber,*) fracpl
         write(fileNumber,*) "fracnfr(MXCOUNT4)"
         write(fileNumber,*) fracnfr

         close(fileNumber)

         deallocate(numRxnsOneActiveReactant)
         deallocate(nallrat)
         deallocate(lastReorderedRxn)
         deallocate(numRxnsThreeActiveReactant)
         deallocate(numRxnsTwoActiveReactant)
         deallocate(nm3bod)
         deallocate(numRxns3rdSpcO2plusN2)
         deallocate(numRxns3rdSpcN2)
         deallocate(numRxns3rdSpcO2)
         deallocate(nmoth)
         deallocate(numKinPhotoRxns)
         deallocate(origSpcToReordSpcMap)
         deallocate(lgasbino)
         deallocate(nreacoth)
         deallocate(lgas3bod)
         deallocate(losinacp)
         deallocate(nreac3b )
         deallocate(nreacair)
         deallocate(nreacn2 )
         deallocate(nreaco2 )
         deallocate(jphotnk)
         deallocate(oldRxnRateNum)
         deallocate(newSpcNumActiveProduct)
         deallocate(numOrigSpcGtrOrEql1PdTerm)
         deallocate(kzthi)
         deallocate(kztlo)

         deallocate(ikztot)

         deallocate(kbh1)
         deallocate(kbh2)
         deallocate(kbh3)
         deallocate( kbh4)
         deallocate( kbh5)
         deallocate(kbl1)
         deallocate( kbl2)
         deallocate(kbl3)
         deallocate( kbl4)
         deallocate( kbl5)

         deallocate(mbh1)
         deallocate( mbh2)
         deallocate(mbh3)
         deallocate( mbh4)
         deallocate( mbh5)
         deallocate(mbl1)
         deallocate( mbl2)
         deallocate(mbl3)
         deallocate( mbl4)
         deallocate( mbl5)

         deallocate(kzeroa)
         deallocate(kzerob)
         deallocate( kzeroc)
         deallocate(kzerod)
         deallocate( kzeroe)

         deallocate(mzeroa)
         deallocate(mzerob)
         deallocate( mzeroc)
         deallocate(mzerod)
         deallocate( mzeroe)
         deallocate(imztot)

         deallocate(ijval)
         deallocate(jzeroa)

         deallocate(idh1)
         deallocate(idh2)
         deallocate(idh3)
         deallocate( idh4)
         deallocate( idh5)
         deallocate(idl1)
         deallocate( idl2)
         deallocate(idl3)
         deallocate( idl4)
         deallocate( idl5)

         deallocate(ikdeca)
         deallocate( ikdecb)
         deallocate(ikdecc)
         deallocate( ikdecd)
         deallocate(ikdece)

         deallocate(kjdeca)
         deallocate(kjdecb)
         deallocate( kjdecc)
         deallocate(kjdecd)
         deallocate(kjdece)

         deallocate(ijthi)
         deallocate(ijtlo)
         deallocate(diagonalTermDecomp)
         deallocate(jhiz1)
         deallocate(jloz1)

         deallocate(sparseMatrixDimension)
         deallocate(npdhi)
         deallocate(npdlo)
         deallocate(nfrhi)
         deallocate(nfrlo)
         deallocate(nplhi)
         deallocate(npllo)

         deallocate(iialpd)
         deallocate(ipospd)
         deallocate(nkpdterm)

         deallocate(jspcnfr)
         deallocate(jspnpl)
         deallocate(nknfr)

         deallocate(reOrderedRxnRateEachLossA)
         deallocate(reOrderedRxnRateEachLossB)
         deallocate(reOrderedRxnRateEachLossC)
         deallocate(reOrderedRxnRateEachLossD)
         deallocate(reOrderedRxnRateEachLossE)

         deallocate(nph1)
         deallocate(nph2)
         deallocate(nph3)
         deallocate(nph4)
         deallocate(nph5)
         deallocate(npl1)
         deallocate(npl2)
         deallocate(npl3)
         deallocate(npl4)
         deallocate(npl5)
         deallocate(nolosp)

         deallocate(newRxnRateNumber)

         deallocate(nknlosp)

         deallocate(nknphotrt)


         deallocate(abst2)
         deallocate(relativeErrorTolerance)
         deallocate(hmaxday)
         deallocate(timeintv)

         deallocate(absoluteErrorTolerance)
         deallocate(enqq1)
         deallocate(enqq2)
         deallocate(enqq3)
         deallocate(conp15)
         deallocate(conpst)

         deallocate(coeffsForSelectingStepAndOrder)

         deallocate(coeffsForIntegrationOrder)
         deallocate(fracpl)
         deallocate(fracnfr)

      end subroutine writeSmv2Chem2Exit
