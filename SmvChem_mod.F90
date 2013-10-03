module SmvChem_mod

   implicit none
   private

#     include "smv2chem_par.h"

   public :: SmvChem_type
   public :: deallocateVariables
   public :: readPhysProc
   public :: writePhysProc

     type SmvChem_type
         integer  :: numQjo, numQks, numQjs, numActive
         ! num lats, longs, and vertical layers
         integer  :: numLat, numLong, numVert
         integer  :: numZones ! ilong * ilat * ivert
         real*8, allocatable :: qqjda(:, :, :, :)
         real*8, allocatable :: qqkda(:, :, :, :)
         real*8, allocatable :: qjGmi(:, :, :, :)
         real*8, allocatable :: qkGmi(:, :, :, :)
         real*8, allocatable :: yda  (:, :, :, :)
         logical, allocatable  :: doCellChem(:)
         real*8, allocatable :: thermalRateConstants(:, :) ! units vary (inverse conc squared, cubed, etc.? )
         real*8, allocatable :: photolysisRateConstants(:, :) ! s^-1
         real*8, allocatable :: surfaceEmissions(:, :) ! molec/cm^3/s
         real*8, allocatable :: speciesConst(:, :) ! molec/cm^3
         integer :: i1, i2, ju1, j2, k1, k2
         ! if pr_qqjk is on, should qqj's & qqk's be determined inside the chemistry solver, or outside?
         logical  :: doQqjkInchem
         ! do surface emissions inside the chemistry solver, or outside?
         logical  :: doSurfEmissInChem
         ! print some diagnostic output to screen?
         logical  :: prDiag
         ! should the periodic qqjk output file be written?
         logical  :: prQqjk
         ! should the SmvgearII output file be written (non-parallel mode only)?
         logical  :: prSmv2
         real*8  :: prNcPeriod ! NetCDF output period
         integer  :: localProc ! local proc #
         real*8  :: timeStep ! model time step (s)

      end type SmvChem_type

      contains

      subroutine readPhysProc (this, fileName)
         implicit none

#     include "smv2chem_par.h"

         ! Arguments
         type(SmvChem_type), intent(inout) :: this
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Reading from: ", fileName
         open(newunit=fileNumber,file=trim(fileName),form="formatted")

         read(fileNumber,'(a)')
         read(fileNumber,*) this%i1, this%i2, this%ju1, this%j2, this%k1, this%k2
         read(fileNumber,'(a)')
         read(fileNumber,*) this%numQjo, this%numQks, this%numQjs, this%numActive

         allocate(this%qjGmi(this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, this%numQjo))
         allocate(this%qkGmi(this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, this%numQks))
         allocate(this%qqjda(this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, this%numQjs))
         allocate(this%qqkda(this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, this%numQks))
         allocate(this%yda  (this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, this%numActive))

         read(fileNumber,'(a)')
         read(fileNumber,*) this%qjGmi(this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, 1:this%numQjo)
         read(fileNumber,'(a)')
         read(fileNumber,*) this%qkGmi(:,:,:,:)
         !!read(fileNumber,*) this%qkGmi(this%i1:i2, ju1:j2, k1:k2, 1:this%numQks) !causing a hang
         read(fileNumber,'(a)')
         read(fileNumber,*) this%qqjda(this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, 1:this%numQjs)
         read(fileNumber,'(a)')
         read(fileNumber,*) this%qqkda(:,:,:,:)
         read(fileNumber,'(a)')
         read(fileNumber,*) this%yda  (this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, 1:this%numActive)
         read(fileNumber,'(a)')
         read(fileNumber,*) this%doQqjkInchem
         read(fileNumber,'(a)')
         read(fileNumber,*) this%doSurfEmissInChem
         read(fileNumber,'(a)')
         read(fileNumber,*) this%prDiag
         read(fileNumber,'(a)')
         read(fileNumber,*) this%prQqjk
         read(fileNumber,'(a)')
         read(fileNumber,*) this%prSmv2
         read(fileNumber,'(a)')
         read(fileNumber,*) this%localProc
         read(fileNumber,'(a)')
         read(fileNumber,*) this%numLat, this%numLong, this%numVert
         read(fileNumber,'(a)')
         read(fileNumber,*) this%numZones
         read(fileNumber,'(a)')
         read(fileNumber,*) this%prNcPeriod
         read(fileNumber,'(a)')
         read(fileNumber,*) this%timeStep

         allocate(this%doCellChem(this%numZones))
         allocate(this%thermalRateConstants(this%numZones, ITHERM))
         allocate(this%photolysisRateConstants(this%numZones, IPHOT))
         allocate(this%surfaceEmissions(this%numLat*this%numLong, IGAS))
         allocate(this%speciesConst(this%numZones, IGAS))

         read(fileNumber,'(a)')
         read (fileNumber,*) this%doCellChem(1:this%numZones)
         read(fileNumber,'(a)')
         !read(fileNumber,*) this%thermalRateConstants(1:this%numZones, 1:ITHERM) ! causing a hang
         read(fileNumber,*) this%thermalRateConstants(:,:)
         read(fileNumber,'(a)')
         read(fileNumber,*) this%photolysisRateConstants(1:this%numZones, 1:IPHOT)
         read(fileNumber,'(a)')
         read(fileNumber,*) this%surfaceEmissions(1:this%numLat*this%numLong, 1:IGAS)
         read(fileNumber,'(a)')
         ! check that species const is being read in correctly

         !read(fileNumber,*) this%speciesConst(:,:)
         read(fileNumber,*) this%speciesConst(1:this%numZones, 1:IGAS)

         close(fileNumber)

      end subroutine readPhysProc

      subroutine writePhysProc (this, fileName)
         implicit none

#     include "smv2chem_par.h"

         ! Arguments
         type(SmvChem_type), intent(inout) :: this
         character(len=100), intent(in) :: fileName

         integer :: fileNumber

         print*, "Writing to: ", fileName
         open(newunit=fileNumber,file=trim(fileName),form="formatted")

         write(fileNumber,*) "i1, i2, ju1, j2, k1, k2"
         write(fileNumber,*) this%i1, this%i2, this%ju1, this%j2, this%k1, this%k2
         write(fileNumber,*) "num_qjo, num_qks, num_qjs, num_active"
         write(fileNumber,*) this%numQjo, this%numQks, this%numQjs, this%numActive
         write(fileNumber,*) "qjgmi(i1:i2, ju1:j2, k1:k2, num_qjo)"
         write(fileNumber,*) this%qjGmi(this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, 1:this%numQjo)
         write(fileNumber,*) "qkgmi(i1:i2, ju1:j2, k1:k2, num_qks)"
         write(fileNumber,*) this%qkGmi(:,:,:,:)
         write(fileNumber,*) "qqjda(i1:i2, ju1:j2, k1:k2, num_qjs)"
         write(fileNumber,*) this%qqjda(this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, 1:this%numQjs)
         write(fileNumber,*) "qqkda(i1:i2, ju1:j2, k1:k2, num_qks)"
         write(fileNumber,*) this%qqkda(:,:,:,:)
         write(fileNumber,*) "yda  (i1:i2, ju1:j2, k1:k2, num_active)"
         write(fileNumber,*) this%yda  (this%i1:this%i2, this%ju1:this%j2, this%k1:this%k2, 1:this%numActive)
         write(fileNumber,*) "do_qqjk_inchem"
         write(fileNumber,*) this%doQqjkInchem
         write(fileNumber,*) "do_semiss_inchem"
         write(fileNumber,*) this%doSurfEmissInChem
         write(fileNumber,*) "pr_diag"
         write(fileNumber,*) this%prDiag
         write(fileNumber,*) "pr_qqjk"
         write(fileNumber,*) this%prQqjk
         write(fileNumber,*) "pr_smv2"
         write(fileNumber,*) this%prSmv2
         write(fileNumber,*) "loc_proc"
         write(fileNumber,*) this%localProc
         write(fileNumber,*) "ilat, ilong, ivert"
         write(fileNumber,*) this%numLat, this%numLong, this%numVert
         write(fileNumber,*) "itloop"
         write(fileNumber,*) this%numZones
         write(fileNumber,*) "pr_nc_period"
         write(fileNumber,*) this%prNcPeriod
         write(fileNumber,*) "tdt"
         write(fileNumber,*) this%timeStep
         write(fileNumber,*) "do_cell_chem(itloop)"
         write(fileNumber,*) this%doCellChem(1:this%numZones)
         write(fileNumber,*) "arate(itloop, ITHERM)"
         !write(fileNumber,*) thermalRateConstants(1:this%numZones, 1:ITHERM)
         write(fileNumber,*) this%thermalRateConstants(:,:)
         write(fileNumber,*) "prate(itloop, IPHOT)"
         write(fileNumber,*) this%photolysisRateConstants(1:this%numZones, 1:IPHOT)
         write(fileNumber,*) "yemis(ilat*ilong, IGAS)"
         write(fileNumber,*) this%surfaceEmissions(1:this%numLat*this%numLong, 1:IGAS)
         write(fileNumber,*) "cx(itloop, IGAS)"
         write(fileNumber,*) this%speciesConst(1:this%numZones, 1:IGAS)

         close(fileNumber)

      end subroutine writePhysProc


      subroutine deallocateVariables (this)
         type (SmvChem_type), intent(inout) :: this

         deallocate(this%qjGmi)
         deallocate(this%qkGmi)
         deallocate(this%qqjda)
         deallocate(this%qqkda)
         deallocate(this%yda)
         deallocate(this%doCellChem)
         deallocate(this%thermalRateConstants)
         deallocate(this%photolysisRateConstants)
         deallocate(this%surfaceEmissions)
         deallocate(this%speciesConst)

      end subroutine deallocateVariables

end module SmvChem_mod
