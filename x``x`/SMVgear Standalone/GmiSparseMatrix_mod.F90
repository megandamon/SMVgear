module GmiSparseMatrix_mod

   implicit none
   private

#     include "smv2chem_par.h"

   public :: calculatePredictor
   public :: LU_Decomp
   public :: BackSub
   public :: SparseMatrix


   type SparseMatrix
       integer :: holder
    end type SparseMatrix

contains

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
!   numGridCellsInBlock : # of grid-cells in a grid-block
!   diurnalGasChemType   : gasChemistryType       => for daytime   gas chemistry
!            gasChemistryType + ICS => for nighttime gas chemistry
!   valuesDecomposedMatrix    : array holding values of decomposed matrix
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!   gloss  : first derivative = sum of prod minus loss rates for a spc
!
!-----------------------------------------------------------------------------

      subroutine Backsub  &
     &  (num1stOEqnsSolve, numGridCellsInBlock, diurnalGasChemType, valuesDecomposedMatrix, vdiag, gloss)

      use Smv2Chem2_mod
      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: num1stOEqnsSolve
      integer, intent(in)  :: numGridCellsInBlock
      integer, intent(in)  :: diurnalGasChemType
      real*8,  intent(in)  :: valuesDecomposedMatrix  (KBLOOP, 0:MXARRAY)
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
      KZTLOOP: do kzt = kztlo(diurnalGasChemType), kzthi(diurnalGasChemType)
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
          do k = 1, numGridCellsInBlock
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (valuesDecomposedMatrix(k,ij0) * gloss(k,j0)) -  &
     &        (valuesDecomposedMatrix(k,ij1) * gloss(k,j1)) -  &
     &        (valuesDecomposedMatrix(k,ij2) * gloss(k,j2)) -  &
     &        (valuesDecomposedMatrix(k,ij3) * gloss(k,j3)) -  &
     &        (valuesDecomposedMatrix(k,ij4) * gloss(k,j4))
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

          do k = 1, numGridCellsInBlock
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (valuesDecomposedMatrix(k,ij0) * gloss(k,j0)) -  &
     &        (valuesDecomposedMatrix(k,ij1) * gloss(k,j1)) -  &
     &        (valuesDecomposedMatrix(k,ij2) * gloss(k,j2)) -  &
     &        (valuesDecomposedMatrix(k,ij3) * gloss(k,j3))
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

          do k = 1, numGridCellsInBlock
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (valuesDecomposedMatrix(k,ij0) * gloss(k,j0)) -  &
     &        (valuesDecomposedMatrix(k,ij1) * gloss(k,j1)) -  &
     &        (valuesDecomposedMatrix(k,ij2) * gloss(k,j2))
          end do

        end do

!       -- Sum 2 terms at a time. --

        do kc = kl2, kh2

          ij0 = ij
          ij1 = ij + 1
          ij  = ij + 2

          j0  = kzeroa(kc)
          j1  = kzerob(kc)

          do k = 1, numGridCellsInBlock
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (valuesDecomposedMatrix(k,ij0) * gloss(k,j0)) -  &
     &        (valuesDecomposedMatrix(k,ij1) * gloss(k,j1))
          end do

        end do

!       -- Sum 1 term at a time. --

        do kc = kl1, kh1

          ij0 = ij
          ij  = ij + 1

          j0  = kzeroa(kc)

          do k = 1, numGridCellsInBlock
            gloss(k,i) =  &
     &        gloss(k,i) -  &
     &        (valuesDecomposedMatrix(k,ij0) * gloss(k,j0))
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

        mzt = imztot(i,diurnalGasChemType)

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

            do k = 1, numGridCellsInBlock
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (valuesDecomposedMatrix(k,ij0) * gloss(k,j0)) -  &
     &          (valuesDecomposedMatrix(k,ij1) * gloss(k,j1)) -  &
     &          (valuesDecomposedMatrix(k,ij2) * gloss(k,j2)) -  &
     &          (valuesDecomposedMatrix(k,ij3) * gloss(k,j3)) -  &
     &          (valuesDecomposedMatrix(k,ij4) * gloss(k,j4))
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

            do k = 1, numGridCellsInBlock
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (valuesDecomposedMatrix(k,ij0) * gloss(k,j0)) -  &
     &          (valuesDecomposedMatrix(k,ij1) * gloss(k,j1)) -  &
     &          (valuesDecomposedMatrix(k,ij2) * gloss(k,j2)) -  &
     &          (valuesDecomposedMatrix(k,ij3) * gloss(k,j3))
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

            do k = 1, numGridCellsInBlock
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (valuesDecomposedMatrix(k,ij0) * gloss(k,j0)) -  &
     &          (valuesDecomposedMatrix(k,ij1) * gloss(k,j1)) -  &
     &          (valuesDecomposedMatrix(k,ij2) * gloss(k,j2))
            end do

          end do

!         -- Sum 2 terms at a time. --

          do mc = ml2, mh2

            ij0 = ij
            ij1 = ij + 1
            ij  = ij + 2

            j0  = mzeroa(mc)
            j1  = mzerob(mc)

            do k = 1, numGridCellsInBlock
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (valuesDecomposedMatrix(k,ij0) * gloss(k,j0)) -  &
     &          (valuesDecomposedMatrix(k,ij1) * gloss(k,j1))
            end do

          end do

!         -- Sum 1 term at a time. --

          do mc = ml1, mh1

            ij0 = ij
            ij  = ij + 1

            j0  = mzeroa(mc)

            do k = 1, numGridCellsInBlock
              gloss(k,i) =  &
     &          gloss(k,i) -  &
     &          (valuesDecomposedMatrix(k,ij0) * gloss(k,j0))
            end do

          end do

!       ============
        end if MZTIF
!       ============

!       -- Adjust gloss with diagonal element. --

        do k = 1, numGridCellsInBlock
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
!   calculatePredictor
!
! DESCRIPTION
!     Calculates the predictor matrix: (P) = I - h * b * J:
!       J = Jacobian matrix of partial derivates
!       I = identity matrix
!       h = time step
!       b = coefficient of method
!       R = h * b = -r1delt
!
! ARGUMENTS
!   nondiag    : # of final matrix positions, excluding diagonal
                             ! terms, filled after all matrix processes
!   nondiag1   : nondiag + 1
!   iarry      : iarray(ncsp)
!   numGridCellsInBlock : # of grid-cells in a grid-block
!   npdh  :
!   npdl  :
!   r1delt : = -aset(nqq,1) * time_step = -coefficient_of_method * dt
!   jacobian : term of Jacobian (J) = partial derivative
!   predictor : array of iarray units holding values of each matrix
!              position actually used;
!              cc2 = P = I - delt * aset(nqq,1) * partial_derivatives
!
! NOTES: MRD: solver should not know number of reactions
!-----------------------------------------------------------------------------

      subroutine calculatePredictor (nondiag, iarry, ktloop, cx, &
      & npdh, npdl, r1delt, Predict, rateConstants)

      use Smv2Chem2_mod
		use ChemTable_mod
		use GmiMechanism_mod

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------
         integer, intent(in) :: nondiag     ! # of final matrix positions, excluding diagonal
                             ! terms, filled after all matrix processes
         integer, intent(in) :: iarry
         integer, intent(in) :: ktloop
			real*8, intent(in)  :: cx(KBLOOP,MXGSAER)
         integer, intent(in) :: npdh, npdl
         real*8,  intent(in)  :: r1delt
         real*8,  intent(inout) :: Predict (KBLOOP, 0:MXARRAY)
         real*8  :: rateConstants (KBLOOP, NMTRATE)


!     ----------------------
!     Variable declarations.
!     ----------------------
         integer :: iar, k, n, nkn, ial, nin, ina, inb
         real*8  :: fracr1


         ! list of non-zero values
         ! MRD: derived type that combines the predictor below,
         ! with this information below (next two loops)
         ! could be called sparseMatrix (stay on the solver type)
         ! is predictor and Jacobian the same size?
         do iar = 1, nondiag
            do k = 1, ktloop 
               Predict(k,iar) = 0.0d0
            end do
         end do

         do iar = nondiag + 1, iarry
            do k = 1, ktloop 
               Predict(k,iar) = 1.0d0
            end do
         end do

         do n = npdl, npdh

            nkn    = nkpdterm(n) ! nkpdterm in common block
            iar    = ipospd  (n) ! ipospd in common block
            ial    = iialpd  (n) ! iialpd in common block
            fracr1 = fracpl  (n) * r1delt ! fracpl in a common block

				nin = GenChem%RxnNumIn(nkn)
				select case (nin)
					case(2) 
						if (ial .eq. 1) then
							ina = GenChem%RxnIn(2,nkn)
						else 
							ina = GenChem%RxnIn(1,nkn)
						end if
						do k=1, ktloop
							Predict(k,iar) = Predict(k,iar) + fracr1* rateConstants(k,nkn) * cx(k,ina)
						end do
					case(1) !Only one var to take deriv wrt
						do k=1, ktloop
							Predict(k,iar) = Predict(k,iar)+ (fracr1*rateConstants(k,nkn))
						end do
					case(3)
						if (ial .eq. 1) then
							ina = GenChem%RxnIn(2,nkn)
							inb = GenChem%RxnIn(3,nkn)
						else if (ial .eq. 2) then
							ina = GenChem%RxnIn(1,nkn)
							inb = GenChem%RxnIn(3,nkn)
						else 
							ina = GenChem%RxnIn(1,nkn)
							inb = GenChem%RxnIn(2,nkn)
						endif
						do k=1, ktloop
							Predict(k,iar) = Predict(k,iar) + fracr1* &
							& rateConstants(k,nkn)*cx(k,ina)*cx(k,inb)
						end do
				end select

         end do

      end subroutine calculatePredictor

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
!   numGridCellsInBlock : # of grid-cells in a grid-block
!   diurnalGasChemType   : gasChemistryType       => for daytime   gas chemistry
!            gasChemistryType + ICS => for nighttime gas chemistry
!   valuesDecomposedMatrix    : array of sparseMatrixDimension units holding values of each matrix
!            position actually used; originally,
!            valuesDecomposedMatrix = P = I - currentTimeStep * coeffsForIntegrationOrder(orderOfIntegrationMethod,1) * partial_derivatives;
!            however, valuesDecomposedMatrix is decomposed here
!   vdiag  : 1 / current diagonal term of the decomposed matrix
!
!-----------------------------------------------------------------------------
      subroutine LU_Decomp  &
     &  (num1stOEqnsSolve, numGridCellsInBlock, diurnalGasChemType, valuesDecomposedMatrix, vdiag)

      use Smv2Chem2_mod
      implicit none

#     include "smv2chem_par.h"


!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: num1stOEqnsSolve
      integer, intent(in)  :: numGridCellsInBlock
      integer, intent(in)  :: diurnalGasChemType

      real*8,  intent(inout) :: valuesDecomposedMatrix  (KBLOOP, 0:MXARRAY)
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
        IJTLOOP: do ijt = ijtlo(j,diurnalGasChemType), ijthi(j,diurnalGasChemType)
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
            do k = 1, numGridCellsInBlock
              valuesDecomposedMatrix(k,ij) =  & !ij is nth location of this matrix
     &          valuesDecomposedMatrix(k,ij) -  &
     &          (valuesDecomposedMatrix(k,ik0) * valuesDecomposedMatrix(k,kj0)) -  &
     &          (valuesDecomposedMatrix(k,ik1) * valuesDecomposedMatrix(k,kj1)) -  &
     &          (valuesDecomposedMatrix(k,ik2) * valuesDecomposedMatrix(k,kj2)) -  &
     &          (valuesDecomposedMatrix(k,ik3) * valuesDecomposedMatrix(k,kj3)) -  &
     &          (valuesDecomposedMatrix(k,ik4) * valuesDecomposedMatrix(k,kj4))
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

            do k = 1, numGridCellsInBlock
              valuesDecomposedMatrix(k,ij) =  &
     &          valuesDecomposedMatrix(k,ij) -  &
     &          (valuesDecomposedMatrix(k,ik0) * valuesDecomposedMatrix(k,kj0)) -  &
     &          (valuesDecomposedMatrix(k,ik1) * valuesDecomposedMatrix(k,kj1)) -  &
     &          (valuesDecomposedMatrix(k,ik2) * valuesDecomposedMatrix(k,kj2)) -  &
     &          (valuesDecomposedMatrix(k,ik3) * valuesDecomposedMatrix(k,kj3))
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

            do k = 1, numGridCellsInBlock
              valuesDecomposedMatrix(k,ij) =  &
     &          valuesDecomposedMatrix(k,ij) -  &
     &          (valuesDecomposedMatrix(k,ik0) * valuesDecomposedMatrix(k,kj0)) -  &
     &          (valuesDecomposedMatrix(k,ik1) * valuesDecomposedMatrix(k,kj1)) -  &
     &          (valuesDecomposedMatrix(k,ik2) * valuesDecomposedMatrix(k,kj2))
            end do

          end do

!         -- Sum 2 terms at a time. --

          do ic = il2, ih2

            ik0 = ikdeca(ic)
            ik1 = ikdecb(ic)

            kj0 = kjdeca(ic)
            kj1 = kjdecb(ic)

            do k = 1, numGridCellsInBlock
              valuesDecomposedMatrix(k,ij) =  &
     &          valuesDecomposedMatrix(k,ij) -  &
     &          (valuesDecomposedMatrix(k,ik0) * valuesDecomposedMatrix(k,kj0)) -  &
     &          (valuesDecomposedMatrix(k,ik1) * valuesDecomposedMatrix(k,kj1))
            end do

          end do

!         -- Sum 1 term  at a time. --

          do ic = il1, ih1

            ik0 = ikdeca(ic)

            kj0 = kjdeca(ic)

            do k = 1, numGridCellsInBlock
              valuesDecomposedMatrix(k,ij) =  &
     &          valuesDecomposedMatrix(k,ij) -  &
     &          (valuesDecomposedMatrix(k,ik0) * valuesDecomposedMatrix(k,kj0))
            end do

          end do

!       ==============
        end do IJTLOOP
!       ==============

        iar = diagonalTermDecomp(j,diurnalGasChemType)

        do k = 1, numGridCellsInBlock
          vdiag(k,j) = 1.0d0 / valuesDecomposedMatrix(k,iar)
        end do

!       ----------------------------
!       Second loop of decompostion.
!       ----------------------------

        jl = jloz1(j,diurnalGasChemType)
        jh = jhiz1(j,diurnalGasChemType)

        do jc = jl, jh

          ija = jzeroa(jc)

          do k = 1, numGridCellsInBlock
            valuesDecomposedMatrix(k,ija) = valuesDecomposedMatrix(k,ija) * vdiag(k,j)
          end do

        end do

!     ============
      end do JLOOP
!     ============


      return

      end subroutine LU_Decomp

end module GmiSparseMatrix_mod
