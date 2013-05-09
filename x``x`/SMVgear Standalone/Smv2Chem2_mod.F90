
!=============================================================================
!
! $Id: Smv2Chem2_mod,v 1.1.1.1 2013-05-09 16:06:36 mrdamon Exp $
!
! CODE DEVELOPER
!   John Tannahill, LLNL (Original code from Mark Jacobson, UCLA/Stanford)
!   jrt@llnl.gov
!
! FILE
!   Smv2Chem2_mod.F90
!
! DESCRIPTION
!   This include file contains the data for SmvgearII for the
!   variables that need to be saved between calls.
!=============================================================================

module Smv2Chem2_mod

      implicit none
      public

      ! # of rxns with one active reactant
      integer, save, allocatable ::numRxnsOneActiveReactant  (:) !ICP
      integer, save, allocatable ::nallrat(:) !ICP
      ! last reordered rxn # prior to sets of two rxns with two reactants
      integer, save, allocatable ::lastReorderedRxn (:) !ICS)
      ! # of rxns with three active reactants
      integer, save, allocatable ::numRxnsThreeActiveReactant  (:) !ICS)
      ! # of rxns with two active reactants
      integer, save, allocatable ::numRxnsTwoActiveReactant  (:) !ICS)
      integer, save, allocatable ::nm3bod (:) !ICS)
      ! # rxns where spc in third position is M = O2 + N2
      integer, save, allocatable ::numRxns3rdSpcO2plusN2  (:) !ICS)
      ! # rxns where spc in third position is N2
      integer, save, allocatable ::numRxns3rdSpcN2 (:) ! (ICS)
      ! # rxns where spc in third position is O2
      integer, save, allocatable ::numRxns3rdSpcO2 (:) ! (ICS)
      ! # occurences of  spc in third position that are not O2, N2,
      ! or M, or of spc in any position that are inactive (?);
      ! # of occurrences where inactive spc appears in rate eqn (?)
      integer, save, allocatable ::nmoth  (:) !ICS)
      ! # of kinetic + photo rxns
      integer, save, allocatable ::numKinPhotoRxns(:) !ICS)
      !maps original spc #s to spc #s reordered for chemistry
      integer, save, allocatable ::origSpcToReordSpcMap   (:,:) !MXGSAER, ICS)
      ! lgasbino = jokd (this is what was in the common block comments)
      integer, save, allocatable ::lgasbino(:,:) !MAXGL2,  ICS)
      ! old rxn rate # corresponding to each reordered rxn
      integer, save, allocatable ::oldRxnRateNum(:,:) !NMTRATE, ICS)
      ! new spc # of each active product in each rxn
      integer, save, allocatable ::newSpcNumActiveProduct(:,:,:) !NMRPROD, NMTRATE, ICS)


      integer, save, allocatable ::nreacoth(:,:) !MAXGL2,  ICS)
      integer, save, allocatable ::lgas3bod(:,:) !MAXGL3,  ICS)
      integer, save, allocatable ::losinacp(:,:) !MAXGL3,  ICS)
      integer, save, allocatable ::nreac3b (:,:) !MAXGL3,  ICS)
      integer, save, allocatable ::nreacair(:,:) !MAXGL3,  ICS)
      integer, save, allocatable ::nreacn2 (:,:) !MAXGL3,  ICS)
      integer, save, allocatable ::nreaco2 (:,:) !MAXGL3,  ICS)
      integer, save, allocatable ::jphotnk (:,:) !NMTRATE, ICS)


      ! # of original nspec spc with >= 1 pd term
      integer, save, allocatable ::numOrigSpcGtrOrEql1PdTerm(:) !ICS)
      integer, save, allocatable ::kzthi(:) !ICP)
      integer, save, allocatable ::kztlo(:) !ICP)

      integer, save, allocatable ::ikztot(:) !MXCOUNT4)

      integer, save, allocatable ::kbh1(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbh2(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbh3(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbh4(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbh5(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbl1(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbl2(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbl3(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbl4(:) ! (MXCOUNT4)
      integer, save, allocatable ::kbl5(:) ! (MXCOUNT4)

      integer, save, allocatable ::mbh1(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbh2(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbh3(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbh4(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbh5(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbl1(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbl2(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbl3(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbl4(:) ! (MXCOUNT4)
      integer, save, allocatable ::mbl5(:) ! (MXCOUNT4)

      !     kzeroa..e : arrays identifying terms in gloss array
      integer, save, allocatable ::kzeroa(:) ! (MXCOUNT4)
      integer, save, allocatable ::kzerob(:) ! (MXCOUNT4)
      integer, save, allocatable ::kzeroc(:) ! (MXCOUNT4)
      integer, save, allocatable ::kzerod(:) ! (MXCOUNT4)
      integer, save, allocatable ::kzeroe(:) ! (MXCOUNT4)

      integer, save, allocatable ::mzeroa(:) ! (MXCOUNT4)
      integer, save, allocatable ::mzerob(:) ! (MXCOUNT4)
      integer, save, allocatable ::mzeroc(:) ! (MXCOUNT4)
      integer, save, allocatable ::mzerod(:) ! (MXCOUNT4)
      integer, save, allocatable ::mzeroe(:) ! (MXCOUNT4)
      integer, save, allocatable ::imztot(:,:) !MXGSAER, ICP)

      integer, save, allocatable ::ijval (:) !MXCOUNT3)

      ! identifies the array position of each jloz1..jhiz1 term
      integer, save, allocatable ::jzeroa(:) !MXCOUNT3)

      integer, save, allocatable ::idh1  (:) ! MXCOUNT3
      integer, save, allocatable ::idh2  (:) ! MXCOUNT3
      integer, save, allocatable ::idh3  (:) ! MXCOUNT3
      integer, save, allocatable ::idh4  (:) ! MXCOUNT3
      integer, save, allocatable ::idh5(:) ! MXCOUNT3
      integer, save, allocatable ::idl1  (:) ! MXCOUNT3
      integer, save, allocatable ::idl2  (:) ! MXCOUNT3
      integer, save, allocatable ::idl3  (:) ! MXCOUNT3
      integer, save, allocatable ::idl4  (:) ! MXCOUNT3
      integer, save, allocatable ::idl5(:) ! MXCOUNT3

      integer, save, allocatable ::ikdeca(:) ! MXCOUNT3
      integer, save, allocatable ::ikdecb(:) ! MXCOUNT3
      integer, save, allocatable ::ikdecc(:) ! MXCOUNT3
      integer, save, allocatable ::ikdecd(:) ! MXCOUNT3
      integer, save, allocatable ::ikdece(:) !MXCOUNT3)

      integer, save, allocatable ::kjdeca(:) !MXCOUNT3)
      integer, save, allocatable ::kjdecb(:) ! MXCOUNT3
      integer, save, allocatable ::kjdecc(:) ! MXCOUNT3
      integer, save, allocatable ::kjdecd(:) ! MXCOUNT3
      integer, save, allocatable ::kjdece(:) ! MXCOUNT3

      integer, save, allocatable ::ijthi   (:, :) ! (MXGSAER, ICP)
      integer, save, allocatable ::ijtlo(:, :) ! (MXGSAER, ICP)
      ! diagonal term of decompostion
      integer, save, allocatable ::diagonalTermDecomp(:,:) !MXGSAER, ICP)
      integer, save, allocatable ::jhiz1   (:, :) ! (MXGSAER, ICP)
      integer, save, allocatable ::jloz1(:, :) ! (MXGSAER, ICP)

      ! length of 1d array holding all sparse matrix points =
      ! sparse matrix dimension (?); = total # of matrix
      ! positions filled after matrix processes (?)
      integer, save, allocatable ::sparseMatrixDimension(:) !ICP)
      integer, save, allocatable ::npdhi (:) !ICP)
      integer, save, allocatable ::npdlo(:) !ICP)
      integer, save, allocatable ::nfrhi(:) !ICP)
      integer, save, allocatable ::nfrlo(:) !ICP)
      integer, save, allocatable ::nplhi(:) !ICP)
      integer, save, allocatable ::npllo(:) !ICP)


      integer, save, allocatable ::iialpd  (:) !MXCOUNT2)
      integer, save, allocatable ::ipospd  (:) !MXCOUNT2)
      integer, save, allocatable ::nkpdterm(:) !MXCOUNT2)

      integer, save, allocatable ::jspcnfr(:) ! (MXCOUNT4)
      integer, save, allocatable ::jspnpl(:) ! (MXCOUNT4)
      integer, save, allocatable ::nknfr  (:) !MXCOUNT4)

      ! reaordered rxn rate #s for each loss (and prod) term
      integer, save, allocatable ::reOrderedRxnRateEachLossA (:) !MXCOUNT4)
      integer, save, allocatable ::reOrderedRxnRateEachLossB (:) ! (MXCOUNT4)
      integer, save, allocatable ::reOrderedRxnRateEachLossC(:) ! (MXCOUNT4)
      integer, save, allocatable ::reOrderedRxnRateEachLossD (:) ! (MXCOUNT4)
      integer, save, allocatable ::reOrderedRxnRateEachLossE(:) ! (MXCOUNT4)

      integer, save, allocatable ::nph1(:) ! (MXCOUNT4)
      integer, save, allocatable ::nph2(:) ! (MXCOUNT4)
      integer, save, allocatable ::nph3(:) ! (MXCOUNT4)
      integer, save, allocatable ::nph4(:) ! (MXCOUNT4)
      integer, save, allocatable ::nph5(:) ! (MXCOUNT4)
      integer, save, allocatable ::npl1(:) ! (MXCOUNT4)
      integer, save, allocatable ::npl2(:) ! (MXCOUNT4)
      integer, save, allocatable ::npl3(:) ! (MXCOUNT4)
      integer, save, allocatable ::npl4(:) ! (MXCOUNT4)
      integer, save, allocatable ::npl5(:) ! (MXCOUNT4)

      integer, save, allocatable ::nolosp(:) !ICP)

      ! new rxn rate # corresponding to each original rate #
      integer, save, allocatable ::newRxnRateNumber(:,:) !NMTRATE*2, ICS)

      integer, save, allocatable ::nknlosp(:,:) !MAXGL3, ICS)

      integer, save, allocatable ::nknphotrt(:,:) !IPHOT,ICS)


      ! 1 / (chemintv * chemintv) (s^-2)
      real*8, save, allocatable ::abst2   (:) !ICS)
      !     errmax   : relative error tolerance; eps should be < 1 for speedy and
      !                reliable results, 10^-3 is reasonable; for many decimal
      !                places of accuracy, decrease eps
      real*8, save, allocatable ::relativeErrorTolerance  (:) !ICS)
      real*8, save, allocatable ::hmaxday (:) !ICS)
      real*8, save, allocatable ::timeintv(:) !ICS)

      ! pre-defined absolute error tolerances; if it is too small,
      ! then integration will take too long; if it is too large,
      ! convergence will be too easy and errors will accumulate, the
      ! time step may be cut too small, and the integration may stop
      ! (delt < hmin or floating point exception in Decomp);
      ! typical gas-phase values of abstol are 10^3 cm^-3 (?),
      ! typical aq -phase values of abstol are
      ! 10^-13 to 10^-15 m l^-1 (?)
      real*8, save, allocatable ::absoluteErrorTolerance(:,:) !6, ICS)
      real*8, save, allocatable ::enqq1 (:) ! (MORDER)
      real*8, save, allocatable ::enqq2 (:) ! (MORDER)
      real*8, save, allocatable ::enqq3 (:) ! (MORDER)
      real*8, save, allocatable ::conp15(:) ! (MORDER)
      real*8, save, allocatable ::conpst(:) ! (MORDER)
      ! coefficients used in selecting step and order (see Ksparse);
      ! coeffsForSelectingStepAndOrder = original pertst^2
      real*8, save, allocatable ::coeffsForSelectingStepAndOrder(:,:) !MORDER, 3)
      !     aset     : coefficients for determining order of integration method and
      !                for calculating matrix P
      real*8, save, allocatable ::coeffsForIntegrationOrder(:,:) !10, 8)

      !     fracpl  : = -1 for all reactants;
      !               = +1 or +fraction for all products
      real*8, save, allocatable ::fracpl (:) ! (MXCOUNT2)
      real*8, save, allocatable ::fracnfr(:) ! (MXCOUNT4)   !end type Smv2Chem2_type



end module Smv2Chem2_mod

