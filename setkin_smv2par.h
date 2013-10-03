!=======================================================================
!
! $Id: setkin_smv2par.h,v 1.1.1.1.22.1 2009-09-16 16:26:33 enielsen Exp $
!
! FILE
!   setkin_smv2par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    10/2006
!  Reaction dictionary:     GMI_Combo_rxns_124species_JPL06.db
!  Setkin files generated:  Thu Feb  7 17:16:01 2008
!
!========1=========2=========3=========4=========5=========6=========7==

      integer &
     &  SK_IGAS &
     & ,SK_IPHOT &
     & ,SK_ITHERM &
     & ,SK_NACT

      parameter (SK_IGAS   = 121)
      parameter (SK_IPHOT  =  81)
      parameter (SK_ITHERM = 321)
      parameter (SK_NACT   = 117)

!                                  --^--

