subroutine testFailOnPurpose()
   use pFUnit

   call assertEqual(1, 2, message = 'Intentional failure for instructional purposes.')

end subroutine testFailOnPurpose

subroutine testSpeciesConst(actual, expected)
   use pFUnit
   real*8 :: actual
   real*8 :: expected

   call assertEqual(actual, expected, message = 'Species concentration not the same.')

end subroutine testSpeciesConst
