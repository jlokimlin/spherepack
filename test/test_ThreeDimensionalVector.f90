program test_ThreeDimensionalVector

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack_wrapper_library, only: &
        Vector => ThreeDimensionalVector

    ! Explicit typing only
    implicit none

    associate( &
        !
        !==> Assigning integers is safe
        !
        A => Vector( 3, 4, 5 ), &
        B => Vector( 4, 3, 5 ), &
        C => Vector( -5, -12, -13 ) &
        )
        write( stdout, '(A)') ''
        write( stdout, '(A)') ' test_ThreeDimensionalVector *** TEST RUN *** '
        write( stdout, '(A)') ''
        write( stdout, '(A,3(1pe15.1))') '          A = ', A
        write( stdout, '(A,3(1pe15.1))') '          B = ', B
        write( stdout, '(A,3(1pe15.1))') '          C = ', C
        write( stdout, '(A,1pe15.1)')    '      A . B = ', A.dot.B
        write( stdout, '(A,3(1pe15.1))') '      A x B = ', A.cross.B
        write( stdout, '(A,1pe15.1)')    'A . (B x C) = ', A.dot.(B.cross.C)
        write( stdout, '(A,3(1pe15.3))') 'A x (B x C) = ', A.cross.(B.cross.C)
        write( stdout, '(A)') ''

    end associate

end program test_ThreeDimensionalVector
