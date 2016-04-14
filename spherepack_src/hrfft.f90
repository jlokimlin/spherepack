!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK version 3.2                   *
!     *                                                               *
!     *       A Package of Fortran77 Subroutines and Programs         *
!     *                                                               *
!     *              for Modeling Geophysical Processes               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *                  John Adams and Paul Swarztrauber             *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! ... file hrfft.f
!
!     this file contains a multiple fft package for spherepack. it
!     includes code and documentation for performing fast fourier
!     transforms (see subroutines hrffti, hrfftf and hrfftb)
!
! **********************************************************************
!
!     subroutine hrffti(n, wsave)
!
!     subroutine hrffti initializes the array wsave which is used in
!     both hrfftf and hrfftb. the prime factorization of n together 
!     with a tabulation of the trigonometric functions are computed and
!     stored in wsave.
!
!     input parameter
!
!     n       the length of the sequence to be transformed.
!
!     output parameter
!
!     wsave   a work array which must be dimensioned at least 2*n+15.
!             the same work array can be used for both hrfftf and 
!             hrfftb as long as n remains unchanged. different wsave 
!             arrays are required for different values of n. the 
!             contents of wsave must not be changed between calls 
!             of hrfftf or hrfftb.
!
! **********************************************************************
!
!     subroutine hrfftf(m, n, r, mdimr, wsave, work)
!
!     subroutine hrfftf computes the fourier coefficients of m real
!     perodic sequences (fourier analysis); i.e. hrfftf computes the
!     real fft of m sequences each with length n. the transform is 
!     defined below at output parameter r.
!
!     input parameters
!
!     m       the number of sequences.
!
!     n       the length of all m sequences.  the method is most
!             efficient when n is a product of small primes. n may
!             change as long as different work arrays are provided
!
!     r       r(m, n) is a two dimensional real array that contains m
!             sequences each with length n.
!
!     mdimr   the first dimension of the r array as it appears
!             in the program that calls hrfftf. mdimr must be
!             greater than or equal to m.
!
!
!     wsave   a work array with at least least 2*n+15 locations
!             in the program that calls hrfftf. the wsave array must be
!             initialized by calling subroutine hrffti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!             the same wsave array can be used by hrfftf and hrfftb.
!
!     work    a real work array with m*n locations.
!
!
!     output parameters
!
!     r      for all j=1, ..., m
!  
!             r(j, 1) = the sum from i=1 to i=n of r(j, i)
!
!             if n is even set l =n/2   , if n is odd set l = (n+1)/2
!
!               then for k = 2, ..., l
!
!                  r(j, 2*k-2) = the sum from i = 1 to i = n of
!
!                       r(j, i)*cos((k-1)*(i-1)*2*pi/n)
!
!                  r(j, 2*k-1) = the sum from i = 1 to i = n of
!
!                      -r(j, i)*sin((k-1)*(i-1)*2*pi/n)
!
!             if n is even
!
!                  r(j, n) = the sum from i = 1 to i = n of
!
!                       (-1)**(i-1)*r(j, i)
!
!      *****  note
!                  this transform is unnormalized since a call of hrfftf
!                  followed by a call of hrfftb will multiply the input
!                  sequence by n.
!
!     wsave   contains results which must not be destroyed between
!             calls of hrfftf or hrfftb.
!
!     work    a real work array with m*n locations that does
!             not have to be saved.
!
! **********************************************************************
!
!     subroutine hrfftb(m, n, r, mdimr, wsave, work)
!
!     subroutine hrfftb computes the real perodic sequence of m
!     sequences from their fourier coefficients (fourier synthesis). 
!     the transform is defined below at output parameter r.
!
!     input parameters
!
!     m       the number of sequences.
!
!     n       the length of all m sequences.  the method is most
!             efficient when n is a product of small primes. n may
!             change as long as different work arrays are provided
!
!     r       r(m, n) is a two dimensional real array that contains
!             the fourier coefficients of m sequences each with 
!             length n.
!
!     mdimr   the first dimension of the r array as it appears
!             in the program that calls hrfftb. mdimr must be
!             greater than or equal to m.
!
!     wsave   a work array which must be dimensioned at least 2*n+15.
!             in the program that calls hrfftb. the wsave array must be
!             initialized by calling subroutine hrffti(n, wsave) and a
!             different wsave array must be used for each different
!             value of n. this initialization does not have to be
!             repeated so long as n remains unchanged thus subsequent
!             transforms can be obtained faster than the first.
!             the same wsave array can be used by hrfftf and hrfftb.
!
!     work    a real work array with m*n locations.
!
!
!     output parameters
!
!     r      for all j=1, ..., m
!  
!             for n even and for i = 1, ..., n
!
!                  r(j, i) = r(j, 1)+(-1)**(i-1)*r(j, n)
!
!                       plus the sum from k=2 to k=n/2 of
!
!                        2.0*r(j, 2*k-2)*cos((k-1)*(i-1)*2*pi/n)
!
!                       -2.0*r(j, 2*k-1)*sin((k-1)*(i-1)*2*pi/n)
!
!             for n odd and for i = 1, ..., n
!
!                  r(j, i) = r(j, 1) plus the sum from k=2 to k=(n+1)/2 of
!
!                       2.0*r(j, 2*k-2)*cos((k-1)*(i-1)*2*pi/n)
!
!                      -2.0*r(j, 2*k-1)*sin((k-1)*(i-1)*2*pi/n)
!
!      *****  note
!                  this transform is unnormalized since a call of hrfftf
!                  followed by a call of hrfftb will multiply the input
!                  sequence by n.
!
!     wsave   contains results which must not be destroyed between
!             calls of hrfftb or hrfftf.
!
!     work    a real work array with m*n locations that does not
!             have to be saved
!
!
!
subroutine hrffti(n, wsave)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: n
    real (wp),    intent (in out) :: wsave(n+15)
    !----------------------------------------------------------------------

    if (n == 1) then
        return                                                    !
    end if

    call hrfti1(n, wsave(1), wsave(n+1))
                                                                 !
end subroutine hrffti



subroutine hrfti1(n, wa, fac)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: n
    real (wp),    intent (in out) :: wa(n)
    real (wp),    intent (in out) :: fac(15)
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), parameter :: ntryh(4) = [4, 2, 3, 5]
    integer (ip)            :: i, j, k1, l1, l2, ib
    integer (ip)            :: ld, ii, nf, ip_rename, nl, is, nq, nr
    real (wp),    parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
    integer (ip)            :: fi, ido, ipm, nfm1, ntry
    real (wp)               :: arg, argh, argld
    !----------------------------------------------------------------------
                     !
    nl = n                                                                  !
    nf = 0
    j = 0
101 j = j+1

    if (j-4 <= 0) then
        goto 102
    else
        goto 103
    end if

102 ntry = ntryh(j)
    go to 104
103 ntry = ntry+2
104 nq = nl/ntry
    nr = nl-ntry*nq

    if (nr < 0) then
        goto 101
    else if (nr == 0) then
        goto 105
    else
        goto 101
    end if

105 nf = nf+1
    fac(nf+2) = real(ntry, kind=wp)
    nl = nq

    if (ntry /= 2) then
        go to 107
    end if

    if (nf == 1) then
        go to 107
    end if

    do i=2, nf
        ib = nf-i+2
        fac(ib+2) = fac(ib+1)
    end do
    fac(3) = 2.0_wp

107 if (nl /= 1) then
        go to 104
    end if

    fac(1) = real(n, kind=wp)
    fac(2) = real(nf, kind=wp)
    argh = TWO_PI/n
    is = 0
    nfm1 = nf-1
    l1 = 1

    if (nfm1 == 0) then
        return
    end if

    do k1=1, nfm1
        ip_rename = int(fac(k1+2), kind=ip)
        ld = 0
        l2 = l1*ip_rename
        ido = n/l2
        ipm = ip_rename-1
        do j=1, ipm
            ld = ld+l1
            i = is
            argld = real(ld, kind=wp)*argh
            fi = 0.0_wp
            do ii=3, ido, 2
                i = i + 2
                fi = fi + 1.0_wp
                arg = fi * argld
                wa(i-1) = cos(arg)
                wa(i) = sin(arg)
            end do
            is = is+ido
        end do
        l1 = l2
    end do

end subroutine hrfti1



subroutine hrfftf(m, n, r, mdimr, whrfft, work)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: m
    integer (ip), intent (in)     :: n
    real (wp),    intent (in out) :: r(mdimr, n)
    integer (ip), intent (in)     :: mdimr
    real (wp),    intent (in out) :: work(1)
    real (wp),    intent (in out) :: whrfft(n+15)
    !----------------------------------------------------------------------

    if (n == 1) then
        return
    end if

    call hrftf1(m, n, r, mdimr, work, whrfft, whrfft(n+1))

end subroutine hrfftf



subroutine hrftf1(m, n, c, mdimc, ch, wa, fac)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: m
    integer (ip), intent (in)     :: n
    real (wp),    intent (in out) :: ch(m, n)
    integer (ip), intent (in)     :: mdimc
    real (wp),    intent (in out) :: c(mdimc, n)
    real (wp),    intent (in out) :: wa(n)
    real (wp),    intent (in out) :: fac(15)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip) :: i, j, k1, l1, l2
    integer (ip) :: na, kh, nf, ip_rename
    integer (ip) :: iw, ix2, ix3, ix4, ido, idl1
    !----------------------------------------------------------------------

    nf = int(fac(2), kind=ip)
    na = 1
    l2 = n
    iw = n

    do k1=1, nf
        kh = nf-k1
        ip_rename = int(fac(kh+3), kind=ip)
        l1 = l2/ip_rename
        ido = n/l2
        idl1 = ido*l1
        iw = iw-(ip_rename-1)*ido
        na = 1-na

        if (ip_rename /= 4) then
            go to 102
        end if

        ix2 = iw+ido
        ix3 = ix2+ido

        if (na /= 0) then
            go to 101
        end if

        call  hradf4(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2), wa(ix3))
        go to 110

101     call hradf4(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2), wa(ix3))
        go to 110

102     if (ip_rename /= 2) then
            go to 104
        end if

        if (na /= 0) then
            go to 103
        end if

        call  hradf2(m, ido, l1, c, mdimc, ch, m, wa(iw))
        go to 110

103     call hradf2(m, ido, l1, ch, m, c, mdimc, wa(iw))
        go to 110

104     if (ip_rename /= 3) then
            go to 106
        end if

        ix2 = iw+ido

        if (na /= 0) then
            go to 105
        end if

        call  hradf3(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2))
        go to 110

105     call hradf3(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2))
        go to 110

106     if (ip_rename /= 5) then
            go to 108
        end if

        ix2 = iw+ido
        ix3 = ix2+ido
        ix4 = ix3+ido

        if (na /= 0) then
            go to 107
        end if

        call hradf5(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2), wa(ix3), wa(ix4))
        go to 110

107     call hradf5(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2), wa(ix3), wa(ix4))
        go to 110

108     if (ido == 1) then
            na = 1-na
        end if

        if (na /= 0) then
            go to 109
        end if

        call  hradfg(m, ido, ip_rename, l1, idl1, c, c, c, mdimc, ch, ch, m, wa(iw))
        na = 1
        go to 110

109     call hradfg(m, ido, ip_rename, l1, idl1, ch, ch, ch, m, c, c, mdimc, wa(iw))
        na = 0

110     l2 = l1
    end do

    if (na == 1) then
        return
    end if

    c(1:m, 1:n) = ch

end subroutine hrftf1



subroutine hradf4(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2, wa3)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: l1
    real (wp),    intent (in out) :: cc(mdimcc, ido, l1, 4)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: ch(mdimch, ido, 4, l1)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: wa1(ido)
    real (wp),    intent (in out) :: wa2(ido)
    real (wp),    intent (in out) :: wa3(ido)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip)         :: i, k, m, ic, idp2
    real (wp), parameter :: HALF_SQRT2 = sqrt(2.0_wp)/2
    !----------------------------------------------------------------------

    do k=1, l1
        do m=1, mp
            ch(m, 1, 1, k) = &
                (cc(m, 1, k, 2)+cc(m, 1, k, 4)) &
                +(cc(m, 1, k, 1) + cc(m, 1, k, 3))
            ch(m, ido, 4, k) = &
                (cc(m, 1, k, 1)+cc(m, 1, k, 3)) &
                -(cc(m, 1, k, 2) + cc(m, 1, k, 4))
            ch(m, ido, 2, k) = &
                cc(m, 1, k, 1)-cc(m, 1, k, 3)
            ch(m, 1, 3, k) = &
                cc(m, 1, k, 4)-cc(m, 1, k, 2)
        end do
    end do

    if (ido-2 < 0) then
        return
    else if (ido-2 /= 0) then
        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                do m=1, mp

                    ch(m, i-1, 1, k) = &
                        ((wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &              !
                        cc(m, i, k, 2))+(wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &
                        cc(m, i, k, 4)))+(cc(m, i-1, k, 1)+(wa2(i-2)*cc(m, i-1, k, 3)+ &          !
                        wa2(i-1)*cc(m, i, k, 3)))

                    ch(m, ic-1, 4, k) = &
                        (cc(m, i-1, k, 1)+(wa2(i-2)*cc(m, i-1, k, 3)+ &        !
                        wa2(i-1)*cc(m, i, k, 3)))-((wa1(i-2)*cc(m, i-1, k, 2)+ &               !
                        wa1(i-1)*cc(m, i, k, 2))+(wa3(i-2)*cc(m, i-1, k, 4)+ &
                        wa3(i-1)*cc(m, i, k, 4)))

                    ch(m, i, 1, k) = &
                        ((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &
                        cc(m, i-1, k, 2))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4)))+(cc(m, i, k, 1)+(wa2(i-2)*cc(m, i, k, 3)- &            !
                        wa2(i-1)*cc(m, i-1, k, 3)))

                    ch(m, ic, 4, k) = &
                        ((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &
                        cc(m, i-1, k, 2))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4)))-(cc(m, i, k, 1)+(wa2(i-2)*cc(m, i, k, 3)- &            !
                        wa2(i-1)*cc(m, i-1, k, 3)))

                    ch(m, i-1, 3, k) = &
                        ((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &
                        cc(m, i-1, k, 2))-(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4)))+(cc(m, i-1, k, 1)-(wa2(i-2)*cc(m, i-1, k, 3)+ &        !
                        wa2(i-1)*cc(m, i, k, 3)))

                    ch(m, ic-1, 2, k) = &
                        (cc(m, i-1, k, 1)-(wa2(i-2)*cc(m, i-1, k, 3)+ &        !
                        wa2(i-1)*cc(m, i, k, 3)))-((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &        !
                        cc(m, i-1, k, 2))-(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                        cc(m, i-1, k, 4)))

                    ch(m, i, 3, k) = &
                        ((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &
                        cc(m, i, k, 4))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                        cc(m, i, k, 2)))+(cc(m, i, k, 1)-(wa2(i-2)*cc(m, i, k, 3)- &              !
                        wa2(i-1)*cc(m, i-1, k, 3)))

                    ch(m, ic, 2, k) = &
                        ((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &               !
                        cc(m, i, k, 4))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                        cc(m, i, k, 2)))-(cc(m, i, k, 1)-(wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &     !
                        cc(m, i-1, k, 3)))
                end do
            end do
        end do

        if(mod(ido, 2) == 1) then
            return
        end if
    end if

    do k=1, l1
        do m=1, mp

            ch(m, ido, 1, k) = &
                (HALF_SQRT2*(cc(m, ido, k, 2)-cc(m, ido, k, 4)))+ &          !
                cc(m, ido, k, 1)

            ch(m, ido, 3, k) = &
                cc(m, ido, k, 1)-(HALF_SQRT2*(cc(m, ido, k, 2)- &            !
                cc(m, ido, k, 4)))

            ch(m, 1, 2, k) = &
                (-HALF_SQRT2*(cc(m, ido, k, 2)+cc(m, ido, k, 4)))- &           !
                cc(m, ido, k, 3)

            ch(m, 1, 4, k) = &
                (-HALF_SQRT2*(cc(m, ido, k, 2)+cc(m, ido, k, 4)))+ &           !
                cc(m, ido, k, 3)
        end do
    end do

end subroutine hradf4



subroutine hradf2(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: l1
    real (wp),    intent (in out) :: ch(mdimch, ido, 2, l1)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: cc(mdimcc, ido, l1, 2)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: wa1(ido)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip) :: i, k, m, ic, idp2
    !----------------------------------------------------------------------

    do k=1, l1
        do m=1, mp
            ch(m, 1, 1, k) = cc(m, 1, k, 1)+cc(m, 1, k, 2)
            ch(m, ido, 2, k) = cc(m, 1, k, 1)-cc(m, 1, k, 2)
        end do
    end do


    if (ido-2 < 0) then
        return
    else if (ido-2 /= 0) then
        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                do m=1, mp

                    ch(m, i, 1, k) = &
                        cc(m, i, k, 1)+(wa1(i-2)*cc(m, i, k, 2)- &
                        wa1(i-1)*cc(m, i-1, k, 2))

                    ch(m, ic, 2, k) = &
                        (wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &
                        cc(m, i-1, k, 2))-cc(m, i, k, 1)

                    ch(m, i-1, 1, k) = &
                        cc(m, i-1, k, 1)+(wa1(i-2)*cc(m, i-1, k, 2)+ &          !
                        wa1(i-1)*cc(m, i, k, 2))

                    ch(m, ic-1, 2, k) = &
                        cc(m, i-1, k, 1)-(wa1(i-2)*cc(m, i-1, k, 2)+ &         !
                        wa1(i-1)*cc(m, i, k, 2))
                end do
            end do
        end do
        if(mod(ido, 2) == 1) then
            return
        end if
    end if

    do k=1, l1
        do m=1, mp
            ch(m, 1, 2, k) = -cc(m, ido, k, 2)
            ch(m, ido, 1, k) = cc(m, ido, k, 1)
        end do
    end do

end subroutine hradf2



subroutine hradf3(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: l1
    real (wp),    intent (in out) :: ch(mdimch, ido, 3, l1)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: cc(mdimcc, ido, l1, 3)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: wa1(ido)
    real (wp),    intent (in out) :: wa2(ido)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip)         :: i, k, m, ic, idp2
    real (wp), parameter :: ARG=2.0_wp * acos(-1.0_wp)/3
    real (wp), parameter :: TAUR=cos(ARG)
    real (wp), parameter :: TAUI=sin(ARG)
    !----------------------------------------------------------------------

    do k=1, l1
        do m=1, mp

            ch(m, 1, 1, k) = &
                cc(m, 1, k, 1)+(cc(m, 1, k, 2)+cc(m, 1, k, 3))

            ch(m, 1, 3, k) = &
                TAUI*(cc(m, 1, k, 3)-cc(m, 1, k, 2))

            ch(m, ido, 2, k) = &
                cc(m, 1, k, 1)+TAUR* &
                (cc(m, 1, k, 2)+cc(m, 1, k, 3))
        end do
    end do

    if (ido == 1) then
        return
    end if

    idp2 = ido+2
    do k=1, l1
        do i=3, ido, 2
            ic = idp2-i
            do m=1, mp

                ch(m, i-1, 1, k) = &
                    cc(m, i-1, k, 1)+((wa1(i-2)*cc(m, i-1, k, 2)+ &         !
                    wa1(i-1)*cc(m, i, k, 2))+(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &        !
                    cc(m, i, k, 3)))

                ch(m, i, 1, k) = &
                    cc(m, i, k, 1)+((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &      !
                    cc(m, i-1, k, 2))+(wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &
                    cc(m, i-1, k, 3)))

                ch(m, i-1, 3, k) = &
                    (cc(m, i-1, k, 1)+TAUR*((wa1(i-2)* &
                    cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2))+(wa2(i-2)* &
                    cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3))))+(TAUI*((wa1(i-2)* &        !
                    cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2))-(wa2(i-2)* &
                    cc(m, i, k, 3)-wa2(i-1)*cc(m, i-1, k, 3))))

                ch(m, ic-1, 2, k) = &
                    (cc(m, i-1, k, 1)+TAUR*((wa1(i-2)* &
                    cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2))+(wa2(i-2)* &
                    cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3))))-(TAUI*((wa1(i-2)* &
                    cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2))-(wa2(i-2)* &
                    cc(m, i, k, 3)-wa2(i-1)*cc(m, i-1, k, 3))))

                ch(m, i, 3, k) = &
                    (cc(m, i, k, 1)+TAUR*((wa1(i-2)*cc(m, i, k, 2)- &
                    wa1(i-1)*cc(m, i-1, k, 2))+(wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &
                    cc(m, i-1, k, 3))))+(TAUI*((wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                    cc(m, i, k, 3))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                    cc(m, i, k, 2))))

                ch(m, ic, 2, k) = &
                    (TAUI*((wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                    cc(m, i, k, 3))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                    cc(m, i, k, 2))))-(cc(m, i, k, 1)+TAUR*((wa1(i-2)*cc(m, i, k, 2)- &
                    wa1(i-1)*cc(m, i-1, k, 2))+(wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &
                    cc(m, i-1, k, 3))))
            end do
        end do
    end do

end subroutine hradf3




subroutine hradf5(mp, ido, l1, cc, mdimcc, ch, mdimch, &
    wa1, wa2, wa3, wa4)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: l1
    real (wp),    intent (in out) :: ch(mdimch, ido, 5, l1)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: cc(mdimcc, ido, l1, 5)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: wa1(ido)
    real (wp),    intent (in out) :: wa2(ido)
    real (wp),    intent (in out) :: wa3(ido)
    real (wp),    intent (in out) :: wa4(ido)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip)         :: i, k, m, ic, idp2
    real (wp), parameter :: ARG = 2.0_wp * acos(-1.0_wp)/5
    real (wp), parameter :: TR11=cos(ARG)
    real (wp), parameter :: TI11=sin(ARG)
    real (wp), parameter :: TR12=cos(2.0*ARG)
    real (wp), parameter :: TI12=sin(2.0*ARG)
    !----------------------------------------------------------------------

    do k=1, l1
        do m=1, mp

            ch(m, 1, 1, k) = &
                cc(m, 1, k, 1)+(cc(m, 1, k, 5)+cc(m, 1, k, 2))+ &
                (cc(m, 1, k, 4)+cc(m, 1, k, 3))

            ch(m, ido, 2, k) = &
                cc(m, 1, k, 1)+TR11*(cc(m, 1, k, 5)+cc(m, 1, k, 2))+ &
                TR12*(cc(m, 1, k, 4)+cc(m, 1, k, 3))

            ch(m, 1, 3, k) = &
                TI11*(cc(m, 1, k, 5)-cc(m, 1, k, 2))+TI12* &
                (cc(m, 1, k, 4)-cc(m, 1, k, 3))

            ch(m, ido, 4, k) = &
                cc(m, 1, k, 1)+TR12*(cc(m, 1, k, 5)+cc(m, 1, k, 2))+ &        !
                TR11*(cc(m, 1, k, 4)+cc(m, 1, k, 3))

            ch(m, 1, 5, k) = &
                TI12*(cc(m, 1, k, 5)-cc(m, 1, k, 2))-TI11* &
                (cc(m, 1, k, 4)-cc(m, 1, k, 3))
        end do
    end do

    if (ido == 1) then
        return
    end if

    idp2 = ido+2
    do k=1, l1
        do i=3, ido, 2
            ic = idp2-i
            do m=1, mp

                ch(m, i-1, 1, k) = &
                    cc(m, i-1, k, 1)+((wa1(i-2)*cc(m, i-1, k, 2)+ &         !
                    wa1(i-1)*cc(m, i, k, 2))+(wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)* &        !
                    cc(m, i, k, 5)))+((wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &               !
                    cc(m, i, k, 3))+(wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4)))      !

                ch(m, i, 1, k) = &
                    cc(m, i, k, 1)+((wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)* &      !
                    cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &
                    cc(m, i-1, k, 5)))+((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &               !
                    cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                    cc(m, i-1, k, 4)))

                ch(m, i-1, 3, k) = &
                    cc(m, i-1, k, 1)+TR11* &
                    ( wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2) &
                    +wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)*cc(m, i, k, 5))+TR12* &            !
                    ( wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3) &
                    +wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4))+TI11* &            !
                    ( wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2) &
                    -(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)*cc(m, i-1, k, 5)))+TI12* &          !
                    ( wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)*cc(m, i-1, k, 3) &
                    -(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)*cc(m, i-1, k, 4)))

                ch(m, ic-1, 2, k) = &
                    cc(m, i-1, k, 1)+TR11* &
                    ( wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2) &
                    +wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)*cc(m, i, k, 5))+TR12* &            !
                    ( wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3) &
                    +wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4))-(TI11* &            !
                    ( wa1(i-2)*cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2) &
                    -(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)*cc(m, i-1, k, 5)))+TI12* &          !
                    ( wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)*cc(m, i-1, k, 3) &
                    -(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)*cc(m, i-1, k, 4))))                 !

                ch(m, i, 3, k) = &
                    (cc(m, i, k, 1)+TR11*((wa1(i-2)*cc(m, i, k, 2)- &         !
                    wa1(i-1)*cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &        !
                    cc(m, i-1, k, 5)))+TR12*((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &          !
                    cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                    cc(m, i-1, k, 4))))+(TI11*((wa4(i-2)*cc(m, i-1, k, 5)+ &               !
                    wa4(i-1)*cc(m, i, k, 5))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &        !
                    cc(m, i, k, 2)))+TI12*((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &          !
                    cc(m, i, k, 4))-(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                    cc(m, i, k, 3))))

                ch(m, ic, 2, k) = &
                    (TI11*((wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)* &         !
                    cc(m, i, k, 5))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                    cc(m, i, k, 2)))+TI12*((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &          !
                    cc(m, i, k, 4))-(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                    cc(m, i, k, 3))))-(cc(m, i, k, 1)+TR11*((wa1(i-2)*cc(m, i, k, 2)- &       !
                    wa1(i-1)*cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &        !
                    cc(m, i-1, k, 5)))+TR12*((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &          !
                    cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                    cc(m, i-1, k, 4))))

                ch(m, i-1, 5, k) = &
                    (cc(m, i-1, k, 1)+TR12*((wa1(i-2)* &
                    cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2))+(wa4(i-2)* &
                    cc(m, i-1, k, 5)+wa4(i-1)*cc(m, i, k, 5)))+TR11*((wa2(i-2)* &          !
                    cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3))+(wa3(i-2)* &
                    cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4))))+(TI12*((wa1(i-2)* &        !
                    cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2))-(wa4(i-2)*cc(m, i, k, 5)- &     !
                    wa4(i-1)*cc(m, i-1, k, 5)))-TI11*((wa2(i-2)*cc(m, i, k, 3)- &          !
                    wa2(i-1)*cc(m, i-1, k, 3))-(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &        !
                    cc(m, i-1, k, 4))))

                ch(m, ic-1, 4, k) = &
                    (cc(m, i-1, k, 1)+TR12*((wa1(i-2)* &
                    cc(m, i-1, k, 2)+wa1(i-1)*cc(m, i, k, 2))+(wa4(i-2)* &
                    cc(m, i-1, k, 5)+wa4(i-1)*cc(m, i, k, 5)))+TR11*((wa2(i-2)* &          !
                    cc(m, i-1, k, 3)+wa2(i-1)*cc(m, i, k, 3))+(wa3(i-2)* &
                    cc(m, i-1, k, 4)+wa3(i-1)*cc(m, i, k, 4))))-(TI12*((wa1(i-2)* &        !
                    cc(m, i, k, 2)-wa1(i-1)*cc(m, i-1, k, 2))-(wa4(i-2)*cc(m, i, k, 5)- &     !
                    wa4(i-1)*cc(m, i-1, k, 5)))-TI11*((wa2(i-2)*cc(m, i, k, 3)- &          !
                    wa2(i-1)*cc(m, i-1, k, 3))-(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &        !
                    cc(m, i-1, k, 4))))

                ch(m, i, 5, k) = &
                    (cc(m, i, k, 1)+TR12*((wa1(i-2)*cc(m, i, k, 2)- &         !
                    wa1(i-1)*cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &        !
                    cc(m, i-1, k, 5)))+TR11*((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &          !
                    cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                    cc(m, i-1, k, 4))))+(TI12*((wa4(i-2)*cc(m, i-1, k, 5)+ &               !
                    wa4(i-1)*cc(m, i, k, 5))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &        !
                    cc(m, i, k, 2)))-TI11*((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &          !
                    cc(m, i, k, 4))-(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                    cc(m, i, k, 3))))

                ch(m, ic, 4, k) = &
                    (TI12*((wa4(i-2)*cc(m, i-1, k, 5)+wa4(i-1)* &         !
                    cc(m, i, k, 5))-(wa1(i-2)*cc(m, i-1, k, 2)+wa1(i-1)* &
                    cc(m, i, k, 2)))-TI11*((wa3(i-2)*cc(m, i-1, k, 4)+wa3(i-1)* &          !
                    cc(m, i, k, 4))-(wa2(i-2)*cc(m, i-1, k, 3)+wa2(i-1)* &
                    cc(m, i, k, 3))))-(cc(m, i, k, 1)+TR12*((wa1(i-2)*cc(m, i, k, 2)- &       !
                    wa1(i-1)*cc(m, i-1, k, 2))+(wa4(i-2)*cc(m, i, k, 5)-wa4(i-1)* &        !
                    cc(m, i-1, k, 5)))+TR11*((wa2(i-2)*cc(m, i, k, 3)-wa2(i-1)* &          !
                    cc(m, i-1, k, 3))+(wa3(i-2)*cc(m, i, k, 4)-wa3(i-1)* &
                    cc(m, i-1, k, 4))))
            end do
        end do
    end do

end subroutine hradf5



subroutine hradfg(mp, ido, ip_rename, l1, idl1, cc, c1, c2, mdimcc, &
    ch, ch2, mdimch, wa)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: ip_rename
    integer (ip), intent (in)     :: l1
    integer (ip), intent (in)     :: idl1
    real (wp),    intent (in out) :: cc(mdimcc, ido, ip_rename, l1)
    real (wp),    intent (in out) :: c1(mdimcc, ido, l1, ip_rename)
    real (wp),    intent (in out) :: c2(mdimcc, idl1, ip_rename)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: ch(mdimch, ido, l1, ip_rename)
    real (wp),    intent (in out) :: ch2(mdimch, idl1, ip_rename)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: wa(ido)
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip)         :: i, j, k, l, m, j2, ic, jc, lc, ik, is, idij
    real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
    real (wp)            :: dc2, ai1, ai2, ar1, ar2, ds2
    real (wp)            :: ar1h, ar2h
    !----------------------------------------------------------------------

    associate( &
        ipph => (ip_rename+1)/2, &
        ipp2 => ip_rename+2, &
        idp2 => ido+2, &
        nbd => (ido-1)/2, &
        arg => TWO_PI/ip_rename &
        )
        associate( &
            dcp => cos(arg), &
            dsp => sin(arg) &
            )

            if (ido /= 1) then
                do ik=1, idl1
                    do m=1, mp
                        ch2(m, ik, 1) = c2(m, ik, 1)
                    end do
                end do
                do j=2, ip_rename
                    do k=1, l1
                        do m=1, mp
                            ch(m, 1, k, j) = c1(m, 1, k, j)
                        end do
                    end do
                end do
                if (nbd <= l1) then
                    is = -ido
                    do j=2, ip_rename
                        is = is+ido
                        idij = is
                        do i=3, ido, 2
                            idij = idij+2
                            do k=1, l1
                                do m=1, mp
                                    ch(m, i-1, k, j) = wa(idij-1)*c1(m, i-1, k, j)+wa(idij) &            !
                                        *c1(m, i, k, j)
                                    ch(m, i, k, j) = wa(idij-1)*c1(m, i, k, j)-wa(idij) &
                                        *c1(m, i-1, k, j)
                                end do
                            end do
                        end do
                    end do
                else
                    is = -ido
                    do j=2, ip_rename
                        is = is+ido
                        do k=1, l1
                            idij = is
                            do i=3, ido, 2
                                idij = idij+2
                                do m=1, mp
                                    ch(m, i-1, k, j) = wa(idij-1)*c1(m, i-1, k, j)+wa(idij) &            !
                                        *c1(m, i, k, j)
                                    ch(m, i, k, j) = wa(idij-1)*c1(m, i, k, j)-wa(idij) &
                                        *c1(m, i-1, k, j)
                                end do
                            end do
                        end do
                    end do
                end if
                if (nbd >= l1) then
                    do j=2, ipph
                        jc = ipp2-j
                        do k=1, l1
                            do i=3, ido, 2
                                do m=1, mp
                                    c1(m, i-1, k, j) = ch(m, i-1, k, j)+ch(m, i-1, k, jc)
                                    c1(m, i-1, k, jc) = ch(m, i, k, j)-ch(m, i, k, jc)
                                    c1(m, i, k, j) = ch(m, i, k, j)+ch(m, i, k, jc)
                                    c1(m, i, k, jc) = ch(m, i-1, k, jc)-ch(m, i-1, k, j)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2, ipph
                        jc = ipp2-j
                        do i=3, ido, 2
                            do k=1, l1
                                do m=1, mp
                                    c1(m, i-1, k, j) = ch(m, i-1, k, j)+ch(m, i-1, k, jc)
                                    c1(m, i-1, k, jc) = ch(m, i, k, j)-ch(m, i, k, jc)
                                    c1(m, i, k, j) = ch(m, i, k, j)+ch(m, i, k, jc)
                                    c1(m, i, k, jc) = ch(m, i-1, k, jc)-ch(m, i-1, k, j)
                                end do
                            end do
                        end do
                    end do
                end if
            else
                do ik=1, idl1
                    do m=1, mp
                        c2(m, ik, 1) = ch2(m, ik, 1)
                    end do
                end do
            end if

            do j=2, ipph
                jc = ipp2-j
                do k=1, l1
                    do m=1, mp
                        c1(m, 1, k, j) = ch(m, 1, k, j)+ch(m, 1, k, jc)
                        c1(m, 1, k, jc) = ch(m, 1, k, jc)-ch(m, 1, k, j)
                    end do
                end do
            end do

            ar1 = 1.0_wp
            ai1 = 0.0_wp
            do l=2, ipph
                lc = ipp2-l
                ar1h = dcp*ar1-dsp*ai1
                ai1 = dcp*ai1+dsp*ar1
                ar1 = ar1h
                do ik=1, idl1
                    do m=1, mp
                        ch2(m, ik, l) = c2(m, ik, 1)+ar1*c2(m, ik, 2)
                        ch2(m, ik, lc) = ai1*c2(m, ik, ip_rename)
                    end do
                end do
                dc2 = ar1
                ds2 = ai1
                ar2 = ar1
                ai2 = ai1
                do j=3, ipph
                    jc = ipp2-j
                    ar2h = dc2*ar2-ds2*ai2
                    ai2 = dc2*ai2+ds2*ar2
                    ar2 = ar2h
                    do ik=1, idl1
                        do m=1, mp
                            ch2(m, ik, l) = ch2(m, ik, l)+ar2*c2(m, ik, j)
                            ch2(m, ik, lc) = ch2(m, ik, lc)+ai2*c2(m, ik, jc)
                        end do
                    end do
                end do
            end do

            do j=2, ipph
                do ik=1, idl1
                    do m=1, mp
                        ch2(m, ik, 1) = ch2(m, ik, 1)+c2(m, ik, j)
                    end do
                end do
            end do

            if (ido >= l1) then
                do k=1, l1
                    do i=1, ido
                        do m=1, mp
                            cc(m, i, 1, k) = ch(m, i, k, 1)
                        end do
                    end do
                end do
            else
                do i=1, ido
                    do k=1, l1
                        do m=1, mp
                            cc(m, i, 1, k) = ch(m, i, k, 1)
                        end do
                    end do
                end do
            end if

            do j=2, ipph
                jc = ipp2-j
                j2 = j+j
                do k=1, l1
                    do m=1, mp
                        cc(m, ido, j2-2, k) = ch(m, 1, k, j)
                        cc(m, 1, j2-1, k) = ch(m, 1, k, jc)
                    end do
                end do
            end do

            if (ido == 1) then
                return
            end if

            if (nbd >= l1) then
                do j=2, ipph
                    jc = ipp2-j
                    j2 = j+j
                    do k=1, l1
                        do i=3, ido, 2
                            ic = idp2-i
                            do m=1, mp
                                cc(m, i-1, j2-1, k) = ch(m, i-1, k, j)+ch(m, i-1, k, jc)                !
                                cc(m, ic-1, j2-2, k) = ch(m, i-1, k, j)-ch(m, i-1, k, jc)               !
                                cc(m, i, j2-1, k) = ch(m, i, k, j)+ch(m, i, k, jc)
                                cc(m, ic, j2-2, k) = ch(m, i, k, jc)-ch(m, i, k, j)
                            end do
                        end do
                    end do
                end do
            else
                do j=2, ipph
                    jc = ipp2-j
                    j2 = j+j
                    do i=3, ido, 2
                        ic = idp2-i
                        do k=1, l1
                            do m=1, mp
                                cc(m, i-1, j2-1, k) = ch(m, i-1, k, j)+ch(m, i-1, k, jc)                !
                                cc(m, ic-1, j2-2, k) = ch(m, i-1, k, j)-ch(m, i-1, k, jc)               !
                                cc(m, i, j2-1, k) = ch(m, i, k, j)+ch(m, i, k, jc)
                                cc(m, ic, j2-2, k) = ch(m, i, k, jc)-ch(m, i, k, j)
                            end do
                        end do
                    end do
                end do
            end if
        end associate
    end associate

end subroutine hradfg



subroutine hrfftb(m, n, r, mdimr, whrfft, work)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: m
    integer (ip), intent (in)     :: n
    real (wp),    intent (in out) :: r(mdimr, n)
    integer (ip), intent (in)     :: mdimr
    real (wp),    intent (in out) :: whrfft(n+15)
    real (wp),    intent (in out) :: work(1)
    !----------------------------------------------------------------------

    if (n == 1) then
        return
    end if

    call hrftb1(m, n, r, mdimr, work, whrfft, whrfft(n+1))

end subroutine hrfftb



subroutine hrftb1(m, n, c, mdimc, ch, wa, fac)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: m
    integer (ip), intent (in)     :: n
    real (wp),    intent (in out) :: c(mdimc, n)
    integer (ip), intent (in)     :: mdimc
    real (wp),    intent (in out) :: ch(m, n)
    real (wp),    intent (in out) :: wa(n)
    real (wp),    intent (in out) :: fac(15)
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip) :: i, j, k1, l1, l2, na
    integer (ip) :: nf, ip_rename, iw, ix2, ix3, ix4, ido, idl1
    !----------------------------------------------------------------------

    nf = int(fac(2), kind=ip)
    na = 0
    l1 = 1
    iw = 1
    do k1=1, nf
        ip_rename = int(fac(k1+2), kind=ip)
        l2 = ip_rename*l1
        ido = n/l2
        idl1 = ido*l1

        if (ip_rename /= 4) then
            go to 103
        end if

        ix2 = iw+ido
        ix3 = ix2+ido

        if (na /= 0) then
            go to 101
        end if

        call  hradb4(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2), wa(ix3))
        go to 102
101     call hradb4(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2), wa(ix3))
102     na = 1-na
        go to 115

103     if (ip_rename /= 2) then
            go to 106
        end if

        if (na /= 0) then
            go to 104
        end if

        call  hradb2(m, ido, l1, c, mdimc, ch, m, wa(iw))
        go to 105
104     call hradb2(m, ido, l1, ch, m, c, mdimc, wa(iw))
105     na = 1-na
        go to 115

106     if (ip_rename /= 3) then
            go to 109
        end if

        ix2 = iw+ido
        if (na /= 0) then
            go to 107
        end if

        call  hradb3(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2))
        go to 108
107     call hradb3(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2))
108     na = 1-na
        go to 115
109     if (ip_rename /= 5) then
            go to 112
        end if

        ix2 = iw+ido
        ix3 = ix2+ido
        ix4 = ix3+ido

        if (na /= 0) then
            go to 110
        end if

        call hradb5(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2), wa(ix3), wa(ix4))
        go to 111
110     call hradb5(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2), wa(ix3), wa(ix4))
111     na = 1-na
        go to 115

112     if (na /= 0) then
            go to 113
        end if

        call  hradbg(m, ido, ip_rename, l1, idl1, c, c, c, mdimc, ch, ch, m, wa(iw))
        go to 114
113     call hradbg(m, ido, ip_rename, l1, idl1, ch, ch, ch, m, c, c, mdimc, wa(iw))
114     if (ido == 1) then
            na = 1-na
        end if

115     l1 = l2
        iw = iw+(ip_rename-1)*ido
    end do

    if (na == 0) then
        return
    end if

    c(1:m, 1:n) = ch

end subroutine hrftb1



subroutine hradbg(mp, ido, ip_rename, l1, idl1, cc, c1, c2, mdimcc, &
    ch, ch2, mdimch, wa)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: ip_rename
    integer (ip), intent (in)     :: l1
    integer (ip), intent (in)     :: idl1
    real (wp),    intent (in out) :: cc(mdimcc, ido, ip_rename, l1)
    real (wp),    intent (in out) :: c1(mdimcc, ido, l1, ip_rename)
    real (wp),    intent (in out) :: c2(mdimcc, idl1, ip_rename)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: ch(mdimch, ido, l1, ip_rename)
    real (wp),    intent (in out) :: ch2(mdimch, idl1, ip_rename)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: wa(ido)
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip)         :: i, j, k, l, m, j2, ic, jc, lc, ik, is, nbd
    integer (ip)         :: idp2, ipp2, idij, ipph
    real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
    real (wp)            :: dc2, ai1, ai2, ar1, ar2, ds2
    real (wp)            :: dcp, arg, dsp, ar1h, ar2h
    !----------------------------------------------------------------------

    arg = TWO_PI/ip_rename
    dcp = cos(arg)
    dsp = sin(arg)
    idp2 = ido+2
    nbd = (ido-1)/2
    ipp2 = ip_rename+2
    ipph = (ip_rename+1)/2

    if (ido < l1) then
        go to 103
    end if

    do k=1, l1
        do i=1, ido
            do m=1, mp
                ch(m, i, k, 1) = cc(m, i, 1, k)
            end do
        end do
    end do
    go to 106
    103 do i=1, ido
        do k=1, l1
            do m=1, mp
                ch(m, i, k, 1) = cc(m, i, 1, k)
            end do
        end do
    end do
    106 do j=2, ipph
        jc = ipp2-j
        j2 = j+j
        do k=1, l1
            do m=1, mp
                ch(m, 1, k, j) = cc(m, ido, j2-2, k)+cc(m, ido, j2-2, k)
                ch(m, 1, k, jc) = cc(m, 1, j2-1, k)+cc(m, 1, j2-1, k)
            end do
        end do
    end do

    if (ido == 1) then
        go to 116
    end if

    if (nbd < l1) then
        go to 112
    end if

    do j=2, ipph
        jc = ipp2-j
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                do m=1, mp
                    ch(m, i-1, k, j) = cc(m, i-1, 2*j-1, k)+cc(m, ic-1, 2*j-2, k)
                    ch(m, i-1, k, jc) = cc(m, i-1, 2*j-1, k)-cc(m, ic-1, 2*j-2, k)
                    ch(m, i, k, j) = cc(m, i, 2*j-1, k)-cc(m, ic, 2*j-2, k)
                    ch(m, i, k, jc) = cc(m, i, 2*j-1, k)+cc(m, ic, 2*j-2, k)
                end do
            end do
        end do
    end do
    go to 116
    112 do j=2, ipph
        jc = ipp2-j
        do i=3, ido, 2
            ic = idp2-i
            do k=1, l1
                do m=1, mp
                    ch(m, i-1, k, j) = cc(m, i-1, 2*j-1, k)+cc(m, ic-1, 2*j-2, k)
                    ch(m, i-1, k, jc) = cc(m, i-1, 2*j-1, k)-cc(m, ic-1, 2*j-2, k)
                    ch(m, i, k, j) = cc(m, i, 2*j-1, k)-cc(m, ic, 2*j-2, k)
                    ch(m, i, k, jc) = cc(m, i, 2*j-1, k)+cc(m, ic, 2*j-2, k)
                end do
            end do
        end do
    end do
116 ar1 = 1.0_wp
    ai1 = 0.0_wp
    do l=2, ipph
        lc = ipp2-l
        ar1h = dcp*ar1-dsp*ai1
        ai1 = dcp*ai1+dsp*ar1
        ar1 = ar1h
        do ik=1, idl1
            do m=1, mp
                c2(m, ik, l) = ch2(m, ik, 1)+ar1*ch2(m, ik, 2)
                c2(m, ik, lc) = ai1*ch2(m, ik, ip_rename)
            end do
        end do
        dc2 = ar1
        ds2 = ai1
        ar2 = ar1
        ai2 = ai1
        do j=3, ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do ik=1, idl1
                do m=1, mp
                    c2(m, ik, l) = c2(m, ik, l)+ar2*ch2(m, ik, j)
                    c2(m, ik, lc) = c2(m, ik, lc)+ai2*ch2(m, ik, jc)
                end do
            end do
        end do
    end do
    do j=2, ipph
        do ik=1, idl1
            do m=1, mp
                ch2(m, ik, 1) = ch2(m, ik, 1)+ch2(m, ik, j)
            end do
        end do
    end do
    do j=2, ipph
        jc = ipp2-j
        do k=1, l1
            do m=1, mp
                ch(m, 1, k, j) = c1(m, 1, k, j)-c1(m, 1, k, jc)
                ch(m, 1, k, jc) = c1(m, 1, k, j)+c1(m, 1, k, jc)
            end do
        end do
    end do

    if (ido /= 1) then
        if (nbd >= l1) then
            do j=2, ipph
                jc = ipp2-j
                do k=1, l1
                    do i=3, ido, 2
                        do m=1, mp
                            ch(m, i-1, k, j) = c1(m, i-1, k, j)-c1(m, i, k, jc)
                            ch(m, i-1, k, jc) = c1(m, i-1, k, j)+c1(m, i, k, jc)
                            ch(m, i, k, j) = c1(m, i, k, j)+c1(m, i-1, k, jc)
                            ch(m, i, k, jc) = c1(m, i, k, j)-c1(m, i-1, k, jc)
                        end do
                    end do
                end do
            end do
        else
            do j=2, ipph
                jc = ipp2-j
                do i=3, ido, 2
                    do k=1, l1
                        do m=1, mp
                            ch(m, i-1, k, j) = c1(m, i-1, k, j)-c1(m, i, k, jc)
                            ch(m, i-1, k, jc) = c1(m, i-1, k, j)+c1(m, i, k, jc)
                            ch(m, i, k, j) = c1(m, i, k, j)+c1(m, i-1, k, jc)
                            ch(m, i, k, jc) = c1(m, i, k, j)-c1(m, i-1, k, jc)
                        end do
                    end do
                end do
            end do
        end if
    end if

    if (ido == 1) then
        return
    end if

    do ik=1, idl1
        do m=1, mp
            c2(m, ik, 1) = ch2(m, ik, 1)
        end do
    end do
    do j=2, ip_rename
        do k=1, l1
            do m=1, mp
                c1(m, 1, k, j) = ch(m, 1, k, j)
            end do
        end do
    end do
    if (nbd <= l1) then
        is = -ido
        do j=2, ip_rename
            is = is+ido
            idij = is
            do i=3, ido, 2
                idij = idij+2
                do k=1, l1
                    do m=1, mp
                        c1(m, i-1, k, j) = wa(idij-1)*ch(m, i-1, k, j)-wa(idij)* &           !
                            ch(m, i, k, j)
                        c1(m, i, k, j) = wa(idij-1)*ch(m, i, k, j)+wa(idij)* &               !
                            ch(m, i-1, k, j)
                    end do
                end do
            end do
        end do
    else
        is = -ido
        do j=2, ip_rename
            is = is+ido
            do k=1, l1
                idij = is
                do i=3, ido, 2
                    idij = idij+2
                    do m=1, mp
                        c1(m, i-1, k, j) = wa(idij-1)*ch(m, i-1, k, j)-wa(idij)* &           !
                            ch(m, i, k, j)
                        c1(m, i, k, j) = wa(idij-1)*ch(m, i, k, j)+wa(idij)* &               !
                            ch(m, i-1, k, j)
                    end do
                end do
            end do
        end do
    end if

end subroutine hradbg



subroutine hradb4(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2, wa3)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: l1
    real (wp),    intent (in out) :: cc(mdimcc, ido, 4, l1)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: ch(mdimch, ido, l1, 4)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: wa1(ido)
    real (wp),    intent (in out) :: wa2(ido)
    real (wp),    intent (in out) :: wa3(ido)
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip)         :: i, k, m, ic, idp2
    real (wp), parameter :: SQRT2 = sqrt(2.0_wp)
    !----------------------------------------------------------------------

    do k=1, l1
        do m=1, mp
            ch(m, 1, k, 3) = (cc(m, 1, 1, k)+cc(m, ido, 4, k)) &
                -(cc(m, ido, 2, k)+cc(m, ido, 2, k))
            ch(m, 1, k, 1) = (cc(m, 1, 1, k)+cc(m, ido, 4, k)) &
                +(cc(m, ido, 2, k)+cc(m, ido, 2, k))
            ch(m, 1, k, 4) = (cc(m, 1, 1, k)-cc(m, ido, 4, k)) &
                +(cc(m, 1, 3, k)+cc(m, 1, 3, k))
            ch(m, 1, k, 2) = (cc(m, 1, 1, k)-cc(m, ido, 4, k)) &
                -(cc(m, 1, 3, k)+cc(m, 1, 3, k))
        end do
    end do

    if (ido-2 < 0) then
        return
    else if (ido-2 /= 0) then
        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                do m=1, mp
                    ch(m, i-1, k, 1) = &
                        (cc(m, i-1, 1, k)+cc(m, ic-1, 4, k)) &
                        +(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k))

                    ch(m, i, k, 1) = &
                        (cc(m, i, 1, k)-cc(m, ic, 4, k)) &
                        +(cc(m, i, 3, k)-cc(m, ic, 2, k))

                    ch(m, i-1, k, 2) = &
                        wa1(i-2)*((cc(m, i-1, 1, k)-cc(m, ic-1, 4, k)) &          !
                        -(cc(m, i, 3, k)+cc(m, ic, 2, k)))-wa1(i-1) &
                        *((cc(m, i, 1, k)+cc(m, ic, 4, k))+(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k)))      !

                    ch(m, i, k, 2)= &
                        wa1(i-2)*((cc(m, i, 1, k)+cc(m, ic, 4, k)) &
                        +(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k)))+wa1(i-1) &
                        *((cc(m, i-1, 1, k)-cc(m, ic-1, 4, k))-(cc(m, i, 3, k)+cc(m, ic, 2, k)))      !

                    ch(m, i-1, k, 3) = &
                        wa2(i-2)*((cc(m, i-1, 1, k)+cc(m, ic-1, 4, k)) &          !
                        -(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)))-wa2(i-1) &
                        *((cc(m, i, 1, k)-cc(m, ic, 4, k))-(cc(m, i, 3, k)-cc(m, ic, 2, k)))          !

                    ch(m, i, k, 3) = &
                        wa2(i-2)*((cc(m, i, 1, k)-cc(m, ic, 4, k)) &
                        -(cc(m, i, 3, k)-cc(m, ic, 2, k)))+wa2(i-1) &
                        *((cc(m, i-1, 1, k)+cc(m, ic-1, 4, k))-(cc(m, i-1, 3, k) &
                        +cc(m, ic-1, 2, k)))

                    ch(m, i-1, k, 4) = &
                        wa3(i-2)*((cc(m, i-1, 1, k)-cc(m, ic-1, 4, k)) &          !
                        +(cc(m, i, 3, k)+cc(m, ic, 2, k)))-wa3(i-1) &
                        *((cc(m, i, 1, k)+cc(m, ic, 4, k))-(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k)))       !

                    ch(m, i, k, 4) = &
                        wa3(i-2)*((cc(m, i, 1, k)+cc(m, ic, 4, k)) &
                        -(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k)))+wa3(i-1) &
                        *((cc(m, i-1, 1, k)-cc(m, ic-1, 4, k))+(cc(m, i, 3, k)+cc(m, ic, 2, k)))      !
                end do
            end do
        end do

        if(mod(ido, 2) == 1) then
            return
        end if
    end if

    do k=1, l1
        do m=1, mp
            ch(m, ido, k, 1) = &
                (cc(m, ido, 1, k)+cc(m, ido, 3, k)) &
                +(cc(m, ido, 1, k)+cc(m, ido, 3, k))

            ch(m, ido, k, 2) = &
                SQRT2*((cc(m, ido, 1, k)-cc(m, ido, 3, k)) &               !
                -(cc(m, 1, 2, k)+cc(m, 1, 4, k)))

            ch(m, ido, k, 3) = &
                (cc(m, 1, 4, k)-cc(m, 1, 2, k)) &
                +(cc(m, 1, 4, k)-cc(m, 1, 2, k))

            ch(m, ido, k, 4) = &
                -SQRT2*((cc(m, ido, 1, k)-cc(m, ido, 3, k)) &              !
                +(cc(m, 1, 2, k)+cc(m, 1, 4, k)))
        end do
    end do

end subroutine hradb4



subroutine hradb2(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: l1
    real (wp),    intent (in out) :: cc(mdimcc, ido, 2, l1)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: ch(mdimch, ido, l1, 2)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: wa1(ido)
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip) :: i, k, m, ic, idp2
    !----------------------------------------------------------------------

    do k=1, l1
        do m=1, mp
            ch(m, 1, k, 1) = cc(m, 1, 1, k)+cc(m, ido, 2, k)
            ch(m, 1, k, 2) = cc(m, 1, 1, k)-cc(m, ido, 2, k)
        end do
    end do

    if (ido-2 < 0) then
        return
    else if (ido-2 /= 0) then
        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                do m=1, mp
                    ch(m, i-1, k, 1) = &
                        cc(m, i-1, 1, k)+cc(m, ic-1, 2, k)

                    ch(m, i, k, 1) = &
                        cc(m, i, 1, k)-cc(m, ic, 2, k)

                    ch(m, i-1, k, 2) = &
                        wa1(i-2)*(cc(m, i-1, 1, k)-cc(m, ic-1, 2, k)) &         !
                        -wa1(i-1)*(cc(m, i, 1, k)+cc(m, ic, 2, k))

                    ch(m, i, k, 2) = &
                        wa1(i-2)*(cc(m, i, 1, k)+cc(m, ic, 2, k))+wa1(i-1) &      !
                        *(cc(m, i-1, 1, k)-cc(m, ic-1, 2, k))
                end do
            end do
        end do
        if(mod(ido, 2) == 1) then
            return
        end if
    end if

    do k=1, l1
        do m=1, mp
            ch(m, ido, k, 1) = cc(m, ido, 1, k)+cc(m, ido, 1, k)
            ch(m, ido, k, 2) = -(cc(m, 1, 2, k)+cc(m, 1, 2, k))
        end do
    end do

end subroutine hradb2



subroutine hradb3(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: l1
    real (wp),    intent (in out) :: cc(mdimcc, ido, 3, l1)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: ch(mdimch, ido, l1, 3)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: wa1(ido)
    real (wp),    intent (in out) :: wa2(ido)
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip)         :: i, k, m, ic, idp2
    real (wp), parameter :: ARG = 2.0_wp * acos(-1.0_wp)/3
    real (wp), parameter :: TAUR = cos(ARG)
    real (wp), parameter :: TAUI = sin(ARG)
    !----------------------------------------------------------------------


    do k=1, l1
        do m=1, mp
            ch(m, 1, k, 1) = &
                cc(m, 1, 1, k)+2.0_wp * cc(m, ido, 2, k)

            ch(m, 1, k, 2) = &
                cc(m, 1, 1, k)+(2.0_wp * TAUR)*cc(m, ido, 2, k) &
                -(2.0_wp *TAUI)*cc(m, 1, 3, k)

            ch(m, 1, k, 3) = &
                cc(m, 1, 1, k)+(2.0_wp * TAUR)*cc(m, ido, 2, k) &
                +2.0_wp *TAUI*cc(m, 1, 3, k)
        end do
    end do

    if (ido == 1) then
        return
    end if

    idp2 = ido+2
    do k=1, l1
        do i=3, ido, 2
            ic = idp2-i
            do m=1, mp

                ch(m, i-1, k, 1) = &
                    cc(m, i-1, 1, k)+(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k))      !

                ch(m, i, k, 1) = &
                    cc(m, i, 1, k)+(cc(m, i, 3, k)-cc(m, ic, 2, k))              !

                ch(m, i-1, k, 2) = &
                    wa1(i-2)* &
                    ((cc(m, i-1, 1, k)+TAUR*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)))- &
                    (TAUI*(cc(m, i, 3, k)+cc(m, ic, 2, k)))) &
                    -wa1(i-1)* &
                    ((cc(m, i, 1, k)+TAUR*(cc(m, i, 3, k)-cc(m, ic, 2, k)))+ &
                    (TAUI*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k))))

                ch(m, i, k, 2) = &
                    wa1(i-2)* &
                    ((cc(m, i, 1, k)+TAUR*(cc(m, i, 3, k)-cc(m, ic, 2, k)))+ &
                    (TAUI*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k)))) &
                    +wa1(i-1)* &
                    ((cc(m, i-1, 1, k)+TAUR*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)))- &
                    (TAUI*(cc(m, i, 3, k)+cc(m, ic, 2, k))))

                ch(m, i-1, k, 3) = &
                    wa2(i-2)* &
                    ((cc(m, i-1, 1, k)+TAUR*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)))+ &
                    (TAUI*(cc(m, i, 3, k)+cc(m, ic, 2, k)))) &
                    -wa2(i-1)* &
                    ((cc(m, i, 1, k)+TAUR*(cc(m, i, 3, k)-cc(m, ic, 2, k)))- &
                    (TAUI*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k))))

                ch(m, i, k, 3) = &
                    wa2(i-2)* &
                    ((cc(m, i, 1, k)+TAUR*(cc(m, i, 3, k)-cc(m, ic, 2, k)))- &
                    (TAUI*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k)))) &
                    +wa2(i-1)* &
                    ((cc(m, i-1, 1, k)+TAUR*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)))+ &
                    (TAUI*(cc(m, i, 3, k)+cc(m, ic, 2, k))))
            end do
        end do
    end do

end subroutine hradb3



subroutine hradb5(mp, ido, l1, cc, mdimcc, ch, mdimch, &
    wa1, wa2, wa3, wa4)
    !
    ! Purpose:
    !
    ! A multiple fft package for spherepack
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)     :: mp
    integer (ip), intent (in)     :: ido
    integer (ip), intent (in)     :: l1
    real (wp),    intent (in out) :: ch(mdimch, ido, l1, 5)
    integer (ip), intent (in)     :: mdimch
    real (wp),    intent (in out) :: cc(mdimcc, ido, 5, l1)
    integer (ip), intent (in)     :: mdimcc
    real (wp),    intent (in out) :: wa1(ido)
    real (wp),    intent (in out) :: wa2(ido)
    real (wp),    intent (in out) :: wa3(ido)
    real (wp),    intent (in out) :: wa4(ido)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip)         :: i, k, m, ic, idp2
    real (wp), parameter :: ARG = 2.0_wp * acos(-1.0_wp)/5
    real (wp), parameter :: TR11=cos(ARG)
    real (wp), parameter :: TI11=sin(ARG)
    real (wp), parameter :: TR12=cos(2.0_wp *ARG)
    real (wp), parameter :: TI12=sin(2.0_wp *ARG)
    !----------------------------------------------------------------------

    do k=1, l1
        do m=1, mp
            ch(m, 1, k, 1) = &
                cc(m, 1, 1, k)+2.0_wp *cc(m, ido, 2, k)+2.0_wp *cc(m, ido, 4, k)          !

            ch(m, 1, k, 2) = &
                (cc(m, 1, 1, k)+TR11*2.0_wp *cc(m, ido, 2, k) &
                +TR12*2.0_wp *cc(m, ido, 4, k))-(TI11*2.0_wp *cc(m, 1, 3, k) &
                +TI12*2.0_wp *cc(m, 1, 5, k))

            ch(m, 1, k, 3) = &
                (cc(m, 1, 1, k)+TR12*2.0_wp *cc(m, ido, 2, k) &
                +TR11*2.0_wp *cc(m, ido, 4, k))-(TI12*2.0_wp *cc(m, 1, 3, k) &
                -TI11*2.0_wp *cc(m, 1, 5, k))

            ch(m, 1, k, 4) = &
                (cc(m, 1, 1, k)+TR12*2.0_wp *cc(m, ido, 2, k) &
                +TR11*2.0_wp *cc(m, ido, 4, k))+(TI12*2.0_wp *cc(m, 1, 3, k) &
                -TI11*2.0_wp *cc(m, 1, 5, k))

            ch(m, 1, k, 5) = &
                (cc(m, 1, 1, k)+TR11*2.0_wp *cc(m, ido, 2, k) &
                +TR12*2.0_wp *cc(m, ido, 4, k))+(TI11*2.0_wp *cc(m, 1, 3, k) &
                +TI12*2.0_wp *cc(m, 1, 5, k))
        end do
    end do

    if (ido == 1) then
        return
    end if

    idp2 = ido+2
    do k=1, l1
        do i=3, ido, 2
            ic = idp2-i
            do m=1, mp

                ch(m, i-1, k, 1) = &
                    cc(m, i-1, 1, k)+(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)) &    !
                    +(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k))

                ch(m, i, k, 1) = &
                    cc(m, i, 1, k)+(cc(m, i, 3, k)-cc(m, ic, 2, k)) &            !
                    +(cc(m, i, 5, k)-cc(m, ic, 4, k))

                ch(m, i-1, k, 2) = &
                    wa1(i-2)*((cc(m, i-1, 1, k)+TR11* &
                    (cc(m, i-1, 3, k)+cc(m, ic-1, 2, k))+TR12 &
                    *(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k)))-(TI11*(cc(m, i, 3, k) &             !
                    +cc(m, ic, 2, k))+TI12*(cc(m, i, 5, k)+cc(m, ic, 4, k)))) &
                    -wa1(i-1)*((cc(m, i, 1, k)+TR11*(cc(m, i, 3, k)-cc(m, ic, 2, k)) &         !
                    +TR12*(cc(m, i, 5, k)-cc(m, ic, 4, k)))+(TI11*(cc(m, i-1, 3, k) &          !
                    -cc(m, ic-1, 2, k))+TI12*(cc(m, i-1, 5, k)-cc(m, ic-1, 4, k))))            !

                ch(m, i, k, 2) = &
                    wa1(i-2)*((cc(m, i, 1, k)+TR11*(cc(m, i, 3, k) &          !
                    -cc(m, ic, 2, k))+TR12*(cc(m, i, 5, k)-cc(m, ic, 4, k))) &
                    +(TI11*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k))+TI12 &
                    *(cc(m, i-1, 5, k)-cc(m, ic-1, 4, k))))+wa1(i-1) &
                    *((cc(m, i-1, 1, k)+TR11*(cc(m, i-1, 3, k) &
                    +cc(m, ic-1, 2, k))+TR12*(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k))) &           !
                    -(TI11*(cc(m, i, 3, k)+cc(m, ic, 2, k))+TI12 &
                    *(cc(m, i, 5, k)+cc(m, ic, 4, k))))

                ch(m, i-1, k, 3) = &
                    wa2(i-2) &
                    *((cc(m, i-1, 1, k)+TR12*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)) &            !
                    +TR11*(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k)))-(TI12*(cc(m, i, 3, k) &        !
                    +cc(m, ic, 2, k))-TI11*(cc(m, i, 5, k)+cc(m, ic, 4, k)))) &
                    -wa2(i-1) &
                    *((cc(m, i, 1, k)+TR12*(cc(m, i, 3, k)- &
                    cc(m, ic, 2, k))+TR11*(cc(m, i, 5, k)-cc(m, ic, 4, k))) &
                    +(TI12*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k))-TI11 &
                    *(cc(m, i-1, 5, k)-cc(m, ic-1, 4, k))))

                ch(m, i, k, 3) = &
                    wa2(i-2) &
                    *((cc(m, i, 1, k)+TR12*(cc(m, i, 3, k)- &
                    cc(m, ic, 2, k))+TR11*(cc(m, i, 5, k)-cc(m, ic, 4, k))) &
                    +(TI12*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k))-TI11 &
                    *(cc(m, i-1, 5, k)-cc(m, ic-1, 4, k)))) &
                    +wa2(i-1) &
                    *((cc(m, i-1, 1, k)+TR12*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)) &            !
                    +TR11*(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k)))-(TI12*(cc(m, i, 3, k) &        !
                    +cc(m, ic, 2, k))-TI11*(cc(m, i, 5, k)+cc(m, ic, 4, k))))

                ch(m, i-1, k, 4) = &
                    wa3(i-2) &
                    *((cc(m, i-1, 1, k)+TR12*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)) &            !
                    +TR11*(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k)))+(TI12*(cc(m, i, 3, k) &        !
                    +cc(m, ic, 2, k))-TI11*(cc(m, i, 5, k)+cc(m, ic, 4, k)))) &
                    -wa3(i-1) &
                    *((cc(m, i, 1, k)+TR12*(cc(m, i, 3, k)- &
                    cc(m, ic, 2, k))+TR11*(cc(m, i, 5, k)-cc(m, ic, 4, k))) &
                    -(TI12*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k))-TI11 &
                    *(cc(m, i-1, 5, k)-cc(m, ic-1, 4, k))))

                ch(m, i, k, 4) = &
                    wa3(i-2) &
                    *((cc(m, i, 1, k)+TR12*(cc(m, i, 3, k)- &
                    cc(m, ic, 2, k))+TR11*(cc(m, i, 5, k)-cc(m, ic, 4, k))) &
                    -(TI12*(cc(m, i-1, 3, k)-cc(m, ic-1, 2, k))-TI11 &
                    *(cc(m, i-1, 5, k)-cc(m, ic-1, 4, k)))) &
                    +wa3(i-1) &
                    *((cc(m, i-1, 1, k)+TR12*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)) &            !
                    +TR11*(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k)))+(TI12*(cc(m, i, 3, k) &        !
                    +cc(m, ic, 2, k))-TI11*(cc(m, i, 5, k)+cc(m, ic, 4, k))))

                ch(m, i-1, k, 5) = &
                    wa4(i-2) &
                    *((cc(m, i-1, 1, k)+TR11*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)) &            !
                    +TR12*(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k)))+(TI11*(cc(m, i, 3, k) &        !
                    +cc(m, ic, 2, k))+TI12*(cc(m, i, 5, k)+cc(m, ic, 4, k)))) &
                    -wa4(i-1) &
                    *((cc(m, i, 1, k)+TR11*(cc(m, i, 3, k)-cc(m, ic, 2, k)) &
                    +TR12*(cc(m, i, 5, k)-cc(m, ic, 4, k)))-(TI11*(cc(m, i-1, 3, k) &          !
                    -cc(m, ic-1, 2, k))+TI12*(cc(m, i-1, 5, k)-cc(m, ic-1, 4, k))))            !

                ch(m, i, k, 5) = &
                    wa4(i-2) &
                    *((cc(m, i, 1, k)+TR11*(cc(m, i, 3, k)-cc(m, ic, 2, k)) &
                    +TR12*(cc(m, i, 5, k)-cc(m, ic, 4, k)))-(TI11*(cc(m, i-1, 3, k) &          !
                    -cc(m, ic-1, 2, k))+TI12*(cc(m, i-1, 5, k)-cc(m, ic-1, 4, k)))) &          !
                    +wa4(i-1) &
                    *((cc(m, i-1, 1, k)+TR11*(cc(m, i-1, 3, k)+cc(m, ic-1, 2, k)) &            !
                    +TR12*(cc(m, i-1, 5, k)+cc(m, ic-1, 4, k)))+(TI11*(cc(m, i, 3, k) &        !
                    +cc(m, ic, 2, k))+TI12*(cc(m, i, 5, k)+cc(m, ic, 4, k))))
            end do
        end do
    end do

end subroutine hradb5

