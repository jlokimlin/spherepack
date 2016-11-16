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
!     *       A Package of Fortran Subroutines and Programs           *
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
! ... file type_HFFTpack.f90
!
!     This file contains a multiple fft package for spherepack. It
!     includes code and documentation for performing fast Fourier
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
module type_HFFTpack

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI, &
        TWO_PI

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: HFFTpack
    public :: hrffti, hrfftf, hrfftb


    ! Declare derived data type
    type, public :: HFFTpack
    contains
        !-------------------------------------------------------
        ! Type-bound procedures
        !-------------------------------------------------------
        procedure, nopass :: initialize => hrffti
        procedure, nopass :: forward => hrfftf
        procedure, nopass :: backward => hrfftb
        !-------------------------------------------------------
    end type HFFTpack



contains



    subroutine hrffti(n, wsave)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: wsave(n+15)
        !----------------------------------------------------------------------

        if (n == 1) return

        call hrfti1(n, wsave(1), wsave(n+1))

    end subroutine hrffti



    subroutine hrfti1(n, wa, fac)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)  :: n
        real(wp),    intent(out) :: wa(n)
        real(wp),    intent(out) :: fac(15)
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer(ip)            :: i, ib, ido, ii, iip, ipm, is
        integer(ip)            :: j, k1, l1, l2, ld
        integer(ip)            :: nf, nfm1, nl, nq, nr, ntry
        integer(ip), parameter :: NTRYH(*) = [4, 2, 3, 5]
        real(wp)               :: arg,  argh, argld, fi
        !--------------------------------------------------------------

        ! Initialize
        ntry = 0
        nl = n
        nf = 0
        j = 0

        factorize_loop: do
            ! Increment j
            j = j+1

            ! Choose ntry
            if (j <= 4) then
                ntry = NTRYH(j)
            else
                ntry = ntry+2
            end if

            inner_loop: do
                nq = nl/ntry
                nr = nl-ntry*nq
                if (nr < 0) then
                    cycle factorize_loop
                else if (nr == 0) then
                    nf = nf+1
                    fac(nf+2) = ntry
                    nl = nq

                    if (ntry == 2 .and. nf /= 1) then
                        do i=2,nf
                            ib = nf-i+2
                            fac(ib+2) = fac(ib+1)
                        end do
                        fac(3) = 2
                    end if

                    if (nl /= 1) then
                        cycle inner_loop
                    end if
                else
                    cycle factorize_loop
                end if
                exit inner_loop
            end do inner_loop
            exit factorize_loop
        end do factorize_loop

        fac(1) = n
        fac(2) = nf
        argh = TWO_PI/n
        is = 0
        nfm1 = nf-1
        l1 = 1

        if (nfm1 /= 0) then
            do k1=1,nfm1
                iip = int(fac(k1+2), kind=ip)
                ld = 0
                l2 = l1*iip
                ido = n/l2
                ipm = iip-1
                do j=1,ipm
                    ld = ld+l1
                    i = is
                    argld = real(ld, kind=wp) * argh
                    fi = 0.0_wp
                    do ii=3,ido,2
                        i = i+2
                        fi = fi + 1.0_wp
                        arg = fi*argld
                        wa(i-1) = cos(arg)
                        wa(i) = sin(arg)
                    end do
                    is = is+ido
                end do
                l1 = l2
            end do
        end if

    end subroutine hrfti1




    subroutine hrfftf(m, n, r, mdimr, whrfft, work)
        !
        ! Purpose:
        !
        ! A multiple fft package for spherepack
        !
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: r(mdimr, n)
        integer(ip), intent(in)     :: mdimr
        real(wp),    intent(in)     :: whrfft(n+15)
        real(wp),    intent(out)    :: work(*)
         !----------------------------------------------------------------------

        if (n == 1) return

        call hrftf1(m, n, r, mdimr, work, whrfft, whrfft(n+1))

    end subroutine hrfftf



    subroutine hrftf1(m, n, c, mdimc, ch, wa, fac)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: c(mdimc, n)
        integer(ip), intent(in)     :: mdimc
        real(wp),    intent(out)    :: ch(m, n)
        real(wp),    intent(in)     :: wa(n)
        real(wp),    intent(in)     :: fac(15)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: k1, l1, l2
        integer(ip) :: na, kh, nf, iip
        integer(ip) :: iw, ix2, ix3, ix4, ido, idl1
        !----------------------------------------------------------------------

        nf = int(fac(2), kind=ip)
        na = 1
        l2 = n
        iw = n

        do k1=1, nf
            kh = nf-k1
            iip = int(fac(kh+3), kind=ip)
            l1 = l2/iip
            ido = n/l2
            idl1 = ido*l1
            iw = iw-(iip-1)*ido
            na = 1-na

            select case (iip)
                case (2)
                    if (na == 0) then
                        call  hradf2(m, ido, l1, c, mdimc, ch, m, wa(iw))
                    else
                        call hradf2(m, ido, l1, ch, m, c, mdimc, wa(iw))
                    end if
                case (3)
                    ix2 = iw+ido
                    if (na == 0) then
                        call  hradf3(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2))
                    else
                        call hradf3(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2))
                    end if
                case(4)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    if (na == 0) then
                        call  hradf4(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2), wa(ix3))
                    else
                        call hradf4(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2), wa(ix3))
                    end if
                case (5)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    ix4 = ix3+ido
                    if (na == 0) then
                        call hradf5(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                    else
                        call hradf5(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                    end if
                case default
                    if (ido == 1) na = 1-na
                    if (na == 0) then
                        call  hradfg(m, ido, iip, l1, idl1, c, c, c, mdimc, ch, ch, m, wa(iw))
                        na = 1
                    else
                        call hradfg(m, ido, iip, l1, idl1, ch, ch, ch, m, c, c, mdimc, wa(iw))
                        na = 0
                    end if
            end select
            l2 = l1
        end do

        if (na /= 1) c(1:m, 1:n) = ch

    end subroutine hrftf1



    subroutine hradf2(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: ch(mdimch, ido, 2, l1)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(inout)  :: cc(mdimcc, ido, l1, 2)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(in)     :: wa1(ido)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: i, k, m, ic, idp2
        !----------------------------------------------------------------------

        ch(1:mp, 1, 1,:) = cc(1:mp, 1,:, 1)+cc(1:mp, 1,:, 2)
        ch(1:mp, ido, 2,:) = cc(1:mp, 1,:, 1)-cc(1:mp, 1,:, 2)


        if (ido < 2) then
            return
        else if (ido /= 2) then
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

            if (mod(ido, 2) == 1) return

        end if

        ch(1:mp, 1, 2, :) = -cc(1:mp, ido,:, 2)
        ch(1:mp, ido, 1, :) = cc(1:mp, ido,:, 1)


    end subroutine hradf2


    subroutine hradf3(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: ch(mdimch, ido, 3, l1)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(inout)  :: cc(mdimcc, ido, l1, 3)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(in)     :: wa1(ido)
        real(wp),    intent(in)     :: wa2(ido)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)         :: i, k, m, ic, idp2
        real(wp), parameter :: ARG=TWO_PI/3
        real(wp), parameter :: TAUR=cos(ARG)! -0.5_wp
        real(wp), parameter :: TAUI=sin(ARG)
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

        if (ido == 1) return

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


    subroutine hradf4(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2, wa3)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, l1, 4)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(inout)  :: ch(mdimch, ido, 4, l1)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa1(ido)
        real(wp),    intent(in)     :: wa2(ido)
        real(wp),    intent(in)     :: wa3(ido)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)         :: i, k, m, ic, idp2
        real(wp), parameter :: HALF_SQRT2 = sqrt(2.0_wp)/2
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

        if (ido < 2) then
            return
        else if (ido /= 2) then
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

            if (mod(ido,2) == 1) return

        end if

        do k=1, l1
            do m=1, mp

                ch(m, ido, 1, k) = &
                    (HALF_SQRT2*(cc(m, ido, k, 2)-cc(m, ido, k, 4)))+ &
                    cc(m, ido, k, 1)

                ch(m, ido, 3, k) = &
                    cc(m, ido, k, 1)-(HALF_SQRT2*(cc(m, ido, k, 2)- &
                    cc(m, ido, k, 4)))

                ch(m, 1, 2, k) = &
                    (-HALF_SQRT2*(cc(m, ido, k, 2)+cc(m, ido, k, 4)))- &
                    cc(m, ido, k, 3)

                ch(m, 1, 4, k) = &
                    (-HALF_SQRT2*(cc(m, ido, k, 2)+cc(m, ido, k, 4)))+ &
                    cc(m, ido, k, 3)
            end do
        end do

    end subroutine hradf4


    subroutine hradf5(mp, ido, l1, cc, mdimcc, ch, mdimch, &
        wa1, wa2, wa3, wa4)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: ch(mdimch, ido, 5, l1)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(inout)  :: cc(mdimcc, ido, l1, 5)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(in)     :: wa1(ido)
        real(wp),    intent(in)     :: wa2(ido)
        real(wp),    intent(in)     :: wa3(ido)
        real(wp),    intent(in)     :: wa4(ido)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)         :: i, k, m, ic, idp2
        real(wp), parameter :: ARG = TWO_PI/5
        real(wp), parameter :: TR11=cos(ARG)
        real(wp), parameter :: TI11=sin(ARG)
        real(wp), parameter :: TR12=cos(2.0_wp*ARG)
        real(wp), parameter :: TI12=sin(2.0_wp*ARG)
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

        if (ido == 1) return

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



    subroutine hradfg(mp, ido, iip, l1, idl1, cc, c1, c2, mdimcc, &
        ch, ch2, mdimch, wa)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: iip
        integer(ip), intent(in)     :: l1
        integer(ip), intent(in)     :: idl1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, iip, l1)
        real(wp),    intent(inout)  :: c1(mdimcc, ido, l1, iip)
        real(wp),    intent(inout)  :: c2(mdimcc, idl1, iip)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(inout)  :: ch(mdimch, ido, l1, iip)
        real(wp),    intent(inout)  :: ch2(mdimch, idl1, iip)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa(ido)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: i, j, k, l, j2, ic, jc, lc, ik, is, idij
        real(wp)    :: dc2, ai1, ai2, ar1, ar2, ds2
        real(wp)    :: ar1h, ar2h
        !----------------------------------------------------------------------

        associate( &
            ipph => (iip+1)/2, &
            ipp2 => iip+2, &
            idp2 => ido+2, &
            nbd => (ido-1)/2, &
            arg => TWO_PI/iip &
            )
            associate( &
                dcp => cos(arg), &
                dsp => sin(arg) &
                )

                if (ido /= 1) then
                    ch2(1: mp, 1: idl1, 1) = c2(1: mp, 1: idl1, 1)
                    ch(1: mp, 1, 1: l1, 2: iip) = c1(1: mp, 1, 1: l1, 2: iip)

                    if (nbd <= l1) then
                        is = -ido
                        do j=2, iip
                            is = is+ido
                            idij = is
                            do i=3, ido, 2
                                ch(1: mp, i-1, 1: l1, j) = wa(idij+1)*c1(1: mp, i-1, 1: l1, j)+wa(idij+2) &            !
                                    *c1(1: mp, i, 1: l1, j)
                                ch(1: mp, i, 1: l1, j) = wa(idij+1)*c1(1: mp, i, 1: l1, j)-wa(idij+2) &
                                    *c1(1: mp, i-1, 1: l1, j)
                            end do
                        end do
                    else
                        is = -ido
                        do j=2, iip
                            is = is+ido
                            do k=1, l1
                                idij = is
                                do i=3, ido, 2
                                    idij = idij+2
                                    ch(1: mp, i-1, k, j) = wa(idij-1)*c1(1: mp, i-1, k, j)+wa(idij) &            !
                                        *c1(1: mp, i, k, j)
                                    ch(1: mp, i, k, j) = wa(idij-1)*c1(1: mp, i, k, j)-wa(idij) &
                                        *c1(1: mp, i-1, k, j)
                                end do
                            end do
                        end do
                    end if
                    if (nbd >= l1) then
                        do j=2, ipph
                            jc = ipp2-j
                            do k=1, l1
                                do i=3, ido, 2
                                    c1(1: mp, i-1, k, j) = ch(1: mp, i-1, k, j)+ch(1: mp, i-1, k, jc)
                                    c1(1: mp, i-1, k, jc) = ch(1: mp, i, k, j)-ch(1: mp, i, k, jc)
                                    c1(1: mp, i, k, j) = ch(1: mp, i, k, j)+ch(1: mp, i, k, jc)
                                    c1(1: mp, i, k, jc) = ch(1: mp, i-1, k, jc)-ch(1: mp, i-1, k, j)
                                end do
                            end do
                        end do
                    else
                        do j=2, ipph
                            jc = ipp2-j
                            do i=3, ido, 2
                                c1(1: mp, i-1, 1: l1, j) = ch(1: mp, i-1, 1: l1, j)+ch(1: mp, i-1, 1: l1, jc)
                                c1(1: mp, i-1, 1: l1, jc) = ch(1: mp, i, 1: l1, j)-ch(1: mp, i, 1: l1, jc)
                                c1(1: mp, i, 1: l1, j) = ch(1: mp, i, 1: l1, j)+ch(1: mp, i, 1: l1, jc)
                                c1(1: mp, i, 1: l1, jc) = ch(1: mp, i-1, 1: l1, jc)-ch(1: mp, i-1, 1: l1, j)
                            end do
                        end do
                    end if
                else
                    c2(1: mp, 1: idl1, 1) = ch2(1: mp, 1: idl1, 1)
                end if

                do j=2, ipph
                    jc = ipp2-j
                    c1(1: mp, 1, 1: l1, j) = ch(1: mp, 1, 1: l1, j)+ch(1: mp, 1, 1: l1, jc)
                    c1(1: mp, 1, 1: l1, jc) = ch(1: mp, 1, 1: l1, jc)-ch(1: mp, 1, 1: l1, j)
                end do

                ar1 = 1.0_wp
                ai1 = 0.0_wp
                do l=2, ipph
                    lc = ipp2-l
                    ar1h = dcp*ar1-dsp*ai1
                    ai1 = dcp*ai1+dsp*ar1
                    ar1 = ar1h
                    ch2(1: mp, 1: idl1, l) = c2(1: mp, 1: idl1, 1)+ar1*c2(1: mp, 1: idl1, 2)
                    ch2(1: mp, 1: idl1, lc) = ai1*c2(1: mp, 1: idl1, iip)
                    dc2 = ar1
                    ds2 = ai1
                    ar2 = ar1
                    ai2 = ai1
                    do j=3, ipph
                        jc = ipp2-j
                        ar2h = dc2*ar2-ds2*ai2
                        ai2 = dc2*ai2+ds2*ar2
                        ar2 = ar2h
                        ch2(1: mp, 1: idl1, l) = ch2(1: mp, 1: idl1, l)+ar2*c2(1: mp, 1: idl1, j)
                        ch2(1: mp, 1: idl1, lc) = ch2(1: mp, 1: idl1, lc)+ai2*c2(1: mp, 1: idl1, jc)
                    end do
                end do

                do j=2, ipph
                    do ik=1, idl1
                        ch2(1: mp, ik, 1) = ch2(1: mp, ik, 1)+c2(1: mp, ik, j)
                    end do
                end do

                if (ido >= l1) then
                    cc(1: mp, 1: ido, 1, 1: l1) = ch(1: mp, 1: ido, 1: l1, 1)
                else
                    cc(1: mp, 1: ido, 1, 1: l1) = ch(1: mp, 1: ido, 1: l1, 1)
                end if

                do j=2, ipph
                    jc = ipp2-j
                    j2 = 2*j
                    cc(1: mp, ido, j2-2, 1: l1) = ch(1: mp, 1, 1: l1, j)
                    cc(1: mp, 1, j2-1, 1: l1) = ch(1: mp, 1, 1: l1, jc)
                end do

                if (ido == 1) then
                    return
                end if

                if (nbd >= l1) then
                    do j=2, ipph
                        jc = ipp2-j
                        j2 = 2*j
                        do k=1, l1
                            do i=3, ido, 2
                                ic = idp2-i
                                cc(1: mp, i-1, j2-1, k) = ch(1: mp, i-1, k, j)+ch(1: mp, i-1, k, jc)                !
                                cc(1: mp, ic-1, j2-2, k) = ch(1: mp, i-1, k, j)-ch(1: mp, i-1, k, jc)               !
                                cc(1: mp, i, j2-1, k) = ch(1: mp, i, k, j)+ch(1: mp, i, k, jc)
                                cc(1: mp, ic, j2-2, k) = ch(1: mp, i, k, jc)-ch(1: mp, i, k, j)
                            end do
                        end do
                    end do
                else
                    do j=2, ipph
                        jc = ipp2-j
                        j2 = 2*j
                        do i=3, ido, 2
                            ic = idp2-i
                            cc(1: mp, i-1, j2-1, 1: l1) = ch(1: mp, i-1, 1: l1, j)+ch(1: mp, i-1, 1: l1, jc)                !
                            cc(1: mp, ic-1, j2-2, 1: l1) = ch(1: mp, i-1, 1: l1, j)-ch(1: mp, i-1, 1: l1, jc)               !
                            cc(1: mp, i, j2-1, 1: l1) = ch(1: mp, i, 1: l1, j)+ch(1: mp, i, 1: l1, jc)
                            cc(1: mp, ic, j2-2, 1: l1) = ch(1: mp, i, 1: l1, jc)-ch(1: mp, i, 1: l1, j)
                        end do
                    end do
                end if
            end associate
        end associate

    end subroutine hradfg



    subroutine hrfftb(m, n, r, mdimr, whrfft, work)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: r(mdimr, n)
        integer(ip), intent(in)     :: mdimr
        real(wp),    intent(in)     :: whrfft(n+15)
        real(wp),    intent(out)    :: work(*)
        !----------------------------------------------------------------------

        if (n == 1) return

        call hrftb1(m, n, r, mdimr, work, whrfft, whrfft(n+1))

    end subroutine hrfftb



    subroutine hrftb1(m, n, c, mdimc, ch, wa, fac)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: m
        integer(ip), intent(in)     :: n
        real(wp),    intent(inout)  :: c(mdimc, n)
        integer(ip), intent(in)     :: mdimc
        real(wp),    intent(out)    :: ch(m, n)
        real(wp),    intent(in)     :: wa(n)
        real(wp),    intent(in)     :: fac(15)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip) :: k1, l1, l2, na
        integer(ip) :: nf, iip, iw, ix2, ix3, ix4, ido, idl1
        !----------------------------------------------------------------------

        nf = int(fac(2), kind=ip)
        na = 0
        l1 = 1
        iw = 1
        do k1=1, nf
            iip = int(fac(k1+2), kind=ip)
            l2 = iip*l1
            ido = n/l2
            idl1 = ido*l1

            select case (iip)
                case (2)
                    if (na == 0) then
                        call hradb2(m, ido, l1, c, mdimc, ch, m, wa(iw))
                    else
                        call hradb2(m, ido, l1, ch, m, c, mdimc, wa(iw))
                    end if
                case (3)
                    ix2 = iw+ido
                    if (na == 0) then
                        call  hradb3(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2))
                    else
                        call hradb3(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2))
                    end if
                case (4)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    if (na == 0) then
                        call hradb4(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2), wa(ix3))
                    else
                        call hradb4(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2), wa(ix3))
                    end if
                case (5)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    ix4 = ix3+ido
                    if (na == 0) then
                        call hradb5(m, ido, l1, c, mdimc, ch, m, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                    else
                        call hradb5(m, ido, l1, ch, m, c, mdimc, wa(iw), wa(ix2), wa(ix3), wa(ix4))
                    end if
                case default
                    if (na == 0) then
                        call hradbg(m, ido, iip, l1, idl1, c, c, c, mdimc, ch, ch, m, wa(iw))
                    else
                        call hradbg(m, ido, iip, l1, idl1, ch, ch, ch, m, c, c, mdimc, wa(iw))
                    end if
                    if (ido /= 1) na = 1-na
            end select
            na = 1-na
            l1 = l2
            iw = iw+(iip-1)*ido
        end do

        if (na /= 0) then
            c(1:m, 1:n) = ch
        end if

    end subroutine hrftb1



    subroutine hradb2(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, 2, l1)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(out)    :: ch(mdimch, ido, l1, 2)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa1(ido)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip) :: i, k, ic, idp2
        !----------------------------------------------------------------------

        ch(1:mp, 1,:, 1) = cc(1:mp, 1, 1, :)+cc(1:mp, ido, 2, :)
        ch(1:mp, 1,:, 2) = cc(1:mp, 1, 1, :)-cc(1:mp, ido, 2, :)

        if (ido < 2) then
            return
        else if (ido /= 2) then
            idp2 = ido+2
            do k=1, l1
                do i=3, ido, 2
                    ic = idp2-i
                    ch(1: mp, i-1, k, 1) = &
                        cc(1: mp, i-1, 1, k)+cc(1: mp, ic-1, 2, k)

                    ch(1: mp, i, k, 1) = &
                        cc(1: mp, i, 1, k)-cc(1: mp, ic, 2, k)

                    ch(1: mp, i-1, k, 2) = &
                        wa1(i-2)*(cc(1: mp, i-1, 1, k)-cc(1: mp, ic-1, 2, k)) &         !
                        -wa1(i-1)*(cc(1: mp, i, 1, k)+cc(1: mp, ic, 2, k))

                    ch(1: mp, i, k, 2) = &
                        wa1(i-2)*(cc(1: mp, i, 1, k)+cc(1: mp, ic, 2, k))+wa1(i-1) &      !
                        *(cc(1: mp, i-1, 1, k)-cc(1: mp, ic-1, 2, k))
                end do
            end do

            if (mod(ido,2) == 1) return

        end if

        ch(1:mp, ido,:, 1) = cc(1:mp, ido, 1,:)+cc(1:mp, ido, 1,:)
        ch(1:mp, ido,:, 2) = -(cc(1:mp, 1, 2,:)+cc(1:mp, 1, 2,:))

    end subroutine hradb2



    subroutine hradb3(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, 3, l1)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(out)    :: ch(mdimch, ido, l1, 3)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa1(ido)
        real(wp),    intent(in)     :: wa2(ido)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip)         :: i, k, ic, idp2
        real(wp), parameter :: ARG = TWO_PI/3
        real(wp), parameter :: TAUR = cos(ARG)
        real(wp), parameter :: TAUI = sin(ARG)
        !----------------------------------------------------------------------

        ch(1: mp, 1, 1: l1, 1) = &
            cc(1: mp, 1, 1, 1: l1)+2.0_wp * cc(1: mp, ido, 2, 1: l1)

        ch(1: mp, 1, 1: l1, 2) = &
            cc(1: mp, 1, 1, 1: l1)+(2.0_wp * TAUR)*cc(1: mp, ido, 2, 1: l1) &
            -(2.0_wp *TAUI)*cc(1: mp, 1, 3, 1: l1)

        ch(1: mp, 1, 1: l1, 3) = &
            cc(1: mp, 1, 1, 1: l1)+(2.0_wp * TAUR)*cc(1: mp, ido, 2, 1: l1) &
            +2.0_wp *TAUI*cc(1: mp, 1, 3, 1: l1)

        if (ido == 1) return

        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                ch(1: mp, i-1, k, 1) = &
                    cc(1: mp, i-1, 1, k)+(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k))      !

                ch(1: mp, i, k, 1) = &
                    cc(1: mp, i, 1, k)+(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k))              !

                ch(1: mp, i-1, k, 2) = &
                    wa1(i-2)* &
                    ((cc(1: mp, i-1, 1, k)+TAUR*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)))- &
                    (TAUI*(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k)))) &
                    -wa1(i-1)* &
                    ((cc(1: mp, i, 1, k)+TAUR*(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)))+ &
                    (TAUI*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k))))

                ch(1: mp, i, k, 2) = &
                    wa1(i-2)* &
                    ((cc(1: mp, i, 1, k)+TAUR*(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)))+ &
                    (TAUI*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k)))) &
                    +wa1(i-1)* &
                    ((cc(1: mp, i-1, 1, k)+TAUR*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)))- &
                    (TAUI*(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k))))

                ch(1: mp, i-1, k, 3) = &
                    wa2(i-2)* &
                    ((cc(1: mp, i-1, 1, k)+TAUR*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)))+ &
                    (TAUI*(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k)))) &
                    -wa2(i-1)* &
                    ((cc(1: mp, i, 1, k)+TAUR*(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)))- &
                    (TAUI*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k))))

                ch(1: mp, i, k, 3) = &
                    wa2(i-2)* &
                    ((cc(1: mp, i, 1, k)+TAUR*(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)))- &
                    (TAUI*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k)))) &
                    +wa2(i-1)* &
                    ((cc(1: mp, i-1, 1, k)+TAUR*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)))+ &
                    (TAUI*(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k))))
            end do
        end do

    end subroutine hradb3



    subroutine hradb4(mp, ido, l1, cc, mdimcc, ch, mdimch, wa1, wa2, wa3)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, 4, l1)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(out)    :: ch(mdimch, ido, l1, 4)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa1(ido)
        real(wp),    intent(in)     :: wa2(ido)
        real(wp),    intent(in)     :: wa3(ido)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip)         :: i, k, ic, idp2
        real(wp), parameter :: SQRT2 = sqrt(2.0_wp)
        !----------------------------------------------------------------------

        ch(1: mp, 1, 1: l1, 3) = (cc(1: mp, 1, 1, 1: l1)+cc(1: mp, ido, 4, 1: l1)) &
            -(cc(1: mp, ido, 2, 1: l1)+cc(1: mp, ido, 2, 1: l1))
        ch(1: mp, 1, 1: l1, 1) = (cc(1: mp, 1, 1, 1: l1)+cc(1: mp, ido, 4, 1: l1)) &
            +(cc(1: mp, ido, 2, 1: l1)+cc(1: mp, ido, 2, 1: l1))
        ch(1: mp, 1, 1: l1, 4) = (cc(1: mp, 1, 1, 1: l1)-cc(1: mp, ido, 4, 1: l1)) &
            +(cc(1: mp, 1, 3, 1: l1)+cc(1: mp, 1, 3, 1: l1))
        ch(1: mp, 1, 1: l1, 2) = (cc(1: mp, 1, 1, 1: l1)-cc(1: mp, ido, 4, 1: l1)) &
            -(cc(1: mp, 1, 3, 1: l1)+cc(1: mp, 1, 3, 1: l1))

        if (ido < 2) then
            return
        else if (ido /= 2) then
            idp2 = ido+2
            do k=1, l1
                do i=3, ido, 2
                    ic = idp2-i
                    ch(1: mp, i-1, k, 1) = &
                        (cc(1: mp, i-1, 1, k)+cc(1: mp, ic-1, 4, k)) &
                        +(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k))

                    ch(1: mp, i, k, 1) = &
                        (cc(1: mp, i, 1, k)-cc(1: mp, ic, 4, k)) &
                        +(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k))

                    ch(1: mp, i-1, k, 2) = &
                        wa1(i-2)*((cc(1: mp, i-1, 1, k)-cc(1: mp, ic-1, 4, k)) &          !
                        -(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k)))-wa1(i-1) &
                        *((cc(1: mp, i, 1, k)+cc(1: mp, ic, 4, k))+(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k)))      !

                    ch(1: mp, i, k, 2)= &
                        wa1(i-2)*((cc(1: mp, i, 1, k)+cc(1: mp, ic, 4, k)) &
                        +(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k)))+wa1(i-1) &
                        *((cc(1: mp, i-1, 1, k)-cc(1: mp, ic-1, 4, k))-(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k)))      !

                    ch(1: mp, i-1, k, 3) = &
                        wa2(i-2)*((cc(1: mp, i-1, 1, k)+cc(1: mp, ic-1, 4, k)) &          !
                        -(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)))-wa2(i-1) &
                        *((cc(1: mp, i, 1, k)-cc(1: mp, ic, 4, k))-(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)))          !

                    ch(1: mp, i, k, 3) = &
                        wa2(i-2)*((cc(1: mp, i, 1, k)-cc(1: mp, ic, 4, k)) &
                        -(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)))+wa2(i-1) &
                        *((cc(1: mp, i-1, 1, k)+cc(1: mp, ic-1, 4, k))-(cc(1: mp, i-1, 3, k) &
                        +cc(1: mp, ic-1, 2, k)))

                    ch(1: mp, i-1, k, 4) = &
                        wa3(i-2)*((cc(1: mp, i-1, 1, k)-cc(1: mp, ic-1, 4, k)) &          !
                        +(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k)))-wa3(i-1) &
                        *((cc(1: mp, i, 1, k)+cc(1: mp, ic, 4, k))-(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k)))       !

                    ch(1: mp, i, k, 4) = &
                        wa3(i-2)*((cc(1: mp, i, 1, k)+cc(1: mp, ic, 4, k)) &
                        -(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k)))+wa3(i-1) &
                        *((cc(1: mp, i-1, 1, k)-cc(1: mp, ic-1, 4, k))+(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k)))      !
                end do
            end do

            if (mod(ido, 2) == 1) return

        end if

        ch(1: mp, ido, 1: l1, 1) = &
            (cc(1: mp, ido, 1, 1: l1)+cc(1: mp, ido, 3, 1: l1)) &
            +(cc(1: mp, ido, 1, 1: l1)+cc(1: mp, ido, 3, 1: l1))

        ch(1: mp, ido, 1: l1, 2) = &
            SQRT2*((cc(1: mp, ido, 1, 1: l1)-cc(1: mp, ido, 3, 1: l1)) &
            -(cc(1: mp, 1, 2, 1: l1)+cc(1: mp, 1, 4, 1: l1)))

        ch(1: mp, ido, 1: l1, 3) = &
            (cc(1: mp, 1, 4, 1: l1)-cc(1: mp, 1, 2, 1: l1)) &
            +(cc(1: mp, 1, 4, 1: l1)-cc(1: mp, 1, 2, 1: l1))

        ch(1: mp, ido, 1: l1, 4) = &
            -SQRT2*((cc(1: mp, ido, 1, 1: l1)-cc(1: mp, ido, 3, 1: l1)) &
            +(cc(1: mp, 1, 2, 1: l1)+cc(1: mp, 1, 4, 1: l1)))

    end subroutine hradb4



    subroutine hradb5(mp, ido, l1, cc, mdimcc, ch, mdimch, &
        wa1, wa2, wa3, wa4)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: l1
        real(wp),    intent(out)    :: ch(mdimch, ido, l1, 5)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(inout)  :: cc(mdimcc, ido, 5, l1)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(in)     :: wa1(ido)
        real(wp),    intent(in)     :: wa2(ido)
        real(wp),    intent(in)     :: wa3(ido)
        real(wp),    intent(in)     :: wa4(ido)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)         :: i, k, ic, idp2
        real(wp), parameter :: ARG = TWO_PI/5
        real(wp), parameter :: TR11=cos(ARG)
        real(wp), parameter :: TI11=sin(ARG)
        real(wp), parameter :: TR12=cos(2.0_wp *ARG)
        real(wp), parameter :: TI12=sin(2.0_wp *ARG)
        !----------------------------------------------------------------------

        do k=1, l1
            ch(1: mp, 1, k, 1) = &
                cc(1: mp, 1, 1, k)+2.0_wp *cc(1: mp, ido, 2, k)+2.0_wp *cc(1: mp, ido, 4, k)          !

            ch(1: mp, 1, k, 2) = &
                (cc(1: mp, 1, 1, k)+TR11*2.0_wp *cc(1: mp, ido, 2, k) &
                +TR12*2.0_wp *cc(1: mp, ido, 4, k))-(TI11*2.0_wp *cc(1: mp, 1, 3, k) &
                +TI12*2.0_wp *cc(1: mp, 1, 5, k))

            ch(1: mp, 1, k, 3) = &
                (cc(1: mp, 1, 1, k)+TR12*2.0_wp *cc(1: mp, ido, 2, k) &
                +TR11*2.0_wp *cc(1: mp, ido, 4, k))-(TI12*2.0_wp *cc(1: mp, 1, 3, k) &
                -TI11*2.0_wp *cc(1: mp, 1, 5, k))

            ch(1: mp, 1, k, 4) = &
                (cc(1: mp, 1, 1, k)+TR12*2.0_wp *cc(1: mp, ido, 2, k) &
                +TR11*2.0_wp *cc(1: mp, ido, 4, k))+(TI12*2.0_wp *cc(1: mp, 1, 3, k) &
                -TI11*2.0_wp *cc(1: mp, 1, 5, k))

            ch(1: mp, 1, k, 5) = &
                (cc(1: mp, 1, 1, k)+TR11*2.0_wp *cc(1: mp, ido, 2, k) &
                +TR12*2.0_wp *cc(1: mp, ido, 4, k))+(TI11*2.0_wp *cc(1: mp, 1, 3, k) &
                +TI12*2.0_wp *cc(1: mp, 1, 5, k))
        end do

        if (ido == 1) return

        idp2 = ido+2
        do k=1, l1
            do i=3, ido, 2
                ic = idp2-i
                ch(1: mp, i-1, k, 1) = &
                    cc(1: mp, i-1, 1, k)+(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)) &    !
                    +(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k))

                ch(1: mp, i, k, 1) = &
                    cc(1: mp, i, 1, k)+(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)) &            !
                    +(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k))

                ch(1: mp, i-1, k, 2) = &
                    wa1(i-2)*((cc(1: mp, i-1, 1, k)+TR11* &
                    (cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k))+TR12 &
                    *(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k)))-(TI11*(cc(1: mp, i, 3, k) &             !
                    +cc(1: mp, ic, 2, k))+TI12*(cc(1: mp, i, 5, k)+cc(1: mp, ic, 4, k)))) &
                    -wa1(i-1)*((cc(1: mp, i, 1, k)+TR11*(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)) &         !
                    +TR12*(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k)))+(TI11*(cc(1: mp, i-1, 3, k) &          !
                    -cc(1: mp, ic-1, 2, k))+TI12*(cc(1: mp, i-1, 5, k)-cc(1: mp, ic-1, 4, k))))            !

                ch(1: mp, i, k, 2) = &
                    wa1(i-2)*((cc(1: mp, i, 1, k)+TR11*(cc(1: mp, i, 3, k) &          !
                    -cc(1: mp, ic, 2, k))+TR12*(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k))) &
                    +(TI11*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k))+TI12 &
                    *(cc(1: mp, i-1, 5, k)-cc(1: mp, ic-1, 4, k))))+wa1(i-1) &
                    *((cc(1: mp, i-1, 1, k)+TR11*(cc(1: mp, i-1, 3, k) &
                    +cc(1: mp, ic-1, 2, k))+TR12*(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k))) &           !
                    -(TI11*(cc(1: mp, i, 3, k)+cc(1: mp, ic, 2, k))+TI12 &
                    *(cc(1: mp, i, 5, k)+cc(1: mp, ic, 4, k))))

                ch(1: mp, i-1, k, 3) = &
                    wa2(i-2) &
                    *((cc(1: mp, i-1, 1, k)+TR12*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)) &            !
                    +TR11*(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k)))-(TI12*(cc(1: mp, i, 3, k) &        !
                    +cc(1: mp, ic, 2, k))-TI11*(cc(1: mp, i, 5, k)+cc(1: mp, ic, 4, k)))) &
                    -wa2(i-1) &
                    *((cc(1: mp, i, 1, k)+TR12*(cc(1: mp, i, 3, k)- &
                    cc(1: mp, ic, 2, k))+TR11*(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k))) &
                    +(TI12*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k))-TI11 &
                    *(cc(1: mp, i-1, 5, k)-cc(1: mp, ic-1, 4, k))))

                ch(1: mp, i, k, 3) = &
                    wa2(i-2) &
                    *((cc(1: mp, i, 1, k)+TR12*(cc(1: mp, i, 3, k)- &
                    cc(1: mp, ic, 2, k))+TR11*(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k))) &
                    +(TI12*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k))-TI11 &
                    *(cc(1: mp, i-1, 5, k)-cc(1: mp, ic-1, 4, k)))) &
                    +wa2(i-1) &
                    *((cc(1: mp, i-1, 1, k)+TR12*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)) &            !
                    +TR11*(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k)))-(TI12*(cc(1: mp, i, 3, k) &        !
                    +cc(1: mp, ic, 2, k))-TI11*(cc(1: mp, i, 5, k)+cc(1: mp, ic, 4, k))))

                ch(1: mp, i-1, k, 4) = &
                    wa3(i-2) &
                    *((cc(1: mp, i-1, 1, k)+TR12*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)) &            !
                    +TR11*(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k)))+(TI12*(cc(1: mp, i, 3, k) &        !
                    +cc(1: mp, ic, 2, k))-TI11*(cc(1: mp, i, 5, k)+cc(1: mp, ic, 4, k)))) &
                    -wa3(i-1) &
                    *((cc(1: mp, i, 1, k)+TR12*(cc(1: mp, i, 3, k)- &
                    cc(1: mp, ic, 2, k))+TR11*(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k))) &
                    -(TI12*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k))-TI11 &
                    *(cc(1: mp, i-1, 5, k)-cc(1: mp, ic-1, 4, k))))

                ch(1: mp, i, k, 4) = &
                    wa3(i-2) &
                    *((cc(1: mp, i, 1, k)+TR12*(cc(1: mp, i, 3, k)- &
                    cc(1: mp, ic, 2, k))+TR11*(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k))) &
                    -(TI12*(cc(1: mp, i-1, 3, k)-cc(1: mp, ic-1, 2, k))-TI11 &
                    *(cc(1: mp, i-1, 5, k)-cc(1: mp, ic-1, 4, k)))) &
                    +wa3(i-1) &
                    *((cc(1: mp, i-1, 1, k)+TR12*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)) &            !
                    +TR11*(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k)))+(TI12*(cc(1: mp, i, 3, k) &        !
                    +cc(1: mp, ic, 2, k))-TI11*(cc(1: mp, i, 5, k)+cc(1: mp, ic, 4, k))))

                ch(1: mp, i-1, k, 5) = &
                    wa4(i-2) &
                    *((cc(1: mp, i-1, 1, k)+TR11*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)) &            !
                    +TR12*(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k)))+(TI11*(cc(1: mp, i, 3, k) &        !
                    +cc(1: mp, ic, 2, k))+TI12*(cc(1: mp, i, 5, k)+cc(1: mp, ic, 4, k)))) &
                    -wa4(i-1) &
                    *((cc(1: mp, i, 1, k)+TR11*(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)) &
                    +TR12*(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k)))-(TI11*(cc(1: mp, i-1, 3, k) &          !
                    -cc(1: mp, ic-1, 2, k))+TI12*(cc(1: mp, i-1, 5, k)-cc(1: mp, ic-1, 4, k))))            !

                ch(1: mp, i, k, 5) = &
                    wa4(i-2) &
                    *((cc(1: mp, i, 1, k)+TR11*(cc(1: mp, i, 3, k)-cc(1: mp, ic, 2, k)) &
                    +TR12*(cc(1: mp, i, 5, k)-cc(1: mp, ic, 4, k)))-(TI11*(cc(1: mp, i-1, 3, k) &          !
                    -cc(1: mp, ic-1, 2, k))+TI12*(cc(1: mp, i-1, 5, k)-cc(1: mp, ic-1, 4, k)))) &          !
                    +wa4(i-1) &
                    *((cc(1: mp, i-1, 1, k)+TR11*(cc(1: mp, i-1, 3, k)+cc(1: mp, ic-1, 2, k)) &            !
                    +TR12*(cc(1: mp, i-1, 5, k)+cc(1: mp, ic-1, 4, k)))+(TI11*(cc(1: mp, i, 3, k) &        !
                    +cc(1: mp, ic, 2, k))+TI12*(cc(1: mp, i, 5, k)+cc(1: mp, ic, 4, k))))
            end do
        end do

    end subroutine hradb5



    subroutine hradbg(mp, ido, iip, l1, idl1, cc, c1, c2, mdimcc, &
        ch, ch2, mdimch, wa)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: mp
        integer(ip), intent(in)     :: ido
        integer(ip), intent(in)     :: iip
        integer(ip), intent(in)     :: l1
        integer(ip), intent(in)     :: idl1
        real(wp),    intent(inout)  :: cc(mdimcc, ido, iip, l1)
        real(wp),    intent(inout)  :: c1(mdimcc, ido, l1, iip)
        real(wp),    intent(inout)  :: c2(mdimcc, idl1, iip)
        integer(ip), intent(in)     :: mdimcc
        real(wp),    intent(inout)  :: ch(mdimch, ido, l1, iip)
        real(wp),    intent(inout)  :: ch2(mdimch, idl1, iip)
        integer(ip), intent(in)     :: mdimch
        real(wp),    intent(in)     :: wa(ido)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip)         :: i, j, k, l, j2, ic, jc, lc, is, nbd
        integer(ip)         :: idp2, ipp2, idij, ipph
        real(wp)            :: dc2, ai1, ai2, ar1, ar2, ds2
        real(wp)            :: dcp, arg, dsp, ar1h, ar2h
        !----------------------------------------------------------------------

        arg = TWO_PI/iip
        dcp = cos(arg)
        dsp = sin(arg)
        idp2 = ido+2
        nbd = (ido-1)/2
        ipp2 = iip+2
        ipph = (iip+1)/2

        ch(1:mp,:, :, 1) = cc(1:mp,:, 1, :)

        do j=2, ipph
            jc = ipp2-j
            j2 = 2*j
            ch(1:mp, 1,:, j) = cc(1:mp, ido, j2-2, :)+cc(1:mp, ido, j2-2, :)
            ch(1:mp, 1,:, jc) = cc(1:mp, 1, j2-1, :)+cc(1:mp, 1, j2-1, :)
        end do

        if (ido /= 1) then
            if (nbd >= l1) then
                do j=2, ipph
                    jc = ipp2-j
                    do k=1, l1
                        do i=3, ido, 2
                            ic = idp2-i
                            ch(1:mp, i-1, k, j) = cc(1:mp, i-1, 2*j-1, k)+cc(1:mp, ic-1, 2*j-2, k)
                            ch(1:mp, i-1, k, jc) = cc(1:mp, i-1, 2*j-1, k)-cc(1:mp, ic-1, 2*j-2, k)
                            ch(1:mp, i, k, j) = cc(1:mp, i, 2*j-1, k)-cc(1:mp, ic, 2*j-2, k)
                            ch(1:mp, i, k, jc) = cc(1:mp, i, 2*j-1, k)+cc(1:mp, ic, 2*j-2, k)
                        end do
                    end do
                end do
            else
                do j=2, ipph
                    jc = ipp2-j
                    do i=3, ido, 2
                        ic = idp2-i
                        ch(1:mp, i-1,:, j) = cc(1:mp, i-1, 2*j-1, :)+cc(1:mp, ic-1, 2*j-2, :)
                        ch(1:mp, i-1,:, jc) = cc(1:mp, i-1, 2*j-1, :)-cc(1:mp, ic-1, 2*j-2, :)
                        ch(1:mp, i,:, j) = cc(1:mp, i, 2*j-1, :)-cc(1:mp, ic, 2*j-2, :)
                        ch(1:mp, i,:, jc) = cc(1:mp, i, 2*j-1, :)+cc(1:mp, ic, 2*j-2, :)
                    end do
                end do
            end if
        end if

        ar1 = 1.0_wp
        ai1 = 0.0_wp
        do l=2, ipph
            lc = ipp2-l
            ar1h = dcp*ar1-dsp*ai1
            ai1 = dcp*ai1+dsp*ar1
            ar1 = ar1h
            c2(1: mp, 1: idl1, l) = ch2(1: mp, 1: idl1, 1)+ar1*ch2(1: mp, 1: idl1, 2)
            c2(1: mp, 1: idl1, lc) = ai1*ch2(1: mp, 1: idl1, iip)
            dc2 = ar1
            ds2 = ai1
            ar2 = ar1
            ai2 = ai1
            do j=3, ipph
                jc = ipp2-j
                ar2h = dc2*ar2-ds2*ai2
                ai2 = dc2*ai2+ds2*ar2
                ar2 = ar2h
                c2(1: mp, 1: idl1, l) = c2(1: mp, 1: idl1, l)+ar2*ch2(1: mp, 1: idl1, j)
                c2(1: mp, 1: idl1, lc) = c2(1: mp, 1: idl1, lc)+ai2*ch2(1: mp, 1: idl1, jc)
            end do
        end do

        do j=2, ipph
            ch2(1: mp, 1: idl1, 1) = ch2(1: mp, 1: idl1, 1)+ch2(1: mp, 1: idl1, j)
        end do

        do j=2, ipph
            jc = ipp2-j
            ch(1: mp, 1, 1: l1, j) = c1(1: mp, 1, 1: l1, j)-c1(1: mp, 1, 1: l1, jc)
            ch(1: mp, 1, 1: l1, jc) = c1(1: mp, 1, 1: l1, j)+c1(1: mp, 1, 1: l1, jc)
        end do

        if (ido /= 1) then
            if (nbd >= l1) then
                do j=2, ipph
                    jc = ipp2-j
                    do k=1, l1
                        do i=3, ido, 2
                            ch(1: mp, i-1, k, j) = c1(1: mp, i-1, k, j)-c1(1: mp, i, k, jc)
                            ch(1: mp, i-1, k, jc) = c1(1: mp, i-1, k, j)+c1(1: mp, i, k, jc)
                            ch(1: mp, i, k, j) = c1(1: mp, i, k, j)+c1(1: mp, i-1, k, jc)
                            ch(1: mp, i, k, jc) = c1(1: mp, i, k, j)-c1(1: mp, i-1, k, jc)
                        end do
                    end do
                end do
            else
                do j=2, ipph
                    jc = ipp2-j
                    do i=3, ido, 2
                        ch(1: mp, i-1, 1: l1, j) = c1(1: mp, i-1, 1: l1, j)-c1(1: mp, i, 1: l1, jc)
                        ch(1: mp, i-1, 1: l1, jc) = c1(1: mp, i-1, 1: l1, j)+c1(1: mp, i, 1: l1, jc)
                        ch(1: mp, i, 1: l1, j) = c1(1: mp, i, 1: l1, j)+c1(1: mp, i-1, 1: l1, jc)
                        ch(1: mp, i, 1: l1, jc) = c1(1: mp, i, 1: l1, j)-c1(1: mp, i-1, 1: l1, jc)
                    end do
                end do
            end if
        end if

        if (ido == 1) return

        c2(1:mp, 1: idl1, 1) = ch2(1:mp, 1: idl1, 1)
        c1(1:mp, 1, 1: l1, 2: iip) = ch(1:mp, 1, 1: l1, 2: iip)

        if (nbd <= l1) then
            is = -ido
            do j=2, iip
                is = is+ido
                idij = is
                do i=3, ido, 2
                    idij = idij+2
                    do k=1, l1
                        c1(1: mp, i-1, k, j) = wa(idij-1)*ch(1: mp, i-1, k, j)-wa(idij)* &           !
                            ch(1: mp, i, k, j)
                        c1(1: mp, i, k, j) = wa(idij-1)*ch(1: mp, i, k, j)+wa(idij)* &               !
                            ch(1: mp, i-1, k, j)
                    end do
                end do
            end do
        else
            is = -ido
            do j=2, iip
                is = is+ido
                do k=1, l1
                    idij = is
                    do i=3, ido, 2
                        idij = idij+2
                        c1(1: mp, i-1, k, j) = wa(idij-1)*ch(1: mp, i-1, k, j)-wa(idij)* &           !
                            ch(1: mp, i, k, j)
                        c1(1: mp, i, k, j) = wa(idij-1)*ch(1: mp, i, k, j)+wa(idij)* &               !
                            ch(1: mp, i-1, k, j)
                    end do
                end do
            end do
        end if

    end subroutine hradbg



end module type_HFFTpack
