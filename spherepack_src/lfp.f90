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
!
! subroutine lfp (init, n, m, l, cp, pb, w)
!
! dimension of           cp((n/2)+1), pb(l), w(5*l+41)
! arguments
!
! purpose                routine lfp uses coefficients computed by
!                        routine alfk to calculate the 64-bit double precision
!                        normalized associated legendre function pbar(n, 
!                        m, theta) at colatitudes theta=(i-1)*pi/(l-1), 
!                        i=1, ..., l. subroutine lfp evaluates pbar
!                        using one of the following trigonometric
!                        expansions
!
!                        1) for n even and m even, pbar(m, n, theta) =
!                           .5*cp(1) plus the sum from k=1 to k=n/2
!                           of cp(k)*cos(2*k*th)
!
!                        2) for n even and m odd, pbar(m, n, theta) =
!                           the sum from k=1 to k=n/2 of
!                           cp(k)*sin(2*k*th)
!
!                        3) for n odd and m even, pbar(m, n, theta) =
!                           the sum from k=1 to k=(n+1)/2 of
!                           cp(k)*cos((2*k-1)*th)
!
!                        4) for n odd and m odd,  pbar(m, n, theta) =
!                           the sum from k=1 to k=(n+1)/2 of
!                           cp(k)*sin((2*k-1)*th)
!
!
! usage                  call lfp(init, n, m, l, cp, pb, w)
!
! arguments
!
! on input               init
!                          = 0 initialization only
!                          = 1 compute pbar(n, m, theta)
!
!                          lfp call with init = 0 initializes array w;
!                          no values of pbar(n, m, theta) are computed.
!                          init=0 should be used on the first call, or
!                          if l or w values differ from those in the
!                          previous call.
!
!                        n
!                          nonnegative integer, less than l, specifying
!                          the degree of pbar(n, m, theta)
!
!                        m
!                          is the order of pbar(n, m, theta). m can be
!                          any integer however pbar(n, m, theta) = 0
!                          if abs(m) is greater than n and
!                          pbar(n, m, theta) = (-1)**m*pbar(n, -m, theta)
!                          for negative m.
!
!                        l
!                          number of colatitudes theta=(i-1)*pi/(l-1)
!                          for i=1, ..., l where l is greater than 1.
!                          l must be an odd integer.
!
!                        cp
!                          64-bit double precision array of length (n/2)+1
!                          containing coefficients computed by routine
!                          alfk
!
!                        w
!                          a 64-bit double precision work array with at
!                          least 5*l+41 locations
!
! on output              pb
!                          64-bit double precision array of length l containing
!                          pbar(n, m, theta), theta=(i-1)*pi/(l-1) for i=1
!                          , ..., l.
!
!                        w
!                          a 64-bit double precision array containing values
!                          which must not be destroyed if the next call
!                          will have the same value of input parameter n
!
! special conditions     calls to routine lfp must be preceded by an
!                        appropriate call to routine alfk.
!
! precision              64-bit double precision
!
! algorithm              the trigonometric series formula used by
!                        routine lfp to calculate pbar(n, m, theta) for
!                        theta=(i-1)*pi/(l-1), i=1, ..., n, depends on
!                        m and n as follows:
!
!                           1) for n even and m even, the formula is
!                              .5*cp(1) plus the sum from k=1 to k=n/2
!                              of cp(k)*cos(2*k*theta)
!                           2) for n even and m odd. the formula is
!                              the sum from k=1 to k=n/2 of
!                              cp(k)*sin(2*k*theta)
!                           3) for n odd and m even, the formula is
!                              the sum from k=1 to k=(n+1)/2 of
!                              cp(k)*cos((2*k-1)*theta)
!                           4) for n odd and m odd, the formula is
!                              the sum from k=1 to k=(n+1)/2 of
!                              cp(k)*sin((2*k-1)*theta)
!
! accuracy               comparison between routines lfp and double
!                        precision dlfp on the cray1 indicates greater
!                        accuracy for smaller values of input parameter
!                        n.  agreement to 12 places was obtained for
!                        n=10 and to 11 places for n=100.
!
! timing                 time per call to routine lfp is dependent on
!                        the input parameters l and n.
!
subroutine lfp (init, n, m, l, cp, pb, w)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)  :: init
    integer (ip), intent (in)  :: n
    integer (ip), intent (in)  :: m
    integer (ip), intent (in)  :: l
    real (wp)                  :: cp(1)
    real (wp)                  :: pb(1)
    real (wp)                  :: w(1)
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip) :: ma, iw1, iw2
    !----------------------------------------------------------------------

    pb(1:l) = 0.0_wp

    ma = abs(m)

    if (ma > n) return

    iw1 = 2*l+12
    iw2 = iw1+3*(l+1)/2+15

    call lfp1(init, n, ma, l, cp, pb, w, w(iw1), w(iw2))

contains

    subroutine lfp1(init, n, m, l, cp, p, wsave1, wsave2, wsave3)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in) :: init
        integer (ip), intent (in) :: n
        integer (ip), intent (in) :: m
        integer (ip), intent (in) :: l
        real (wp)                 :: cp(*)
        real (wp)                 :: p(*)
        real (wp)                 :: wsave1(*)
        real (wp)                 :: wsave2(*)
        real (wp)                 :: wsave3(*)
        !----------------------------------------------------------------------
        ! Dictionary: calling arguments
        !----------------------------------------------------------------------
        integer (ip), save   :: lc, lq, ls
        integer (ip)         :: i
        integer (ip)         :: lm1, np1, ls2, kdp, lmi
        real (wp)            :: dt
        real (wp), parameter :: PI = acos(-1.0_wp)
        real (wp), parameter :: ONE_OVER_SQRT2 = 1.0_wp/sqrt(2.0_wp)
        !----------------------------------------------------------------------

        select case (init)
            case (0)

                lc=(l+1)/2
                ls=lc-2
                lq=lc-1
                call sinti(ls, wsave1)
                call costi(lc, wsave2)
                call cosqi(lq, wsave3)

            case default

                if ((n <= 0) .and. (m <= 0)) then
                    p(1:l) = ONE_OVER_SQRT2
                else
                    ls2 = (l+1)/2
                    lm1 = l-1
                    np1 = n+1

                    dt = PI/lm1

                    if (mod(n, 2) <= 0) then
                        if (mod(m, 2) <= 0) then
                            kdp = n/2+1
                            p(1:kdp)=0.5_wp*cp(1:kdp)
                            p(lc)= 2.0_wp*p(lc)

                            call cost(lc, p, wsave2)

                            do i=1, lc
                                lmi=l-i
                                p(lmi+1)=p(i)
                            end do
                        else
                            kdp=n/2
                            p(2:kdp+1)=0.5_wp*cp(1:kdp)
                            p(ls+2)=0.0_wp

                            call sint(ls, p(2), wsave1)

                            do i=1, ls
                                lmi=l-i
                                p(lmi)=-p(i+1)
                            end do

                            p(l)=0.0_wp
                        end if
                    else
                        kdp=(n+1)/2

                        if (mod(m, 2) <= 0) then

                            p(1:kdp)=0.25_wp*cp(1:kdp)

                            call cosqb(lq, p, wsave3)

                            do i=1, lq
                                lmi=l-i
                                p(lmi+1)=-p(i)
                            end do

                        else
                            p(2:kdp+1)=0.25_wp*cp(1:kdp)

                            call sinqb(lq, p(2), wsave3)

                            do i=1, lq
                                lmi=l-i
                                p(lmi)=p(i+1)
                            end do

                            p(l)=0.0_wp
                        end if
                    end if
                end if
        end select

    end subroutine lfp1

end subroutine lfp
