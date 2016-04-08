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
!
! subroutine lfim (init, theta, l, n, nm, pb, id, wlfim)
!
! dimension of           theta(l),  pb(id, nm+1),  wlfim(4*l*(nm+1))
! arguments
!
! purpose                given n and l, routine lfim calculates
!                        the normalized associated legendre functions
!                        pbar(n, m, theta) for m=0, ..., n and theta(i)
!                        for i=1, ..., l where
!
!                        pbar(m, n, theta) = sqrt((2*n+1)*factorial(n-m)
!                        /(2*factorial(n+m)))*sin(theta)**m/(2**n*
!                        factorial(n)) times the (n+m)th derivative of
!                        (x**2-1)**n with respect to x=cos(theta)
!
! usage                  call lfim (init, theta, l, n, nm, pb, id, wlfim)
!
! arguments
! on input               init
!                        = 0
!                            initialization only - using parameters
!                            l, nm and array theta, subroutine lfim
!                            initializes array wlfim for subsequent
!                            use in the computation of the associated
!                            legendre functions pb. initialization
!                            does not have to be repeated unless
!                            l, nm, or array theta are changed.
!                        = 1
!                            subroutine lfim uses the array wlfim that
!                            was computed with init = 0 to compute pb.
!
!                        theta
!                          an array that contains the colatitudes
!                          at which the associated legendre functions
!                          will be computed. the colatitudes must be
!                          specified in radians.
!
!                        l
!                          the length of the theta array. lfim is
!                          vectorized with vector length l.
!
!                        n
!                          nonnegative integer, less than nm, specifying
!                          degree of pbar(n, m, theta). subroutine lfim
!                          must be called starting with n=0. n must be
!                          incremented by one in subsequent calls and
!                          must not exceed nm.
!
!                        nm
!                          the maximum value of n and m
!
!                        id
!                          the first dimension of the two dimensional
!                          array pb as it appears in the program that
!                          calls lfim. (see output parameter pb)
!
!                        wlfim
!                          an array with length 4*l*(nm+1) which
!                          must be initialized by calling lfim
!                          with init=0 (see parameter init)  it
!                          must not be altered between calls to
!                          lfim.
!
!
! on output              pb
!                          a two dimensional array with first
!                          dimension id in the program that calls
!                          lfim. the second dimension of pb must
!                          be at least nm+1. starting with n=0
!                          lfim is called repeatedly with n being
!                          increased by one between calls. on each
!                          call, subroutine lfim computes
!                          = pbar(m, n, theta(i)) for m=0, ..., n and
!                          i=1, ...l.
!
!                        wlfim
!                          array containing values which must not
!                          be altered unless l, nm or the array theta
!                          are changed in which case lfim must be
!                          called with init=0 to reinitialize the
!                          wlfim array.
!
! special conditions     n must be increased by one between calls
!                        of lfim in which n is not zero.
!
! precision              single
!
!
! algorithm              routine lfim calculates pbar(n, m, theta) using
!                        a four term recurrence relation. (unpublished
!                        notes by paul n. swarztrauber)
!
subroutine lfim(init, theta, l, n, nm, pb, id, wlfim)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: init
    real,    intent (in)     :: theta(*)
    integer, intent (in)     :: l
    integer, intent (in)     :: n
    integer, intent (in)     :: nm
    real,    intent (in out) :: pb(1)
    integer, intent (in)     :: id
    real,    intent (in out) :: wlfim(1)
    !----------------------------------------------------------------------

    !
    !     total length of wlfim is 4*l*(nm+1)
    !
    associate(lnx => l*(nm+1) )
        associate( iw1=> lnx+1)
            associate( iw2 => iw1+lnx)
                associate( iw3 => iw2+lnx)
                    call lfim1(init, theta, l, n, nm, id, pb, wlfim, wlfim(iw1), &
                        wlfim(iw2), wlfim(iw3), wlfim(iw2))
                end associate
            end associate
        end associate
    end associate

end subroutine lfim


subroutine lfim1(init, theta, l, n, nm, id, p3, phz, ph1, p1, p2, cp)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: init
    real,    intent (in)     :: theta(*)
    integer, intent (in)     :: l
    integer, intent (in)     :: n
    integer, intent (in)     :: nm
    integer, intent (in)     :: id
    real,    intent (in out) :: p3(id,*)
    real,    intent (in out) :: phz(l,*)
    real,    intent (in out) :: ph1(l,*)
    real,    intent (in out) :: p1(l,*)
    real,    intent (in out) :: p2(l,*)
    real,    intent (in out) :: cp(*)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer         :: i, m, nm1, nh, mp1, np1, mp3, nmp1
    real            :: cc, dd, ee, cn, fm, fn, fnmm, fnpm
    real            :: tm, tn, temp
    real, parameter :: SQRT2 = sqrt(2.0)
    real, parameter :: SQRT5 = sqrt(5.0)
    real, parameter :: SQRT6 = sqrt(6.0)
    real, parameter :: ONE_OVER_SQRT2 = 1.0/SQRT2
    real, parameter :: ONE_OVER_SQRT6 = 1.0/SQRT6
    real, parameter :: SQRT5_OVER_SQRT6 = SQRT5/SQRT6
    !----------------------------------------------------------------------

    nmp1 = nm+1

    if (init /= 0) go to 5

    phz(:, 1) = ONE_OVER_SQRT2

    do np1=2, nmp1
        nh = np1-1
        call alfk(nh, 0, cp)
        do i=1, l
            call lfpt(nh, 0, theta(i), cp, phz(i, np1))
        end do
        call alfk(nh, 1, cp)
        do  i=1, l
            call lfpt(nh, 1, theta(i), cp, ph1(i, np1))
        end do
    end do


    return
5   if (n > 2) go to 60
    if (n-1< 0) then
        goto 25
    else if (n-1 == 0) then 
        goto 30
    else 
        goto 35
    end if

25  p3(:, 1) = phz(:, 1)
    return

30  p3(:, 1) = phz(:, 2)
    p3(:, 2) = ph1(:, 2)
    return

35  p3(:, 1) = phz(:, 3)
    p3(:, 2) = ph1(:, 3)
    p3(:, 3) = SQRT5_OVER_SQRT6 * phz(:, 1) - ONE_OVER_SQRT6 * p3(:, 1)
    p1(:, 1) = phz(:, 2)
    p1(:, 2) = ph1(:, 2)
    p2(:, 1) = phz(:, 3)
    p2(:, 2) = ph1(:, 3)
    p2(:, 3) = p3(:, 3)
    return
60  nm1 = n-1
    np1 = n+1
    fn = real(n)
    tn = fn+fn
    cn = (tn+1.0)/(tn-3.0)
    p3(:, 1) = phz(:, np1)
    p3(:, 2) = ph1(:, np1)

    if (nm1 < 3) go to 71

    do mp1=3, nm1
        m = mp1-1
        fm = real(m)
        fnpm = fn+fm
        fnmm = fn-fm
        temp = fnpm*(fnpm-1.0)
        cc = sqrt(cn*(fnpm-3.0)*(fnpm-2.0)/temp)
        dd = sqrt(cn*fnmm*(fnmm-1.0)/temp)
        ee = sqrt((fnmm+1.0)*(fnmm+2.0)/temp)
        p3(:, mp1) = cc*p1(i, mp1-2)+dd*p1(:, mp1)-ee*p3(:, mp1-2)
    end do

71  fnpm = fn+fn-1.
    temp = fnpm*(fnpm-1.0)
    cc = sqrt(cn*(fnpm-3.0)*(fnpm-2.0)/temp)
    ee = sqrt(6.0/temp)

    p3(:, n) = cc*p1(:, n-2)-ee*p3(:, n-2)

    fnpm = fn+fn
    temp = fnpm*(fnpm-1.0)
    cc = sqrt(cn*(fnpm-3.0)*(fnpm-2.0)/temp)
    ee = sqrt(2.0/temp)
    p3(:, n+1) = cc*p1(:, n-1)-ee*p3(:, n-1)


    do mp1=1, np1
        p1(:, mp1) = p2(:, mp1)
        p2(:, mp1) = p3(:, mp1)
    end do

end subroutine lfim1
