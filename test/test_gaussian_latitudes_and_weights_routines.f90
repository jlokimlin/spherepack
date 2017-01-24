!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK                               *
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
!     this program tests subroutine compute_gaussian_latitudes_and_weights for computing 
!     the Gauss Legendre points and weights in file compute_gaussian_latitudes_and_weights.f
!     It tests only the april 2002 version and not the
!     older version in file compute_gaussian_latitudes_and_weights.old.f
!     gauss points and weights are computed using newtons method
!     with equally spaced points as first guess. Points are
!     computed as theta_i where x_i = cos(theta_i)
!     the method is discussed in "On computing the Gauss_Legendre
!     points and weights" by Paul N Swarztrauber, accepted for 
!     publication in the SIAM journal of scientific computing.
!                         April 2002
!
program test_gaussian_latitudes_and_weights_routines

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT, &
        wp => REAL64, &
        sp => REAL32

    use spherepack_library, only: &
        compute_gaussian_latitudes_and_weights, pi

    ! Explicit typing only
    implicit none

    real :: dmax
    real :: hold
    integer :: i
    integer :: ido
    integer :: ierror
    integer :: irm
    integer :: ldw
    integer :: lwork
    integer, parameter :: NLAT = 63
    real :: pmax
    real :: rerr
    real :: rmax

    real :: sums

    real :: tdoub
    real :: tsing
    real :: wmax

    real(wp) :: theta(NLAT), wts(NLAT), work(NLAT+2)
    real(wp) :: dtheta(NLAT), dwts(NLAT)
    real(wp) :: dwmx, tmax, dwmax, dtmax
    real     :: theta_sp(NLAT), wts_sp(NLAT), work_sp(NLAT+1)
    real     :: diff(NLAT)
    real(sp) :: t1(2),t2(2)
    real     :: sumw
    
    ! Description
    write( stdout, '(/a/)') &
        '     gaussian_latitudes_and_weights_routines *** TEST RUN *** '

    lwork = NLAT+1
    hold = etime(t1)
    hold = t1(1)
    call gsqd(NLAT,theta,wts,work,lwork,ierror)
    tdoub = etime(t1)
    tdoub = t1(1)-hold

    ! Check error flag
    if (ierror /= 0) write( stdout, '(a,i5)') ' ierror = ', ierror

    write( stdout, '(a,1pe15.6)') ' tdoub ', tdoub

    ! Compare double with double
    ldw = NLAT+2
    hold = etime(t1)
    hold = t1(1)
    call compute_gaussian_latitudes_and_weights(NLAT,dtheta,dwts,ierror)
    tdoub = etime(t1)
    tdoub = t1(1)-hold

    if (ierror /= 0) write( stdout,'(a,i5)') ' ierror = ', ierror

    write( stdout, '(a,1pe15.6)') &
        ' tdoub compute_gaussian_latitudes_and_weights', tdoub

    dwmx = 0.0
    tmax = 0.0
    dwmax = 0.0
    dtmax = 0.0
    ido = (NLAT+1)/2

    do i=1,ido
        dtmax = max(dtmax, abs(theta(i)-dtheta(i)))
        dwmax = max(dwmax, abs(wts(i)-dwts(i)))
        tmax = max(tmax, abs(theta(i)))
        dwmx = max(dwmx, abs(wts(i)))
    end do

    dtmax = dtmax/tmax
    dwmax = dwmax/dwmx
    write( stdout, 33) NLAT,dtmax,dwmax
33  format(' nlat',i6,'  points ',1pd15.6,' weights ',d15.6)

    sumw = sum(wts)

    write( stdout, 638) sumw
638 format('  sumw ',1pd39.30)

    hold = etime(t2)
    hold = t2(1)
    lwork = NLAT+2
    call sgaqd(NLAT,theta_sp,wts_sp,work_sp,lwork,ierror)
    tsing = etime(t2)
    tsing = t2(1)-hold
    if (ierror /= 0) write( stdout, 5) ierror
5   format(' iserror=',i5)

    sums =sum(wts_sp)

    write( stdout, 636) sums
636 format('  sums ',1pe24.15)
    dmax = 0.0
    wmax = 0.0
    rmax = 0.0
    ido = (NLAT+1)/2
    do i=1,ido
        diff(i) = wts(i)-wts_sp(i)
        dmax = max(dmax,abs(diff(i)))
        wmax = max(wmax,abs(wts_sp(i)))
        rerr = abs(diff(i)/wts_sp(i))
        if (rerr>rmax) then
            rmax = rerr
            irm = i
        end if
6       format(' wts ',1pd15.6,' swts ',e15.6,' rpter ',e15.6)
    end do

    !      write( stdout, 7) (diff(i),i=nlat-25,nlat)
7   format(' diff in weights'/(1p8e10.03))
    dmax = dmax/wmax
    write( stdout, 9) NLAT,irm,dmax,rmax
9   format(' weights: nlat ',i6,' irele ',i6, &
        ' dmax ',1pe15.6,' rmax ',1pe15.6)
    dmax = 0.0
    pmax = 0.0
    rmax = 0.0
    do i=1,ido
        diff(i) = theta(i)-theta_sp(i)
        dmax = max(dmax,abs(diff(i)))
        pmax = max(pmax,abs(theta_sp(i)))
        rerr = abs(diff(i)/theta_sp(i))
        if (rerr>rmax) then
            rmax = rerr
            irm = i
        end if
    end do
    !      write( stdout, 11) (diff(i),i=nlat-25,nlat)
11  format(' diff in points'/(1p8e10.03))
    dmax = dmax/pmax
    write( stdout, 12) NLAT,irm,dmax,rmax
12  format(' points:  nlat ',i6,' irele ',i6, &
        ' dmax ',1pe15.6,' rmax ',1pe15.6)

    dmax = 0.0
    do i=1,NLAT
        diff(i) = cos(theta(i))-cos(theta_sp(i))
        dmax = max(dmax,abs(diff(i)))
    end do

    ! format(' diff in points'/(1p8e10.03))
    write( stdout, 112) dmax
112 format(' max difference in mu',1pe15.6)
    write( stdout, 1) NLAT,tsing,tdoub
1   format(' nlat',i6,' tsing',1pe15.6,' tdoub',e15.6)

contains
    !
    !     subroutine gsqd is a single precision version of compute_gaussian_latitudes_and_weights.
    !     gauss points and weights are computed using newtons method
    !     with equally spaced points as first guess. Points are
    !     computed as theta_i where x_i = cos(theta_i)
    !     the method is discussed in "On computing the Gauss_Legendre
    !     points and weights" by Paul N Swarztrauber, accepted for 
    !     publication in the SIAM journal of scientific computing.
    !                         April 2002
    !
    subroutine gsqd(nlat,theta,wts,dwork,ldwork,ierror)

        integer :: i
        integer :: ierror
        integer :: ldwork
        integer :: nlat
        integer :: ns2
        !
        !
        !     subroutine gsqd computes the nlat gaussian colatitudes and weights
        !     in real. The colatitudes are in radians and lie in the
        !     in the interval (0,pi).
        !
        !     input parameters
        !
        !     nlat    the number of gaussian colatitudes in the interval (0,pi)
        !             (between the two poles).  nlat must be greater than zero.
        !
        !     dwork   a real temporary work space.
        !
        !     ldwork  the length of the work space  in the routine calling gsqd
        !             ldwork must be at least nlat+1.
        !
        !     output parameters
        !
        !     theta   a real vector of length nlat containing the
        !             nlat gaussian colatitudes on the sphere in increasing radians
        !             in the interval (0,pi).
        !
        !     wts     a real vector of length nlat containing the
        !             nlat gaussian weights.
        !
        !     ierror = 0 no errors
        !            = 1 if ldwork.lt.nlat+1
        !            = 2 if nlat.le.0
        !
        !  *****************************************************************
        !
        dimension dwork(nlat+1),theta(nlat),wts(nlat)
        real HALF_PI,x,theta,wts,dwork
        ierror = 1
        !
        !     check work space length
        !
        if (ldwork<nlat+1) return
        ierror = 2
        if (nlat<=0) return
        ierror = 0
        !
        !     compute weights and points analytically when nlat=1,2
        !
        if (nlat==1) then
            theta(1) = acos(0.0)
            wts(1) = 2.0
            return
        end if
        if (nlat==2) then
            x = sqrt(1.0/3.0)
            theta(1) = acos(x)
            theta(2) = acos(-x)
            wts(1) = 1.0
            wts(2) = 1.0
            return
        end if
        !
        !     compute points
        !
        call gsqd_lower_routine(nlat,theta,dwork)
        !
        !     compute weights
        !
        call egwts(nlat,theta,wts,dwork)
        !
        !     extend points and weights via symmetries
        !
        HALF_PI = acos(-1.0)
        ns2 = nlat/2
        do i=1,ns2
            wts(nlat-i+1) = wts(i)
            theta(nlat-i+1) = HALF_PI-theta(i)
        end do
        return
    end subroutine gsqd
    !
    subroutine gsqd_lower_routine(nlat,theta,cp)

        real :: eps
        real :: summation
        
        integer :: it
        integer :: nhalf
        integer :: nix
        integer :: nlat
        integer :: ns2
        real  theta((nlat+1)/2),cp(nlat/2+1)
        real pi,HALF_PI,dtheta,dthalf, &
            cmax,dcor,pb,dpb,sgnd,zero,zlast
        !
        eps = sqrt(epsilon(1.0))
        eps = eps*sqrt(eps)
        PI = acos(-1.0)
        HALF_PI = PI/2
        theta(1) = HALF_PI
        if (nlat==1) goto 30
        ns2 = nlat/2
        nhalf = (nlat+1)/2
        !
        call dlfcz (nlat,cp)
        !
        !      check fourier-legendre coefficients
        !
        summation = sum(cp(1:ns2+1))/ns2

        !format(' check on dble f-l coefficients ',1pe15.6)
        !
        dtheta = HALF_PI/nhalf
        dthalf = dtheta/2.0
        cmax = .2*dtheta
        if (mod(nlat,2)/=0) then
            theta(nhalf) = HALF_PI
            zero = HALF_PI-dtheta
            nix = nhalf-1
        else
            zero = HALF_PI-dthalf
            nix = nhalf
        end if
9       it = 0
10      it = it+1
        zlast = zero
        !
        !     newton iterations to convergence
        !
        call lft (0,nlat,zero,cp,pb)
        call dlft (0,nlat,zero,cp,dpb)
        dcor = pb/dpb
        sgnd = 1.0
        if (dcor /= 0.0) sgnd = dcor/abs(dcor)
        !      write( stdout, 2) nix,zero,theta(nix),dcor,cmax,sgnd
2       format(i7,1p5d15.6)
        dcor = sgnd*min(abs(dcor),cmax)
        zero = zero-dcor
        if (abs(zero-zlast)>eps*abs(zero)) goto 10
        theta(nix) = zero
        nix = nix-1
        if (nix==0) goto 30
        if (nix==nhalf-1)  zero = 3.0*zero-pi
        if (nix<nhalf-1)  zero = zero+zero-theta(nix+2)
        goto 9
30      return
    end subroutine gsqd_lower_routine

    subroutine egwts(n,theta,wts,work)

        integer :: n
        !
        !     computes gauss weights as described in swarztrauber
        !     and spotz, generalized harmonic transforms
        !
        real theta(n),wts(n),work(n+1)
        !
        call egwts1(n,theta,wts,work,work(n/2+2))

    end subroutine egwts

    subroutine egwts1(n,theta,wts,dcp,cp)

        integer :: i
        integer :: n
        integer :: nhalf
        real theta((n+1)/2),wts((n+1)/2),cp((n-1)/2+1), &
            dcp(n/2+1),fn,sqnn,pb,dpb
        fn = n
        sqnn = sqrt((fn+fn-1)*(fn+fn+1))
        call dlfcz (n-1,cp)
        call dlfcz (n,dcp)
        nhalf = (n+1)/2
        do i=1,nhalf
            call lft(0,n-1,theta(i),cp,pb)
            call dlft(0,n,theta(i),dcp,dpb)
            wts(i) = -sqnn*sin(theta(i))/(fn*pb*dpb)
        end do

    end subroutine egwts1

    subroutine lfc (m,n,cp)

        integer :: i
        integer :: l
        integer :: m
        integer :: ma
        integer :: n
        integer :: nex
        integer :: nmms2
        !
        real cp,fnum,fden,fnmh,a1,b1,c1,cp2,fnnp1,fnmsq,fk, &
            t1,t2,pm1,sc10,sc20,sc40
        dimension       cp(n/2+1)
        parameter (sc10=1024.0)
        parameter (sc20=sc10*sc10)
        parameter (sc40=sc20*sc20)
        !
        cp(1) = 0.0
        ma = iabs(m)
        if (ma > n) return
        if (n-1< 0) then
            goto 2
        else if (n-1 == 0) then
            goto 3
        else
            goto 5
        end if
2       cp(1) = sqrt(2.0)
        return
3       if (ma /= 0) goto 4
        cp(1) = sqrt(1.5)
        return
4       cp(1) = sqrt(.75)
        if (m == -1) cp(1) = -cp(1)
        return
5       if (mod(n+ma,2) /= 0) goto 10
        nmms2 = (n-ma)/2
        fnum = n+ma+1
        fnmh = n-ma+1
        pm1 = 1.0
        goto 15
10      nmms2 = (n-ma-1)/2
        fnum = n+ma+2
        fnmh = n-ma+2
        pm1 = -1.0
        !      t1 = 1.
        !      t1 = 2.0**(n-1)
        !      t1 = 1.0/t1
15      t1 = 1.0/sc20
        nex = 20
        fden = 2.0
        if (nmms2 < 1) goto 20
        do 18 i=1,nmms2
            t1 = fnum*t1/fden
            if (t1 > sc20) then
                t1 = t1/sc40
                nex = nex+40
            end if
            fnum = fnum+2.
            fden = fden+2.
18      continue
20      t1 = t1/2.0**(n-1-nex)
        if (mod(ma/2,2) /= 0) t1 = -t1
        t2 = 1.
        if (ma == 0) goto 26
        do 25 i=1,ma
            t2 = fnmh*t2/(fnmh+pm1)
            fnmh = fnmh+2.
25      continue
26      cp2 = t1*sqrt((n+.5)*t2)
        fnnp1 = n*(n+1)
        fnmsq = fnnp1-2.0*ma*ma
        l = (n+1)/2
        if (mod(n,2) == 0 .and. mod(ma,2) == 0) l = l+1
        cp(l) = cp2
        if (m >= 0) goto 29
        if (mod(ma,2) /= 0) cp(l) = -cp(l)
29      if (l <= 1) return
        fk = n
        a1 = (fk-2.)*(fk-1.)-fnnp1
        b1 = 2.*(fk*fk-fnmsq)
        cp(l-1) = b1*cp(l)/a1
30      l = l-1
        if (l <= 1) return
        fk = fk-2.
        a1 = (fk-2.)*(fk-1.)-fnnp1
        b1 = -2.*(fk*fk-fnmsq)
        c1 = (fk+1.)*(fk+2.)-fnnp1
        cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
        goto 30
    end subroutine lfc
    subroutine lft (m,n,theta,cp,pb)

        integer :: k
        integer :: kdo
        integer :: m
        integer :: mmod
        integer :: n
        integer :: nmod
        real cp(*),pb,theta,cdt,sdt,cth,sth,chh
        cdt = cos(2.0*theta)
        sdt = sin(2.0*theta)
        nmod=mod(n,2)
        mmod=mod(m,2)
        if (nmod< 0) then
            goto 1
        else if (nmod == 0) then
            goto 1
        else
            goto 2
        end if
1       if (mmod< 0) then
            goto 3
        else if (mmod == 0) then
            goto 3
        else
            goto 4
        end if
        !
        !     n even, m even
        !
3       kdo=n/2
        pb = .5*cp(1)
        if (n == 0) return
        cth = cdt
        sth = sdt
        do 170 k=1,kdo
            !     pb = pb+cp(k+1)*cos(2*k*theta)
            pb = pb+cp(k+1)*cth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
170     continue
        return
        !
        !     n even, m odd
        !
4       kdo = n/2
        pb = 0.0
        cth = cdt
        sth = sdt
        do 180 k=1,kdo
            !     pb = pb+cp(k)*sin(2*k*theta)
            pb = pb+cp(k)*sth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
180     continue
        return
2       if (mmod< 0) then
            goto 13
        else if (mmod == 0) then
            goto 13
        else
            goto 14
        end if
        !
        !     n odd, m even
        !
13      kdo = (n+1)/2
        pb = 0.0
        cth = cos(theta)
        sth = sin(theta)
        do 190 k=1,kdo
            !     pb = pb+cp(k)*cos((2*k-1)*theta)
            pb = pb+cp(k)*cth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
190     continue
        return
        !
        !     n odd, m odd
        !
14      kdo = (n+1)/2
        pb = 0.0
        cth = cos(theta)
        sth = sin(theta)
        do 200 k=1,kdo
            !     pb = pb+cp(k)*sin((2*k-1)*theta)
            pb = pb+cp(k)*sth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
200     continue
        return
    end subroutine lft
    subroutine dlft (m,n,theta,cp,pb)

        integer :: k
        integer :: kdo
        integer :: m
        integer :: mmod
        integer :: n
        integer :: nmod
        !
        !     computes the derivative of pmn(theta) with respect to theta
        !
        dimension cp(1)
        real cp,pb,theta,cdt,sdt,cth,sth,chh
        cdt = cos(2.0*theta)
        sdt = sin(2.0*theta)
        nmod=mod(n,2)
        mmod=mod(abs(m),2)
        if (nmod< 0) then
            goto 1
        else if (nmod == 0) then
            goto 1
        else
            goto 2
        end if
1       if (mmod< 0) then
            goto 3
        else if (mmod == 0) then
            goto 3
        else
            goto 4
        end if
        !
        !     n even, m even
        !
3       kdo=n/2
        pb = 0.0
        if (n == 0) return
        cth = cdt
        sth = sdt
        do 170 k=1,kdo
            !     pb = pb+cp(k+1)*cos(2*k*theta)
            pb = pb-2.0*k*cp(k+1)*sth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
170     continue
        return
        !
        !     n even, m odd
        !
4       kdo = n/2
        pb = 0.0
        cth = cdt
        sth = sdt
        do 180 k=1,kdo
            !     pb = pb+cp(k)*sin(2*k*theta)
            pb = pb+2.0*k*cp(k)*cth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
180     continue
        return
2       if (mmod< 0) then
            goto 13
        else if (mmod == 0) then
            goto 13
        else
            goto 14
        end if
        !
        !     n odd, m even
        !
13      kdo = (n+1)/2
        pb = 0.0
        cth = cos(theta)
        sth = sin(theta)
        do 190 k=1,kdo
            !     pb = pb+cp(k)*cos((2*k-1)*theta)
            pb = pb-(2.0*k-1)*cp(k)*sth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
190     continue
        return
        !
        !     n odd, m odd
        !
14      kdo = (n+1)/2
        pb = 0.0
        cth = cos(theta)
        sth = sin(theta)
        do 200 k=1,kdo
            !     pb = pb+cp(k)*sin((2*k-1)*theta)
            pb = pb+(2.0*k-1)*cp(k)*cth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
200     continue

    end subroutine dlft

    subroutine dlfcz(n,cp)

        integer :: i
        integer :: ic
        integer :: j
        integer :: n
        integer :: ncp
        !
        !     computes the fourier coefficients of the legendre
        !     polynomials. n is the degree and integer(n/2+1)
        !     coefficients are returned in cp
        !
        real cp(n/2+1)
        real cn,t1,t2,t3,t4,coef,fi
        !
        cn = 2.0
        write( stdout, 9) cn
9       format(' check1 on dble cn ',1pd20.011)
        if (n>0) then
            ic = 0
            fi = 0.0
            do i=1,n
                fi = fi+2.0
                cn = (1.0-1.0/fi**2)*cn
                if (abs(cn) > 5.0 .and. ic == 0) then
                    ic = 1
                    write( stdout, 7) i,cn
7                   format('  i ',i7,' check3 on cn',1pd15.6)
                end if
            end do
        end if
        write( stdout, 8) cn
8       format(' check2 on dble cn ',1pd20.011)
        cn = sqrt(cn)
        ncp = n/2+1
        t1 = -1.0
        t2 = n+1.0
        t3 = 0.0
        t4 = n+n+1.0
        cp(ncp) = cn
        coef = 1.0
        write( stdout, 11) cn
11      format(' check on dble cn ',1pd20.011)
        !      do j = ncp-1,1,-1
        j = ncp
10      j = j-1
        t1 = t1+2.0
        t2 = t2-1.0
        t3 = t3+1.0
        t4 = t4-2.0
        coef = (t1*t2)/(t3*t4)*coef
        cp(j) = coef*cn
        if (j>1) goto 10
        !      end do
        return
    end subroutine dlfcz
    !
    subroutine sgaqd(nlat,theta,wts,w,lwork,ierror)

        real :: cmax
        real :: cz
        real :: dcor
        real :: dpb
        real :: dthalf
        real :: dtheta
        real :: eps
        integer :: i
        integer :: idx
        integer :: ierror
        integer :: it
        integer :: itmax
        integer :: lwork
        integer :: mnlat
        integer :: nhalf
        integer :: nix
        integer :: nlat
        integer :: ns2
        real :: pb
        real :: pi
        real :: HALF_PI
        real :: sgnd
        real :: theta
        real :: w
        real :: wts
        real :: x
        real :: zero
        real :: zhold
        real :: zlast
        real :: zprev
        !
        !                             February 2002
        !
        !     gauss points and weights are computed using the fourier-newton
        !     described in "on computing the points and weights for 
        !     gauss-legendre quadrature", paul n. swarztrauber, siam journal 
        !     on scientific computing that has been accepted for publication.
        !     This routine is faster and more accurate than older program
        !     with the same name.
        !
        !     subroutine sgaqd computes the nlat gaussian colatitudes and
        !     weights in single precision. the colatitudes are in radians
        !     and lie in the interval (0,pi).
        !
        !     input parameters
        !
        !     nlat    the number of gaussian colatitudes in the interval (0,pi)
        !             (between the two poles).  nlat must be greater than zero.
        !
        !     w       unused variable that permits a simple exchange with the
        !             old routine with the same name in spherepack.
        !
        !     lwork   unused variable that permits a simple exchange with the
        !             old routine with the same name in spherepack.
        !
        !     output parameters
        !
        !     theta   a vector of length nlat containing the nlat gaussian
        !             colatitudes on the sphere in increasing radians
        !             in the interval (0,pi).
        !
        !     wts     a vector of length nlat containing the
        !                   nlat gaussian weights.
        !
        !     ierror = 0 no errors
        !            = 1 if nlat.le.0
        !
        !  *****************************************************************
        !
        dimension theta(nlat),wts(nlat), w(*)
        real :: summation
        !
        !     check work space length
        !
        ierror = 1
        if (nlat<=0) return
        ierror = 0
        !
        !     compute weights and points analytically when nlat=1,2
        !
        if (nlat==1) then
            theta(1) = acos(0.0)
            wts(1) = 2.0
            return
        end if
        if (nlat==2) then
            x = sqrt(1.0/3.0)
            theta(1) = acos(x)
            theta(2) = acos(-x)
            wts(1) = 1.0
            wts(2) = 1.0
            return
        end if
        eps = sqrt(epsilon(1.0))
        eps = eps*sqrt(eps)
        HALF_PI = acos(0.0)
        pi = HALF_PI+HALF_PI
        mnlat = mod(nlat,2)
        ns2 = nlat/2
        nhalf = (nlat+1)/2
        idx = ns2+2
        !
        call lfcz (nlat,cz,theta(ns2+1),wts(ns2+1))
        !
        dtheta = HALF_PI/nhalf
        dthalf = dtheta/2.0
        cmax = .2*dtheta
        !
        !     estimate first point next to theta = pi/2
        !
        if (mnlat/=0) then
            zprev = HALF_PI
            zero = HALF_PI-dtheta
            nix = nhalf-1
        else
            zero = HALF_PI-dthalf
            nix = nhalf
        end if
        itmax = 0
9       it = 0
10      it = it+1
        zlast = zero
        !
        !     newton iterations
        !
        call slpdp (nlat,zero,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
        dcor = pb/dpb
        sgnd = 1.0
        if (dcor /= 0.0) sgnd = dcor/abs(dcor)
        dcor = sgnd*min(abs(dcor),cmax)
        zero = zero-dcor
        if (abs(zero-zlast)>eps*abs(zero)) goto 10
        theta(nix) = zero
        zhold = zero
        !      wts(nix) = (nlat+nlat+1)/(dpb*dpb)
        !
        !     yakimiw's formula permits using old pb and dpb
        !
        wts(nix) = (nlat+nlat+1)/(dpb+pb*cos(zlast)/sin(zlast))**2
        nix = nix-1
        if (nix==0) goto 30
        if (nix==nhalf-1)  zero = 3.0*zero-pi
        if (nix<nhalf-1)  zero = zero+zero-zprev
        zprev = zhold
        goto 9
        !
        !     extend points and weights via symmetries
        !
30      if (mnlat/=0) then
            theta(nhalf) = HALF_PI
            call slpdp (nlat,HALF_PI,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
            wts(nhalf) = (nlat+nlat+1)/(dpb*dpb)
        end if
        do i=1,ns2
            wts(nlat-i+1) = wts(i)
            theta(nlat-i+1) = pi-theta(i)
        end do
        summation = 0.0
        do i=1,nlat
            summation = summation+wts(i)
        end do
        do i=1,nlat
            wts(i) = 2.0*wts(i)/summation
        end do
        return
    end subroutine sgaqd
    subroutine lfcz(n,cz,cp,dcp)

        real :: cp
        real :: cz
        real :: dcp
        integer :: j
        integer :: n
        integer :: ncp
        real :: t1
        real :: t2
        real :: t3
        real :: t4
        !
        !     computes the fourier coefficients of the legendre
        !     polynomial p_n^0 and its derivative.
        !     n is the degree and n/2 or (n+1)/2
        !     coefficients are returned in cp depending on whether
        !     n is even or odd. The same number of coefficients
        !     are returned in dcp. For n even the constant
        !     coefficient is returned in cz.
        !
        dimension cp((n+1)/2),dcp((n+1)/2)
        ncp = (n+1)/2
        t1 = -1.0
        t2 = n+1.0
        t3 = 0.0
        t4 = n+n+1.0
        if (mod(n,2)==0) then
            cp(ncp) = 1.0
            do j = ncp,2,-1
                t1 = t1+2.0
                t2 = t2-1.0
                t3 = t3+1.0
                t4 = t4-2.0
                cp(j-1) = (t1*t2)/(t3*t4)*cp(j)
            end do
            t1 = t1+2.0
            t2 = t2-1.0
            t3 = t3+1.0
            t4 = t4-2.0
            cz = (t1*t2)/(t3*t4)*cp(1)
            do j=1,ncp
                dcp(j) = (j+j)*cp(j)
            end do
        else
            cp(ncp) = 1.0
            do j = ncp-1,1,-1
                t1 = t1+2.0
                t2 = t2-1.0
                t3 = t3+1.0
                t4 = t4-2.0
                cp(j) = (t1*t2)/(t3*t4)*cp(j+1)
            end do
            do j=1,ncp
                dcp(j) = (j+j-1)*cp(j)
            end do
        end if

    end subroutine lfcz

    subroutine slpdp (n,theta,cz,cp,dcp,pb,dpb)

        real :: cdt
        real :: chh
        real :: cp
        real :: cost
        real :: cz
        real :: dcp
        real :: dpb
        real :: fn
        integer :: k
        integer :: kdo
        integer :: n
        real :: pb
        real :: sdt
        real :: sint
        real :: theta
        !
        !     computes pn(theta) and its derivative dpb(theta) with
        !                          respect to theta
        !
        dimension cp(n/2+1),dcp(n/2+1)
        !
        fn = real(n)
        cdt = cos(2.0*theta)
        sdt = sin(2.0*theta)
        if (mod(n,2) ==0) then
            !
            !     n even
            !
            kdo = n/2
            pb = .5*cz
            dpb = 0.0
            if (n > 0) then
                cost = cdt
                sint = sdt
                do k=1,kdo
                    !      pb = pb+cp(k)*cos(2*k*theta)
                    pb = pb+cp(k)*cost
                    !      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta)
                    dpb = dpb-dcp(k)*sint
                    chh = cdt*cost-sdt*sint
                    sint = sdt*cost+cdt*sint
                    cost = chh
                end do
            end if
        else
            !
            !     n odd
            !
            kdo = (n+1)/2
            pb = 0.0
            dpb = 0.0
            cost = cos(theta)
            sint = sin(theta)
            do k=1,kdo
                !      pb = pb+cp(k)*cos((2*k-1)*theta)
                pb = pb+cp(k)*cost
                !      dpb = dpb-(k+k-1)*cp(k)*sin((2*k-1)*theta)
                dpb = dpb-dcp(k)*sint
                chh = cdt*cost-sdt*sint
                sint = sdt*cost+cdt*sint
                cost = chh
            end do
        end if

    end subroutine slpdp

end program test_gaussian_latitudes_and_weights_routines
