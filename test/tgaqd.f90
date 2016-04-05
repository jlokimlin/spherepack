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
!     this program tests subroutine gaqd for computing 
!     the Gauss Legendre points and weights in file gaqd.f
!     It tests only the april 2002 version and not the
!     older version in file gaqd.old.f
!     gauss points and weights are computed using newtons method
!     with equally spaced points as first guess. Points are
!     computed as theta_i where x_i = cos(theta_i)
!     the method is discussed in "On computing the Gauss_Legendre
!     points and weights" by Paul N Swarztrauber, accepted for 
!     publication in the SIAM journal of scientific computing.
!                         April 2002
!
program testint
    parameter (nlat=63)
    real theta(nlat),wts(nlat),work(nlat+2), &
        dtheta(nlat),dwts(nlat),dw(nlat+1)
    real dwmx,tmax,dwmax,dtmax
    dimension stheta(nlat),swts(nlat),swork(nlat+1)
    dimension diff(nlat)
    real*4 t1(2),t2(2)
    real sumw
    !
    lwork = nlat+1
    hold = etime(t1)
    hold = t1(1)
    call gsqd(nlat,theta,wts,work,lwork,ierror)
    tdoub = etime(t1)
    tdoub = t1(1)-hold
    if(ierror /= 0) write(6,4) ierror
4   format(' ierror=',i5)
    write(*,739) tdoub
739 format(' tdoub ',1pe15.6)
    !
    !     compare double with double
    !
    ldw = nlat+2
    hold = etime(t1)
    hold = t1(1)
    call gaqd(nlat,dtheta,dwts,dw,ldw,ierror)
    tdoub = etime(t1)
    tdoub = t1(1)-hold
    if(ierror /= 0) write(6,30) ierror
30  format(' ierror=',i5)
    write(*,31) tdoub
31  format(' tdoub gaqd',1pe15.6)
    !
    dwmx = 0.0
    tmax = 0.0
    dwmax = 0.0
    dtmax = 0.0
    ido = (nlat+1)/2
    do 32 i=1,ido
        dtmax = max(dtmax,abs(theta(i)-dtheta(i)))
        dwmax = max(dwmax,abs(wts(i)-dwts(i)))
        tmax = max(tmax,abs(theta(i)))
        dwmx = max(dwmx,abs(wts(i)))
32  continue
    dtmax = dtmax/tmax
    dwmax = dwmax/dwmx
    write(*,33) nlat,dtmax,dwmax
33  format(' nlat',i6,'  points ',1pd15.6,' weights ',d15.6)
    !
    sumw = 0.
    do i=1,nlat
        sumw = sumw+wts(i)
    end do
    write(*,638) sumw
638 format('  sumw ',1pd39.30)
    !      if(0.eq.0) go to 671
    hold = etime(t2)
    hold = t2(1)
    lwork = nlat+2
    call sgaqd(nlat,stheta,swts,swork,lwork,ierror)
    tsing = etime(t2)
    tsing = t2(1)-hold
    if(ierror /= 0) write(6,5) ierror
5   format(' iserror=',i5)
    sums = 0.
    do i=1,nlat
        sums = sums+swts(i)
    end do
    write(*,636) sums
636 format('  sums ',1pe24.15)
    dmax = 0.
    wmax = 0.
    rmax = 0.
    ido = (nlat+1)/2
    do 6 i=1,ido
        diff(i) = wts(i)-swts(i)
        dmax = max(dmax,abs(diff(i)))
        wmax = max(wmax,abs(swts(i)))
        rerr = abs(diff(i)/swts(i))
        if(rerr>rmax) then
            rmax = rerr
            irm = i
        end if
        !      if(i.lt.121) write(*,27) wts(i),swts(i),rerr
27      format(' wts ',1pd15.6,' swts ',e15.6,' rpter ',e15.6)
6   continue
    !      write(*,7) (diff(i),i=nlat-25,nlat)
7   format(' diff in weights'/(1p8e10.3))
    dmax = dmax/wmax
    write(*,9) nlat,irm,dmax,rmax
9   format(' weights: nlat ',i6,' irele ',i6, &
        ' dmax ',1pe15.6,' rmax ',1pe15.6)
    dmax = 0.
    pmax = 0.
    rmax = 0.
    do 10 i=1,ido
        diff(i) = theta(i)-stheta(i)
        dmax = max(dmax,abs(diff(i)))
        pmax = max(pmax,abs(stheta(i)))
        rerr = abs(diff(i)/stheta(i))
        if(rerr>rmax) then
            rmax = rerr
            irm = i
        end if
10  continue
    !      write(*,11) (diff(i),i=nlat-25,nlat)
11  format(' diff in points'/(1p8e10.3))
    dmax = dmax/pmax
    write(*,12) nlat,irm,dmax,rmax
12  format(' points:  nlat ',i6,' irele ',i6, &
        ' dmax ',1pe15.6,' rmax ',1pe15.6)
    dmax = 0.
    do 110 i=1,nlat
        diff(i) = cos(theta(i))-cos(stheta(i))
        dmax = max(dmax,abs(diff(i)))
110 continue
    !      write(*,111) (diff(i),i=nlat-25,nlat)
111 format(' diff in points'/(1p8e10.3))
    write(*,112) dmax
112 format(' max difference in mu',1pe15.6)
    write(*,1) nlat,tsing,tdoub
1   format(' nlat',i6,' tsing',1pe15.6,' tdoub',e15.6)
    671 end program testint
    !
    !     subroutine gsqd is a single precision version of gaqd.
    !     gauss points and weights are computed using newtons method
    !     with equally spaced points as first guess. Points are
    !     computed as theta_i where x_i = cos(theta_i)
    !     the method is discussed in "On computing the Gauss_Legendre
    !     points and weights" by Paul N Swarztrauber, accepted for 
    !     publication in the SIAM journal of scientific computing.
    !                         April 2002
    !
    subroutine gsqd(nlat,theta,wts,dwork,ldwork,ierror)
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
        real pis2,x,theta,wts,dwork
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
            wts(1) = 2.d0
            return
        end if
        if (nlat==2) then
            x = sqrt(1.0/3.0)
            theta(1) = acos(x)
            theta(2) = acos(-x)
            wts(1) = 1.d0
            wts(2) = 1.d0
            return
        end if
        !
        !     compute points
        !
        call gsqd1(nlat,theta,dwork)
        !
        !     compute weights
        !
        call egwts(nlat,theta,wts,dwork)
        !
        !     extend points and weights via symmetries
        !
        pis2 = 4.0*atan(1.0)
        ns2 = nlat/2
        do i=1,ns2
            wts(nlat-i+1) = wts(i)
            theta(nlat-i+1) = pis2-theta(i)
        end do
        return
    end subroutine gsqd
    !
    subroutine gsqd1(nlat,theta,cp)
        real  theta((nlat+1)/2),cp(nlat/2+1)
        real pi,pis2,dtheta,dthalf, &
            cmax,dcor,pb,dpb,sgnd,zero,zlast
        !
        eps = sqrt(epsilon(1.0))
        eps = eps*sqrt(eps)
        pis2 = 2.0*atan(1.0)
        pi = pis2+pis2
        theta(1) = pis2
        if(nlat==1) go to 30
        ns2 = nlat/2
        nhalf = (nlat+1)/2
        !
        call dlfcz (nlat,cp)
        !
        !      check fourier-legendre coefficients
        !
        sum = 0.
        do i=1,ns2+1
            sum = sum+cp(i)
        end do
        sum = sum/ns2
        !      write(*,689)  sum
689     format(' check on dble f-l coefficients ',1pe15.6)
        !
        dtheta = pis2/nhalf
        dthalf = dtheta/2.0
        cmax = .2d0*dtheta
        if(mod(nlat,2)/=0) then
            theta(nhalf) = pis2
            zero = pis2-dtheta
            nix = nhalf-1
        else
            zero = pis2-dthalf
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
        if(dcor /= 0.0) sgnd = dcor/abs(dcor)
        !      write(*,2) nix,zero,theta(nix),dcor,cmax,sgnd
2       format(i7,1p5d15.6)
        dcor = sgnd*min(abs(dcor),cmax)
        zero = zero-dcor
        if(abs(zero-zlast)>eps*abs(zero)) go to 10
        theta(nix) = zero
        nix = nix-1
        if(nix==0) go to 30
        if(nix==nhalf-1)  zero = 3.0*zero-pi
        if(nix<nhalf-1)  zero = zero+zero-theta(nix+2)
        go to 9
30      return
    end subroutine gsqd1
    subroutine egwts(n,theta,wts,work)
        !
        !     computes gauss weights as described in swarztrauber
        !     and spotz, generalized harmonic transforms
        !
        real theta(n),wts(n),work(n+1)
        !
        call egwts1(n,theta,wts,work,work(n/2+2))
        return
    end subroutine egwts
    subroutine egwts1(n,theta,wts,dcp,cp)
        real theta((n+1)/2),wts((n+1)/2),cp((n-1)/2+1), &
            dcp(n/2+1),fn,sqnn,pb,dpb
        fn = n
        sqnn = sqrt((fn+fn-1)*(fn+fn+1))
        call dlfcz (n-1,cp)
        call dlfcz (n,dcp)
        nhalf = (n+1)/2
        do i=1,nhalf
            call lft (0,n-1,theta(i),cp,pb)
            call dlft (0,n,theta(i),dcp,dpb)
            wts(i) = -sqnn*sin(theta(i))/(fn*pb*dpb)
        end do
        return
    end subroutine egwts1
    subroutine lfc (m,n,cp)
        !
        real cp,fnum,fden,fnmh,a1,b1,c1,cp2,fnnp1,fnmsq,fk, &
            t1,t2,pm1,sc10,sc20,sc40
        dimension       cp(n/2+1)
        parameter (sc10=1024.d0)
        parameter (sc20=sc10*sc10)
        parameter (sc40=sc20*sc20)
        !
        cp(1) = 0.
        ma = iabs(m)
        if(ma > n) return
        if(n-1< 0) then
            goto 2
        else if(n-1 == 0) then 
            goto 3
        else 
            goto 5
        end if
2       cp(1) = sqrt(2.d0)
        return
3       if(ma /= 0) go to 4
        cp(1) = sqrt(1.5d0)
        return
4       cp(1) = sqrt(.75d0)
        if(m == -1) cp(1) = -cp(1)
        return
5       if(mod(n+ma,2) /= 0) go to 10
        nmms2 = (n-ma)/2
        fnum = n+ma+1
        fnmh = n-ma+1
        pm1 = 1.d0
        go to 15
10      nmms2 = (n-ma-1)/2
        fnum = n+ma+2
        fnmh = n-ma+2
        pm1 = -1.d0
        !      t1 = 1.
        !      t1 = 2.d0**(n-1)
        !      t1 = 1.d0/t1
15      t1 = 1.d0/sc20
        nex = 20
        fden = 2.d0
        if(nmms2 < 1) go to 20
        do 18 i=1,nmms2
            t1 = fnum*t1/fden
            if(t1 > sc20) then
                t1 = t1/sc40
                nex = nex+40
            end if
            fnum = fnum+2.
            fden = fden+2.
18      continue
20      t1 = t1/2.d0**(n-1-nex)
        if(mod(ma/2,2) /= 0) t1 = -t1
        t2 = 1.
        if(ma == 0) go to 26
        do 25 i=1,ma
            t2 = fnmh*t2/(fnmh+pm1)
            fnmh = fnmh+2.
25      continue
26      cp2 = t1*sqrt((n+.5d0)*t2)
        fnnp1 = n*(n+1)
        fnmsq = fnnp1-2.d0*ma*ma
        l = (n+1)/2
        if(mod(n,2) == 0 .and. mod(ma,2) == 0) l = l+1
        cp(l) = cp2
        if(m >= 0) go to 29
        if(mod(ma,2) /= 0) cp(l) = -cp(l)
29      if(l <= 1) return
        fk = n
        a1 = (fk-2.)*(fk-1.)-fnnp1
        b1 = 2.*(fk*fk-fnmsq)
        cp(l-1) = b1*cp(l)/a1
30      l = l-1
        if(l <= 1) return
        fk = fk-2.
        a1 = (fk-2.)*(fk-1.)-fnnp1
        b1 = -2.*(fk*fk-fnmsq)
        c1 = (fk+1.)*(fk+2.)-fnnp1
        cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
        go to 30
    end subroutine lfc
    subroutine lft (m,n,theta,cp,pb)
        real cp(*),pb,theta,cdt,sdt,cth,sth,chh
        cdt = cos(theta+theta)
        sdt = sin(theta+theta)
        nmod=mod(n,2)
        mmod=mod(m,2)
        if(nmod< 0) then
            goto 1
        else if(nmod == 0) then 
            goto 1
        else 
            goto 2
        end if
1       if(mmod< 0) then
            goto 3
        else if(mmod == 0) then 
            goto 3
        else 
            goto 4
        end if
        !
        !     n even, m even
        !
3       kdo=n/2
        pb = .5*cp(1)
        if(n == 0) return
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
        pb = 0.
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
2       if(mmod< 0) then
            goto 13
        else if(mmod == 0) then 
            goto 13
        else 
            goto 14
        end if
        !
        !     n odd, m even
        !
13      kdo = (n+1)/2
        pb = 0.
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
        pb = 0.
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
        !
        !     computes the derivative of pmn(theta) with respect to theta
        !
        dimension cp(1)
        real cp,pb,theta,cdt,sdt,cth,sth,chh
        cdt = cos(theta+theta)
        sdt = sin(theta+theta)
        nmod=mod(n,2)
        mmod=mod(abs(m),2)
        if(nmod< 0) then
            goto 1
        else if(nmod == 0) then 
            goto 1
        else 
            goto 2
        end if
1       if(mmod< 0) then
            goto 3
        else if(mmod == 0) then 
            goto 3
        else 
            goto 4
        end if
        !
        !     n even, m even
        !
3       kdo=n/2
        pb = 0.d0
        if(n == 0) return
        cth = cdt
        sth = sdt
        do 170 k=1,kdo
            !     pb = pb+cp(k+1)*cos(2*k*theta)
            pb = pb-2.d0*k*cp(k+1)*sth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
170     continue
        return
        !
        !     n even, m odd
        !
4       kdo = n/2
        pb = 0.
        cth = cdt
        sth = sdt
        do 180 k=1,kdo
            !     pb = pb+cp(k)*sin(2*k*theta)
            pb = pb+2.d0*k*cp(k)*cth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
180     continue
        return
2       if(mmod< 0) then
            goto 13
        else if(mmod == 0) then 
            goto 13
        else 
            goto 14
        end if
        !
        !     n odd, m even
        !
13      kdo = (n+1)/2
        pb = 0.
        cth = cos(theta)
        sth = sin(theta)
        do 190 k=1,kdo
            !     pb = pb+cp(k)*cos((2*k-1)*theta)
            pb = pb-(2.d0*k-1)*cp(k)*sth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
190     continue
        return
        !
        !     n odd, m odd
        !
14      kdo = (n+1)/2
        pb = 0.
        cth = cos(theta)
        sth = sin(theta)
        do 200 k=1,kdo
            !     pb = pb+cp(k)*sin((2*k-1)*theta)
            pb = pb+(2.d0*k-1)*cp(k)*cth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
200     continue
        return
    end subroutine dlft
    subroutine dlfcz(n,cp)
        !
        !     computes the fourier coefficients of the legendre
        !     polynomials. n is the degree and integer(n/2+1)
        !     coefficients are returned in cp
        !
        real cp(n/2+1)
        real cn,t1,t2,t3,t4,coef,fi
        !
        cn = 2.0
        write(*,9) cn
9       format(' check1 on dble cn ',1pd20.11)
        if(n>0) then
            ic = 0
            fi = 0.0
            do i=1,n
                fi = fi+2.0
                cn = (1.0-1.0/fi**2)*cn
                if(abs(cn)> 5.0.and.ic==0) then
                    ic = 1
                    write(*,7) i,cn
7                   format('  i ',i7,' check3 on cn',1pd15.6)
                end if
            end do
        end if
        write(*,8) cn
8       format(' check2 on dble cn ',1pd20.11)
        cn = sqrt(cn)
        ncp = n/2+1
        t1 = -1.0
        t2 = n+1.0
        t3 = 0.
        t4 = n+n+1.0
        cp(ncp) = cn
        coef = 1.0
        write(*,11) cn
11      format(' check on dble cn ',1pd20.11)
        !      do j = ncp-1,1,-1
        j = ncp
10      j = j-1
        t1 = t1+2.0
        t2 = t2-1.0
        t3 = t3+1.0
        t4 = t4-2.0
        coef = (t1*t2)/(t3*t4)*coef
        cp(j) = coef*cn
        if(j>1) go to 10
        !      end do
        return
    end subroutine dlfcz
    !
    subroutine sgaqd(nlat,theta,wts,w,lwork,ierror)
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
        dimension theta(nlat),wts(nlat)
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
        eps = sqrt(zeps(1.0))
        eps = eps*sqrt(eps)
        pis2 = 2.0*atan(1.0)
        pi = pis2+pis2
        mnlat = mod(nlat,2)
        ns2 = nlat/2
        nhalf = (nlat+1)/2
        idx = ns2+2
        !
        call lfcz (nlat,cz,theta(ns2+1),wts(ns2+1))
        !
        dtheta = pis2/nhalf
        dthalf = dtheta/2.0
        cmax = .2*dtheta
        !
        !     estimate first point next to theta = pi/2
        !
        if(mnlat/=0) then
            zprev = pis2
            zero = pis2-dtheta
            nix = nhalf-1
        else
            zero = pis2-dthalf
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
        if(dcor /= 0.0) sgnd = dcor/abs(dcor)
        dcor = sgnd*min(abs(dcor),cmax)
        zero = zero-dcor
        if(abs(zero-zlast)>eps*abs(zero)) go to 10
        theta(nix) = zero
        zhold = zero
        !      wts(nix) = (nlat+nlat+1)/(dpb*dpb)
        !
        !     yakimiw's formula permits using old pb and dpb
        !
        wts(nix) = (nlat+nlat+1)/(dpb+pb*cos(zlast)/sin(zlast))**2
        nix = nix-1
        if(nix==0) go to 30
        if(nix==nhalf-1)  zero = 3.0*zero-pi
        if(nix<nhalf-1)  zero = zero+zero-zprev
        zprev = zhold
        go to 9
        !
        !     extend points and weights via symmetries
        !
30      if(mnlat/=0) then
            theta(nhalf) = pis2
            call slpdp (nlat,pis2,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
            wts(nhalf) = (nlat+nlat+1)/(dpb*dpb)
        end if
        do i=1,ns2
            wts(nlat-i+1) = wts(i)
            theta(nlat-i+1) = pi-theta(i)
        end do
        sum = 0.
        do i=1,nlat
            sum = sum+wts(i)
        end do
        do i=1,nlat
            wts(i) = 2.0*wts(i)/sum
        end do
        return
    end subroutine sgaqd
    subroutine lfcz(n,cz,cp,dcp)
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
        t3 = 0.
        t4 = n+n+1.0
        if(mod(n,2)==0) then
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
        !      write(*,23) (cp(j),dcp(j),j=1,ncp)
23      format(' coefficients '/(1p4e15.6))
        return
    end subroutine lfcz
    subroutine slpdp (n,theta,cz,cp,dcp,pb,dpb)
        !
        !     computes pn(theta) and its derivative dpb(theta) with
        !                          respect to theta
        !
        dimension cp(n/2+1),dcp(n/2+1)
        !
        fn = real(n)
        cdt = cos(theta+theta)
        sdt = sin(theta+theta)
        if(mod(n,2) ==0) then
            !
            !     n even
            !
            kdo = n/2
            pb = .5*cz
            dpb = 0.0
            if(n > 0) then
                cth = cdt
                sth = sdt
                do 170 k=1,kdo
                    !      pb = pb+cp(k)*cos(2*k*theta)
                    pb = pb+cp(k)*cth
                    !      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta)
                    dpb = dpb-dcp(k)*sth
                    chh = cdt*cth-sdt*sth
                    sth = sdt*cth+cdt*sth
                    cth = chh
170             continue
            end if
        else
            !
            !     n odd
            !
            kdo = (n+1)/2
            pb = 0.
            dpb = 0.
            cth = cos(theta)
            sth = sin(theta)
            do 190 k=1,kdo
                !      pb = pb+cp(k)*cos((2*k-1)*theta)
                pb = pb+cp(k)*cth
                !      dpb = dpb-(k+k-1)*cp(k)*sin((2*k-1)*theta)
                dpb = dpb-dcp(k)*sth
                chh = cdt*cth-sdt*sth
                sth = sdt*cth+cdt*sth
                cth = chh
190         continue
        end if
        return
    end subroutine slpdp
    real function zeps (x)
        !
        !     estimate unit roundoff in quantities of size x.
        !
        !
        !     this program should function properly on all systems
        !     satisfying the following two assumptions,
        !        1.  the base used in representing floating point
        !            numbers is not a power of three.
        !        2.  the quantity  a  in statement 10 is represented to
        !            the accuracy used in floating point variables
        !            that are stored in memory.
        !     the statement number 10 and the go to 10 are intended to
        !     force optimizing compilers to generate code satisfying
        !     assumption 2.
        !     under these assumptions, it should be true that,
        !            a  is not exactly equal to four-thirds,
        !            b  has a zero for its last bit or digit,
        !            c  is not exactly equal to one,
        !            eps  measures the separation of 1.0 from
        !                 the next larger floating point number.
        !     the developers of eispack would appreciate being informed
        !     about any systems where these assumptions do not hold.
        !
        !     this version dated 4/6/83.
        !
        a = 4.0/3.0
10      b = a - 1.0
        c = b + b + b
        eps = abs(c-1.0)
        if (eps == 0.0) go to 10
        zeps = eps*abs(x)
        return
    end function zeps
!
