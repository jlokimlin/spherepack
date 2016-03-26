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
!     This version of gaqd implements the method presented in:
!     P. N. swarztrauber, Computing the points and weights for 
!     Gauss-Legendre quadrature, SIAM J. Sci. Comput.,
!     24(2002) pp. 945-954.
!     
!     The w and lwork arrays are dummy and included only to
!     permit a simple pluggable exchange with the
!     old gaqd in previous versions of SPHEREPACK
!
!
!
subroutine gaqd(nlat,theta,wts,w,lwork,ierror)
    !  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    !
    !     gauss points and weights are computed using the fourier-newton
    !     described in "on computing the points and weights for 
    !     gauss-legendre quadrature", paul n. swarztrauber, siam journal 
    !     on scientific computing that has been accepted for publication.
    !     This routine is faster and more accurate than older program
    !     with the same name.
    !
    !     subroutine gaqd computes the nlat gaussian colatitudes and weights
    !     in real. the colatitudes are in radians and lie in the
    !     in the interval (0,pi).
    !
    !     input parameters
    !
    !     nlat    the number of gaussian colatitudes in the interval (0,pi)
    !             (between the two poles).  nlat must be greater than zero.
    !
    !     w       unused real variable that permits a simple
    !             exchange with the old routine with the same name
    !             in spherepack.
    !
    !     lwork   unused variable that permits a simple exchange with the
    !             old routine with the same name in spherepack.
    !
    !     output parameters
    !
    !     theta   a real array with length nlat
    !             containing the gaussian colatitudes in
    !             increasing radians on the interval (0,pi).
    !
    !     wts     a real array with lenght nlat
    !             containing the gaussian weights.
    !
    !     ierror = 0 no errors
    !            = 1 if nlat.le.0
    !
    !  *****************************************************************
    !
    real theta(nlat),wts(nlat),w, &
        x,pi,pis2,dtheta,dthalf,cmax,zprev,zlast,zero, &
        zhold,pb,dpb,dcor,summation,cz
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
    pis2 = 2.0*atan(1.0)
    pi = pis2+pis2
    mnlat = mod(nlat,2)
    ns2 = nlat/2
    nhalf = (nlat+1)/2
    idx = ns2+2
    !
    call cpdp (nlat,cz,theta(ns2+1),wts(ns2+1))
    !
    dtheta = pis2/nhalf
    dthalf = dtheta/2.0
    cmax = 0.2*dtheta
    !
    !     estimate first point next to theta = pi/2
    !
    if(mnlat/=0) then
        zero = pis2-dtheta
        zprev = pis2
        nix = nhalf-1
    else
        zero = pis2-dthalf
        nix = nhalf
    end if
9   it = 0
10  it = it+1
    zlast = zero
    !
    !     newton iterations
    !
    call tpdp (nlat,zero,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
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
    if(nix==nhalf-1)  then
        zero = 3.0*zero-pi
    end if
    if(nix<nhalf-1)  zero = zero+zero-zprev
    zprev = zhold
    go to 9
    !
    !     extend points and weights via symmetries
    !
30  if(mnlat/=0) then
        theta(nhalf) = pis2
        call tpdp (nlat,pis2,cz,theta(ns2+1),wts(ns2+1),pb,dpb)
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

end subroutine gaqd

subroutine cpdp(n,cz,cp,dcp)
    !
    !     computes the fourier coefficients of the legendre
    !     polynomial p_n^0 and its derivative.
    !     n is the degree and n/2 or (n+1)/2
    !     coefficients are returned in cp depending on whether
    !     n is even or odd. The same number of coefficients
    !     are returned in dcp. For n even the constant
    !     coefficient is returned in cz.
    !
    real cp(n/2+1),dcp(n/2+1), &
        t1,t2,t3,t4,cz
    ncp = (n+1)/2
    t1 = -1.0
    t2 = n+1.0
    t3 = 0.0
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

end subroutine cpdp


subroutine tpdp (n,theta,cz,cp,dcp,pb,dpb)
    !
    !     computes pn(theta) and its derivative dpb(theta) with
    !     respect to theta
    !
    real cp(n/2+1),dcp(n/2+1),cz, &
        pb,dpb,fn,theta,cdt,sdt,cth,sth,chh
    !
    fn = n
    cdt = cos(theta+theta)
    sdt = sin(theta+theta)
    if(mod(n,2) ==0) then
        !
        !     n even
        !
        kdo = n/2
        pb = .5d0*cz
        dpb = 0.0
        if(n > 0) then
            cth = cdt
            sth = sdt
            do k=1,kdo
                !      pb = pb+cp(k)*cos(2*k*theta)
                pb = pb+cp(k)*cth
                !      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta)
                dpb = dpb-dcp(k)*sth
                chh = cdt*cth-sdt*sth
                sth = sdt*cth+cdt*sth
                cth = chh
            end do
        end if
    else
        !
        !     n odd
        !
        kdo = (n+1)/2
        pb = 0.0
        dpb = 0.0
        cth = cos(theta)
        sth = sin(theta)
        do k=1,kdo
            pb = pb+cp(k)*cth
            dpb = dpb-dcp(k)*sth
            chh = cdt*cth-sdt*sth
            sth = sdt*cth+cdt*sth
            cth = chh
        end do
    end if

end subroutine tpdp

