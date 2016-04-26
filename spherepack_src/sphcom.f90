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
! ... file sphcom.f
!
!     this file must be loaded with all main program files
!     in spherepack.  it includes undocumented subroutines
!     called by some or all of main programs
!
pure subroutine dnlfk(m, n, cp)
    !
    !     cp requires n/2+1 real locations
    !
    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)  :: m
    integer (ip), intent (in)  :: n
    real (wp),    intent (out) :: cp(1)
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip)         :: i, l, ma, nex,  nmms2
    real (wp)            :: a1, b1, c1, t1, t2
    real (wp)            :: fk, cp2, pm1
    real (wp)            :: fden, fnmh, fnum, fnnp1, fnmsq
    real (wp), parameter :: SC10=1024.0_wp
    real (wp), parameter :: SC20=SC10**2
    real (wp), parameter :: SC40=SC20**2
    !----------------------------------------------------------------------

    ma = abs(m)

    if (ma > n) then
        cp(1) = 0.0_wp
        return
    end if

    if (n < 1) then
        cp(1) = sqrt(2.0_wp)
    else if (n == 1) then
        if (ma /= 0) then
            cp(1) = sqrt(0.75_wp)
            if (m == -1) then
                cp(1) = -cp(1)
            end if
        else
            cp(1) = sqrt(1.5_wp)
        end if
    else
        if (mod(n+ma, 2) /= 0) then
            nmms2 = (n-ma-1)/2
            fnum = real(n+ma+2, kind=wp)
            fnmh = real(n-ma+2, kind=wp)
            pm1 = -1.0_wp
        else
            nmms2 = (n-ma)/2
            fnum = real(n+ma+1, kind=wp)
            fnmh = real(n-ma+1, kind=wp)
            pm1 = 1.0_wp
        end if

        t1 = 1.0_wp/SC20
        nex = 20
        fden = 2.0_wp
        if (.not.(nmms2 < 1)) then
            do i=1, nmms2
                t1 = fnum*t1/fden
                if (t1 > SC20) then
                    t1 = t1/SC40
                    nex = nex+40
                end if
                fnum = fnum+2.0_wp
                fden = fden+2.0_wp
            end do
        end if

        t1 = t1/2.0**(n-1-nex)

        if (mod(ma/2, 2) /= 0) then
            t1 = -t1
        end if

        t2 = 1.0_wp

        if (.not.(ma == 0)) then
            do i=1, ma
                t2 = fnmh*t2/(fnmh+pm1)
                fnmh = fnmh+2.0_wp
            end do
        end if

        cp2 = t1*sqrt((real(n, kind=wp)+0.5_wp)*t2)
        fnnp1 = real(n*(n+1), kind=wp)
        fnmsq = fnnp1 - 2.0_wp * real(ma**2, kind=wp)
        l = (n+1)/2

        if (mod(n, 2) == 0 .and. mod(ma, 2) == 0) then
            l = l+1
        end if

        cp(l) = cp2

        if (.not.(m >= 0)) then
            if (mod(ma, 2) /= 0) then
                cp(l) = -cp(l)
            end if
        end if

        if (l <= 1) then
            return
        end if

        fk = real(n, kind=wp)
        a1 = (fk-2.0_wp)*(fk-1.0_wp)-fnnp1
        b1 = 2.0_wp*(fk*fk-fnmsq)
        cp(l-1) = b1*cp(l)/a1

        l = l - 1

        do while (.not.(l <= 1))
            fk = fk-2.0_wp
            a1 = (fk-2.0_wp)*(fk-1.0_wp)-fnnp1
            b1 = -2.0_wp*(fk*fk-fnmsq)
            c1 = (fk+1.0_wp)*(fk+2.0_wp)-fnnp1
            cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
            l=l-1
        end do
    end if

end subroutine dnlfk



pure subroutine dnlft(m, n, theta, cp, pb)

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32

    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer (ip), intent (in)  :: m
    integer (ip), intent (in)  :: n
    real (wp),    intent (in)  :: theta
    real (wp),    intent (out) :: cp(*)
    real (wp),    intent (out) :: pb
    !----------------------------------------------------------------------
    ! Dictionary: local variables
    !----------------------------------------------------------------------
    integer (ip) ::  k, kdo, mmod, nmod
    real (wp)    :: chh, cdt, cth, sdt, sth
    !----------------------------------------------------------------------

    cdt = cos(2.0_wp * theta)
    sdt = sin(2.0_wp * theta)

    nmod = mod(n, 2)
    mmod = mod(m, 2)

    if (nmod <= 0) then
        if (mmod <= 0) then
            !
            !==>  n even, m even
            !
            kdo = n/2
            pb = 0.5_wp * cp(1)

            if (n == 0) then
                return
            end if

            cth = cdt
            sth = sdt

            do k=1, kdo
                pb = pb+cp(k+1)*cth
                chh = cdt*cth-sdt*sth
                sth = sdt*cth+cdt*sth
                cth = chh
            end do
        else
            !
            !==> n even, m odd
            !
            kdo = n/2
            pb = 0.0_wp
            cth = cdt
            sth = sdt

            do k=1, kdo
                pb = pb+cp(k)*sth
                chh = cdt*cth-sdt*sth
                sth = sdt*cth+cdt*sth
                cth = chh
            end do
        end if
    else
        if (mmod <= 0) then
            !
            !==> n odd, m even
            !
            kdo = (n+1)/2
            pb = 0.0_wp
            cth = cos(theta)
            sth = sin(theta)

            do k=1, kdo
                pb = pb+cp(k)*cth
                chh = cdt*cth-sdt*sth
                sth = sdt*cth+cdt*sth
                cth = chh
            end do
        else
            !
            !==>  n odd, m odd
            !
            kdo = (n+1)/2
            pb = 0.0_wp
            cth = cos(theta)
            sth = sin(theta)

            do k=1, kdo
                pb = pb+cp(k)*sth
                chh = cdt*cth-sdt*sth
                sth = sdt*cth+cdt*sth
                cth = chh
            end do
        end if
    end if

end subroutine dnlft



subroutine dnlftd (m, n, theta, cp, pb)
!
!     computes the derivative of pmn(theta) with respect to theta
!
dimension cp(1)
real cp, pb, theta, cdt, sdt, cth, sth, chh
cdt = cos(2.0*theta)
sdt = sin(2.0*theta)
nmod=mod(n, 2)
mmod=mod(abs(m), 2)
if (nmod< 0) then
    goto 1
else if (nmod == 0) then
    goto 1
else
    goto 2
end if
1 if (mmod< 0) then
    goto 3
else if (mmod == 0) then
    goto 3
else
    goto 4
end if
!
!     n even, m even
!
3 kdo=n/2
pb = 0.0
if (n == 0) return
cth = cdt
sth = sdt
do 170 k=1, kdo
!     pb = pb+cp(k+1)*cos(2*k*theta)
pb = pb-2.0*k*cp(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
170 continue
return
!
!     n even, m odd
!
4 kdo = n/2
pb = 0.
cth = cdt
sth = sdt
do 180 k=1, kdo
!     pb = pb+cp(k)*sin(2*k*theta)
pb = pb+2.0*k*cp(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
180 continue
return
2 if (mmod< 0) then
    goto 13
else if (mmod == 0) then
    goto 13
else
    goto 14
end if
!
!     n odd, m even
!
13 kdo = (n+1)/2
pb = 0.
cth = cos(theta)
sth = sin(theta)
do 190 k=1, kdo
!     pb = pb+cp(k)*cos((2*k-1)*theta)
pb = pb-(2.0*k-1)*cp(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
190 continue
return
!
!     n odd, m odd
!
14 kdo = (n+1)/2
pb = 0.
cth = cos(theta)
sth = sin(theta)
do 200 k=1, kdo
!     pb = pb+cp(k)*sin((2*k-1)*theta)
pb = pb+(2.0*k-1)*cp(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
200 continue

end subroutine dnlftd



subroutine legin(mode, l, nlat, m, w, pmn, km)
!     this subroutine computes legendre polynomials for n=m, ..., l-1
!     and  i=1, ..., late (late=((nlat+mod(nlat, 2))/2)gaussian grid
!     in pmn(n+1, i, km) using swarztrauber's recursion formula.
!     the vector w contains quantities precomputed in shigc.
!     legin must be called in the order m=0, 1, ..., l-1
!     (e.g., if m=10 is sought it must be preceded by calls with
!     m=0, 1, 2, ..., 9 in that order)
dimension w(1), pmn(1)
!     set size of pole to equator gaussian grid
late = (nlat+mod(nlat, 2))/2
!     partition w (set pointers for p0n, p1n, abel, bbel, cbel, pmn)
i1 = 1+nlat
i2 = i1+nlat*late
i3 = i2+nlat*late
i4 = i3+(2*nlat-l)*(l-1)/2
i5 = i4+(2*nlat-l)*(l-1)/2
call legin1(mode, l, nlat, late, m, w(i1), w(i2), w(i3), w(i4), &
            w(i5), pmn, km)

end subroutine legin



subroutine legin1(mode, l, nlat, late, m, p0n, p1n, abel, bbel, cbel, &
                  pmn, km)
dimension p0n(nlat, late), p1n(nlat, late)
dimension abel(1), bbel(1), cbel(1), pmn(nlat, late, 3)
data km0, km1, km2/ 1, 2, 3/
save km0, km1, km2
!     define index function used in storing triangular
!     arrays for recursion coefficients (functions of (m, n))
!     for 2.le.m.le.n-1 and 2.le.n.le.l-1
indx(m, n) = (n-1)*(n-2)/2+m-1
!     for l.le.n.le.nlat and 2.le.m.le.l
imndx(m, n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1

!     set do loop indices for full or half sphere
ms = m+1
ninc = 1
select case (mode)
	case (1)
		!     only compute pmn for n-m odd
		ms = m+2
		ninc = 2
	case (2)
		!     only compute pmn for n-m even
		ms = m+1
		ninc = 2
end select


if (m>1) then
do 100 np1=ms, nlat, ninc
n = np1-1
imn = indx(m, n)
if (n >= l) imn = imndx(m, n)
do 100 i=1, late
pmn(np1, i, km0) = abel(imn)*pmn(n-1, i, km2) &
            +bbel(imn)*pmn(n-1, i, km0) &
            -cbel(imn)*pmn(np1, i, km2)
100 continue

else if (m==0) then
do 101 np1=ms, nlat, ninc
do 101 i=1, late
pmn(np1, i, km0) = p0n(np1, i)
101 continue

else if (m==1) then
do 102 np1=ms, nlat, ninc
do 102 i=1, late
pmn(np1, i, km0) = p1n(np1, i)
102 continue
end if

!     permute column indices
!     km0, km1, km2 store m, m-1, m-2 columns
kmt = km0
km0 = km2
km2 = km1
km1 = kmt
!     set current m index in output param km
km = kmt

end subroutine legin1



subroutine zfin (isym, nlat, nlon, m, z, i3, wzfin)
    implicit none
    !----------------------------------------------------------------------
    ! Dictionary: calling arguments
    !----------------------------------------------------------------------
    integer, intent (in)     :: isym
    integer, intent (in)     :: nlat
    integer, intent (in)     :: nlon
    integer, intent (in)     :: m
    real,    intent (inout)  :: z(1)
    integer, intent (in)     :: i3
    real,    intent (in out) :: wzfin(1)
    !----------------------------------------------------------------------

    associate( imid => (nlat+1)/2 )
        associate( mmax => min(nlat, nlon/2+1) )
            associate( lim => nlat*imid )
                associate( labc => ((mmax-2)*(nlat+nlat-mmax-1))/2 )
                    associate( iw1 => lim+1 )
                        associate( iw2 => iw1+lim )
                            associate( iw3 => iw2+labc )
                                associate( iw4 => iw3+labc )
                                    !
                                    !     the length of wzfin is 2*lim+3*labc
                                    !
                                    call zfin1(isym, nlat, m, z, imid, i3, &
                                        wzfin, wzfin(iw1), wzfin(iw2), &
                                        wzfin(iw3), wzfin(iw4))
                                end associate
                            end associate
                        end associate
                    end associate
                end associate
            end associate
        end associate
    end associate

end subroutine zfin

subroutine zfin1 (isym, nlat, m, z, imid, i3, zz, z1, a, b, c)
dimension       z(imid, nlat, 3), zz(imid, 1), z1(imid, 1), &
                a(1), b(1), c(1)
save i1, i2
ihold = i1
i1 = i2
i2 = i3
i3 = ihold
if (m-1< 0) then
    goto 2225
else if (m-1 == 0) then
    goto 30
else
    goto 35
end if
2225 i1 = 1
i2 = 2
i3 = 3
do 45 np1=1, nlat
do 45 i=1, imid
z(i, np1, i3) = zz(i, np1)
45 continue
return
30 do 50 np1=2, nlat
do 50 i=1, imid
z(i, np1, i3) = z1(i, np1)
50 continue
return
35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
if (isym == 1) go to 36
do 85 i=1, imid
z(i, m+1, i3) = a(ns)*z(i, m-1, i1)-c(ns)*z(i, m+1, i1)
85 continue
36 if (m == nlat-1) return
if (isym == 2) go to 71
ns = ns+1
do 70 i=1, imid
z(i, m+2, i3) = a(ns)*z(i, m, i1)-c(ns)*z(i, m+2, i1)
70 continue
71 nstrt = m+3
if (isym == 1) nstrt = m+4
if (nstrt > nlat) go to 80
nstp = 2
if (isym == 0) nstp = 1
do 75 np1=nstrt, nlat, nstp
ns = ns+nstp
do 75 i=1, imid
z(i, np1, i3) = a(ns)*z(i, np1-2, i1)+b(ns)*z(i, np1-2, i3) &
                              -c(ns)*z(i, np1, i1)
75 continue
80 return

end subroutine zfin1



subroutine zfinit (nlat, nlon, wzfin, dwork)
dimension       wzfin(*)
real dwork(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     the length of wzfin is 3*((l-3)*l+2)/2 + 2*l*imid
!     the length of dwork is nlat+2
!
call zfini1 (nlat, nlon, imid, wzfin, wzfin(iw1), dwork, &
                                       dwork(nlat/2+1))

end subroutine zfinit



subroutine zfini1 (nlat, nlon, imid, z, abc, cz, work)
!
!     abc must have 3*((mmax-2)*(nlat+nlat-mmax-1))/2 locations
!     where mmax = min(nlat, nlon/2+1)
!     cz and work must each have nlat+1 locations
!
dimension z(imid, nlat, 2), abc(1)
real pi, dt, th, zh, cz(*), work(*)
pi = acos(-1.0)
dt = pi/(nlat-1)
do mp1=1, 2
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dnzfk(nlat, m, n, cz, work)
do i=1, imid
th = (i-1)*dt
call dnzft(nlat, m, n, th, cz, zh)
z(i, np1, mp1) = zh
end do
z(1, np1, mp1) = .5*z(1, np1, mp1)
end do
end do

call rabcp(nlat, nlon, abc)

end subroutine zfini1


subroutine dnzfk(nlat, m, n, cz, work)
!
!     dnzfk computes the coefficients in the trigonometric
!     expansion of the z functions that are used in spherical
!     harmonic analysis.
!
dimension  cz(1), work(1)
!
!     cz and work must both have nlat/2+1 locations
!
real sum, sc1, t1, t2, work, cz
lc = (nlat+1)/2
sc1 = 2.0/real(nlat-1)
call dnlfk(m, n, work)
nmod = mod(n, 2)
mmod = mod(m, 2)
if (nmod< 0) then
    goto 1
else if (nmod == 0) then
    goto 1
else
    goto 2
end if
1 if (mmod< 0) then
    goto 3
else if (mmod == 0) then
    goto 3
else
    goto 4
end if
!
!     n even, m even
!
3 kdo = n/2+1
do 5 idx=1, lc
i = idx+idx-2
sum = work(1)/(1.0-real(i**2))
if (kdo<2) go to 29
do 6 kp1=2, kdo
k = kp1-1
t1 = 1.0-(k+k+i)**2
t2 = 1.0-(k+k-i)**2
8 sum = sum+work(kp1)*(t1+t2)/(t1*t2)
6 continue
29 cz(idx) = sc1*sum
5 continue
return
!
!     n even, m odd
!
4 kdo = n/2
do 9 idx=1, lc
i = idx+idx-2
sum = 0.
do 101 k=1, kdo
t1 = 1.0-real(2*k+i)**2
t2 = 1.0-real(2*k-i)**2
12 sum=sum+work(k)*(t1-t2)/(t1*t2)
101 continue
cz(idx) = sc1*sum
9 continue
return
2 if (mmod< 0) then
    goto 13
else if (mmod == 0) then
    goto 13
else
    goto 14
end if
!
!     n odd, m even
!
13 kdo = (n+1)/2
do 15 idx=1, lc
i = idx+idx-1
sum = 0.
do 16 k=1, kdo
t1 = 1.0-(k+k-1+i)**2
t2 = 1.0-(k+k-1-i)**2
18 sum=sum+work(k)*(t1+t2)/(t1*t2)
16 continue
cz(idx)=sc1*sum
15 continue
return
!
!     n odd, m odd
!
14 kdo = (n+1)/2
do 19 idx=1, lc
i = idx+idx-3
sum=0.
do 20 k=1, kdo
t1 = 1.0-(k+k-1+i)**2
t2 = 1.0-(k+k-1-i)**2
22 sum=sum+work(k)*(t1-t2)/(t1*t2)
20 continue
cz(idx)=sc1*sum
19 continue

end subroutine dnzfk



subroutine dnzft(nlat, m, n, th, cz, zh)
dimension cz(1)
real cz, zh, th, cdt, sdt, cth, sth, chh
zh = 0.
cdt = cos(th+th)
sdt = sin(th+th)
lmod = mod(nlat, 2)
mmod = mod(m, 2)
nmod = mod(n, 2)
if (lmod< 0) then
    goto 20
else if (lmod == 0) then
    goto 20
else
    goto 10
end if
10 lc = (nlat+1)/2
lq = lc-1
ls = lc-2
if (nmod< 0) then
    goto 1
else if (nmod == 0) then
    goto 1
else
    goto 2
end if
1 if (mmod< 0) then
    goto 3
else if (mmod == 0) then
    goto 3
else
    goto 4
end if
!
!     nlat odd n even m even
!
3 zh = .5*(cz(1)+cz(lc)*cos(2*lq*th))
cth = cdt
sth = sdt
do 201 k=2, lq
!     zh = zh+cz(k)*cos(2*(k-1)*th)
zh = zh+cz(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
201 continue
return
!
!     nlat odd n even m odd
!
4 cth = cdt
sth = sdt
do 202 k=1, ls
!     zh = zh+cz(k+1)*sin(2*k*th)
zh = zh+cz(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
202 continue
return
!
!     nlat odd n odd, m even
!
2 if (mmod< 0) then
    goto 5
else if (mmod == 0) then
    goto 5
else
    goto 6
end if
5 cth = cos(th)
sth = sin(th)
do 203 k=1, lq
!     zh = zh+cz(k)*cos((2*k-1)*th)
zh = zh+cz(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
203 continue
return
!
!     nlat odd n odd m odd
!
6 cth = cos(th)
sth = sin(th)
do 204 k=1, lq
!     zh = zh+cz(k+1)*sin((2*k-1)*th)
zh = zh+cz(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
204 continue
return
20 lc = nlat/2
lq = lc-1
if (nmod< 0) then
    goto 30
else if (nmod == 0) then
    goto 30
else
    goto 80
end if
30 if (mmod< 0) then
    goto 40
else if (mmod == 0) then
    goto 40
else
    goto 60
end if
!
!     nlat even n even m even
!
40 zh = .5*cz(1)
cth = cdt
sth = sdt
do 50 k=2, lc
!     zh = zh+cz(k)*cos(2*(k-1)*th)
zh = zh+cz(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
50 continue
return
!
!     nlat even n even m odd
!
60 cth = cdt
sth = sdt
do 70 k=1, lq
!     zh = zh+cz(k+1)*sin(2*k*th)
zh = zh+cz(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
70 continue
return
!
!     nlat even n odd m even
!
80 if (mmod< 0) then
    goto 90
else if (mmod == 0) then
    goto 90
else
    goto 110
end if
90 zh = .5*cz(lc)*cos((nlat-1)*th)
cth = cos(th)
sth = sin(th)
do 100 k=1, lq
!     zh = zh+cz(k)*cos((2*k-1)*th)
zh = zh+cz(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
100 continue
return
!
!     nlat even n odd m odd
!
110 cth = cos(th)
sth = sin(th)
do 120 k=1, lq
!     zh = zh+cz(k+1)*sin((2*k-1)*th)
zh = zh+cz(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
120 continue

end subroutine dnzft



subroutine alin (isym, nlat, nlon, m, p, i3, walin)
    dimension       p(1)        , walin(1)
    imid = (nlat+1)/2
    lim = nlat*imid
    mmax = min(nlat, nlon/2+1)
    labc = ((mmax-2)*(nlat+nlat-mmax-1))/2
    iw1 = lim+1
    iw2 = iw1+lim
    iw3 = iw2+labc
    iw4 = iw3+labc
    !
    !     the length of walin is ((5*l-7)*l+6)/2
    !
    call alin1 (isym, nlat, m, p, imid, i3, walin, walin(iw1), walin(iw2), &
        walin(iw3), walin(iw4))

end subroutine alin



subroutine alin1(isym, nlat, m, p, imid, i3, pz, p1, a, b, c)
dimension       p(imid, nlat, 3), pz(imid, 1), p1(imid, 1), &
                a(1), b(1), c(1)
save i1, i2
ihold = i1
i1 = i2
i2 = i3
i3 = ihold
if (m-1< 0) then
    goto 25
else if (m-1 == 0) then
    goto 30
else
    goto 35
end if
25 i1 = 1
i2 = 2
i3 = 3
do 45 np1=1, nlat
do 45 i=1, imid
p(i, np1, i3) = pz(i, np1)
45 continue
return
30 do 50 np1=2, nlat
do 50 i=1, imid
p(i, np1, i3) = p1(i, np1)
50 continue
return
35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
if (isym == 1) go to 36
do 85 i=1, imid
p(i, m+1, i3) = a(ns)*p(i, m-1, i1)-c(ns)*p(i, m+1, i1)
85 continue
36 if (m == nlat-1) return
if (isym == 2) go to 71
ns = ns+1
do 70 i=1, imid
p(i, m+2, i3) = a(ns)*p(i, m, i1)-c(ns)*p(i, m+2, i1)
70 continue
71 nstrt = m+3
if (isym == 1) nstrt = m+4
if (nstrt > nlat) go to 80
nstp = 2
if (isym == 0) nstp = 1
do 75 np1=nstrt, nlat, nstp
ns = ns+nstp
do 75 i=1, imid
p(i, np1, i3) = a(ns)*p(i, np1-2, i1)+b(ns)*p(i, np1-2, i3) &
                              -c(ns)*p(i, np1, i1)
75 continue
80 return

end subroutine alin1


subroutine alinit (nlat, nlon, walin, dwork)
dimension       walin(*)
real dwork(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     the length of walin is 3*((l-3)*l+2)/2 + 2*l*imid
!     the length of work is nlat+1
!
call alini1 (nlat, nlon, imid, walin, walin(iw1), dwork)

end subroutine alinit



subroutine alini1 (nlat, nlon, imid, p, abc, cp)
dimension p(imid, nlat, 2), abc(1), cp(1)
real pi, dt, th, cp, ph
pi = acos(-1.0)
dt = pi/(nlat-1)
do 160 mp1=1, 2
m = mp1-1
do 160 np1=mp1, nlat
n = np1-1
call dnlfk (m, n, cp)
do 160 i=1, imid
th = (i-1)*dt
call dnlft (m, n, th, cp, ph)
p(i, np1, mp1) = ph
160 continue
call rabcp(nlat, nlon, abc)

end subroutine alini1



subroutine rabcp(nlat, nlon, abc)
!
!     subroutine rabcp computes the coefficients in the recurrence
!     relation for the associated legendre fuctions. array abc
!     must have 3*((mmax-2)*(nlat+nlat-mmax-1))/2 locations.
!
dimension abc(1)
mmax = min(nlat, nlon/2+1)
labc = ((mmax-2)*(nlat+nlat-mmax-1))/2
iw1 = labc+1
iw2 = iw1+labc
call rabcp1(nlat, nlon, abc, abc(iw1), abc(iw2))

end subroutine rabcp



subroutine rabcp1(nlat, nlon, a, b, c)
!
!     coefficients a, b, and c for computing pbar(m, n, theta) are
!     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
!
dimension a(1), b(1), c(1)
mmax = min(nlat, nlon/2+1)
do 215 mp1=3, mmax
m = mp1-1
ns = ((m-2)*(nlat+nlat-m-1))/2+1
fm = real(m)
tm = fm+fm
temp = tm*(tm-1.)
a(ns) = sqrt((tm+1.)*(tm-2.)/temp)
c(ns) = sqrt(2./temp)
if (m == nlat-1) go to 215
ns = ns+1
temp = tm*(tm+1.)
a(ns) = sqrt((tm+3.)*(tm-2.)/temp)
c(ns) = sqrt(6./temp)
mp3 = m+3
if (mp3 > nlat) go to 215
do 210 np1=mp3, nlat
n = np1-1
ns = ns+1
fn = real(n)
tn = fn+fn
cn = (tn+1.)/(tn-3.)
fnpm = fn+fm
fnmm = fn-fm
temp = fnpm*(fnpm-1.)
a(ns) = sqrt(cn*(fnpm-3.)*(fnpm-2.)/temp)
b(ns) = sqrt(cn*fnmm*(fnmm-1.)/temp)
c(ns) = sqrt((fnmm+1.)*(fnmm+2.)/temp)
210 continue
215 continue

end subroutine rabcp1



subroutine sea1(nlat, nlon, imid, z, idz, zin, wzfin, dwork)
dimension z(idz, *), zin(imid, nlat, 3), wzfin(*)
real dwork(*)
call zfinit(nlat, nlon, wzfin, dwork)
mmax = min(nlat, nlon/2+1)
do mp1=1, mmax
m = mp1-1
call zfin (0, nlat, nlon, m, zin, i3, wzfin)
do np1=mp1, nlat
mn = m*(nlat-1)-(m*(m-1))/2+np1
do i=1, imid
z(mn, i) = zin(i, np1, i3)
end do
end do
end do

end subroutine sea1



subroutine ses1(nlat, nlon, imid, p, pin, walin, dwork)
dimension p(imid, *), pin(imid, nlat, 3), walin(*)
real dwork(*)
call alinit (nlat, nlon, walin, dwork)
mmax = min(nlat, nlon/2+1)
do 10 mp1=1, mmax
m = mp1-1
call alin(0, nlat, nlon, m, pin, i3, walin)
do 10 np1=mp1, nlat
mn = m*(nlat-1)-(m*(m-1))/2+np1
do 10 i=1, imid
p(i, mn) = pin(i, np1, i3)
10 continue

end subroutine ses1



subroutine zvinit (nlat, nlon, wzvin, dwork)
dimension       wzvin(1)
real dwork(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     the length of wzvin is
!         2*nlat*imid +3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     the length of dwork is nlat+2
!
call zvini1 (nlat, nlon, imid, wzvin, wzvin(iw1), dwork, &
                                    dwork(nlat/2+2))

end subroutine zvinit



subroutine zvini1 (nlat, nlon, imid, zv, abc, czv, work)
!
!     abc must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     locations where mmax = min(nlat, (nlon+1)/2)
!     czv and work must each have nlat/2+1  locations
!
dimension zv(imid, nlat, 2), abc(1)
real dt, czv(1), zvh, th, work(1)
real, parameter :: pi = acos(-1.0)
dt = pi/(nlat-1)
mdo = min(2, nlat, (nlon+1)/2)
do mp1=1, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dzvk(nlat, m, n, czv, work)
do i=1, imid
th = (i-1)*dt
call dzvt(nlat, m, n, th, czv, zvh)
zv(i, np1, mp1) = zvh
end do
zv(1, np1, mp1) = .5*zv(1, np1, mp1)
end do
end do

call rabcv(nlat, nlon, abc)

end subroutine zvini1



subroutine zwinit (nlat, nlon, wzwin, dwork)
dimension       wzwin(1)
real dwork(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     the length of wzvin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of dwork is nlat+2
!
call zwini1 (nlat, nlon, imid, wzwin, wzwin(iw1), dwork, &
                                        dwork(nlat/2+2))

end subroutine zwinit



subroutine zwini1 (nlat, nlon, imid, zw, abc, czw, work)
!
!     abc must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     locations where mmax = min(nlat, (nlon+1)/2)
!     czw and work must each have nlat+1 locations
!
dimension zw(imid, nlat, 2), abc(1)
real  dt, czw(1), zwh, th, work(1)
real, parameter :: pi = acos(-1.0)
dt = pi/(nlat-1)
mdo = min(3, nlat, (nlon+1)/2)

if (mdo < 2) return

do mp1=2, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dzwk(nlat, m, n, czw, work)
do i=1, imid
th = (i-1)*dt
call dzwt(nlat, m, n, th, czw, zwh)
zw(i, np1, m) = zwh
end do
zw(1, np1, m) = .5*zw(1, np1, m)
end do
end do

call rabcw(nlat, nlon, abc)

end subroutine zwini1



subroutine zvin (ityp, nlat, nlon, m, zv, i3, wzvin)
dimension       zv(1)        , wzvin(1)
imid = (nlat+1)/2
lim = nlat*imid
mmax = min(nlat, (nlon+1)/2)
labc = (max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
iw1 = lim+1
iw2 = iw1+lim
iw3 = iw2+labc
iw4 = iw3+labc
!
!     the length of wzvin is 2*lim+3*labc
!
call zvin1 (ityp, nlat, m, zv, imid, i3, wzvin, wzvin(iw1), wzvin(iw2), &
            wzvin(iw3), wzvin(iw4))

end subroutine zvin



subroutine zvin1 (ityp, nlat, m, zv, imid, i3, zvz, zv1, a, b, c)
dimension       zv(imid, nlat, 3), zvz(imid, 1), zv1(imid, 1), &
                a(1), b(1), c(1)
save i1, i2
ihold = i1
i1 = i2
i2 = i3
i3 = ihold
if (m-1< 0) then
    goto 25
else if (m-1 == 0) then
    goto 30
else
    goto 35
end if
25 i1 = 1
i2 = 2
i3 = 3
do 45 np1=1, nlat
do 45 i=1, imid
zv(i, np1, i3) = zvz(i, np1)
45 continue
return
30 do 50 np1=2, nlat
do 50 i=1, imid
zv(i, np1, i3) = zv1(i, np1)
50 continue
return
35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
if (ityp == 1) go to 36
do 85 i=1, imid
zv(i, m+1, i3) = a(ns)*zv(i, m-1, i1)-c(ns)*zv(i, m+1, i1)
85 continue
36 if (m == nlat-1) return
if (ityp == 2) go to 71
ns = ns+1
do 70 i=1, imid
zv(i, m+2, i3) = a(ns)*zv(i, m, i1)-c(ns)*zv(i, m+2, i1)
70 continue
71 nstrt = m+3
if (ityp == 1) nstrt = m+4
if (nstrt > nlat) go to 80
nstp = 2
if (ityp == 0) nstp = 1
do 75 np1=nstrt, nlat, nstp
ns = ns+nstp
do 75 i=1, imid
zv(i, np1, i3) = a(ns)*zv(i, np1-2, i1)+b(ns)*zv(i, np1-2, i3) &
                              -c(ns)*zv(i, np1, i1)
75 continue
80 return

end subroutine zvin1



subroutine zwin (ityp, nlat, nlon, m, zw, i3, wzwin)
dimension       zw(1)        , wzwin(1)
imid = (nlat+1)/2
lim = nlat*imid
mmax = min(nlat, (nlon+1)/2)
labc = (max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
iw1 = lim+1
iw2 = iw1+lim
iw3 = iw2+labc
iw4 = iw3+labc
!
!     the length of wzwin is 2*lim+3*labc
!
call zwin1 (ityp, nlat, m, zw, imid, i3, wzwin, wzwin(iw1), wzwin(iw2), &
            wzwin(iw3), wzwin(iw4))

end subroutine zwin



subroutine zwin1 (ityp, nlat, m, zw, imid, i3, zw1, zw2, a, b, c)
dimension       zw(imid, nlat, 3), zw1(imid, 1), zw2(imid, 1), &
                a(1), b(1), c(1)
save i1, i2
ihold = i1
i1 = i2
i2 = i3
i3 = ihold
if (m-2< 0) then
    goto 25
else if (m-2 == 0) then
    goto 30
else
    goto 35
end if
25 i1 = 1
i2 = 2
i3 = 3
do 45 np1=2, nlat
do 45 i=1, imid
zw(i, np1, i3) = zw1(i, np1)
45 continue
return
30 do 50 np1=3, nlat
do 50 i=1, imid
zw(i, np1, i3) = zw2(i, np1)
50 continue
return
35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
if (ityp == 1) go to 36
do 85 i=1, imid
zw(i, m+1, i3) = a(ns)*zw(i, m-1, i1)-c(ns)*zw(i, m+1, i1)
85 continue
36 if (m == nlat-1) return
if (ityp == 2) go to 71
ns = ns+1
do 70 i=1, imid
zw(i, m+2, i3) = a(ns)*zw(i, m, i1)-c(ns)*zw(i, m+2, i1)
70 continue
71 nstrt = m+3
if (ityp == 1) nstrt = m+4
if (nstrt > nlat) go to 80
nstp = 2
if (ityp == 0) nstp = 1
do 75 np1=nstrt, nlat, nstp
ns = ns+nstp
do 75 i=1, imid
zw(i, np1, i3) = a(ns)*zw(i, np1-2, i1)+b(ns)*zw(i, np1-2, i3) &
                              -c(ns)*zw(i, np1, i1)
75 continue
80 return

end subroutine zwin1



subroutine vbinit (nlat, nlon, wvbin, dwork)
dimension wvbin(1)
real dwork(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of dwork is nlat+2
!
call vbini1 (nlat, nlon, imid, wvbin, wvbin(iw1), dwork, &
                                       dwork(nlat/2+2))

end subroutine vbinit



subroutine vbini1 (nlat, nlon, imid, vb, abc, cvb, work)
!
!     abc must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     locations where mmax = min(nlat, (nlon+1)/2)
!     cvb and work must each have nlat+1 locations
!
dimension vb(imid, nlat, 2), abc(1)
real dt, cvb(1), th, vbh, work(1)
real, parameter :: pi = acos(-1.0)
dt = pi/(nlat-1)
mdo = min(2, nlat, (nlon+1)/2)
do mp1=1, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dvbk(m, n, cvb, work)
do  i=1, imid
th = (i-1)*dt
call dvbt(m, n, th, cvb, vbh)
vb(i, np1, mp1) = vbh
end do
end do
end do

call rabcv(nlat, nlon, abc)

end subroutine vbini1



subroutine wbinit (nlat, nlon, wwbin, dwork)
dimension       wwbin(1)
real dwork(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of dwork is nlat+2
!
call wbini1 (nlat, nlon, imid, wwbin, wwbin(iw1), dwork, &
                                        dwork(nlat/2+2))

end subroutine wbinit



subroutine wbini1 (nlat, nlon, imid, wb, abc, cwb, work)
!
!     abc must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     locations where mmax = min(nlat, (nlon+1)/2)
!     cwb and work must each have nlat/2+1 locations
!
dimension wb(imid, nlat, 2), abc(1)
real dt, cwb(1), wbh, th, work(1)
real, parameter :: pi = acos(-1.0)

dt = pi/(nlat-1)
mdo = min(3, nlat, (nlon+1)/2)
if (mdo < 2) return
do mp1=2, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dwbk(m, n, cwb, work)
do i=1, imid
th = (i-1)*dt
call dwbt(m, n, th, cwb, wbh)
wb(i, np1, m) = wbh
end do
end do
end do

call rabcw(nlat, nlon, abc)

end subroutine wbini1



subroutine vbin (ityp, nlat, nlon, m, vb, i3, wvbin)
dimension       vb(1)        , wvbin(1)
imid = (nlat+1)/2
lim = nlat*imid
mmax = min(nlat, (nlon+1)/2)
labc = (max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
iw1 = lim+1
iw2 = iw1+lim
iw3 = iw2+labc
iw4 = iw3+labc
!
!     the length of wvbin is 2*lim+3*labc
!
call vbin1 (ityp, nlat, m, vb, imid, i3, wvbin, wvbin(iw1), wvbin(iw2), &
            wvbin(iw3), wvbin(iw4))

end subroutine vbin



subroutine vbin1 (ityp, nlat, m, vb, imid, i3, vbz, vb1, a, b, c)
dimension       vb(imid, nlat, 3), vbz(imid, 1), vb1(imid, 1), &
                a(1), b(1), c(1)
save i1, i2
ihold = i1
i1 = i2
i2 = i3
i3 = ihold
if (m-1< 0) then
    goto 25
else if (m-1 == 0) then
    goto 30
else
    goto 35
end if
25 i1 = 1
i2 = 2
i3 = 3
do 45 np1=1, nlat
do 45 i=1, imid
vb(i, np1, i3) = vbz(i, np1)
45 continue
return
30 do 50 np1=2, nlat
do 50 i=1, imid
vb(i, np1, i3) = vb1(i, np1)
50 continue
return
35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
if (ityp == 1) go to 36
do 85 i=1, imid
vb(i, m+1, i3) = a(ns)*vb(i, m-1, i1)-c(ns)*vb(i, m+1, i1)
85 continue
36 if (m == nlat-1) return
if (ityp == 2) go to 71
ns = ns+1
do 70 i=1, imid
vb(i, m+2, i3) = a(ns)*vb(i, m, i1)-c(ns)*vb(i, m+2, i1)
70 continue
71 nstrt = m+3
if (ityp == 1) nstrt = m+4
if (nstrt > nlat) go to 80
nstp = 2
if (ityp == 0) nstp = 1
do 75 np1=nstrt, nlat, nstp
ns = ns+nstp
do 75 i=1, imid
vb(i, np1, i3) = a(ns)*vb(i, np1-2, i1)+b(ns)*vb(i, np1-2, i3) &
                              -c(ns)*vb(i, np1, i1)
75 continue
80 return

end subroutine vbin1



subroutine wbin (ityp, nlat, nlon, m, wb, i3, wwbin)
dimension       wb(1)        , wwbin(1)
imid = (nlat+1)/2
lim = nlat*imid
mmax = min(nlat, (nlon+1)/2)
labc = (max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
iw1 = lim+1
iw2 = iw1+lim
iw3 = iw2+labc
iw4 = iw3+labc
!
!     the length of wwbin is 2*lim+3*labc
!
call wbin1 (ityp, nlat, m, wb, imid, i3, wwbin, wwbin(iw1), wwbin(iw2), &
            wwbin(iw3), wwbin(iw4))

end subroutine wbin



subroutine wbin1 (ityp, nlat, m, wb, imid, i3, wb1, wb2, a, b, c)
dimension       wb(imid, nlat, 3), wb1(imid, 1), wb2(imid, 1), &
                a(1), b(1), c(1)
save i1, i2
ihold = i1
i1 = i2
i2 = i3
i3 = ihold
if (m-2< 0) then
    goto 25
else if (m-2 == 0) then
    goto 30
else
    goto 35
end if
25 i1 = 1
i2 = 2
i3 = 3
do 45 np1=2, nlat
do 45 i=1, imid
wb(i, np1, i3) = wb1(i, np1)
45 continue
return
30 do 50 np1=3, nlat
do 50 i=1, imid
wb(i, np1, i3) = wb2(i, np1)
50 continue
return
35 ns = ((m-2)*(nlat+nlat-m-1))/2+1
if (ityp == 1) go to 36
do 85 i=1, imid
wb(i, m+1, i3) = a(ns)*wb(i, m-1, i1)-c(ns)*wb(i, m+1, i1)
85 continue
36 if (m == nlat-1) return
if (ityp == 2) go to 71
ns = ns+1
do 70 i=1, imid
wb(i, m+2, i3) = a(ns)*wb(i, m, i1)-c(ns)*wb(i, m+2, i1)
70 continue
71 nstrt = m+3
if (ityp == 1) nstrt = m+4
if (nstrt > nlat) go to 80
nstp = 2
if (ityp == 0) nstp = 1
do 75 np1=nstrt, nlat, nstp
ns = ns+nstp
do 75 i=1, imid
wb(i, np1, i3) = a(ns)*wb(i, np1-2, i1)+b(ns)*wb(i, np1-2, i3) &
                              -c(ns)*wb(i, np1, i1)
75 continue
80 return


end subroutine wbin1



subroutine dzvk(nlat, m, n, czv, work)
!
!     subroutine dzvk computes the coefficients in the trigonometric
!     expansion of the quadrature function zvbar(n, m, theta)
!
!     input parameters
!
!     nlat      the number of colatitudes including the poles.
!
!     n      the degree (subscript) of wbarv(n, m, theta)
!
!     m      the order (superscript) of wbarv(n, m, theta)
!
!     work   a work array with at least nlat/2+1 locations
!
!     output parameter
!
!     czv     the fourier coefficients of zvbar(n, m, theta).
!
dimension czv(1), work(1)
real czv, sc1, sum, work, t1, t2
if (n <= 0) return
lc = (nlat+1)/2
sc1 = 2.0/(nlat-1)
call dvbk(m, n, work, czv)
nmod = mod(n, 2)
mmod = mod(m, 2)
if (nmod /= 0) go to 1
if (mmod /= 0) go to 2
!
!     n even, m even
!
kdo = n/2
do 9 id=1, lc
i = id+id-2
sum = 0.
do 10 k=1, kdo
t1 = 1.0-(k+k+i)**2
t2 = 1.0-(k+k-i)**2
sum = sum+work(k)*(t1-t2)/(t1*t2)
10 continue
czv(id) = sc1*sum
9 continue
return
!
!     n even, m odd
!
2 kdo = n/2
do 5 id=1, lc
i = id+id-2
sum = 0.
do 6 k=1, kdo
t1 = 1.0-(k+k+i)**2
t2 = 1.0-(k+k-i)**2
sum = sum+work(k)*(t1+t2)/(t1*t2)
6 continue
czv(id) = sc1*sum
5 continue
return
1 if (mmod /= 0) go to 3
!
!     n odd, m even
!
kdo = (n+1)/2
do 19 id=1, lc
i = id+id-3
sum = 0.
do 20 k=1, kdo
t1 = 1.0-(k+k-1+i)**2
t2 = 1.0-(k+k-1-i)**2
sum = sum+work(k)*(t1-t2)/(t1*t2)
20 continue
czv(id) = sc1*sum
19 continue
return
!
!     n odd, m odd
!
3 kdo = (n+1)/2
do 15 id=1, lc
i = id+id-1
sum = 0.
do 16 k=1, kdo
t1 = 1.0-(k+k-1+i)**2
t2 = 1.0-(k+k-1-i)**2
sum = sum+work(k)*(t1+t2)/(t1*t2)
16 continue
czv(id) = sc1*sum
15 continue

end subroutine dzvk



subroutine dzvt(nlat, m, n, th, czv, zvh)
!
!     subroutine dzvt tabulates the function zvbar(n, m, theta)
!     at theta = th in real
!
!     input parameters
!
!     nlat      the number of colatitudes including the poles.
!
!     n      the degree (subscript) of zvbar(n, m, theta)
!
!     m      the order (superscript) of zvbar(n, m, theta)
!
!     czv     the fourier coefficients of zvbar(n, m, theta)
!             as computed by subroutine zwk.
!
!     output parameter
!
!     zvh     zvbar(m, n, theta) evaluated at theta = th
!
dimension czv(1)
real th, czv, zvh, cth, sth, cdt, sdt, chh
zvh = 0.
if (n <= 0) return
lc = (nlat+1)/2
lq = lc-1
ls = lc-2
cth = cos(th)
sth = sin(th)
cdt = cth*cth-sth*sth
sdt = 2.*sth*cth
lmod = mod(nlat, 2)
mmod = mod(m, 2)
nmod = mod(n, 2)
if (lmod == 0) go to 50
if (nmod /= 0) go to 1
cth = cdt
sth = sdt
if (mmod /= 0) go to 2
!
!     nlat odd  n even  m even
!
do 10 k=1, ls
zvh = zvh+czv(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
10 continue
return
!
!     nlat odd  n even  m odd
!
2 zvh = .5*czv(1)
do 20 k=2, lq
zvh = zvh+czv(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
20 continue
zvh = zvh+.5*czv(lc)*cos((nlat-1)*th)
return
1 if (mmod /= 0) go to 3
!
!     nlat odd  n odd  m even
!
do 30 k=1, lq
zvh = zvh+czv(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
30 continue
return
!
!     nlat odd  n odd  m odd
!
3 do 40 k=1, lq
zvh = zvh+czv(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
40 continue
return
50 if (nmod /= 0) go to 51
cth = cdt
sth = sdt
if (mmod /= 0) go to 52
!
!     nlat even  n even  m even
!
do 55 k=1, lq
zvh = zvh+czv(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
55 continue
return
!
!     nlat even  n even  m odd
!
52 zvh = .5*czv(1)
do 57 k=2, lc
zvh = zvh+czv(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
57 continue
return
51 if (mmod /= 0) go to 53
!
!     nlat even  n odd  m even
!
do 58 k=1, lq
zvh = zvh+czv(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
58 continue
return
!
!     nlat even  n odd  m odd
!
53 zvh = .5*czv(lc)*cos((nlat-1)*th)
do 60 k=1, lq
zvh = zvh+czv(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
60 continue

end subroutine dzvt



subroutine dzwk(nlat, m, n, czw, work)
!
!     subroutine dzwk computes the coefficients in the trigonometric
!     expansion of the quadrature function zwbar(n, m, theta)
!
!     input parameters
!
!     nlat      the number of colatitudes including the poles.
!
!     n      the degree (subscript) of zwbar(n, m, theta)
!
!     m      the order (superscript) of zwbar(n, m, theta)
!
!     work   a work array with at least nlat/2+1 locations
!
!     output parameter
!
!     czw     the fourier coefficients of zwbar(n, m, theta).
!
dimension czw(1), work(1)
real czw, work, sc1, sum, t1, t2
if (n <= 0) return
lc = (nlat+1)/2
sc1 = 2.0/(nlat-1)
call dwbk(m, n, work, czw)
nmod = mod(n, 2)
mmod = mod(m, 2)
if (nmod /= 0) go to 1
if (mmod /= 0) go to 2
!
!     n even, m even
!
kdo = n/2
do 19 id=1, lc
i = id+id-3
sum = 0.
do 20 k=1, kdo
t1 = 1.0-(k+k-1+i)**2
t2 = 1.0-(k+k-1-i)**2
sum = sum+work(k)*(t1-t2)/(t1*t2)
20 continue
czw(id) = sc1*sum
19 continue
return
!
!     n even, m odd
!
2 kdo = n/2
do 15 id=1, lc
i = id+id-1
sum = 0.
do 16 k=1, kdo
t1 = 1.0-(k+k-1+i)**2
t2 = 1.0-(k+k-1-i)**2
sum = sum+work(k)*(t1+t2)/(t1*t2)
16 continue
czw(id) = sc1*sum
15 continue
return
1 if (mmod /= 0) go to 3
!
!     n odd, m even
!
kdo = (n-1)/2
do 9 id=1, lc
i = id+id-2
sum = 0.
do 10 k=1, kdo
t1 = 1.0-(k+k+i)**2
t2 = 1.0-(k+k-i)**2
sum = sum+work(k)*(t1-t2)/(t1*t2)
10 continue
czw(id) = sc1*sum
9 continue
return
!
!     n odd, m odd
!
3 kdo = (n+1)/2
do 5 id=1, lc
i = id+id-2
sum = work(1)/(1.0-i*i)
if (kdo < 2) go to 29
do 6 kp1=2, kdo
k = kp1-1
t1 = 1.0-(k+k+i)**2
t2 = 1.0-(k+k-i)**2
sum = sum+work(kp1)*(t1+t2)/(t1*t2)
6 continue
29 czw(id) = sc1*sum
5 continue

end subroutine dzwk



subroutine dzwt(nlat, m, n, th, czw, zwh)
!
!     subroutine dzwt tabulates the function zwbar(n, m, theta)
!     at theta = th in real
!
!     input parameters
!
!     nlat      the number of colatitudes including the poles.
!            nlat must be an odd integer
!
!     n      the degree (subscript) of zwbar(n, m, theta)
!
!     m      the order (superscript) of zwbar(n, m, theta)
!
!     czw     the fourier coefficients of zwbar(n, m, theta)
!             as computed by subroutine zwk.
!
!     output parameter
!
!     zwh     zwbar(m, n, theta) evaluated at theta = th
!
dimension czw(1)
real czw, zwh, th, cth, sth, cdt, sdt, chh
zwh = 0.
if (n <= 0) return
lc = (nlat+1)/2
lq = lc-1
ls = lc-2
cth = cos(th)
sth = sin(th)
cdt = cth*cth-sth*sth
sdt = 2.*sth*cth
lmod = mod(nlat, 2)
mmod = mod(m, 2)
nmod = mod(n, 2)
if (lmod == 0) go to 50
if (nmod /= 0) go to 1
if (mmod /= 0) go to 2
!
!     nlat odd  n even  m even
!
do 30 k=1, lq
zwh = zwh+czw(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
30 continue
return
!
!     nlat odd  n even  m odd
!
2 do 40 k=1, lq
zwh = zwh+czw(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
40 continue
return
1 cth = cdt
sth = sdt
if (mmod /= 0) go to 3
!
!     nlat odd  n odd  m even
!
do 10 k=1, ls
zwh = zwh+czw(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
10 continue
return
!
!     nlat odd  n odd  m odd
!
3 zwh = .5*czw(1)
do 20 k=2, lq
zwh = zwh+czw(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
20 continue
zwh = zwh+.5*czw(lc)*cos((nlat-1)*th)
return
50 if (nmod /= 0) go to 51
if (mmod /= 0) go to 52
!
!     nlat even  n even  m even
!
do 55 k=1, lq
zwh = zwh+czw(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
55 continue
return
!
!     nlat even  n even  m odd
!
52 zwh = .5*czw(lc)*cos((nlat-1)*th)
do 60 k=1, lq
zwh = zwh+czw(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
60 continue
return
51 cth = cdt
sth = sdt
if (mmod /= 0) go to 53
!
!     nlat even  n odd  m even
!
do 65 k=1, lq
zwh = zwh+czw(k+1)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
65 continue
return
!
!     nlat even  n odd  m odd
!
53 zwh = .5*czw(1)
do 70 k=2, lc
zwh = zwh+czw(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
70 continue

end subroutine dzwt



subroutine dvbk(m, n, cv, work)
real cv(1), work(1), fn, fk, cf
cv(1) = 0.
if (n <= 0) return
fn = n
srnp1 = sqrt(fn * (fn + 1.0))
cf = 2.*m/srnp1
modn = mod(n, 2)
modm = mod(m, 2)
call dnlfk(m, n, work)
if (modn /= 0) go to 70
ncv = n/2
if (ncv == 0) return
fk = 0.
if (modm /= 0) go to 60
!
!     n even m even
!
do 55 l=1, ncv
fk = fk+2.
cv(l) = -fk*work(l+1)/srnp1
55 continue
return
!
!     n even m odd
!
60 do 65 l=1, ncv
fk = fk+2.
cv(l) = fk*work(l)/srnp1
65 continue
return
70 ncv = (n+1)/2
fk = -1.
if (modm /= 0) go to 80
!
!     n odd m even
!
do l=1, ncv
fk = fk+2.
cv(l) = -fk*work(l)/srnp1
end do
return
!
!     n odd m odd
!
80 do l=1, ncv
fk = fk+2.
cv(l) = fk*work(l)/srnp1
end do

end subroutine dvbk



subroutine dwbk(m, n, cw, work)
real cw(1), work(1), fn, cf, srnp1
cw(1) = 0.
if (n<=0 .or. m<=0) return
fn = n
srnp1 = sqrt(fn * (fn + 1.0))
cf = 2.*m/srnp1
modn = mod(n, 2)
modm = mod(m, 2)
call dnlfk(m, n, work)
if (m == 0) go to 50
if (modn /= 0) go to 30
l = n/2
if (l == 0) go to 50
if (modm /= 0) go to 20
!
!     n even m even
!
cw(l) = -cf*work(l+1)
10 l = l-1
if (l <= 0) go to 50
cw(l) = cw(l+1)-cf*work(l+1)
go to 10
!
!     n even m odd
!
20 cw(l) = cf*work(l)
25 l = l-1
if (l <= 0) go to 50
cw(l) = cw(l+1)+cf*work(l)
go to 25
30 if (modm /= 0) go to 40
l = (n-1)/2
if (l == 0) go to 50
!
!     n odd m even
!
cw(l) = -cf*work(l+1)
35 l = l-1
if (l <= 0) go to 50
cw(l) = cw(l+1)-cf*work(l+1)
go to 35
!
!     n odd m odd
!
40 l = (n+1)/2
cw(l) = cf*work(l)
45 l = l-1
if (l <= 0) go to 50
cw(l) = cw(l+1)+cf*work(l)
go to 45
50 return

end subroutine dwbk



subroutine dvbt(m, n, theta, cv, vh)
dimension cv(1)
real cv, vh, theta, cth, sth, cdt, sdt, chh
vh = 0.
if (n==0) return
cth = cos(theta)
sth = sin(theta)
cdt = cth*cth-sth*sth
sdt = 2.*sth*cth
mmod = mod(m, 2)
nmod = mod(n, 2)
if (nmod /= 0) go to 1
cth = cdt
sth = sdt
if (mmod /= 0) go to 2
!
!     n even  m even
!
ncv = n/2
do k=1, ncv
vh = vh+cv(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
end do
return
!
!     n even  m odd
!
2 ncv = n/2
do k=1, ncv
vh = vh+cv(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
end do

return
1 if (mmod /= 0) go to 3
!
!     n odd m even
!
ncv = (n+1)/2
do k=1, ncv
vh = vh+cv(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
end do
return
!
! case m odd and n odd
!
3 ncv = (n+1)/2
do k=1, ncv
vh = vh+cv(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
end do

end subroutine dvbt



subroutine dwbt(m, n, theta, cw, wh)
dimension cw(1)
real theta, cw, wh, cth, sth, cdt, sdt, chh
wh = 0.
if (n<=0 .or. m<=0) return
cth = cos(theta)
sth = sin(theta)
cdt = cth*cth-sth*sth
sdt = 2.*sth*cth
mmod=mod(m, 2)
nmod=mod(n, 2)
if (nmod /= 0) go to 1
if (mmod /= 0) go to 2
!
!     n even  m even
!
ncw = n/2
do 10 k=1, ncw
wh = wh+cw(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
10 continue
return
!
!     n even  m odd
!
2 ncw = n/2
do 8 k=1, ncw
wh = wh+cw(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
8 continue
return
1 cth = cdt
sth = sdt
if (mmod /= 0) go to 3
!
!     n odd m even
!
ncw = (n-1)/2
do 20 k=1, ncw
wh = wh+cw(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
20 continue
return
!
! case m odd and n odd
!
3 ncw = (n+1)/2
wh = .5*cw(1)
if (ncw<2) return
do 25 k=2, ncw
wh = wh+cw(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
25 continue

end subroutine dwbt



subroutine rabcv(nlat, nlon, abc)
!
!     subroutine rabcp computes the coefficients in the recurrence
!     relation for the functions vbar(m, n, theta). array abc
!     must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2 locations.
!
dimension abc(1)
mmax = min(nlat, (nlon+1)/2)
labc = (max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
iw1 = labc+1
iw2 = iw1+labc
call rabcv1(nlat, nlon, abc, abc(iw1), abc(iw2))

end subroutine rabcv



subroutine rabcv1(nlat, nlon, a, b, c)
!
!     coefficients a, b, and c for computing vbar(m, n, theta) are
!     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
!
dimension a(1), b(1), c(1)
mmax = min(nlat, (nlon+1)/2)
if (mmax < 3) return
do 215 mp1=3, mmax
m = mp1-1
ns = ((m-2)*(nlat+nlat-m-1))/2+1
fm = real(m)
tm = fm+fm
temp = tm*(tm-1.)
tpn = (fm-2.)*(fm-1.)/(fm*(fm+1.))
a(ns) = sqrt(tpn*(tm+1.)*(tm-2.)/temp)
c(ns) = sqrt(2./temp)
if (m == nlat-1) go to 215
ns = ns+1
temp = tm*(tm+1.)
tpn = (fm-1.)*fm/((fm+1.)*(fm+2.))
a(ns) = sqrt(tpn*(tm+3.)*(tm-2.)/temp)
c(ns) = sqrt(6./temp)
mp3 = m+3
if (mp3 > nlat) go to 215
do 210 np1=mp3, nlat
n = np1-1
ns = ns+1
fn = real(n)
tn = fn+fn
cn = (tn+1.)/(tn-3.)
tpn = (fn-2.)*(fn-1.)/(fn*(fn + 1.0))
fnpm = fn+fm
fnmm = fn-fm
temp = fnpm*(fnpm-1.)
a(ns) = sqrt(tpn*cn*(fnpm-3.)*(fnpm-2.)/temp)
b(ns) = sqrt(tpn*cn*fnmm*(fnmm-1.)/temp)
c(ns) = sqrt((fnmm+1.)*(fnmm+2.)/temp)
210 continue
215 continue

end subroutine rabcv1



subroutine rabcw(nlat, nlon, abc)
!
!     subroutine rabcw computes the coefficients in the recurrence
!     relation for the functions wbar(m, n, theta). array abc
!     must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2 locations.
!
dimension abc(1)
mmax = min(nlat, (nlon+1)/2)
labc = (max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
iw1 = labc+1
iw2 = iw1+labc
call rabcw1(nlat, nlon, abc, abc(iw1), abc(iw2))

end subroutine rabcw



subroutine rabcw1(nlat, nlon, a, b, c)
!
!     coefficients a, b, and c for computing wbar(m, n, theta) are
!     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1
!
dimension a(1), b(1), c(1)
mmax = min(nlat, (nlon+1)/2)
if (mmax < 4) return
do 215 mp1=4, mmax
m = mp1-1
ns = ((m-2)*(nlat+nlat-m-1))/2+1
fm = real(m)
tm = fm+fm
temp = tm*(tm-1.)
tpn = (fm-2.)*(fm-1.)/(fm*(fm+1.))
tph = fm/(fm-2.)
a(ns) = tph*sqrt(tpn*(tm+1.)*(tm-2.)/temp)
c(ns) = tph*sqrt(2./temp)
if (m == nlat-1) go to 215
ns = ns+1
temp = tm*(tm+1.)
tpn = (fm-1.)*fm/((fm+1.)*(fm+2.))
tph = fm/(fm-2.)
a(ns) = tph*sqrt(tpn*(tm+3.)*(tm-2.)/temp)
c(ns) = tph*sqrt(6./temp)
mp3 = m+3
if (mp3 > nlat) go to 215
do 210 np1=mp3, nlat
n = np1-1
ns = ns+1
fn = real(n)
tn = fn+fn
cn = (tn+1.)/(tn-3.)
fnpm = fn+fm
fnmm = fn-fm
temp = fnpm*(fnpm-1.)
tpn = (fn-2.)*(fn-1.)/(fn*(fn + 1.0))
tph = fm/(fm-2.)
a(ns) = tph*sqrt(tpn*cn*(fnpm-3.)*(fnpm-2.)/temp)
b(ns) = sqrt(tpn*cn*fnmm*(fnmm-1.)/temp)
c(ns) = tph*sqrt((fnmm+1.)*(fnmm+2.)/temp)
210 continue
215 continue

end subroutine rabcw1



subroutine vtinit (nlat, nlon, wvbin, dwork)
dimension       wvbin(*)
real dwork(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of dwork is nlat+2
!
call vtini1 (nlat, nlon, imid, wvbin, wvbin(iw1), dwork, &
                                       dwork(nlat/2+2))

end subroutine vtinit



subroutine vtini1 (nlat, nlon, imid, vb, abc, cvb, work)
!
!     abc must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     locations where mmax = min(nlat, (nlon+1)/2)
!     cvb and work must each have nlat/2+1 locations
!
dimension vb(imid, nlat, 2), abc(1), cvb(1)
real dt, cvb, th, vbh, work(*)
real, parameter :: pi = acos(-1.0)
dt = pi/(nlat-1)
mdo = min(2, nlat, (nlon+1)/2)
do mp1=1, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dvtk(m, n, cvb, work)
do i=1, imid
th = (i-1)*dt
call dvtt(m, n, th, cvb, vbh)
vb(i, np1, mp1) = vbh
end do
end do
end do
call rabcv(nlat, nlon, abc)

end subroutine vtini1



subroutine wtinit (nlat, nlon, wwbin, dwork)
dimension       wwbin(1)
real dwork(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of dwork is nlat+2
!
call wtini1 (nlat, nlon, imid, wwbin, wwbin(iw1), dwork, &
                                       dwork(nlat/2+2))

end subroutine wtinit



subroutine wtini1 (nlat, nlon, imid, wb, abc, cwb, work)
!
!     abc must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     locations where mmax = min(nlat, (nlon+1)/2)
!     cwb and work must each have nlat/2+1 locations
!
dimension wb(imid, nlat, 2), abc(1)
real dt, cwb(*), wbh, th, work(*)
real, parameter :: pi = acos(-1.0)
dt = pi/(nlat-1)
mdo = min(3, nlat, (nlon+1)/2)
if (mdo < 2) return
do mp1=2, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dwtk(m, n, cwb, work)
do i=1, imid
th = (i-1)*dt
call dwtt(m, n, th, cwb, wbh)
wb(i, np1, m) = wbh
end do
end do
end do

call rabcw(nlat, nlon, abc)

end subroutine wtini1



subroutine vtgint (nlat, nlon, theta, wvbin, work)
dimension       wvbin(*)
real theta(*), work(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     theta is a real array with (nlat+1)/2 locations
!     nlat is the maximum value of n+1
!     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of work is nlat+2
!
call vtgit1 (nlat, nlon, imid, theta, wvbin, wvbin(iw1), &
                        work, work(nlat/2+2))

end subroutine vtgint



subroutine vtgit1 (nlat, nlon, imid, theta, vb, abc, cvb, work)
!
!     abc must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     locations where mmax = min(nlat, (nlon+1)/2)
!     cvb and work must each have nlat/2+1   locations
!
dimension vb(imid, nlat, 2), abc(*)
real theta(*), cvb(*), work(*), vbh
mdo = min(2, nlat, (nlon+1)/2)
do mp1=1, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dvtk(m, n, cvb, work)
do i=1, imid
call dvtt(m, n, theta(i), cvb, vbh)
vb(i, np1, mp1) = vbh
end do
end do
end do

call rabcv(nlat, nlon, abc)

end subroutine vtgit1



subroutine wtgint (nlat, nlon, theta, wwbin, work)
dimension       wwbin(*)
real theta(*), work(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     theta is a real array with (nlat+1)/2 locations
!     nlat is the maximum value of n+1
!     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of work is nlat+2
!
call wtgit1 (nlat, nlon, imid, theta, wwbin, wwbin(iw1), &
                        work, work(nlat/2+2))

end subroutine wtgint



subroutine wtgit1 (nlat, nlon, imid, theta, wb, abc, cwb, work)
!
!     abc must have 3*((nlat-3)*nlat+2)/2 locations
!     cwb and work must each have nlat/2+1 locations
!
dimension wb(imid, nlat, 2), abc(1)
real theta(*), cwb(*), work(*), wbh
mdo = min(3, nlat, (nlon+1)/2)
if (mdo < 2) return
do mp1=2, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dwtk(m, n, cwb, work)
do i=1, imid
call dwtt(m, n, theta(i), cwb, wbh)
wb(i, np1, m) = wbh
end do
end do
end do

call rabcw(nlat, nlon, abc)

end subroutine wtgit1



subroutine dvtk(m, n, cv, work)
real cv(*), work(*), fn, fk, cf, srnp1
cv(1) = 0.
if (n <= 0) return
fn = n
srnp1 = sqrt(fn * (fn + 1.0))
cf = 2.*m/srnp1
modn = mod(n, 2)
modm = mod(m, 2)
call dnlfk(m, n, work)
if (modn /= 0) go to 70
ncv = n/2
if (ncv == 0) return
fk = 0.
if (modm /= 0) go to 60
!
!     n even m even
!
do 55 l=1, ncv
fk = fk+2.
cv(l) = -fk*fk*work(l+1)/srnp1
55 continue
return
!
!     n even m odd
!
60 do 65 l=1, ncv
fk = fk+2.
cv(l) = -fk*fk*work(l)/srnp1
65 continue
return
70 ncv = (n+1)/2
fk = -1.
if (modm /= 0) go to 80
!
!     n odd m even
!
do 75 l=1, ncv
fk = fk+2.
cv(l) = -fk*fk*work(l)/srnp1
75 continue
return
!
!     n odd m odd
!
80 do 85 l=1, ncv
fk = fk+2.
cv(l) = -fk*fk*work(l)/srnp1
85 continue

end subroutine dvtk



subroutine dwtk(m, n, cw, work)
real cw(*), work(*), fn, cf, srnp1
cw(1) = 0.
if (n<=0 .or. m<=0) return
fn = n
srnp1 = sqrt(fn * (fn + 1.0))
cf = 2.*m/srnp1
modn = mod(n, 2)
modm = mod(m, 2)
call dnlfk(m, n, work)
if (m == 0) go to 50
if (modn /= 0) go to 30
l = n/2
if (l == 0) go to 50
if (modm /= 0) go to 20
!
!     n even m even
!
cw(l) = -cf*work(l+1)
10 l = l-1
if (l <= 0) go to 50
cw(l) = cw(l+1)-cf*work(l+1)
cw(l+1) = (l+l+1)*cw(l+1)
go to 10
!
!     n even m odd
!
20 cw(l) = cf*work(l)
25 l = l-1
   if (l< 0) then
       goto 50
   else if (l == 0) then
       goto 27
   else
       goto 26
   end if
26 cw(l) = cw(l+1)+cf*work(l)
27 cw(l+1) = -(l+l+1)*cw(l+1)
go to 25
30 if (modm /= 0) go to 40
l = (n-1)/2
if (l == 0) go to 50
!
!     n odd m even
!
cw(l) = -cf*work(l+1)
35 l = l-1
   if (l< 0) then
       goto 50
   else if (l == 0) then
       goto 37
   else
       goto 36
   end if
36 cw(l) = cw(l+1)-cf*work(l+1)
37 cw(l+1) = (l+l+2)*cw(l+1)
go to 35
!
!     n odd m odd
!
40 l = (n+1)/2
cw(l) = cf*work(l)
45 l = l-1
   if (l< 0) then
       goto 50
   else if (l == 0) then
       goto 47
   else
       goto 46
   end if
46 cw(l) = cw(l+1)+cf*work(l)
47 cw(l+1) = -(l+l)*cw(l+1)
go to 45
50 return

end subroutine dwtk



subroutine dvtt(m, n, theta, cv, vh)
dimension cv(1)
real cv, vh, theta, cth, sth, cdt, sdt, chh
vh = 0.
if (n==0) return
cth = cos(theta)
sth = sin(theta)
cdt = cth*cth-sth*sth
sdt = 2.*sth*cth
mmod = mod(m, 2)
nmod = mod(n, 2)
if (nmod /= 0) go to 1
cth = cdt
sth = sdt
if (mmod /= 0) go to 2
!
!     n even  m even
!
ncv = n/2
do 10 k=1, ncv
vh = vh+cv(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
10 continue
return
!
!     n even  m odd
!
2 ncv = n/2
do 15 k=1, ncv
vh = vh+cv(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
15 continue
return
1 if (mmod /= 0) go to 3
!
!     n odd m even
!
ncv = (n+1)/2
do 20 k=1, ncv
vh = vh+cv(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
20 continue
return
!
! case m odd and n odd
!
3 ncv = (n+1)/2
do 25 k=1, ncv
vh = vh+cv(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
25 continue

end subroutine dvtt



subroutine dwtt(m, n, theta, cw, wh)
dimension cw(1)
real theta, cw, wh, cth, sth, cdt, sdt, chh
wh = 0.
if (n<=0 .or. m<=0) return
cth = cos(theta)
sth = sin(theta)
cdt = cth*cth-sth*sth
sdt = 2.*sth*cth
mmod=mod(m, 2)
nmod=mod(n, 2)
if (nmod /= 0) go to 1
if (mmod /= 0) go to 2
!
!     n even  m even
!
ncw = n/2
do 10 k=1, ncw
wh = wh+cw(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
10 continue
return
!
!     n even  m odd
!
2 ncw = n/2
do 8 k=1, ncw
wh = wh+cw(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
8 continue
return
1 cth = cdt
sth = sdt
if (mmod /= 0) go to 3
!
!     n odd m even
!
ncw = (n-1)/2
do 20 k=1, ncw
wh = wh+cw(k)*cth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
20 continue
return
!
! case m odd and n odd
!
3 ncw = (n+1)/2
wh = 0.
if (ncw<2) return
do 25 k=2, ncw
wh = wh+cw(k)*sth
chh = cdt*cth-sdt*sth
sth = sdt*cth+cdt*sth
cth = chh
25 continue

end subroutine dwtt



subroutine vbgint (nlat, nlon, theta, wvbin, work)
dimension       wvbin(1)
real theta(*), work(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     theta is a real array with (nlat+1)/2 locations
!     nlat is the maximum value of n+1
!     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of work is nlat+2
!
call vbgit1 (nlat, nlon, imid, theta, wvbin, wvbin(iw1), &
                        work, work(nlat/2+2))

end subroutine vbgint



subroutine vbgit1 (nlat, nlon, imid, theta, vb, abc, cvb, work)
!
!     abc must have 3*(max(mmax-2, 0)*(nlat+nlat-mmax-1))/2
!     locations where mmax = min(nlat, (nlon+1)/2)
!     cvb and work must each have nlat/2+1 locations
!
dimension vb(imid, nlat, 2), abc(1)
real cvb(1), theta(1), vbh, work(1)
mdo = min(2, nlat, (nlon+1)/2)
do mp1=1, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dvbk(m, n, cvb, work)
do i=1, imid
call dvbt(m, n, theta(i), cvb, vbh)
vb(i, np1, mp1) = vbh
end do
end do
end do

call rabcv(nlat, nlon, abc)

end subroutine vbgit1



subroutine wbgint (nlat, nlon, theta, wwbin, work)
dimension       wwbin(1)
real work(*), theta(*)
imid = (nlat+1)/2
iw1 = 2*nlat*imid+1
!
!     theta is a real array with (nlat+1)/2 locations
!     nlat is the maximum value of n+1
!     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2
!     the length of work is nlat+2
!
call wbgit1 (nlat, nlon, imid, theta, wwbin, wwbin(iw1), &
                        work, work(nlat/2+2))

end subroutine wbgint



subroutine wbgit1 (nlat, nlon, imid, theta, wb, abc, cwb, work)
!
!     abc must have 3*((nlat-3)*nlat+2)/2 locations
!     cwb and work must each have nlat/2+1 locations
!
dimension wb(imid, nlat, 2), abc(1)
real cwb(1), theta(1), wbh, work(1)
mdo = min(3, nlat, (nlon+1)/2)
if (mdo < 2) return
do mp1=2, mdo
m = mp1-1
do np1=mp1, nlat
n = np1-1
call dwbk(m, n, cwb, work)
do i=1, imid
call dwbt(m, n, theta(i), cwb, wbh)
wb(i, np1, m) = wbh
end do
end do
end do

call rabcw(nlat, nlon, abc)

end subroutine wbgit1
