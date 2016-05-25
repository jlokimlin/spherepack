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
! ... file visgau.f
!
!     contains documentation and code for subroutine visgau
!
subroutine visequ (nlat, nlon, h, len, eyer, eyelat, eyelon, &
                                    wk, lwk, iwk, liwk, ierror)
!
!     subroutine visequ will display a function on the sphere
!     as a solid. ie. as a "lumpy" sphere. visequ calls subroutine
!     vsurf to produce the visible surface rendering.
!
!     requires routines visequ1 interp sptc diag stride triang vsurf
!                       vsurf1 prjct box icvmg projct
!
!     visgeo uses the ncar graphics package.
!     compile with: ncargf77 (all programs above)
!
!     execute with:  a.out
!
!     on screen display with:  ctrans -d x11 gmeta
!                          
!     print with:  ctrans -d ps.color gmeta > gmeta.ps
!                  lpr -p(your printer) gmeta.ps 
!
!     input parameters
!
!     nlat   the number of colatitudes on the full sphere including the
!            poles. for example, nlat = 37 for a five degree grid.
!            nlat determines the grid increment in colatitude as
!            pi/(nlat-1).  if nlat is odd the equator is located at
!            grid point i=(nlat+1)/2. if nlat is even the equator is
!            located half way between points i=nlat/2 and i=nlat/2+1.
!            nlat must be at least 3. note: on the half sphere, the
!            number of grid points in the colatitudinal direction is
!            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!
!     h      a two dimensional array that contains the discrete
!            function to be displayed. h(i, j) is the distance from the
!            center of the sphere to the surface at the colatitude 
!            point theta(i) = (i-1)*pi/(nlat-1) and longitude point 
!            phi(j) = (j-1)*2*pi/nlon.
!
!     len    the first dimension of the array h as it appears in the
!            program that calls visequ.
!
!     eyer   the distance from the center of the sphere to the eye.
!
!     eyelat the colatitudinal coordinate of the eye (in degrees).
!
!     eyelon the longitudinal  coordinate of the eye (in degrees).
!
!     wk     a real work array 
!
!     lwk    the dimension of the array wk as it appears in the
!            program that calls visequ. lwk must be at least 
!                       46*nlat*(nlon+1).
!
!     iwk    an integer work array
!
!     liwk   the dimension of the array iwk as it appears in the
!            program that calls visequ. liwk must be at least 
!                       14*nlat*(nlon+1).
!
!     ierror = 0    no error
!            = 1    the eye is positioned inside the sphere
!            = 2    lwk  is less than 46*nlat*(nlon+1)
!            = 3    liwk is less than 14*nlat*(nlon+1)
!
!  
dimension h(len, *), wk(*)
integer iwk(*)
n = nlat
m = nlon+1
mn = m*n
ierror = 2
if (lwk < 46*mn) return
ierror = 3
if (liwk < 14*mn) return
ierror = 1
do 10 j=1, nlon 
do 10 i=1, nlat
if (eyer <= h(i, j)) return
10 continue
ierror = 0
pi = acos(-1.0)
dtr = pi/180.
!
!     ****     set up pointers to sub work arrays in wk and iwk
!
ntri = mn+mn
nw1 = 1
nw2 = nw1+mn
nclat = 1
nslat = nclat+n
nxp = 1
nyp = nxp+mn
nx1 = 1
ny1 = nx1+ntri
nz1 = ny1+ntri
nx2 = nz1+ntri
ny2 = nx2+ntri
nz2 = ny2+ntri
nx3 = nz2+ntri
ny3 = nx3+ntri
nz3 = ny3+ntri
nx  = nz3+ntri
ny = nx+mn
nz = ny+mn
nwrk = nx
nitype = 1
niflag = ntri+1
nmst = niflag+mn
nmfac = nmst+n
!     total iwk is 7*ntri
!     total wk is 58*nlat*(nlon+1)
!     ****     mid-cell interpolation, calculation of polar values
call interp(h, len, m, n, wk(nw1), wk(nw2), iwk(niflag))
!     ****     transform grid points to cartesian coordinates
call sptc(h, len, m, n, wk(nclat), wk(nslat), wk(nx), wk(ny), wk(nz))
!     ****     transform eye position to cartesian coordinates
xeye=eyer*sin(dtr*eyelat)
yeye=xeye*sin(dtr*eyelon)
xeye=xeye*cos(dtr*eyelon)
zeye=eyer*cos(dtr*eyelat)
!     ****     project grid points
call projct(m, n, xeye, yeye, zeye, wk(nx), wk(ny), wk(nz), wk(nxp), &
            wk(nyp)) 
!     ****     check for visibility of cell boundaries
call diag(m, n, wk(nxp), wk(nyp), iwk(niflag))
!     ****     compute longitude stride as a function of latitude
call stride(m, n, iwk(nmst), iwk(nmfac))
!     ****     perform triangulation
call triang(m, n, wk(nx), wk(ny), wk(nz), itri, wk(nx1), wk(ny1), &
wk(nz1), wk(nx2), wk(ny2), wk(nz2), wk(nx3), wk(ny3), wk(nz3), &
iwk(nitype), iwk(niflag), iwk(nmst))
!     ****     call surface plotting routine
call vsurf(xeye, yeye, zeye, itri, wk(nx1), wk(ny1), wk(nz1), wk(nx2), &
wk(ny2), wk(nz2), wk(nx3), wk(ny3), wk(nz3), iwk(nitype), wk(nwrk), &
iwk(niflag))
return
end subroutine visequ
subroutine interp(h, len, m, n, w1, w2, iflag)
!     ****     interpolates to mid points of grid cells using second
!     ****     order formula
dimension h(len, 1), w1(n, m), w2(n+2, m+2), iflag(n, m), sten(4, 4)
data sten/.015625, 2*-.078125, .015625, -.078125, 2*.390625, &
2*-.078125, 2*.390625, -.078125, .015625, 2*-.078125, .015625/
!     ****     copy h to w2
mm1 = m-1
do 1 i=1, mm1
do 1 j=1, n
w2(j, i+1)=h(j, i)
1 continue
!     ****     add periodic points
do 2 j=1, n
w2(j, 1)=w2(j, m)
w2(j, m+1)=w2(j, 2)
w2(j, m+2)=w2(j, 3)
2 continue
n1=2
n2=n-2
!     ****     perform interpolation
!     ****     set w1 to zero

do 7 i=1, m
do 7 j=1, n
w1(j, i)=0.
7 continue
!     ****     interpolate
do 8 k=1, 4
do 8 l=1, 4
do 8 i=1, m-1
do 8 j=n1, n2
w1(j, i)=w1(j, i)+w2(j+l-2, i+k-1)*sten(k, l)
8 continue
!     ****     set up iflag array
!     ****     iflag(j, i)=0  if diagonal is (j, i) to (j+1, i+1)
!     ****     iflag(j, i)=16 if diagonal is (j+1, i), (j, i+1)
do 9 i=1, m-1
do 9 j=n1, n2
iflag(j, i)=icvmg(16, 0, abs(.5*(w2(j, i+1)+w2(j+1, i+2))-w1(j, i))- &
abs(.5*(w2(j, i+2)+w2(j+1, i+1))-w1(j, i)))
9 continue
return
end subroutine interp
subroutine sptc(r, len, m, n, clat, slat, x, y, z)
!     ****     transforms from spherical to cartesian coordinates
dimension r(len, 1), clat(n), slat(n), x(n, m), y(n, m), z(n, m)
pi = acos(-1.0)
dt = pi/(n-1)
dp = (pi+pi)/(m-1)
do 10 j=1, n
clat(j) = cos((j-1)*dt)
slat(j) = sin((j-1)*dt)
10 continue
do 20 i=1, m-1
clon = cos((i-1)*dp)
slon = sin((i-1)*dp)
do 20 j=1, n
x(j, i)=r(j, i)*slat(j)
y(j, i)=x(j, i)*slon
x(j, i)=x(j, i)*clon
z(j, i)=r(j, i)*clat(j)
20 continue
do 30 j=1, n
x(j, m)=x(j, 1)
y(j, m)=y(j, 1)
z(j, m)=z(j, 1)
30 continue
return
end subroutine sptc
subroutine diag(m, n, xp, yp, iflag)
!
!     ****     label visibility of cell sides
!
!     north side corresponds to j
!     south side corresponds to j+1
!     west  side corresponds to i
!     east  side corresponds to i+1
!
!     let iflag = b4 b3 b2 b1 b0 (in binary) then b0 through b3 are 
!     either o or 1 depending on whether the east, south, north
!     or west side is either invisible or visible, respectively.
!
!     b4 is o if the diagonal is from (i, j) to (i+1, j+1) and 1 if
!     the diagonal is from (i, j+1) to (i+1, j).
!
dimension xp(n, m), yp(n, m), iflag(n, m)
!     ****     arithmetic statement function
cp(j1, i1, j2, i2, j3, i3)=((xp(j1, i1)-xp(j2, i2))*(yp(j3, i3)-yp(j2, i2)) &
-(xp(j3, i3)-xp(j2, i2))*(yp(j1, i1)-yp(j2, i2)))
do 100 j=2, n-2
do 100 i=1, m-1
if (iflag(j, i) >= 16) go to 20
if (cp(j+1, i+1, j+1, i, j, i) <= 0) go to 10 
!     west and south are visible
iflag(j, i) = iflag(j, i)+10  
10 if (cp(j, i, j, i+1, j+1, i+1) <= 0) go to 100 
!     east and north are visible
iflag(j, i) = iflag(j, i)+5  
go to 100
20 if (cp(j+1, i, j, i, j, i+1) <= 0) go to 30 
!     west and north are visible
iflag(j, i) = iflag(j, i)+12  
30 if (cp(j, i+1, j+1, i+1, j+1, i) <= 0) go to 100 
!     east and south are visible
iflag(j, i) = iflag(j, i)+3  
100 continue
! 
!     classify the poles
!
do 200 i=1, m-1
iflag(1, i) = 0
if (cp(2, i+1, 2, i, 1, i) > 0) iflag(1, i) = 15
iflag(n-1, i) = 0
if (cp(n, i, n-1, i, n-1, i+1) > 0) iflag(n-1, i) = 31
200 continue
do 250 j=1, n-1
iflag(j, m) = iflag(j, 1)
250 continue 
return
end subroutine diag
subroutine stride(m, n, mst, mfac)
dimension mfac(*), mtryh(3), mst(n), icl(8)
data mtryh(1), mtryh(2), mtryh(3)/2, 3, 5/                        
data icl(1), icl(2), icl(3), icl(4), icl(5), icl(6), icl(7), icl(8) &
     /0, 1, 2, 12, 3, 13, 23, 123/
!
!     find prime factors of m-1
!
ml = m-1
nf = 0           
j = 0            
101 j = j+1          
    if (j-3< 0) then                
        goto 102
    else if (j-3 == 0) then 
        goto 102
    else 
        goto 103
    end if
102 mtry = mtryh(j)                     
go to 104        
103 mtry = mtry+2    
104 mq = ml/mtry     
mr = ml-mtry*mq                     
if (mr< 0) then                 
    goto 101
else if (mr == 0) then 
    goto 105
else 
    goto 101
end if
105 nf = nf+1        
mfac(nf) = mtry                    
ml = mq          
if (ml /= 1) go to 104            
if (mfac(nf) > 2) go to 106
nf = nf-1
mfac(nf) = 4
106 tphi = .707/real(m-1)
ns2 = n/2
mf1 = mfac(nf)
mst(1) = (m-1)/mf1
pi = acos(-1.0)
dt = pi/real(n - 1)
jf = nf-1
do 110 jdo=2, ns2 
j = jdo
theta = (j-1)*dt
st = sin(theta)
mf2 = mf1*mfac(jf)
if (abs(st/mf1-tphi) > abs(st/mf2-tphi)) go to 115
mst(j) = mst(j-1)
go to 110
115 mst(j) = (m-1)/mf2
mf1 = mf2
jf = jf-1
if (jf == 0) go to 120
110 continue
120 do 125 jdo=j, ns2
mst(jdo) = 1
125 continue
do 130 jdo=1, ns2
mst(n-jdo) = mst(jdo)
130 continue
!      write (6, 135) (mst(j), j=1, n)
135 format(' colatitude strides'/(15i5))
!
return
end subroutine stride
subroutine triang(m, n, x, y, z, itri, x1, y1, z1, x2, y2, z2, x3, y3, z3, &
                  ityp, iflag, mst)
!     ****     performs triangulation
dimension x(n, m), y(n, m), z(n, m), x1(1), y1(1), z1(1), &
x2(1), y2(1), z2(1), x3(1), y3(1), z3(1), ityp(1), iflag(n, m), &
mst(n), icl(8) 
data icl(1), icl(2), icl(3), icl(4), icl(5), icl(6), icl(7), icl(8) &
     /0, 1, 2, 12, 3, 13, 23, 123/
itri = 0
n1=2
n2=n-2
do 100 j=n1, n2
do 100 i=1, m-1
if (iflag(j, i) >= 16) go to 50
if (mod(iflag(j, i), 16) < 8) go to 70
itri = itri+1
x1(itri) = x(j, i)
y1(itri) = y(j, i)
z1(itri) = z(j, i)
x2(itri) = x(j+1, i)
y2(itri) = y(j+1, i)
z2(itri) = z(j+1, i)
x3(itri) = x(j+1, i+1)
y3(itri) = y(j+1, i+1)
z3(itri) = z(j+1, i+1)
ityph = 3
if (mod(i-1, mst(j)) == 0) go to 60
if (mod(iflag(j, i-1), 2) == 0) go to 60
ityph = ityph-1
60 if (mod(iflag(j, i), 2) == 0) ityph = ityph+4
ityp(itri) = icl(ityph+1)
70 if (mod(iflag(j, i), 2) == 0) go to 100
itri = itri+1
x1(itri) = x(j, i)
y1(itri) = y(j, i)
z1(itri) = z(j, i)
x2(itri) = x(j+1, i+1)
y2(itri) = y(j+1, i+1)
z2(itri) = z(j+1, i+1)
x3(itri) = x(j, i+1)
y3(itri) = y(j, i+1)
z3(itri) = z(j, i+1)
ityph = 0
if (mod(iflag(j, i), 16) < 8) ityph = ityph+1
if (mod(iflag(j, i+1), 16) < 8) ityph = ityph+2
if (mod(iflag(j-1, i), 4) < 2) ityph = ityph+4
ityp(itri) = icl(ityph+1)
go to 100
50 if (mod(iflag(j, i), 16) < 8) go to 20
itri = itri+1
x1(itri) = x(j, i)
y1(itri) = y(j, i)
z1(itri) = z(j, i)
x2(itri) = x(j+1, i)
y2(itri) = y(j+1, i)
z2(itri) = z(j+1, i)
x3(itri) = x(j, i+1)
y3(itri) = y(j, i+1)
z3(itri) = z(j, i+1)
ityph = 1
if (mod(i-1, mst(j)) == 0) go to 10
if (mod(iflag(j, i-1), 2) == 0) go to 10
ityph = 0
10 if (mod(iflag(j, i), 2) == 0) ityph = ityph+2
if (mod(iflag(j-1, i), 4) < 2) ityph = ityph+4
ityp(itri) = icl(ityph+1)
20 if (mod(iflag(j, i), 2) == 0) go to 100
itri = itri+1
x1(itri) = x(j+1, i)
y1(itri) = y(j+1, i)
z1(itri) = z(j+1, i)
x2(itri) = x(j+1, i+1)
y2(itri) = y(j+1, i+1)
z2(itri) = z(j+1, i+1)
x3(itri) = x(j, i+1)
y3(itri) = y(j, i+1)
z3(itri) = z(j, i+1)
ityph = 1
if (mod(iflag(j, i+1), 16) < 8) ityph = ityph+2
if (mod(iflag(j, i), 16) < 8) ityph = ityph+4
ityp(itri) = icl(ityph+1)
100 continue
!
!     ****     triangles around north and south poles
!
do 200 i=1, m-1
if (mod(iflag(1, i), 16) < 8) go to 250
itri = itri+1
x1(itri) = x(1, i)
y1(itri) = y(1, i)
z1(itri) = z(1, i)
x2(itri) = x(2, i)
y2(itri) = y(2, i)
z2(itri) = z(2, i)
x3(itri) = x(2, i+1)
y3(itri) = y(2, i+1)
z3(itri) = z(2, i+1)
ityph = 3
if (mod(i-1, mst(1)) == 0) go to 260
if (mod(iflag(1, i-1), 2) == 0) go to 260
ityph = ityph-1
260 if (mod(iflag(1, i+1), 16) < 8) ityph = ityph+4
ityp(itri) = icl(ityph+1)
250 if (mod(iflag(n-1, i), 16) < 8) go to 200 
itri = itri+1
x1(itri)=x(n-1, i)
y1(itri)=y(n-1, i)
z1(itri)=z(n-1, i)
x2(itri)=x(n, i)
y2(itri)=y(n, i)
z2(itri)=z(n, i)
x3(itri)=x(n-1, i+1)
y3(itri)=y(n-1, i+1)
z3(itri)=z(n-1, i+1)
ityph = 1
if (mod(i-1, mst(n-1)) == 0) go to 210
if (mod(iflag(n-1, i-1), 2) == 0) go to 210
ityph = 0
210 if (mod(iflag(n-1, i+1), 16) < 8) ityph = ityph+2
if (mod(iflag(n-2, i), 4) < 2) ityph = ityph+4
ityp(itri) = icl(ityph+1)
200 continue
return
end subroutine triang
subroutine VSURF(xeye, yeye, zeye, ntri, x1, y1, z1, x2, y2, z2, &
                 x3, y3, z3, itype, work, iwork)
!
!    subroutine vsurf is like subroutine hidel except the triangles
!    are categorized. vsurf is also like solid except triangles rather
!    than lines are covered.
!
!     written by paul n. swarztrauber, national center for atmospheric
!     research, p.o. box 3000, boulder, colorado, 80307  
!
!    this program plots visible lines for the surface defined
!    by the input 3-d triangles with corners at (x1, y1, z1), (x2, y2, z2)
!    and (x3, y3, z3). the sides of these these triangles may or
!    may not be plotted depending on itype. if itype is 1 then the
!    side between points (x1, y1, z1) and (x2, y2, z2) is plotted if it
!    is visible. if itype is 2 then the side between (x2, y2, z2)
!    and (x3, y3, z3) is plotted. if itype is 3 then the visible portion
!    of the side between (x3, y3, z3) and (x1, y1, z1) is plotted.
!    any combination is possible by specifying itype to be one
!    of the following values: 0, 1, 2, 3, 12, 13, 23, 123.
!
!    the length of real    array  work must be at least 14*ntri
!
!    the length of integer array iwork must be at least  6*ntri
!
!
!    the vertices of the triangles are renumbered by vsurf so that
!    their projections are orientated counterclockwise. the user need
!    only be aware that the vertices may be renumbered by vsurf.
!
dimension x1(ntri), y1(ntri), z1(ntri), x2(ntri), y2(ntri), z2(ntri), &
          x3(ntri), y3(ntri), z3(ntri), itype(ntri), work(14*ntri)
integer iwork(6*ntri)
!
call vsurf1(xeye, yeye, zeye, ntri, x1, y1, z1, x2, y2, z2, x3, y3, z3, &
 itype, work, work(ntri+1), work(2*ntri+1), work(3*ntri+1), &
 work(4*ntri+1), work(5*ntri+1), work(6*ntri+1), work(7*ntri+1), &
 work(8*ntri+1), work(9*ntri+1), work(10*ntri+1), work(11*ntri+1), &
 work(12*ntri+1), work(13*ntri+1), iwork, IWORK(ntri+1), &
 IWORK(2*ntri+1), IWORK(4*ntri+1))
return
end subroutine VSURF
subroutine vsurf1(xeye, yeye, zeye, ntri, x1, y1, z1, x2, y2, z2, x3, y3, z3, &
 itype, px1, py1, px2, py2, px3, py3, vx1, vy1, vx2, vy2, vx3, vy3, tl, tr, kh, &
 next, istart, ifinal)
!
dimension x1(ntri), y1(ntri), z1(ntri), x2(ntri), y2(ntri), z2(ntri), &
          x3(ntri), y3(ntri), z3(ntri), itype(ntri), &
          px1(ntri), py1(ntri), px2(ntri), py2(ntri), &
          px3(ntri), py3(ntri), vx1(ntri), vy1(ntri), &
          vx2(ntri), vy2(ntri), vx3(ntri), vy3(ntri), &
          tl(ntri), tr(ntri), next(ntri), kh(ntri), &
          istart(2*ntri), ifinal(2*ntri), ltp(3), &
          ird(11), ip2(11), nct(11), ncv(11), last(11)
!
real l2e
real le2
!
!     compute projections of 3-d points
!
le2 = log(2.0) !.6931471805599453094172321
l2e = 1.0/le2
fntri = ntri
irmax = .5*l2e*log(fntri)
irmax = min(irmax, 10)
irmp1 = irmax+1
do 4 icv=1, 11
ncv(icv) = 0
4 continue
nct(1) = 0
ip2(1) = 1
ird(1) = 0
isize = 4
do 7 irp1=2, irmp1
ir = irp1-1
nct(irp1) = 0
ip2(irp1) = 2**ir
ird(irp1) = ird(ir)+isize
isize = (ip2(irp1)+1)**2
7 continue 
isxm = ird(irmp1)+isize+1
do 8 isx=1, isxm
istart(isx) = 0
ifinal(isx) = 0
8 continue
do 6 i=1, ntri
next(i) = 0
6 continue 
call prjct(0, xeye, yeye, zeye, x, y, z, dum1, dum2)
!      write(6, 127) ntri
127 format(' ntri in hidel', i5)
do 86 k=1, ntri
call prjct(1, xeye, yeye, zeye, x1(k), y1(k), z1(k), px1(k), py1(k))
call prjct(1, xeye, yeye, zeye, x2(k), y2(k), z2(k), px2(k), py2(k))
call prjct(1, xeye, yeye, zeye, x3(k), y3(k), z3(k), px3(k), py3(k))
  if (k < 3) then
!          write(6, 333) xeye, yeye, zeye, x1(k), y1(k), z1(k), px1(k), py1(k)
333     format(' xeye, etc.', 8e8.1)
  endif
86 continue
!
!     orientate triangles counter clockwise
!
do 70 k=1, ntri
cprod = (px2(k)-px1(k))*(py3(k)-py1(k))-(py2(k)-py1(k)) &
       *(px3(k)-px1(k))
!      if (cprod.eq.0.) write(6, 79) k, px1(k), px2(k), px3(k), 
!     -                              py1(k), py2(k), py3(k)
79 format('  cprod=0 at k=', i5, 6e9.2)
if (cprod>=0.) go to 70
px1h = px1(k)
py1h = py1(k)
px1(k) = px2(k)
py1(k) = py2(k)
px2(k) = px1h
py2(k) = py1h
x1hold = x1(k)
y1hold = y1(k)
z1hold = z1(k)
x1(k) = x2(k)
y1(k) = y2(k)
z1(k) = z2(k)
x2(k) = x1hold
y2(k) = y1hold
z2(k) = z1hold
ityp = itype(k)
if (ityp==2) itype(k) = 3
if (ityp==3) itype(k) = 2
if (ityp==12) itype(k) = 13
if (ityp==13) itype(k) = 12
70 continue
!
!     set screen limits
!
pmax = px1(1)
pmin = px1(1)
do 87 k=1, ntri
pmin = amin1(pmin, px1(k), py1(k), px2(k), py2(k), px3(k), py3(k))
pmax = amax1(pmax, px1(k), py1(k), px2(k), py2(k), px3(k), py3(k))
87 continue
pmin = 1.1*pmin
pmax = 1.1*pmax
call set(0., 1., 0., 1., pmin, pmax, pmin, pmax, 1)
xmin = amin1(px1(1), px2(1), px3(1)) 
xmax = amax1(px1(1), px2(1), px3(1)) 
ymin = amin1(py1(1), py2(1), py3(1)) 
ymax = amax1(py1(1), py2(1), py3(1)) 
do 1 i=2, ntri
xmin = amin1(xmin, px1(i), px2(i), px3(i)) 
xmax = amax1(xmax, px1(i), px2(i), px3(i)) 
ymin = amin1(ymin, py1(i), py2(i), py3(i)) 
ymax = amax1(ymax, py1(i), py2(i), py3(i)) 
1 continue
dmx = xmax-xmin
dmy = ymax-ymin
if (dmx > dmy) go to 2 
c = ymin
d = ymax
xmid = .5*(xmin+xmax)
hdy = .5*dmy
a = xmid-hdy
b = xmid+hdy
go to 3
2 a = xmin
b = xmax
ymid = .5*(ymin+ymax)
hdx = .5*dmx
c = ymid-hdx
d = ymid+hdx
3 hgr = b-a
!
!     categorize triangles
!
do 100 i=1, ntri
xmin = amin1(px1(i), px2(i), px3(i))
xmax = amax1(px1(i), px2(i), px3(i))
ymin = amin1(py1(i), py2(i), py3(i))
ymax = amax1(py1(i), py2(i), py3(i))
dxt = amax1(xmax-xmin, ymax-ymin)
if (dxt > 0.) go to 10
ir = irmax
go to 20
10 ir = l2e*log(hgr/dxt)  
ir = min(ir, irmax)
20 irp1 = ir+1
nct(irp1) = nct(irp1)+1
hr = hgr/ip2(irp1)
xmid = .5*(xmin+xmax)
id = (xmid-a)/hr+1.5
ymid = .5*(ymin+ymax)
jd = (ymid-c)/hr+1.5
ijd = ip2(irp1)+1
isx = id+(jd-1)*ijd+ird(irp1)
ifx = ifinal(isx)
if (ifx > 0) go to 50
istart(isx) = i
go to 60
50 next(ifx) = i 
60 ifinal(isx) = i
100 continue
!      write(6, 106) tcat, (irp1, nct(irp1), irp1=1, irmp1)
106 format(' time to categorize   ', e15.6/(' ir+1', i3, ' ntri', i7))
!
!     sort triangles into boxes
!     
l = 0
do 30 irp1=1, irmp1
if (nct(irp1) == 0) go to 30
ist = ird(irp1)+1    
isd = ip2(irp1)+1
call box(isd, istart(ist), next, l, ifinal)
last(irp1) = l+1
30 continue
do 35 irp1=1, irmp1
il = ird(irp1)+(ip2(irp1)+1)**2+1
if (istart(il) == 0) istart(il) = last(irp1)
35 continue
!      write(6, 31) tsort, l, ntri
31 format(' time to sort  ', e15.6, '   l', i8, '   ntri', i8)
do 90 k=1, ntri
vx1(k) = px2(k)-px1(k)
vy1(k) = py2(k)-py1(k)
vx2(k) = px3(k)-px2(k)
vy2(k) = py3(k)-py2(k)
vx3(k) = px1(k)-px3(k)
vy3(k) = py1(k)-py3(k)
90 continue
tl1 = 0.
tl2 = 0.
maxs = 0
do 500 ir2=1, irmp1
if (nct(ir2) == 0) go to 500
ist = ird(ir2)    
isd = ip2(ir2)+1
do 490 j2=1, isd
do 480 i2=1, isd
ist = ist+1
ls = istart(ist)
lf = istart(ist+1)-1
if (lf < ls) go to 480
!
!     define coverings
!
kcv = 0
i2m = i2-1
j2m = j2-1
do 300 ir1=1, irmp1
if (nct(ir1) == 0) go to 300
if (ir1 >= ir2) go to 260
irdp = 2**(ir2-ir1)
i1s = (i2m-1)/irdp
i1f = (i2m+1)/irdp
if = i2m+1-i1f*irdp
if (if > 0) i1f = i1f+1
j1s = (j2m-1)/irdp
j1f = (j2m+1)/irdp
jf = j2m+1-j1f*irdp
if (jf > 0) j1f = j1f+1
go to 270
260 irdp = 2**(ir1-ir2)
i1s = irdp*(i2m-1)
i1f = irdp*(i2m+1)
j1s = irdp*(j2m-1)
j1f = irdp*(j2m+1)
270 ijd = ip2(ir1)+1
i1s = max(i1s+1, 1)
i1f = min(i1f+1, ijd)
j1s = max(j1s+1, 1)
j1f = min(j1f+1, ijd)
ixh = (j1s-2)*ijd+ird(ir1)
ixs = i1s+ixh
ixf = i1f+ixh
do 290 j1=j1s, j1f
ixs = ixs+ijd
kds = istart(ixs)
ixf = ixf+ijd
kdf = istart(ixf+1)-1
if (kdf < kds) go to 290
do 280 kd=kds, kdf
kcv = kcv+1
kh(kcv) = ifinal(kd) 
280 continue
290 continue
300 continue
do 310 icv=1, 10    
if (kcv <= ncv(icv)) go to 310
ncv(icv) = kcv 
go to 320
310 continue
!
!
320 do 470 ldo=ls, lf
l = ifinal(ldo)
ith = itype(l)
if (ith == 0) go to 470
ltp(1) = 0
ltp(2) = 0
ltp(3) = 0
id1 = ith/100 
ith = ith-100*id1
id2 = ith/10 
id3 = ith-10*id2
if (id1 /= 0) ltp(id1) = 1
if (id2 /= 0) ltp(id2) = 1
if (id3 /= 0) ltp(id3) = 1
!     if ((ith.eq.123) .or. (ith.eq.12) .or.(ith.eq.13)) ltp(1) = 1
!     if ((ith.eq.123) .or. (ith.eq.23) .or.(ith.eq.12)) ltp(2) = 1
!     if ((ith.eq.123) .or. (ith.eq.13) .or.(ith.eq.23)) ltp(3) = 1
do 460 ns=1, 3
go to (101, 102, 103), ns
101 if (ltp(ns) == 0) go to 460
px4 = px1(l)
py4 = py1(l)
px5 = px2(l)
py5 = py2(l)
x4 = x1(l)
y4 = y1(l)
z4 = z1(l)
x5 = x2(l)
y5 = y2(l)
z5 = z2(l)
go to 105
102 if (ltp(ns) == 0) go to 460
px4 = px2(l)
py4 = py2(l)
px5 = px3(l)
py5 = py3(l)
x4 = x2(l)
y4 = y2(l)
z4 = z2(l)
x5 = x3(l)
y5 = y3(l)
z5 = z3(l)
go to 105
103 if (ltp(ns) == 0) go to 460
px4 = px1(l)
py4 = py1(l)
px5 = px3(l)
py5 = py3(l)
x4 = x1(l)
y4 = y1(l)
z4 = z1(l)
x5 = x3(l)
y5 = y3(l)
z5 = z3(l)
105 x54 = px5-px4
y54 = py5-py4
nseg = 0
do 440 kd=1, kcv
k = kh(kd) 
c17 = vx1(k)*y54-vy1(k)*x54
c27 = vx2(k)*y54-vy2(k)*x54
c37 = vx3(k)*y54-vy3(k)*x54
c14 = vy1(k)*(px4-px1(k))-vx1(k)*(py4-py1(k))
c25 = vy2(k)*(px4-px2(k))-vx2(k)*(py4-py2(k))
c36 = vy3(k)*(px4-px3(k))-vx3(k)*(py4-py3(k))
tmin = 0.
tmax = 1.
if (c17< 0) then
    goto 151
else if (c17 == 0) then 
    goto 152
else 
    goto 153
end if
151 tmax = amin1(c14/c17, tmax)   
go to 154
152 if (c14< 0) then
    goto 154
else if (c14 == 0) then 
    goto 440
else 
    goto 440
end if
153 tmin = amax1(c14/c17, tmin)
154 if (c27< 0) then
        goto 155
    else if (c27 == 0) then 
        goto 156
    else 
        goto 157
    end if
155 tmax = amin1(c25/c27, tmax)   
go to 158
156 if (c25< 0) then
    goto 158
else if (c25 == 0) then 
    goto 440
else 
    goto 440
end if
157 tmin = amax1(c25/c27, tmin)
158 if (c37< 0) then
        goto 159
    else if (c37 == 0) then 
        goto 160
    else 
        goto 161
    end if
159 tmax = amin1(c36/c37, tmax)   
go to 162
160 if (c36< 0) then
    goto 162
else if (c36 == 0) then 
    goto 440
else 
    goto 440
end if
161 tmin = amax1(c36/c37, tmin)
162 if (tmax-tmin < .00001) go to 440
xpl = x4+tmin*(x5-x4)
ypl = y4+tmin*(y5-y4)
zpl = z4+tmin*(z5-z4)
xpr = x4+tmax*(x5-x4)
ypr = y4+tmax*(y5-y4)
zpr = z4+tmax*(z5-z4)
!
!     the projections of line and plane intersect
!     now determine if plane covers line
!
vx1t = x2(k)-x1(k)
vy1t = y2(k)-y1(k)
vz1t = z2(k)-z1(k)
vx2t = x3(k)-x1(k)
vy2t = y3(k)-y1(k)
vz2t = z3(k)-z1(k)
apl = vy1t*vz2t-vy2t*vz1t
bpl = vx2t*vz1t-vx1t*vz2t
cpl = vx1t*vy2t-vx2t*vy1t
dpl = apl*x1(k)+bpl*y1(k)+cpl*z1(k)
vx3t = xpl-xeye
vy3t = ypl-yeye
vz3t = zpl-zeye
den = apl*vx3t+bpl*vy3t+cpl*vz3t
til = 0.
if (den == 0.) go to 410
til = (dpl-apl*xeye-bpl*yeye-cpl*zeye)/den
410 vx3t = xpr-xeye
vy3t = ypr-yeye
vz3t = zpr-zeye
den = apl*vx3t+bpl*vy3t+cpl*vz3t
tir = 0.
if (den == 0.) go to 412
tir = (dpl-apl*xeye-bpl*yeye-cpl*zeye)/den
412 if (til>=.99999.and.tir>=.99999) go to 440
if (til<1..and.tir<1.) go to 164
vx3t = xpr-xpl
vy3t = ypr-ypl
vz3t = zpr-zpl
den = apl*vx3t+bpl*vy3t+cpl*vz3t
tim = 0.
if (den == 0.) go to 414
tim = (dpl-apl*xpl-bpl*ypl-cpl*zpl)/den
414 thold = tmin+tim*(tmax-tmin)
if (til>=1.) go to 163
tmax = thold
go to 164
163 tmin = thold
164 nseg = nseg+1
tl(nseg) = tmin
tr(nseg) = tmax
440 continue
maxs = max(maxs, nseg)
if (nseg-1< 0) then
    goto 171
else if (nseg-1 == 0) then 
    goto 180
else 
    goto 172
end if
171 call line(px4, py4, px5, py5)
go to 460
!
!     order the segments according to left end point tl(k)
!
172 do 173 k=2, nseg
do 173 i=k, nseg
if (tl(k-1)<=tl(i)) go to 173
tlh = tl(k-1)
trh = tr(k-1)
tl(k-1) = tl(i)
tr(k-1) = tr(i)
tl(i) = tlh
tr(i) = trh
173 continue
!
!     eliminate segment overlap
!
k1 = 1
k2 = 1
174 k2 = k2+1
if (k2>nseg) go to 176
if (tr(k1)<tl(k2)) go to 175
tr(k1) = amax1(tr(k1), tr(k2))
go to 174
175 k1 = k1+1
tl(k1) = tl(k2)
tr(k1) = tr(k2)
go to 174
176 nseg = k1
!
!     plot all segments of the line
!
180 do 181 ks =1, nseg
kb = nseg-ks+1
tl(kb+1) = tr(kb)
tr(kb) = tl(kb)
181 continue
tl(1) = 0.
tr(nseg+1) = 1.
nsegp = nseg+1
do 450 k=1, nsegp
if (abs(tr(k)-tl(k))<.000001) go to 450
xa = px4+tl(k)*(px5-px4)
ya = py4+tl(k)*(py5-py4)
xb = px4+tr(k)*(px5-px4)
yb = py4+tr(k)*(py5-py4)
call line(xa, ya, xb, yb)
450 continue
460 continue
470 continue
480 continue
490 continue
500 continue
!      write(6, 903) tl1, tl2
903 format(' time to cover', e15.6/ &
       ' time to test ', e15.6)
!      write(6, 904) maxs
904 format(' maximum number of segments', i5)
!      write(6, 250) (ncv(icv), icv=1, 10)
250 format('  the ten largest coverings'/(10i5))
call frame
end subroutine vsurf1
subroutine prjct(init, xeye, yeye, zeye, x, y, z, px, py)
!
!     subroutine prjct projects the point x, y, z onto a plane through
!     the origin that is perpendicular to a line between the origin
!     and the eye. the projection is along the line between the eye
!     and the point x, y, z. px and py are the coordinates of the
!     projection in the plane.
!     (version 2 , 12-10-82)
!
save
if (init/=0) go to 1
rads1 = xeye**2+yeye**2
rads2 = rads1+zeye**2
d1 = sqrt(rads1)
d2 = sqrt(rads2)
cx1 = -yeye/d1
cy1 = xeye/d1
cx2 = -xeye*zeye/(d1*d2)
cy2 = -yeye*zeye/(d1*d2)
cz2 = d1/d2
cx3 = xeye/d2
cy3 = yeye/d2
cz3 = zeye/d2
return
1 x1 = cx1*x+cy1*y
y1 = cx2*x+cy2*y+cz2*z
z1 = cx3*x+cy3*y+cz3*z
ratio = d2/(d2-z1)
px = ratio*x1
py = ratio*y1
return
end subroutine prjct
subroutine box(isd, istart, next, l, list)
dimension istart(isd, isd), next(1), list(1)
do 30 jd=1, isd
do 10 id=1, isd
idx = istart(id, jd)
istart(id, jd) = l+1
if (idx == 0) go to 10
20 l = l+1
list(l) = idx
if (next(idx) == 0) go to 10
idx = next(idx)
go to 20
10 continue
30 continue
return
end subroutine box
integer function icvmg(i1, i2, r)
integer i1, i2
real    r 
!
!     returns i1 if i3.ge.0 and returns i2 if i3.lt.0 .
!
icvmg = i1
if (r < 0.) icvmg = i2
return
end function icvmg
subroutine projct(m, n, xeye, yeye, zeye, x, y, z, px, py)
!     ****     projects point (x, y, z) onto plane thru origin and perp
!     ****     to line joining origin and eye
dimension x(n, m), y(n, m), z(n, m), px(n, m), py(n, m)
call prjct(0, xeye, yeye, zeye, rdum1, rdum2, rdum3, rdum4, rdum5)
do 100 i=1, m
do 100 j=1, n
call prjct(1, xeye, yeye, zeye, x(j, i), y(j, i), z(j, i), px(j, i), py(j, i))
100 continue
return
end subroutine projct
