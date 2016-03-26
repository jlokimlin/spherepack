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
! ... file shags.f
!
!     this file contains code and documentation for subroutines
!     shags and shagsi
!
! ... files which must be loaded with shags.f
!
!     sphcom.f, hrfft.f, gaqd.f
!
!     subroutine shags(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
!    1                    wshags,lshags,work,lwork,ierror)
!
!     subroutine shags performs the spherical harmonic analysis
!     on the array g and stores the result in the arrays a and b.
!     the analysis is performed on a gaussian grid in colatitude
!     and an equally spaced grid in longitude.  the associated
!     legendre functions are stored rather than recomputed as they
!     are in subroutine shagc.  the analysis is described below
!     at output parameters a,b.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0,pi) and are compu
!            in radians in theta(1),...,theta(nlat) by subroutine gaqd.
!            if nlat is odd the equator will be included as the grid poi
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!     isym   = 0  no symmetries exist about the equator. the analysis
!                 is performed on the entire sphere.  i.e. on the
!                 array g(i,j) for i=1,...,nlat and j=1,...,nlon.
!                 (see description of g below)
!
!            = 1  g is antisymmetric about the equator. the analysis
!                 is performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the analysis is performed on the
!                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
!                 if nlat is even the analysis is performed on the
!                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!
!
!            = 2  g is symmetric about the equator. the analysis is
!                 performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the analysis is performed on the
!                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
!                 if nlat is even the analysis is performed on the
!                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!
!     nt     the number of analyses.  in the program that calls shags,
!            the arrays g,a and b can be three dimensional in which
!            case multiple analyses will be performed.  the third
!            index is the analysis index which assumes the values
!            k=1,...,nt.  for a single analysis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that the arrays g,a and b
!            have only two dimensions.
!
!     g      a two or three dimensional array (see input parameter
!            nt) that contains the discrete function to be analyzed.
!            g(i,j) contains the value of the function at the gaussian
!            point theta(i) and longitude point phi(j) = (j-1)*2*pi/nlon
!            the index ranges are defined above at the input parameter
!            isym.
!
!     idg    the first dimension of the array g as it appears in the
!            program that calls shags. if isym equals zero then idg
!            must be at least nlat.  if isym is nonzero then idg must
!            be at least nlat/2 if nlat is even or at least (nlat+1)/2
!            if nlat is odd.
!
!     jdg    the second dimension of the array g as it appears in the
!            program that calls shags. jdg must be at least nlon.
!
!     mdab   the first dimension of the arrays a and b as it appears
!            in the program that calls shags. mdab must be at least
!            min((nlon+2)/2,nlat) if nlon is even or at least
!            min((nlon+1)/2,nlat) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears
!            in the program that calls shags. ndab must be at least nlat
!
!     wshags an array which must be initialized by subroutine shagsi.
!            once initialized, wshags can be used repeatedly by shags
!            as long as nlat and nlon remain unchanged.  wshags must
!            not be altered between calls of shags.
!
!     lshags the dimension of the array wshags as it appears in the
!            program that calls shags. define
!
!               l1 = min(nlat,(nlon+2)/2) if nlon is even or
!               l1 = min(nlat,(nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshags must be at least
!
!            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
!
!     work   a real work space which need not be saved
!
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shags. define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!
!            if isym is zero then lwork must be at least
!
!                  nlat*nlon*(nt+1)
!
!            if isym is nonzero then lwork must be at least
!
!                  l2*nlon*(nt+1)
!
!     **************************************************************
!
!     output parameters
!
!     a,b    both a,b are two or three dimensional arrays (see input
!            parameter nt) that contain the spherical harmonic
!            coefficients in the representation of g(i,j) given in the
!            discription of subroutine shags. for isym=0, a(m,n) and
!            b(m,n) are given by the equations listed below. symmetric
!            versions are used when isym is greater than zero.
!
!     definitions
!
!     1. the normalized associated legendre functions
!
!     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                       *sin(theta)**m/(2**n*factorial(n)) times the
!                       (n+m)th derivative of (x**2-1)**n with respect
!                       to x=cos(theta).
!
!     2. the fourier transform of g(i,j).
!
!     c(m,i)          = 2/nlon times the sum from j=1 to j=nlon of
!                       g(i,j)*cos((m-1)*(j-1)*2*pi/nlon)
!                       (the first and last terms in this sum
!                       are divided by 2)
!
!     s(m,i)          = 2/nlon times the sum from j=2 to j=nlon of
!                       g(i,j)*sin((m-1)*(j-1)*2*pi/nlon)
!
!
!     3. the gaussian points and weights on the sphere
!        (computed by subroutine gaqd).
!
!        theta(1),...,theta(nlat) (gaussian pts in radians)
!        wts(1),...,wts(nlat) (corresponding gaussian weights)
!
!
!     4. the maximum (plus one) longitudinal wave number
!
!            mmax = min(nlat,(nlon+2)/2) if nlon is even or
!            mmax = min(nlat,(nlon+1)/2) if nlon is odd.
!
!
!     then for m=0,...,mmax-1 and n=m,...,nlat-1 the arrays a,b
!     are given by
!
!     a(m+1,n+1)     =  the sum from i=1 to i=nlat of
!                       c(m+1,i)*wts(i)*pbar(m,n,theta(i))
!
!     b(m+1,n+1)      = the sum from i=1 to nlat of
!                       s(m+1,i)*wts(i)*pbar(m,n,theta(i))
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of isym
!            = 4  error in the specification of nt
!            = 5  error in the specification of idg
!            = 6  error in the specification of jdg
!            = 7  error in the specification of mdab
!            = 8  error in the specification of ndab
!            = 9  error in the specification of lshags
!            = 10 error in the specification of lwork
!
!
! ****************************************************************
!
!     subroutine shagsi(nlat,nlon,wshags,lshags,work,lwork,dwork,ldwork,
!    +                  ierror)
!
!     subroutine shagsi initializes the array wshags which can then
!     be used repeatedly by subroutines shags. it precomputes
!     and stores in wshags quantities such as gaussian weights,
!     legendre polynomial coefficients, and fft trigonometric tables.
!
!     input parameters
!
!     nlat   the number of points in the gaussian colatitude grid on the
!            full sphere. these lie in the interval (0,pi) and are compu
!            in radians in theta(1),...,theta(nlat) by subroutine gaqd.
!            if nlat is odd the equator will be included as the grid poi
!            theta((nlat+1)/2).  if nlat is even the equator will be
!            excluded as a grid point and will lie half way between
!            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3.
!            note: on the half sphere, the number of grid points in the
!            colatitudinal direction is nlat/2 if nlat is even or
!            (nlat+1)/2 if nlat is odd.
!
!     nlon   the number of distinct londitude points.  nlon determines
!            the grid increment in longitude as 2*pi/nlon. for example
!            nlon = 72 for a five degree grid. nlon must be greater
!            than or equal to 4. the efficiency of the computation is
!            improved when nlon is a product of small prime numbers.
!
!     wshags an array which must be initialized by subroutine shagsi.
!            once initialized, wshags can be used repeatedly by shags
!            as long as nlat and nlon remain unchanged.  wshags must
!            not be altered between calls of shags.
!
!     lshags the dimension of the array wshags as it appears in the
!            program that calls shags. define
!
!               l1 = min(nlat,(nlon+2)/2) if nlon is even or
!               l1 = min(nlat,(nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshags must be at least
!
!            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
!
!     work   a real work space which need not be saved
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shagsi. lwork must be at least
!            4*nlat*(nlat+2)+2 in the routine calling shagsi
!
!     dwork   a real work array that does not have to be saved.
!
!     ldwork  the length of dwork in the calling routine.  ldwork must
!             be at least nlat*(nlat+4)
!
!     output parameter
!
!     wshags an array which must be initialized before calling shags or
!            once initialized, wshags can be used repeatedly by shags or
!            as long as nlat and nlon remain unchanged.  wshags must not
!            altered between calls of shasc.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshags
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!            = 6  failure in gaqd to compute gaussian points
!                 (due to failure in eigenvalue routine)
!
!
! ****************************************************************
subroutine shags(nlat,nlon,mode,nt,g,idg,jdg,a,b,mdab,ndab, &
    wshags,lshags,work,lwork,ierror)
    !     subroutine shags performs the spherical harmonic analysis on
    !     a gaussian grid on the array(s) in g and returns the coefficients
    !     in array(s) a,b. the necessary legendre polynomials are fully
    !     stored in this version.
    !
    dimension g(idg,jdg,1),a(mdab,ndab,1),b(mdab,ndab,1), &
        wshags(lshags),work(lwork)
    !     check input parameters
    ierror = 1
    if (nlat<3) return
    ierror = 2
    if (nlon<4) return
    ierror = 3
    if (mode<0 .or.mode>2) return
    !     set m limit for pmn
    l = min((nlon+2)/2,nlat)
    !     set gaussian point nearest equator pointer
    late = (nlat+mod(nlat,2))/2
    !     set number of grid points for analysis/synthesis
    lat = nlat
    if (mode/=0) lat = late
    ierror = 4
    if (nt<1) return
    ierror = 5
    if (idg<lat) return
    ierror = 6
    if (jdg<nlon) return
    ierror = 7
    if(mdab < l) return
    ierror = 8
    if(ndab < nlat) return
    l1 = l
    l2 = late
    ierror = 9
    !     check permanent work space length
    !
    lp= nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    if(lshags<lp) return
    ierror = 10
    !     check temporary work space length
    if (mode==0 .and. lwork<nlat*nlon*(nt+1)) return
    if (mode/=0 .and. lwork<l2*nlon*(nt+1)) return
    ierror = 0
    !     set starting address for gaussian wts ,fft values,
    !     and fully stored legendre polys in wshags
    iwts = 1
    ifft = nlat+2*nlat*late+3*(l*(l-1)/2+(nlat-l)*(l-1))+1
    ipmn = ifft+nlon+15
    !     set pointer for internal storage of g
    iw = lat*nlon*nt+1

    call shags1(nlat,nlon,l,lat,mode,g,idg,jdg,nt,a,b,mdab,ndab, &
        wshags(iwts),wshags(ifft),wshags(ipmn),late,work,work(iw))

end subroutine shags



subroutine shags1(nlat,nlon,l,lat,mode,gs,idg,jdg,nt,a,b,mdab, &
    ndab,wts,wfft,pmn,late,g,work)
    dimension gs(idg,jdg,nt),a(mdab,ndab,nt), &
        b(mdab,ndab,nt),g(lat,nlon,nt)
    dimension wfft(1),pmn(late,1),wts(nlat),work(1)

    !     set gs array internally in shags1
    do k=1,nt
        do j=1,nlon
            do i=1,lat
                g(i,j,k) = gs(i,j,k)
            end do
        end do
    end do

    !     do fourier transform
    do k=1,nt
        call hrfftf(lat,nlon,g(1,1,k),lat,wfft,work)
    end do

    !     scale result
    sfn = 2.0/nlon
    do k=1,nt
        do j=1,nlon
            do i=1,lat
                g(i,j,k) = sfn*g(i,j,k)
            end do
        end do
    end do

    !     compute using gaussian quadrature
    !     a(n,m) = s (ga(theta,m)*pnm(theta)*sin(theta)*dtheta)
    !     b(n,m) = s (gb(theta,m)*pnm(theta)*sin(theta)*dtheta)
    !     here ga,gb are the cos(phi),sin(phi) coefficients of
    !     the fourier expansion of g(theta,phi) in phi.  as a result
    !     of the above fourier transform they are stored in array
    !     g as follows:
    !     for each theta(i) and k= l-1
    !     ga(0),ga(1),gb(1),ga(2),gb(2),...,ga(k-1),gb(k-1),ga(k)
    !     correspond to
    !     g(i,1),g(i,2),g(i,3),g(i,4),g(i,5),...,g(i,2l-4),g(i,2l-3),g(i,2l-2)
    !     whenever 2*l-2 = nlon exactly
    !     initialize coefficients to zero
    do k=1,nt
        do np1=1,nlat
            do mp1=1,l
                a(mp1,np1,k) = 0.0
                b(mp1,np1,k) = 0.0
            end do
        end do
    end do

    !     set mp1 limit on b(mp1) calculation
    lm1 = l
    if (nlon == l+l-2) then
        lm1 = l-1
    end if

    if (mode==0) then
        !     for full sphere (mode=0) and even/odd reduction:
        !     overwrite g(i) with (g(i)+g(nlat-i+1))*wts(i)
        !     overwrite g(nlat-i+1) with (g(i)-g(nlat-i+1))*wts(i)
        nl2 = nlat/2
        do k=1,nt
            do j=1,nlon
                do i=1,nl2
                    is = nlat-i+1
                    t1 = g(i,j,k)
                    t2 = g(is,j,k)
                    g(i,j,k) = wts(i)*(t1+t2)
                    g(is,j,k) = wts(i)*(t1-t2)
                end do
                !     adjust equator if necessary(nlat odd)
                if (mod(nlat,2)/=0) g(late,j,k) = wts(late)*g(late,j,k)
            end do
        end do
        !     set m = 0 coefficients first
        mp1 = 1
        m = 0
        mml1 = m*(2*nlat-m-1)/2
        do k=1,nt
            do i=1,late
                is = nlat-i+1
                do np1=1,nlat,2
                    !     n even
                    a(1,np1,k) = a(1,np1,k)+g(i,1,k)*pmn(i,mml1+np1)
                end do
                do np1=2,nlat,2
                    !     n odd
                    a(1,np1,k) = a(1,np1,k)+g(is,1,k)*pmn(i,mml1+np1)
                end do
            end do
        end do
        !     compute m.ge.1  coefficients next
        do mp1=2,lm1
            m = mp1-1
            mml1 = m*(2*nlat-m-1)/2
            mp2 = mp1+1
            do k=1,nt
                do i=1,late
                    is = nlat-i+1
                    !     n-m even
                    do np1=mp1,nlat,2
                        a(mp1,np1,k) = a(mp1,np1,k)+g(i,2*m,k)*pmn(i,mml1+np1)
                        b(mp1,np1,k) = b(mp1,np1,k)+g(i,2*m+1,k)*pmn(i,mml1+np1)
                    end do
                    !     n-m odd
                    do np1=mp2,nlat,2
                        a(mp1,np1,k) = a(mp1,np1,k)+g(is,2*m,k)*pmn(i,mml1+np1)
                        b(mp1,np1,k) = b(mp1,np1,k)+g(is,2*m+1,k)*pmn(i,mml1+np1)
                    end do
                end do
            end do
        end do

        if (nlon == l+l-2) then
            !     compute m=l-1, n=l-1,l,...,nlat-1 coefficients
            m = l-1
            mml1 = m*(2*nlat-m-1)/2
            do k=1,nt
                do i=1,late
                    is = nlat-i+1
                    do np1=l,nlat,2
                        mn = mml1+np1
                        a(l,np1,k) = a(l,np1,k)+0.5*g(i,nlon,k)*pmn(i,mn)
                    end do
                    !     n-m  odd
                    lp1 = l+1
                    do np1=lp1,nlat,2
                        mn = mml1+np1
                        a(l,np1,k) = a(l,np1,k)+0.5*g(is,nlon,k)*pmn(i,mn)
                    end do
                end do
            end do
        end if
    else
        !     half sphere
        !     overwrite g(i) with wts(i)*(g(i)+g(i)) for i=1,...,nlate/2
        nl2 = nlat/2
        do  k=1,nt
            do j=1,nlon
                do i=1,nl2
                    g(i,j,k) = wts(i)*(g(i,j,k)+g(i,j,k))
                end do
                !     adjust equator separately if a grid point
                if (nl2<late) g(late,j,k) = wts(late)*g(late,j,k)
            end do
        end do
        !     set m = 0 coefficients first
        mp1 = 1
        m = 0
        mml1 = m*(2*nlat-m-1)/2
        ms = 1
        if (mode==1) ms = 2
        do k=1,nt
            do i=1,late
                do np1=ms,nlat,2
                    a(1,np1,k) = a(1,np1,k)+g(i,1,k)*pmn(i,mml1+np1)
                end do
            end do
        end do
        !     compute m.ge.1  coefficients next
        do mp1=2,lm1
            m = mp1-1
            mml1 = m*(2*nlat-m-1)/2
            ms = mp1
            if (mode==1) ms = mp1+1
            do k=1,nt
                do  i=1,late
                    do np1=ms,nlat,2
                        a(mp1,np1,k) = a(mp1,np1,k)+g(i,2*m,k)*pmn(i,mml1+np1)
                        b(mp1,np1,k) = b(mp1,np1,k)+g(i,2*m+1,k)*pmn(i,mml1+np1)
                    end do
                end do
            end do
        end do

        if (nlon==l+l-2) then
            !     compute n=m=l-1 coefficients last
            m = l-1
            mml1 = m*(2*nlat-m-1)/2
            !     set starting n for mode even
            ns = l
            !     set starting n for mode odd
            if (mode==1) ns = l+1
            do k=1,nt
                do i=1,late
                    do np1=ns,nlat,2
                        mn = mml1+np1
                        a(l,np1,k) = a(l,np1,k)+0.5*g(i,nlon,k)*pmn(i,mn)
                    end do
                end do
            end do
        end if
    end if

end subroutine shags1



subroutine shagsi(nlat,nlon,wshags,lshags,work,lwork,dwork,ldwork, &
    ierror)
    !
    !     this subroutine must be called before calling shags or shsgs with
    !     fixed nlat,nlon. it precomputes the gaussian weights, points
    !     and all necessary legendre polys and stores them in wshags.
    !     these quantities must be preserved when calling shags or shsgs
    !     repeatedly with fixed nlat,nlon.  dwork must be of length at
    !     least nlat*(nlat+4) in the routine calling shagsi.  This is
    !     not checked.  undetectable errors will result if dwork is
    !     smaller than nlat*(nlat+4).
    !
    dimension wshags(lshags),work(lwork)
    real dwork(ldwork)
    ierror = 1
    if (nlat<3) return
    ierror = 2
    if (nlon<4) return
    !     set triangular truncation limit for spherical harmonic basis
    l = min((nlon+2)/2,nlat)
    !     set equator or nearest point (if excluded) pointer
    late = (nlat+1)/2
    l1 = l
    l2 = late
    !     check permanent work space length
    ierror = 3
    lp=nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    if(lshags<lp) return
    ierror = 4
    !     check temporary work space
    if (lwork<4*nlat*(nlat+2)+2) return
    ierror = 5
    !     check real temporary space
    if (ldwork < nlat*(nlat+4)) return
    ierror = 0
    !     set preliminary quantites needed to compute and store legendre polys
    ldw = nlat*(nlat+4)
    call shagsp(nlat,nlon,wshags,lshags,dwork,ldwork,ierror)
    if (ierror/=0) return
    !     set legendre poly pointer in wshags
    ipmnf = nlat+2*nlat*late+3*(l*(l-1)/2+(nlat-l)*(l-1))+nlon+16
    call shagss1(nlat,l,late,wshags,work,wshags(ipmnf))

end subroutine shagsi


subroutine shagss1(nlat,l,late,w,pmn,pmnf)
    dimension w(1),pmn(nlat,late,3),pmnf(late,1)
    !     compute and store legendre polys for i=1,...,late,m=0,...,l-1
    do i=1,nlat
        do j=1,late
            do k=1,3
                pmn(i,j,k) = 0.0
            end do
        end do
    end do
    do mp1=1,l
        m = mp1-1
        mml1 = m*(2*nlat-m-1)/2
        !     compute pmn for n=m,...,nlat-1 and i=1,...,(l+1)/2
        mode = 0
        call legin(mode,l,nlat,m,w,pmn,km)
        !     store above in pmnf
        do np1=mp1,nlat
            mn = mml1+np1
            do i=1,late
                pmnf(i,mn) = pmn(np1,i,km)

            end do
        end do
    end do

end subroutine shagss1


subroutine shagsp(nlat,nlon,wshags,lshags,dwork,ldwork,ierror)
    dimension wshags(lshags)
    real dwork(ldwork)

    ierror = 1
    if (nlat<3) return
    ierror = 2
    if (nlon<4) return
    !     set triangular truncation limit for spherical harmonic basis
    l = min((nlon+2)/2,nlat)
    !     set equator or nearest point (if excluded) pointer
    late = (nlat+mod(nlat,2))/2
    l1 = l
    l2 = late
    ierror = 3
    !     check permanent work space length
    if (lshags < nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15)return
    ierror = 4
    !     if (lwork.lt.4*nlat*(nlat+2)+2) return
    if (ldwork<nlat*(nlat+4))return
    ierror = 0
    !     set pointers
    i1 = 1
    i2 = i1+nlat
    i3 = i2+nlat*late
    i4 = i3+nlat*late
    i5 = i4+l*(l-1)/2 +(nlat-l)*(l-1)
    i6 = i5+l*(l-1)/2 +(nlat-l)*(l-1)
    i7 = i6+l*(l-1)/2 +(nlat-l)*(l-1)
    !     set indices in temp work for real gaussian wts and pts
    idth = 1
    !     idwts = idth+2*nlat
    !     iw = idwts+2*nlat
    idwts = idth+nlat
    iw = idwts+nlat
    call shagsp1(nlat,nlon,l,late,wshags(i1),wshags(i2),wshags(i3), &
        wshags(i4),wshags(i5),wshags(i6),wshags(i7),dwork(idth), &
        dwork(idwts),dwork(iw),ierror)
    if (ierror/=0) ierror = 6

end subroutine shagsp



subroutine shagsp1(nlat,nlon,l,late,wts,p0n,p1n,abel,bbel,cbel, &
    wfft,dtheta,dwts,work,ier)
    dimension wts(nlat),p0n(nlat,late),p1n(nlat,late),abel(1),bbel(1), &
        cbel(1),wfft(1),dtheta(nlat),dwts(nlat)
    real pb,dtheta,dwts,work(*)

    indx(m,n) = (n-1)*(n-2)/2+m-1
    imndx(m,n) = l*(l-1)/2+(n-l-1)*(l-1)+m-1

    call hrffti(nlon,wfft)

    !     compute real gaussian points and weights
    !     lw = 4*nlat*(nlat+2)
    lw = nlat*(nlat+2)

    call gaqd(nlat,dtheta,dwts,work,lw,ier)

    if (ier/=0) return

    !     store gaussian weights single precision to save computation
    !     in inner loops in analysis
    do i=1,nlat
        wts(i) = dwts(i)
    end do
    !     initialize p0n,p1n using real dnlfk,dnlft
    do np1=1,nlat
        do i=1,late
            p0n(np1,i) = 0.0
            p1n(np1,i) = 0.0
        end do
    end do
    !     compute m=n=0 legendre polynomials for all theta(i)
    np1 = 1
    n = 0
    m = 0
    call dnlfk(m,n,work)
    do  i=1,late
        call dnlft(m,n,dtheta(i),work,pb)
        p0n(1,i) = pb
    end do
    !     compute p0n,p1n for all theta(i) when n.gt.0
    do np1=2,nlat
        n = np1-1
        m = 0
        call dnlfk(m,n,work)
        do i=1,late
            call dnlft(m,n,dtheta(i),work,pb)
            p0n(np1,i) = pb
        end do
        !     compute m=1 legendre polynomials for all n and theta(i)
        m = 1
        call dnlfk(m,n,work)
        do i=1,late
            call dnlft(m,n,dtheta(i),work,pb)
            p1n(np1,i) = pb
        end do
    end do
    !
    !     compute and store swarztrauber recursion coefficients
    !     for 2.le.m.le.n and 2.le.n.le.nlat in abel,bbel,cbel
    do n=2,nlat
        mlim = min(n,l)
        do  m=2,mlim
            imn = indx(m,n)
            if (n>=l) imn = imndx(m,n)
            abel(imn)=sqrt(real((2*n+1)*(m+n-2)*(m+n-3))/ &
                real(((2*n-3)*(m+n-1)*(m+n))))
            bbel(imn)=sqrt(real((2*n+1)*(n-m-1)*(n-m))/ &
                real(((2*n-3)*(m+n-1)*(m+n))))
            cbel(imn)=sqrt(real((n-m+1)*(n-m+2))/ &
                real(((n+m-1)*(n+m))))
        end do
    end do

end subroutine shagsp1
