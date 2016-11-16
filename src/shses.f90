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
! ... file shses.f90
!
!     this file contains code and documentation for subroutines
!     shses and shsesi
!
! ... files which must be loaded with shses.f90
!
!     type_SpherepackAux.f90, type_HFFTpack.f90
!
!     subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab,
!    +                 wshses,lshses,work,lwork,ierror)
!
!     subroutine shses performs the spherical harmonic synthesis
!     on the arrays a and b and stores the result in the array g.
!     the synthesis is performed on an equally spaced grid.  the
!     associated legendre functions are stored rather than recomputed
!     as they are in subroutine shsec.  the synthesis is described
!     below at output parameter g.
!
! *** required files from spherepack2
!
!     type_SpherepackAux.f90, type_HFFTpack.f90
!
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
!     isym   = 0  no symmetries exist about the equator. the synthesis
!                 is performed on the entire sphere.  i.e. on the
!                 array g(i,j) for i=1,...,nlat and j=1,...,nlon.
!                 (see description of g below)
!
!            = 1  g is antisymmetric about the equator. the synthesis
!                 is performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the synthesis is performed on the
!                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
!                 if nlat is even the synthesis is performed on the
!                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!
!
!            = 2  g is symmetric about the equator. the synthesis is
!                 performed on the northern hemisphere only.  i.e.
!                 if nlat is odd the synthesis is performed on the
!                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.
!                 if nlat is even the synthesis is performed on the
!                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon.
!
!     nt     the number of syntheses.  in the program that calls shses,
!            the arrays g,a and b can be three dimensional in which
!            case multiple syntheses will be performed.  the third
!            index is the synthesis index which assumes the values
!            k=1,...,nt.  for a single synthesis set nt=1. the
!            discription of the remaining parameters is simplified
!            by assuming that nt=1 or that the arrays g,a and b
!            have only two dimensions.
!
!     idg    the first dimension of the array g as it appears in the
!            program that calls shses.  if isym equals zero then idg
!            must be at least nlat.  if isym is nonzero then idg
!            must be at least nlat/2 if nlat is even or at least
!            (nlat+1)/2 if nlat is odd.
!
!     jdg    the second dimension of the array g as it appears in the
!            program that calls shses.  jdg must be at least nlon.
!
!     a,b    two or three dimensional arrays (see the input parameter
!            nt) that contain the coefficients in the spherical harmonic
!            expansion of g(i,j) given below at the definition of the
!            output parameter g.  a(m,n) and b(m,n) are defined for
!            indices m=1,...,mmax and n=m,...,nlat where mmax is the
!            maximum (plus one) longitudinal wave number given by
!            mmax = min(nlat,(nlon+2)/2) if nlon is even or
!            mmax = min(nlat,(nlon+1)/2) if nlon is odd.
!
!     mdab   the first dimension of the arrays a and b as it appears
!            in the program that calls shses. mdab must be at least
!            min(nlat,(nlon+2)/2) if nlon is even or at least
!            min(nlat,(nlon+1)/2) if nlon is odd.
!
!     ndab   the second dimension of the arrays a and b as it appears
!            in the program that calls shses. ndab must be at least nlat
!
!     wshses an array which must be initialized by subroutine shsesi.
!            once initialized, wshses can be used repeatedly by shses
!            as long as nlon and nlat remain unchanged.  wshses must
!            not be altered between calls of shses.
!
!     lshses the dimension of the array wshses as it appears in the
!            program that calls shses. define
!
!               l1 = min(nlat,(nlon+2)/2) if nlon is even or
!               l1 = min(nlat,(nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshses must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in the
!            program that calls shses.  define
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            if isym is zero then lwork must be at least
!
!               (nt+1)*nlat*nlon
!
!            if isym is nonzero lwork must be at least
!
!               (nt+1)*l2*nlon.
!
!     **************************************************************
!
!     output parameters
!
!     g      a two or three dimensional array (see input parameter
!            nt) that contains the spherical harmonic synthesis of
!            the arrays a and b at the colatitude point theta(i) =
!            (i-1)*pi/(nlat-1) and longitude point phi(j) =
!            (j-1)*2*pi/nlon. the index ranges are defined above at
!            at the input parameter isym.  for isym=0, g(i,j) is
!            given by the the equations listed below.  symmetric
!            versions are used when isym is greater than zero.
!
!     the normalized associated legendre functions are given by
!
!     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))
!                       *sin(theta)**m/(2**n*factorial(n)) times the
!                       (n+m)th derivative of (x**2-1)**n with respect
!                       to x=cos(theta)
!
!     define the maximum (plus one) longitudinal wave number
!     as   mmax = min(nlat,(nlon+2)/2) if nlon is even or
!          mmax = min(nlat,(nlon+1)/2) if nlon is odd.
!
!     then g(i,j) = the sum from n=0 to n=nlat-1 of
!
!                   .5*pbar(0,n,theta(i))*a(1,n+1)
!
!              plus the sum from m=1 to m=mmax-1 of
!
!                   the sum from n=m to n=nlat-1 of
!
!              pbar(m,n,theta(i))*(a(m+1,n+1)*cos(m*phi(j))
!                                    -b(m+1,n+1)*sin(m*phi(j)))
!
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
!            = 9  error in the specification of lshses
!            = 10 error in the specification of lwork
!
!
! ****************************************************************
!     subroutine shsesi(nlat,nlon,wshses,lshses,work,lwork,dwork,
!    +                  ldwork,ierror)
!
!     subroutine shsesi initializes the array wshses which can then
!     be used repeatedly by subroutine shses.
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
!     lshses the dimension of the array wshses as it appears in the
!            program that calls shsesi. define
!
!               l1 = min(nlat,(nlon+2)/2) if nlon is even or
!               l1 = min(nlat,(nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lshses must be at least
!
!               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
!
!     work   a real   work array that does not have to be saved.
!
!     lwork  the dimension of the array work as it appears in
!            the program that calls shsesi.  define
!
!               l1 = min(nlat,(nlon+2)/2) if nlon is even or
!               l1 = min(nlat,(nlon+1)/2) if nlon is odd
!
!            and
!
!               l2 = nlat/2        if nlat is even or
!               l2 = (nlat+1)/2    if nlat is odd
!
!            then lwork must be at least
!
!               5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
!
!
!     dwork  a real work array that does not have to be saved.
!
!     ldwork the dimension of the array dwork as it appears in the
!            program that calls shsesi.  ldwork must be at least nlat+1
!
!
!     output parameters
!
!     wshses an array which is initialized for use by subroutine shses.
!            once initialized, wshses can be used repeatedly by shses
!            as long as nlon and nlat remain unchanged.  wshses must
!            not be altered between calls of shses.
!
!     ierror = 0  no errors
!            = 1  error in the specification of nlat
!            = 2  error in the specification of nlon
!            = 3  error in the specification of lshses
!            = 4  error in the specification of lwork
!            = 5  error in the specification of ldwork
!
submodule (scalar_synthesis_routines) scalar_synthesis_shses

contains

    module subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab, &
        wshses,lshses,work,lwork,ierror)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: nlon
        integer(ip), intent(in)     :: isym
        integer(ip), intent(in)     :: nt
        real(wp),    intent(out)    :: g(idg,jdg,nt)
        integer(ip), intent(in)     :: idg
        integer(ip), intent(in)     :: jdg
        real(wp),    intent(inout)  :: a(mdab,ndab,nt)
        real(wp),    intent(inout)  :: b(mdab,ndab,nt)
        integer(ip), intent(in)     :: mdab
        integer(ip), intent(in)     :: ndab
        real(wp),    intent(inout)  :: wshses(lshses)
        integer(ip), intent(in)     :: lshses
        real(wp),    intent(inout)  :: work(lwork)
        integer(ip), intent(in)     :: lwork
        integer(ip), intent(out)    :: ierror
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip) :: imid, ist, lpimn, ls, mmax, nln
        !----------------------------------------------------------------------

        mmax = min(nlat,nlon/2+1)
        imid = (nlat+1)/2
        lpimn = (imid*mmax*(nlat+nlat-mmax+1))/2

        if (isym > 0) then
            ls = imid
        else
            ls = nlat
        end if

        select case (isym)
            case (0)
                ist = imid
            case default
                ist = 0
        end select

        nln = nt*ls*nlon

        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 4) then
            ierror = 2
            return
        else if (isym < 0 .or. isym > 2) then
            ierror = 3
            return
        else if (nt < 0) then
            ierror = 4
            return
        else if (&
            (isym == 0 .and. idg < nlat) &
            .or. &
            (isym /= 0 .and. idg < (nlat+1)/2)&
            ) then
            ierror = 5
            return
        else if (jdg < nlon) then
            ierror = 6
            return
        else if (mdab < mmax) then
            ierror = 7
            return
        else if (ndab < nlat) then
            ierror = 8
            return
        else if (lshses < lpimn+nlon+15) then
            ierror = 9
            return
        else if (lwork < nln+ls*nlon) then
            ierror = 10
            return
        else
            ierror = 0
        end if

        associate( &
            iw1 => ist+1, &
            iw2 => nln+1, &
            iw3 => lpimn+1 &
            )
            call shses1(nlat,isym,nt,g,idg,jdg,a,b,mdab,ndab,wshses,imid, &
                ls,nlon,work,work(iw1),work(iw2),wshses(iw3))
        end associate

    end subroutine shses


    module subroutine shsesi(nlat,nlon,wshses,lshses,work,lwork,dwork, &
        ldwork,ierror)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)  :: nlat
        integer(ip), intent(in)  :: nlon
        real(wp),    intent(out) :: wshses(lshses)
        integer(ip), intent(in)  :: lshses
        real(wp),    intent(out) :: work(lwork)
        integer(ip), intent(in)  :: lwork
        real(wp),    intent(out) :: dwork(ldwork)
        integer(ip), intent(in)  :: ldwork
        integer(ip), intent(out) :: ierror
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)         :: imid, labc, lpimn, mmax
        type(HFFTpack)      :: hfft
        type(SpherepackAux) :: sphere_aux
        !----------------------------------------------------------------------

        mmax = min(nlat,nlon/2+1)
        imid = (nlat+1)/2
        lpimn = (imid*mmax*(2*nlat-mmax+1))/2
        labc = 3*((mmax-2)*(2*nlat-mmax-1))/2

        !
        !==> Check validity of input arguments
        !
        if (nlat < 3) then
            ierror = 1
            return
        else if (nlon < 4) then
            ierror = 2
            return
        else if (lshses < lpimn+nlon+15) then
            ierror = 3
            return
        else if (lwork < 5*nlat*imid + labc) then
            ierror = 4
            return
        else if (ldwork < nlat+1) then
            ierror = 5
            return
        else
            ierror = 0
        end if

        associate( &
            iw1 => 3*nlat*imid+1, &
            iw2 => lpimn+1 &
            )
            call sphere_aux%ses1(nlat, nlon, imid, wshses, work, work(iw1), dwork)
            call hfft%initialize(nlon, wshses(iw2))
        end associate

    end subroutine shsesi


    subroutine shses1(nlat,isym,nt,g,idgs,jdgs,a,b,mdab,ndab,p,imid, &
        idg,jdg,ge,go,work,whrfft)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer(ip), intent(in)     :: nlat
        integer(ip), intent(in)     :: isym
        integer(ip), intent(in)     :: nt
        real(wp),    intent(inout)  :: g(idgs, jdgs, nt)
        integer(ip), intent(in)     :: idgs
        integer(ip), intent(in)     :: jdgs
        real(wp),    intent(in)     :: a(mdab, ndab, nt)
        real(wp),    intent(in)     :: b(mdab, ndab, nt)
        integer(ip), intent(in)     :: mdab
        integer(ip), intent(in)     :: ndab
        real(wp),    intent(in)     :: p(imid,*)
        integer(ip), intent(in)     :: idg
        integer(ip), intent(in)     :: jdg
        real(wp),    intent(inout)  :: ge(idg,jdg,*)
        real(wp),    intent(inout)  :: go(idg,jdg,*)
        real(wp),    intent(inout)  :: work(*)
        real(wp),    intent(inout)  :: whrfft(*)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer(ip)    :: i, j, imid, imm1, k, ls
        integer(ip)    :: m, mb, mdo, mmax, mn, modl, nlon
        integer(ip)    :: mp1, mp2,  ndo, nlp1, np1
        type(HFFTpack) :: hfft
        !----------------------------------------------------------------------

        ls = idg
        nlon = jdg
        mmax = min(nlat,nlon/2+1)

        if (2*mmax-1 > nlon) then
            mdo = mmax-1
        else
            mdo = mmax
        end if

        nlp1 = nlat+1
        modl = mod(nlat,2)

        if (modl /= 0) then
            imm1 = imid-1
        else
            imm1 = imid
        end if

        ge(1:ls,1:nlon,1:nt) = 0.0_wp

        block_construct: block

            if (isym /= 1) then

                do k=1,nt
                    do np1=1,nlat,2
                        do i=1,imid
                            ge(i,1,k)=ge(i,1,k)+a(1,np1,k)*p(i,np1)
                        end do
                    end do
                end do

                select case (mod(nlat,2))
                    case (0)
                        ndo = nlat-1
                    case default
                        ndo = nlat
                end select

                do mp1=2,mdo
                    m = mp1-1
                    mb = m*(nlat-1)-(m*(m-1))/2
                    do np1=mp1,ndo,2
                        mn = mb+np1
                        do k=1,nt
                            do i=1,imid
                                ge(i,2*mp1-2,k) = ge(i,2*mp1-2,k)+a(mp1,np1,k)*p(i,mn)
                                ge(i,2*mp1-1,k) = ge(i,2*mp1-1,k)+b(mp1,np1,k)*p(i,mn)
                            end do
                        end do
                    end do
                end do

                if (.not.(mdo == mmax .or. mmax > ndo)) then
                    mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
                    do np1=mmax,ndo,2
                        mn = mb+np1
                        do k=1,nt
                            do i=1,imid
                                ge(i,2*mmax-2,k) = ge(i,2*mmax-2,k)+a(mmax,np1,k)*p(i,mn)
                            end do
                        end do
                    end do
                end if

                if(isym == 2) exit block_construct

            end if

            do k=1,nt
                do np1=2,nlat,2
                    do i=1,imm1
                        go(i,1,k)=go(i,1,k)+a(1,np1,k)*p(i,np1)
                    end do
                end do
            end do

            if(mod(nlat,2) /= 0) then
                ndo = nlat-1
            else
                ndo = nlat
            end if

            do mp1=2,mdo
                mp2 = mp1+1
                m = mp1-1
                mb = m*(nlat-1)-(m*(m-1))/2
                do np1=mp2,ndo,2
                    mn = mb+np1
                    do k=1,nt
                        do i=1,imm1
                            go(i,2*mp1-2,k) = go(i,2*mp1-2,k)+a(mp1,np1,k)*p(i,mn)
                            go(i,2*mp1-1,k) = go(i,2*mp1-1,k)+b(mp1,np1,k)*p(i,mn)
                        end do
                    end do
                end do
            end do

            mp2 = mmax+1

            if (.not.(mdo == mmax .or. mp2 > ndo)) then
                mb = mdo*(nlat-1)-(mdo*(mdo-1))/2
                do np1=mp2,ndo,2
                    mn = mb+np1
                    do k=1,nt
                        do i=1,imm1
                            go(i,2*mmax-2,k) = go(i,2*mmax-2,k)+a(mmax,np1,k)*p(i,mn)
                        end do
                    end do
                end do
            end if

        end block block_construct

        do k=1,nt
            if(mod(nlon,2) == 0) ge(1:ls,nlon,k) = 2.0_wp*ge(1:ls,nlon,k)
            call hfft%backward(ls,nlon,ge(1,1,k),ls,whrfft,work)
        end do

        select case (isym)
            case (0)
                do k=1,nt
                    do j=1,nlon
                        do i=1,imm1
                            g(i,j,k) = 0.5_wp*(ge(i,j,k)+go(i,j,k))
                            g(nlp1-i,j,k) = 0.5_wp*(ge(i,j,k)-go(i,j,k))
                        end do

                        if (modl /= 0) g(imid,j,k) = 0.5_wp*ge(imid,j,k)

                    end do
                end do
            case default
                do k=1,nt
                    g(1:imid,1:nlon,k) = 0.5_wp*ge(1:imid,1:nlon,k)
                end do
        end select

    end subroutine shses1


end submodule scalar_synthesis_shses
