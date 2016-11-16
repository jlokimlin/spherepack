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
! ... file trssph.f
!
!     contains documentation and code for subroutine trssph
!
! ... required files
!
!     type_SpherepackAux.f, type_HFFTpack.f, compute_gaussian_latitudes_and_weights.f, shaec.f, shsec.f, shagc.f, shsgc.f
!
!
!     subroutine trssph(intl, igrida, nlona, nlata, da, igridb, nlonb, nlatb, 
!    +db, wsave, lsave, lsvmin, work, lwork, lwkmin, dwork, ldwork, ier)
!
! *** purpose
!
!     subroutine trssph transfers data given in array da on a grid on the
!     full sphere to data in array db on a grid on the full sphere.  the
!     grids on which da is given and db is generated can be specified
!     independently of each other (see description below and the arguments
!     igrida, igridb).  for transferring vector data on the sphere, use
!     subroutine trvsph.

!     notice that scalar and vector quantities are fundamentally different
!     on the sphere.  for example, vectors are discontinuous and multiple
!     valued at the poles.  scalars are continuous and single valued at the
!     poles. erroneous results would be produced if one attempted to transfer
!     vector fields between grids with subroutine trssph applied to each
!     component of the vector.
!
!
! *** underlying grid assumptions and a description
!
!     discussions with the ncar scd data support group and others indicate
!     there is no standard grid for storing observational or model generated
!     data on the sphere.  subroutine trssph was designed to handle most
!     cases likely to be encountered when moving data from one grid format
!     to another.
!
!     the grid on which da is given must be equally spaced in longitude
!     and either equally spaced or gaussian in latitude (or colatitude).
!     longitude, which can be either the first or second dimension of da, 
!     subdivides [0, 2pi) excluding the periodic point 2pi.  (co)latitude, 
!     which can be the second or first dimension of da, has south
!     to north or north to south orientation with increasing subscript
!     value in da (see the argument igrida).
!
!     the grid on which db is generated must be equally spaced in longitude
!     and either equally spaced or gaussian in latitude (or colatitude).
!     longitude, which can be either the first or second dimension of db, 
!     subdivides [0, 2pi) excluding the periodic point 2pi.  (co)latitude, 
!     which can be the second or first dimension of db, has south
!     to north or north to south orientation with increasing subscript
!     value in db (see the argument igridb).
!
!     let nlon be either nlona or nlonb (the number of grid points in
!     longitude.  the longitude grid subdivides [0, 2pi) into nlon spaced
!     points
!
!          (j-1)*2.*pi/nlon  (j=1, ..., nlon).
!
!     it is not necessary to communicate to subroutine trssph whether the
!     underlying grids are in latitude or colatitude.  it is only necessary
!     to communicate whether they run south to north or north to south with
!     increasing subscripts.  a brief discussion of latitude and colatitude
!     follows.  equally spaced latitude grids are assumed to subdivide
!     [-pi/2, pi/2] with the south pole at -pi/2 and north pole at pi/2.
!     equally spaced colatitude grids subdivide [0, pi] with the north pole
!     at 0 and south pole at pi.  equally spaced partitions on the sphere
!     include both poles.  gaussian latitude grids subdivide (-pi/2, pi/2)
!     and gaussian colatitude grids subdivide (0, pi).  gaussian grids do not
!     include the poles.  the gaussian grid points are uniquely determined by
!     the size of the partition.  they can be computed in colatitude in
!     (0, pi) (north to south) in real by the spherepack subroutine
!     compute_gaussian_latitudes_and_weights.  let nlat be nlata or nlatb if either the da or db grid is
!     gaussian.  let
!
!        north pole                             south pole
!        ----------                             ----------
!           0.0    <  cth(1) < ... < cth(nlat)  <   pi
!
!
!     be nlat gaussian colatitude points in the interval (0, pi) and let
!
!        south pole                        north pole
!        ----------                        ----------
!           -pi/2  < th(1) < ... < th(nlat) < pi/2
!
!     be nlat gaussian latitude points in the open interval (-pi/2, pi/2).
!     these are related by
!
!          th(i) = -pi/2 + cth(i)  (i=1, ..., nlat)
!
!     if the da or db grid is equally spaced in (co)latitude then
!
!          ctht(i) = (i-1)*pi/(nlat-1)
!                                               (i=1, ..., nlat)
!          tht(i) = -pi/2 + (i-1)*pi/(nlat-1)
!
!     define the equally spaced (north to south) colatitude and (south to
!     north) latitude grids.
!
!
! *** method (simplified description)
!
!     for simplicity, assume da is a nlat by nlon data tabulation and da(i, j)
!     is the value at latitude theta(i) and longitude phi(j).  then
!     coefficients a(m, n) and b(m, n) can be determined so that da(i, j) is
!     approximated by the sum
!
!         l-1  n
!     (a) sum sum pbar(m, n, theta(i))*(a(m, n)*cos(m*phi(j)+b(m, n)*sin(m*phi(j))
!         n=0 m=0
!
!     here pbar(n, m, theta) are the normalized associated legendre functions
!     and l = min(nlat, (nlon+2)/2).  the determination of a(m, n) and b(m, n)
!     is called spherical harmonic analysis. a sum of this form can then be
!     used to regenerate the data in db on the new grid with the known
!     a(m, n) and b(m, n).  this is referred to spherical harmonic synthesis.
!     analysis and synthesis subroutines from the software package spherepack, 
!     are used for these purposes.
!
!     if da or db is not in mathematical spherical coordinates then array
!     transposition and/or subscript reordering is used prior to harmonic
!     analysis and after harmonic synthesis.
!
! *** advantages
!
!     the use of surface spherical harmonics to transfer spherical grid data
!     has advantages over pointwise grid interpolation schemes on the sphere.
!     it is highly accurate.  if p(x, y, z) is any polynomial of degree n or
!     less in x, y, z cartesian coordinates which is restricted to the surface
!     of the sphere, then p is exactly represented by sums of the form (a)
!     whenever n = mino(nlat, nlon/2) (i.e., transfers with spherical harmonics
!     have n(th) order accuracy.  by way of contrast, bilinear interpolation
!     schemes are exact for polynomials of degree one.  bicubic interpolation
!     is exact only for polynomials of degree three or less.  the method
!     also produces a weighted least squares fit to the data in which waves
!     are resolved uniformly on the full sphere.  high frequencies, induced
!     by closeness of grid points near the poles (due to computational
!     or observational errors) are smoothed.  finally, the method is
!     consistent with methods used to generate data in numerical spectral
!     models based on spherical harmonics.  for more discussion of these and
!     related issues,  see the article: "on the spectral approximation of
!     discrete scalar and vector functions on the sphere, " siam j. numer.
!     anal., vol 16. dec 1979, pp. 934-949, by paul swarztrauber.
!
!
! *** comment
!
!     on a nlon by nlat or nlat by nlon grid (gaussian or equally spaced)
!     spherical harmonic analysis generates and synthesis utilizes
!     min(nlat, (nlon+2)/2)) by nlat coefficients.  consequently, for
!     da and db,  if either
!
!             min(nlatb, (nlonb+2)/2) < min(nlata, (nlona+2)/2)
!
!     or if
!
!             nlatb < nlata
!
!     then all the coefficients generated by an analysis of da cannot be used
!     in the synthesis which generates db.  in this case "information" can be
!     lost in generating db.  more precisely, information will be lost if the
!     analysis of da yields nonzero coefficients which are outside the bounds
!     determined by the db grid.  nevertheless, transference of values with
!     spherical harmonics will yield results consistent with grid resolution
!     and is highly accurate.
!
!
! *** input arguments
!
! ... intl
!
!     an initialization argument which should be zero on an initial call to
!     trssph.  intl should be one if trssph is being recalled and
!
!          igrida, nlona, nlata, igridb, nlonb, nlatb
!
!     have not changed from the previous call.  if any of these arguments
!     have changed, intl=0 must be used to avoid undetectable errors.  calls
!     with intl=1 bypass redundant computation and save time.  it can be used
!     when transferring multiple data sets with the same underlying grids.
!
!
! ... igrida
!
!     an integer vector dimensioned two which identifies the underlying grid
!     on the full sphere for the given data array da as follows:
!
!     igrida(1)
!
!     = -1
!     if the latitude (or colatitude) grid for da is an equally spaced
!     partition of [-pi/2, pi/2] ( or [0, pi] ) including the poles which
!     runs north to south
!
!     = +1
!     if the latitude (or colatitude) grid for da is an equally spaced
!     partition of [-pi/2, pi/2] ( or [0, pi] ) including the poles which
!     runs south to north
!
!     = -2
!     if the latitude (or colatitude) grid for da is a gaussian partition
!     of (-pi/2, pi/2) ( or (0, pi) ) excluding the poles which runs north
!     to south
!
!     = +2
!     if the latitude (or colatitude) grid for da is a gaussian partition
!     of (-pi/2, pi/2) ( or (0, pi) ) excluding the poles which runs south
!     north
!
!     igrida(2)
!
!     = 0 if the underlying grid for da is a nlona by nlata
!
!     = 1 if the underlying grid for da is a nlata by nlona
!
!
! ... nlona
!
!     the number of longitude points on the uniform grid which partitions
!     [0, 2pi) for the given data array da.  nlona is also the first or second
!     dimension of da (see igrida(2)) in the program which calls trssph.
!     nlona determines the grid increment in longitude as 2*pi/nlona. for
!     example nlona = 72 for a five degree grid.  nlona must be greater than
!     or equal to 4.  the efficiency of the computation is improved when
!     nlona is a product of small prime numbers
!
! ... nlata
!
!     the number of points in the latitude (or colatitude) grid
!     for the given data array da.  nlata is also the first or second
!     dimension of da (see igrida(2)) in the program which calls trssph.
!     if nlata is odd then the equator will be located at the (nlata+1)/2
!     gaussian grid point.  if nlata is even then the equator will be
!     located half way between the nlata/2 and nlata/2+1 grid points.
!
! *** note:
!     igrida(1)=-1 or igrida(1)=-2 and igrida(2)=1 corresponds to
!     the "usual" mathematical spherical coordinate system required
!     by most of the drivers in spherepack2.  igrida(1)=1 or igrida(1)=2
!     and igrida(2)=0 corresponds to the "usual" geophysical spherical
!     coordinate system.
!
! ... da
!
!     a two dimensional array that contains the data to be transferred.
!     da must be dimensioned nlona by nlata in the program calling trssph if
!     igrida(2) = 0.  da must be dimensioned nlata by nlona in the program
!     calling trssph if igrida(2) = 1.  if da is not properly dimensioned
!     and if the latitude (colatitude) values do not run south to north or
!     north to south as flagged by igrida(1) (self cannot be checked!) then
!     incorrect results will be produced.
!
! ... igridb
!
!     an integer vector dimensioned two which identifies the underlying grid
!     on the full sphere for the transformed data array db as follows:
!
!     igridb(1)
!
!     = -1
!     if the latitude (or colatitude) grid for db is an equally spaced
!     partition of [-pi/2, pi/2] ( or [0, pi] ) including the poles which
!     north to south
!
!     = +1
!     if the latitude (or colatitude) grid for db is an equally spaced
!     partition of [-pi/2, pi/2] ( or [0, pi] ) including the poles which
!     south to north
!
!     = -2
!     if the latitude (or colatitude) grid for db is a gaussian partition
!     of (-pi/2, pi/2) ( or (0, pi) ) excluding the poles which runs north to
!     south
!
!     = +2
!     if the latitude (or colatitude) grid for db is a gaussian partition
!     of (-pi/2, pi/2) ( or (0, pi) ) excluding the poles which runs south to
!     north
!
!
!     igridb(2)
!
!     = 0 if the underlying grid for db is a nlonb by nlatb
!
!     = 1 if the underlying grid for db is a nlatb by nlonb
!
!
! ... nlonb
!
!     the number of longitude points on the uniform grid which partitions
!     [0, 2pi) for the transformed data array db.  nlonb is also the first or
!     second dimension of db (see igridb(2)) in the program which calls
!     trssph.  nlonb determines the grid increment in longitude as 2*pi/nlonb.
!     for example nlonb = 72 for a five degree grid.  nlonb must be greater
!     than or equal to 4.  the efficiency of the computation is improved when
!     nlonb is a product of small prime numbers
!
! ... nlatb
!
!     the number of points in the latitude (or colatitude) grid
!     for the transformed data array db.  nlatb is also the first or second
!     dimension of db (see igridb(2)) in the program which calls trssph.
!     if nlatb is odd then the equator will be located at the (nlatb+1)/2
!     gaussian grid point.  if nlatb is even then the equator will be
!     located half way between the nlatb/2 and nlatb/2+1 grid points.
!
! ... wsave
!
!     a saved work space array that can be utilized repeatedly by trssph
!     as long as the arguments nlata, nlona, nlatb, nlonb remain unchanged.
!     wsave is set by a intl=0 call to trssph.  wsave must not be altered
!     when trssph is being recalled with intl=1.
!
! ... lsave
!
!     the dimension of the work space wsave as it appears in the program
!     that calls trssph.  the minimum required value of lsave for the
!     current set of input arguments is set in the output argument lsvmin.
!     it can be determined by calling trssph with lsave=0 and printing lsvmin.
!     let
!
!          lwa =  2*nlata*la2+3*((la1-2)*(nlata+nlata-la1-1))/2+nlona+15
!
!     if the grid for da is equally spaced in (co)latitude.  let
!
!          lwa = nlata*(2*la2+3*la1-2)+3*la1*(1-la1)/2+nlona+15
!
!     if the grid for da is gaussian in (co)latitude.
!     let
!
!          lwb = nlatb*(2*lb2+3*lb1-2)+3*lb1*(1-lb1)/2+nlonb+15
!
!     if the grid for db is gaussian in (co)latitude.  let
!
!          lwb = 2*nlatb*lb2+3*((lb1-2)*(nlatb+nlatb-lb1-1))/2+nlonb+15
!
!     if the grid for db is equally spaced in (co)latitude.  then
!     the quantity
!
!          lwa + lwb
!
!     is the minimum required length of wsave.  this value is returned
!     in the output argument lsvmin even if lsave is to small (ierror=10)
!
! ... work
!
!     a real work array that does not have to be preserved
!
! ... lwork
!
!     the dimension of the array work as it appears in the program
!     calling trssph. the minimum required value of lwork for the current
!     set of input arguments is set in the output argument lwkmin.
!     it can be determined by calling trssph with lwork=0 and printing
!     lwkmin.  an estimate for lwork follows.  let nlat, nlon, l1, l2 be
!     defined by
!
!       nlat = max(nlata, nlatb), nlon = nax0(nlona, nlonb), 
!       l1 = min(nlat, (nlon+2)/2), l2 = (nlat+1)/2
!
!     then the quantity
!
!          nlat*(4*l1+nlon+2*nlat+4)+3*((l1-2)*2*(2*nlat-l1-1))/2
!
!     will suffice as a length for the unsaved work space.
!
!  *  both of the formulas above for lsave and lwork may overestimate the
!     required minimum values.  they can be predetermined by calling trssph
!     with lsave=lwork=0 and printout of lsvmin and lwkmin.
!
! ... dwork
!
!     a real work array that does not have to be preserved.
!
! ... ldwork
!
!     The length of dwork in the routine calling trssph.
!     Let
!
!       nlat = max(nlata, nlatb)
!
!     ldwork must be at least nlat*(nlat+4)
!
! *** output arguments
!
!
! ... db
!
!     a two dimensional array that contains the transformed data.  db
!     must be dimensioned nlonb by nlatb in the program calling trssph if
!     igridb(2) = 0 or 1.  db must be dimensioned nlatb by nlonb in the
!     program calling trssph if igridb(2) = 1.  if db is not properly
!     dimensioned and if the latitude (colatitude) values do not run south
!     north or north to south as flagged by igrdb(1) (self cannot be checked!)
!     then incorrect results will be produced.
!
! ... lsvmin
!
!     the minimum length of the saved work space in wsave.
!     lsvmin is computed even if lsave < lsvmin (ier = 10).
!
! ... lwkmin
!
!     the minimum length of the unsaved work space in work.
!     lwkmin is computed even if lwork < lwkmin (ier = 11).
!
! *** error argument
!
! ... ier = 0  if no errors are detected
!
!         = 1  if intl is not 0 or 1
!
!         = 2  if igrida(1) is not -1 or +1 or -2 or +2
!
!         = 3  if igrida(2) is not 0 or 1
!
!         = 4  if nlona is less than 4
!
!         = 5  if nlata is less than 3
!
!         = 6  if igridb(1) is not -1 or +1 or -2 or +2
!
!         = 7  if igridb(2) is not 0 or 1
!
!         = 8  if nlonb is less than 4
!
!         = 9  if nlatb is less than 3

!         =10  if there is insufficient saved work space (lsave < lsvmin)
!
!         =11  if there is insufficient unsaved work space (lwork < lwkmin)
!
!         =12  indicates failure in an eigenvalue routine which computes
!              gaussian weights and points
!
!         =13  if ldwork is too small (insufficient unsaved real
!              work space)
!
! *****************************************************
! *****************************************************
!
!     end of argument description ... code follows
!
! *****************************************************
! *****************************************************
!
module module_trssph

    use spherepack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use module_shaec, only: &
        shaec, shaeci

    use module_shsec, only: &
        shsec, shseci

    use module_shagc, only: &
        shagc, shagci

    use module_shsgc, only: &
        shsgc, shsgci

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    public :: trssph

contains

    subroutine trssph(intl, igrida, nlona, nlata, da, igridb, nlonb, nlatb, &
        db, wsave, lsave, lsvmin, work, lwork, lwkmin, dwork, ldwork, ier)

        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer, intent(in)     :: intl
        integer, intent(in)     :: igrida(2)
        integer, intent(in)     :: nlona
        integer, intent(in)     :: nlata
        integer, intent(in)     :: igridb(2)
        integer, intent(in)     :: nlonb
        integer, intent(in)     :: nlatb
        integer, intent(in)     :: lsave
        integer, intent(inout)  :: lsvmin
        integer, intent(in)     :: lwork
        integer, intent(inout)  :: lwkmin
        integer, intent(in)     :: ldwork
        integer, intent(out)    :: ier
        real,    intent(inout)  :: da(*)
        real,    intent(out)    :: db(*)
        real,    intent(inout)  :: wsave(*)
        real,    intent(inout)  :: work(*)
        real,    intent(inout)  :: dwork(*)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer :: ig, igrda, igrdb, la1, la2, lb1, lb2, lwa, lwb, iaa, iab, iba, ibb
        integer :: lwk3, lwk4, lw, iw, jb, nt, isym, nlat
        !----------------------------------------------------------------------
        !
        !     include a save statement to ensure local variables in trssph, set during
        !     an intl=0 call, are preserved if trssph is recalled with intl=1
        !
        save
        !
        !     check input arguments
        !
        ier = 1
        if (intl*(intl-1)/=0) return
        ier = 2
        ig = igrida(1)
        if ((ig-1)*(ig+1)*(ig-2)*(ig+2)/=0) return
        ier = 3
        ig = igrida(2)
        if (ig*(ig-1)/=0) return
        ier = 4
        if (nlona < 4) return
        ier = 5
        if (nlata <3) return
        ier = 6
        ig = igridb(1)
        if ((ig-1)*(ig+1)*(ig-2)*(ig+2)/=0) return
        ier = 7
        ig = igridb(2)
        if (ig*(ig-1)/=0) return
        ier = 8
        if (nlonb <4) return
        ier = 9
        if (nlatb <3) return
        ier = 0

        igrda = abs(igrida(1))
        igrdb = abs(igridb(1))
        if (intl==0) then
            la1 = min(nlata, (nlona+2)/2)
            la2 = (nlata+1)/2
            lb1 = min(nlatb, (nlonb+2)/2)
            lb2 = (nlatb+1)/2
            !
            !     set saved work space length for analysis
            !
            if (igrda == 1) then
                !
                !     saved space for analysis on  equally spaced grid
                !
                lwa =  2*nlata*la2+3*((la1-2)*(nlata+nlata-la1-1))/2+nlona+15
            else
                !
                !     saved space for analysis on gaussian grid
                !
                lwa = nlata*(2*la2+3*la1-2)+3*la1*(1-la1)/2+nlona+15
            end if
            !
            !     set wsave pointer
            !
            jb = 1+lwa
            !
            !     set pointers for spherical harmonic coefs
            !
            iaa = 1
            iba = iaa+la1*nlata
            iab = iba+la1*nlata
            if (igrdb == 2) then
                !
                !     set saved work space length for gaussian synthesis
                !
                lwb = nlatb*(2*lb2+3*lb1-2)+3*lb1*(1-lb1)/2+nlonb+15
            else
                !
                !     set saved work space length for equally spaced synthesis
                !
                lwb = 2*nlatb*lb2+3*((lb1-2)*(nlatb+nlatb-lb1-1))/2+nlonb+15
            end if
            !
            !     set minimum saved work space length
            !
            lsvmin = lwa + lwb
            !
            !     set remaining harmonic pointer
            !
            ibb = iab+lb1*nlatb
            !
            !     set pointers for remaining work
            !
            iw = ibb+lb1*nlatb
            !
            !     set remaining work space length in lw
            !
            lw = lwork - iw
            lwk3 = nlata*nlona*2
            lwk4 = nlatb*nlonb*2
            !
            !     set minimum unsaved work space required by trssph
            !
            lwkmin = iw + max(lwk3, lwk4)
            !
            !     set error flags if saved or unsaved work spaces are insufficient
            !
            ier = 10
            if (lsave < lsvmin) return
            ier = 11
            if (lwork < lwkmin) return
            ier = 13
            nlat = max(nlata, nlatb)
            if (ldwork < nlat*(nlat+4)) return
            ier = 0
            if (igrda == 1) then
                !
                !     initialize wsave for equally spaced analysis
                !
                call shaeci(nlata, nlona, wsave, lwa, dwork, ldwork, ier)
            else
                !
                !     initialize wsave for gaussian analysis
                !
                call shagci(nlata, nlona, wsave, lwa, dwork, ldwork, ier)
                if (ier/=0) then
                    !
                    !     flag failure in spherepack gaussian software
                    !
                    ier = 12
                    return
                end if
            end if

            if (igrdb == 2) then
                !
                !     initialize wsave for gaussian synthesis
                !
                call shsgci(nlatb, nlonb, wsave(jb), lwb, dwork, ldwork, ier)
                if (ier/=0) then
                    !
                    !     flag failure in spherepack gaussian software
                    !
                    ier = 12
                    return
                end if
            else
                !
                !     initialize wsave for equally spaced synthesis
                !
                call shseci(nlatb, nlonb, wsave(jb), lwb, dwork, ldwork, ier)
            end if
        !
        !     end of initialization (intl=0) call
        !
        end if
        !
        !     transpose and/or reorder (co)latitude if necessary for da
        !     (arrays must have latitude (colatitude) as the first dimension
        !     and run north to south for spherepack software)
        !
        if (igrida(2) == 0) call trsplat(nlona, nlata, da, work)
        if (igrida(1) > 0) call convlat(nlata, nlona, da)

        nt = 1
        isym = 0
        if (igrda == 2) then
            !
            !     do spherical harmonic analysis of "adjusted" da on gaussian grid
            !
            call shagc(nlata, nlona, isym, nt, da, nlata, nlona, work(iaa), &
                work(iba), la1, nlata, wsave, lwa, work(iw), lw, ier)
        else
            !
            !     do spherical harmonic analysis of "adjusted" da on equally spaced grid
            !
            call shaec(nlata, nlona, isym, nt, da, nlata, nlona, work(iaa), &
                work(iba), la1, nlata, wsave, lwa, work(iw), lw, ier)
        end if
        !
        !     transfer da grid coefficients to db grid coefficients
        !     truncating to zero as necessary
        !
        call trab(la1, nlata, work(iaa), work(iba), lb1, nlatb, work(iab), &
            work(ibb))

        if (igrdb == 1) then
            !
            !     do spherical harmonic synthesis on nlatb by nlonb equally spaced grid
            !
            call shsec(nlatb, nlonb, isym, nt, db, nlatb, nlonb, work(iab), &
                work(ibb), lb1, nlatb, wsave(jb), lwb, work(iw), lw, ier)
        else
            !
            !     do spherical harmonic synthesis on nlatb by nlonb gaussian grid
            !
            call shsgc(nlatb, nlonb, isym, nt, db, nlatb, nlonb, work(iab), &
                work(ibb), lb1, nlatb, wsave(jb), lwb, work(iw), lw, ier)
        end if
        !
        !     both da, db are currently latitude by longitude north to south arrays
        !     restore da and set db to agree with flags in igrida and igridb
        !
        if (igrida(1) > 0) then
            call convlat(nlata, nlona, da)
        end if

        if (igridb(1) > 0) then
            call convlat(nlatb, nlonb, db)
        end if

        if (igrida(2) == 0) then
            call trsplat(nlata, nlona, da, work)
        end if

        if (igridb(2) == 0) then
            call trsplat(nlatb, nlonb, db, work)
        end if

    end subroutine trssph



    subroutine trab(ma, na, aa, ba, mb, nb, ab, bb)
        !
        ! Purpose:
        !
        !     set coefficients for b grid from coefficients for a grid
        !

        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer, intent(in)  :: ma
        integer, intent(in)  :: na
        integer, intent(in)  :: mb
        integer, intent(in)  :: nb
        real,    intent(in)  :: aa(ma, na)
        real,    intent(in)  :: ba(ma, na)
        real,    intent(out) :: ab(mb, nb)
        real,    intent(out) :: bb(mb, nb)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer :: i, j, m, n !! Counters
        !----------------------------------------------------------------------

        m = min(ma, mb)
        n = min(na, nb)


        do j=1, n
            do i=1, m
                ab(i, j) = aa(i, j)
                bb(i, j) = ba(i, j)
            end do
        end do
        !
        !     set coefs outside triangle to zero
        !
        do i=m+1, mb
            do j=1, nb
                ab(i, j) = 0.0
                bb(i, j) = 0.0
            end do
        end do

        do j=n+1, nb
            do i=1, mb
                ab(i, j) = 0.0
                bb(i, j) = 0.0
            end do
        end do

    end subroutine trab



    subroutine trsplat(n, m, data, work)
        !
        !     transpose the n by m array data to a m by n array data
        !     work must be at least n*m words long
        !

        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer, intent(in)     :: n
        integer, intent(in)     :: m
        real,    intent(inout)  :: data(*)
        real,    intent(inout)  :: work(*)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer :: i, j, ij, ji !! Counters
        !----------------------------------------------------------------------

        do j=1, m
            do i=1, n
                ij = (j-1)*n+i
                work(ij) = data(ij)
            end do
        end do

        do i=1, n
            do j=1, m
                ji = (i-1)*m+j
                ij = (j-1)*n+i
                data(ji) = work(ij)
            end do
        end do

    end subroutine trsplat



    subroutine convlat(nlat, nlon, data)
        !
        ! Purpose:
        !
        ! Reverse order of latitude (colatitude) grids
        !

        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer, intent(in)     :: nlat
        integer, intent(in)     :: nlon
        real,    intent(inout)  :: data(nlat,nlon)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer :: i, j, half_nlat, ib !! Counters
        real    :: temp
         !----------------------------------------------------------------------

        half_nlat = nlat/2
        do i=1, half_nlat
            ib = nlat-i+1
            do j=1, nlon
                temp = data(i, j)
                data(i, j) = data(ib, j)
                data(ib, j) = temp
            end do
        end do

    end subroutine convlat

end module module_trssph
