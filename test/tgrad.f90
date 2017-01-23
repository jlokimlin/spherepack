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
!
!     12/96
!
!     a program for testing all gradient and inverse gradient routines
!
!     (1) first a scalar field is set in st by restricting a poly in x,y,z
!         to the sphere surface
!
!     (2) a scalar analysis is used to compute the coefs a,b of st
!
!     (3) a,b are input to the various gradient routines to compute a vector field
!         (v,w)
!
!     (4) the vector field (v,w) is compared with the gradient of st obtained
!         analytically
!
!     (5) the inverse gradient of (v,w) is computed and compared with (1)
!
program tgrad

    use spherepack_library

    implicit none

    real :: a
    real :: b
    real :: bi
    real :: br
    real :: ci
    real :: cosp
    real :: cost
    real :: cr
    real :: dlat
    real :: dphi
    real :: dsfdp
    real :: dsfdt
    real :: dxdp
    real :: dxdt
    real :: dydp
    real :: dydt
    real :: dzdp
    real :: dzdt
    real :: err2s
    real :: err2v
    real :: err2w
    integer :: i
    integer :: icase
    integer :: ier
    integer :: ierror
    integer :: isym
    integer :: ityp
    integer :: j
    integer :: k
    integer :: ldwork
    integer :: lldwork
    integer :: lleng
    integer :: llsav
    integer :: lsave
    integer :: lwork
    integer :: mdab
    integer :: mdb
    integer :: mmdab
    integer :: mmdb
    integer :: nlat
    integer :: nlon
    integer :: nnlat
    integer :: nnlon
    integer :: nnt
    integer :: nt
    real :: phi

    real :: se
    real :: sf
    real :: sinp
    real :: sint
    real :: theta
    real :: gaussian_latitudes
    real :: v
    real :: ve
    real :: w
    real :: we
    real :: work
    real :: wsave
    real :: x
    real :: y
    real :: z
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat= 33,nnlon= 18, nnt = 4)
    parameter  (mmdb=(nnlon+1)/2, mmdab=(nnlon+2)/2)
    parameter (lleng= 5*nnlat*nnlat*nnlon,llsav= 5*nnlat*nnlat*nnlon)
    parameter (lldwork = 4*nnlat*nnlat)
    real dwork(lldwork)
    dimension work(lleng),wsave(llsav)
    dimension br(mmdb,nnlat,nnt),bi(mmdb,nnlat,nnt)
    dimension cr(mmdb,nnlat,nnt),ci(mmdb,nnlat,nnt)
    dimension a(mmdab,nnlat,nnt),b(mmdab,nnlat,nnt)
    dimension sf(nnlat,nnlon,nnt)
    dimension gaussian_latitudes(nnlat),dtheta(nnlat),dwts(nnlat)
    dimension v(nnlat,nnlon,nnt),w(nnlat,nnlon,nnt)
    real dtheta, dwts
    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    lwork = lleng
    lsave = llsav
    mdb = mmdb
    mdab = mmdab
    nt = nnt
    call iout(nlat,"nlat")
    call iout(nlon,"nlon")
    call iout(nt,"  nt")
    isym = 0
    ityp = 0
    !
    !     set equally spaced colatitude and longitude increments
    !
    dphi = TWO_PI/nlon
    dlat = PI/(nlat-1)
    !
    !     compute nlat gaussian points in thetag
    !
    ldwork = lldwork
    call compute_gaussian_latitudes_and_weights(nlat,dtheta,dwts,ier)
    do  i=1,nlat
        gaussian_latitudes(i) = dtheta(i)
    end do
    call name("gaqd")
    call iout(ier," ier")
    call vecout(gaussian_latitudes,"thtg",nlat)
    !
    !     test all analysis and synthesis subroutines
    !
    do icase=1,4
        !
        !     icase=1 test gradec,igradec
        !     icase=2 test grades,igrades
        !     icase=3 test gradgc,igradgc
        !     icase=4 test gradgs,igradgs
        !
        call name("****")
        call name("****")
        call iout(icase,"icas")
        !
        !
        !     set scalar field as polys in x,y,z and set vector field (v,w) using
        !     v = dsfdt, w = 1/sint*dsfdp
        !
        do k=1,nt
            do j=1,nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1,nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=gaussian_latitudes(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    dxdp = -sint*sinp
                    dydt = cost*sinp
                    dydp = sint*cosp
                    dzdt = -sint
                    dzdp = 0.0
                    select case (k)
                        case (1)
                            sf(i,j,k) = x*y
                            dsfdt = x*dydt+y*dxdt
                            dsfdp = x*dydp+y*dxdp
                            v(i,j,k) = dsfdt
                            w(i,j,k) = (cosp*dydp+sinp*dxdp)
                        case (2)
                            sf(i,j,k) = x*z
                            dsfdp = x*dzdp+z*dxdp
                            dsfdt = x*dzdt+z*dxdt
                            v(i,j,k) = dsfdt
                            w(i,j,k) = cosp*dzdp-z*sinp
                        case (3)
                            sf(i,j,k) = y*z
                            dsfdt = y*dzdt + z*dydt
                            dsfdp = y*dzdp+ z*dydp
                            v(i,j,k) = dsfdt
                            w(i,j,k) = sinp*dzdp + z*cosp
                        case (4)
                            sf(i,j,k) = x*y*z
                            dsfdt = x*y*dzdt + x*z*dydt + y*z*dxdt
                            dsfdp = x*y*dzdp + x*z*dydp + y*z*dxdp
                            v(i,j,k) = dsfdt
                            w(i,j,k) = cosp*y*dzdp+cosp*z*dydp+sinp*z*dxdp
                    end select
                end do
            end do
        end do

        !     do kk=1,nt
        !     call iout(kk,"**kk")
        !     call aout(sf(1,1,kk),"  sf",nlat,nlon)
        !     call aout(v(1,1,kk),"   v",nlat,nlon)
        !     call aout(w(1,1,kk),"   w",nlat,nlon)
        !     end do

        select case(icase)
            case(1)

                call name("**ec")
                !
                !     analyze scalar field st
                !
                call shaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
                call name("shai")
                call iout(ierror,"ierr")

                call shaec(nlat,nlon,isym,nt,sf,nlat,nlon,a,b,mdab,nlat,wsave, &
                    lsave,work,lwork,ierror)
                call name("sha ")
                call iout(ierror,"ierr")

                !     do kk=1,nt
                !     call iout(kk,"**kk")
                !     call aout(a(1,1,kk),"   a",nlat,nlat)
                !     call aout(b(1,1,kk),"   b",nlat,nlat)
                !     end do

                !
                !     compute gradient of st in v,w
                !

                call vhseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
                call name("vhci")
                call iout(ierror,"ierr")

                call gradec(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,mdab,nlat,wsave, &
                    lsave,work,lwork,ierror)
                call name("grad")
                call iout(ierror,"ierr")

            case(2)

                call name("**es")
                !
                !     analyze scalar field st
                !
                call shaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
                call name("shai")
                call iout(ierror,"ierr")

                call shaes(nlat,nlon,isym,nt,sf,nlat,nlon,a,b,mdab,nlat,wsave, &
                    lsave,work,lwork,ierror)
                call name("sha ")
                call iout(ierror,"ierr")

                !     do kk=1,nt
                !     call iout(kk,"**kk")
                !     call aout(a(1,1,kk),"   a",nlat,nlat)
                !     call aout(b(1,1,kk),"   b",nlat,nlat)
                !     end do

                !
                !     compute gradient of st in v,w
                !

                call vhsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
                call name("vhsi")
                call iout(ierror,"ierr")

                call grades(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,mdab,nlat,wsave, &
                    lsave,work,lwork,ierror)
                call name("grad")
                call iout(ierror,"ierr")

            case(3)


                call name("**gc")
                call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
                call name("shai")
                call iout(ierror,"ierr")

                call shagc(nlat,nlon,isym,nt,sf,nlat,nlon,a,b,mdab,nlat,wsave, &
                    lsave,work,lwork,ierror)
                call name("sha ")
                call iout(ierror,"ierr")

                !     do kk=1,nt
                !     call iout(kk,"**kk")
                !     call aout(a(1,1,kk),"   a",nlat,nlat)
                !     call aout(b(1,1,kk),"   b",nlat,nlat)
                !     end do

                !
                !     compute gradient of st in v,w
                !

                call vhsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
                call name("vhgc")
                call iout(ierror,"ierr")

                call gradgc(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,mdab,nlat,wsave, &
                    lsave,work,lwork,ierror)
                call name("grad")
                call iout(ierror,"ierr")


            case(4)

                call name("**gs")

                call shagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
                call name("shai")
                call iout(ierror,"ierr")

                call shags(nlat,nlon,isym,nt,sf,nlat,nlon,a,b,mdab,nlat,wsave, &
                    lsave,work,lwork,ierror)
                call name("sha ")
                call iout(ierror,"ierr")

                !     do kk=1,nt
                !     call iout(kk,"**kk")
                !     call aout(a(1,1,kk),"   a",nlat,nlat)
                !     call aout(b(1,1,kk),"   b",nlat,nlat)
                !     end do

                !
                !     compute gradient of st in v,w
                !

                call vhsgsi(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
                call name("vhgs")
                call iout(ierror,"ierr")

                call gradgs(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,mdab,nlat,wsave, &
                    lsave,work,lwork,ierror)
                call name("grad")
                call iout(ierror,"ierr")
        end select


        !     do kk=1,nt
        !     call iout(kk,"**kk")
        !     call aout(v(1,1,kk),"   v",nlat,nlon)
        !     call aout(w(1,1,kk),"   w",nlat,nlon)
        !     end do

        !
        !     compute "error" in v,w
        !
        err2v = 0.0
        err2w = 0.0
        do k=1,nt
            do j=1,nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1,nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=gaussian_latitudes(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    dxdp = -sint*sinp
                    dydt = cost*sinp
                    dydp = sint*cosp
                    dzdt = -sint
                    dzdp = 0.0
                    select case (k)
                        case (1)
                            dsfdt = x*dydt+y*dxdt
                            dsfdp = x*dydp+y*dxdp
                            ve = dsfdt
                            we = (cosp*dydp+sinp*dxdp)
                        case (2)
                            dsfdp = x*dzdp+z*dxdp
                            dsfdt = x*dzdt+z*dxdt
                            ve = dsfdt
                            we = cosp*dzdp-z*sinp
                        case (3)
                            dsfdt = y*dzdt + z*dydt
                            dsfdp = y*dzdp+ z*dydp
                            ve = dsfdt
                            we = sinp*dzdp + z*cosp
                        case (4)
                            dsfdt = x*y*dzdt + x*z*dydt + y*z*dxdt
                            dsfdp = x*y*dzdp + x*z*dydp + y*z*dxdp
                            ve = dsfdt
                            we = cosp*y*dzdp+cosp*z*dydp+sinp*z*dxdp
                    end select
                    err2v = err2v + (v(i,j,k)-ve)**2
                    err2w = err2w + (w(i,j,k)-we)**2
                end do
            end do
        end do
        !
        !     set and print least squares error in v,w
        !
        err2v = sqrt(err2v/(nt*nlat*nlon))
        err2w = sqrt(err2w/(nt*nlat*nlon))
        call vout(err2v,"errv")
        call vout(err2w,"errw")
        !
        !     now recompute sf by inverting (v,w) using igrad(ec,es,gc,gs)
        !
        select case (icase)
            case(1)

                call vhaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
                call name("vhai")
                call iout(ierror,"ierr")

                call vhaec(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci, &
                    mdb,nlat,wsave,lsave,work,lwork,ierror)
                call name("vha ")
                call iout(ierror,"ierr")

                call shseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
                call name("shec")
                call iout(ierror,"ierr")

                call igradec(nlat,nlon,isym,nt,sf,nlat,nlon,br,bi, &
                    mdb,nlat,wsave,lsave,work,lwork,ierror)
                call name("igra")
                call iout(ierror,"ierr")

            case(2)
                !
                !     analyze vector field (v,w)
                !
                call vhaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
                call name("vhai")
                call iout(ierror,"ierr")

                call vhaes(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci, &
                    mdb,nlat,wsave,lsave,work,lwork,ierror)
                call name("vha ")
                call iout(ierror,"ierr")

                !     do kk=1,nt
                !     call iout(kk,"**kk")
                !     call aout(br(1,1,kk),"  br",nlat,nlat)
                !     call aout(bi(1,1,kk),"  bi",nlat,nlat)
                !     call aout(cr(1,1,kk),"  cr",nlat,nlat)
                !     call aout(ci(1,1,kk),"  ci",nlat,nlat)
                !     end do

                call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
                call name("shes")
                call iout(ierror,"ierr")

                call igrades(nlat,nlon,isym,nt,sf,nlat,nlon,br,bi, &
                    mdb,nlat,wsave,lsave,work,lwork,ierror)
                call name("igra")
                call iout(ierror,"ierr")

            case(3)

                call vhagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
                call name("vhai")
                call iout(ierror,"ierr")

                call vhagc(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci, &
                    mdb,nlat,wsave,lsave,work,lwork,ierror)
                call name("vha ")
                call iout(ierror,"ierr")

                call shsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
                call name("shgc")
                call iout(ierror,"ierr")

                call igradgc(nlat,nlon,isym,nt,sf,nlat,nlon,br,bi, &
                    mdb,nlat,wsave,lsave,work,lwork,ierror)
                call name("igra")
                call iout(ierror,"ierr")

            case(4)
                call vhagsi(nlat, nlon, wsave, lsave, dwork, ldwork, ierror)
                call name("vhai")
                call iout(ierror,"ierr")

                call vhags(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci, &
                    mdb,nlat,wsave,lsave,work,lwork,ierror)
                call name("vha ")
                call iout(ierror,"ierr")

                call shsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
                call name("shgs")
                call iout(ierror,"ierr")

                call igradgs(nlat,nlon,isym,nt,sf,nlat,nlon,br,bi, &
                    mdb,nlat,wsave,lsave,work,lwork,ierror)
                call name("igra")
                call iout(ierror,"ierr")

        end select

        !     do kk=1,nt
        !     call iout(kk,"**kk")
        !     call aout(sf(1,1,kk),"  sf",nlat,nlon)
        !     end do

        !
        !     compare this sf with original
        !
        err2s = 0.0
        do k=1,nt
            do j=1,nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1,nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=gaussian_latitudes(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    select case (k)
                        case (1)
                            se = x*y
                        case (2)
                            se = x*z
                        case (3)
                            se = y*z
                        case (4)
                            se = x*y*z
                    end select
                    err2s = err2s + (sf(i,j,k)-se)**2
                end do
            end do
        end do
        err2s = sqrt(err2s/(nlat*nlon*nt))
        call vout(err2s,"errs")
    !
    !     end of icase loop
    !
    end do
end program tgrad
!
subroutine iout(ivar,nam)
    implicit none
    integer :: ivar
    character(len=*), intent(in) :: nam
    write(6,10) nam , ivar
10  format(1h a4, 3h = ,i8)
    return
end subroutine iout
!
subroutine vout(var,nam)
    implicit none
    real :: var
    character(len=*), intent(in) :: nam
    write(6,10) nam , var
10  format(1h a4,3h = ,e12.5)
    return
end subroutine vout
!
subroutine name(nam)
    implicit none
    character(len=*), intent(in) :: nam
    write(6,100) nam
100 format(1h a8)
    return
end subroutine name

