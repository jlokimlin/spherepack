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
!
!     4/97
!
!     a program for testing all theta derivative subroutines
!     vtses,vtsec,vtsgs,vtsgc
!
!
!     (1) first set a valid vector field (v,w) in terms of x,y,z
!         cartesian coordinates
!
!     (2) analytically compute (vt,wt) from (1)
!
!     (3) compute (vt,wt) using vtses,vtsec,vtsgs,vtsgc and compare with (2)
!
program tvts

    use, intrinsic :: ISO_Fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
        wp, & ! working precision
        PI, &
        compute_gaussian_latitudes_and_weights, vhaeci, vhaec, vtseci, vtsec, vhaesi, &
        vhaes, vtsesi, vtses, vhagci, vhagc, vtsgci, vtsgc, &
        vhagsi, vhags, vtsgsi, vtsgs

    ! Explicit typing only
    implicit none

    real :: bi
    real :: br
    real :: ci
    real :: cosp
    real :: cost
    real :: cr
    real :: dlat
    real :: dphi
    real :: dxdt
    real :: dydt
    real :: dzdt
    real :: emz
    real :: err2v
    real :: err2w
    real :: ex
    real :: ey
    real :: ez
    integer :: i
    integer :: icase
    integer :: ier
    integer :: ierror
    integer :: ityp
    integer :: j
    integer :: k
    integer :: ldwork
    integer :: lldwork
    integer :: lleng
    integer :: llsav
    integer :: lsave
    integer :: lwork
    integer :: nlat
    integer :: nlon
    integer :: nnlat
    integer :: nnlon
    integer :: nnt
    integer :: nt
    real :: phi
    real :: sinp
    real :: sint
    real :: theta
    real :: thetag
    real :: v
    real :: vt
    real :: vtsav
    real :: w
    real :: work
    real :: wsave
    real :: wt
    real :: wtsav
    real :: x
    real :: y
    real :: z
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat= 25,nnlon= 19, nnt = 3)
    parameter (lleng= 5*nnlat*nnlat*nnlon,llsav= 5*nnlat*nnlat*nnlon)
    dimension work(lleng),wsave(llsav)
    parameter (lldwork = 4*nnlat*nnlat )
    real dwork(lldwork)
    dimension br(nnlat,nnlat,nnt),bi(nnlat,nnlat,nnt)
    dimension cr(nnlat,nnlat,nnt),ci(nnlat,nnlat,nnt)
    dimension thetag(nnlat),dtheta(nnlat),dwts(nnlat)
    dimension v(nnlat,nnlon,nnt),w(nnlat,nnlon,nnt)
    dimension vt(nnlat,nnlon,nnt),wt(nnlat,nnlon,nnt)
    dimension vtsav(nnlat,nnlon,nnt),wtsav(nnlat,nnlon,nnt)
    real dtheta, dwts
    real dummy_variable

    write( stdout, '(/a/)') '     tvts *** TEST RUN *** '

    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    lwork = lleng
    lsave = llsav
    nt = nnt
    call iout(nlat,"nlat")
    call iout(nlon,"nlon")
    call iout(nt,"  nt")
    ityp = 0
    !
    !     set equally spaced colatitude and longitude increments
    !
    dphi = (pi+pi)/nlon
    dlat = pi/(nlat-1)
    !
    !     compute nlat gaussian points in thetag
    !
    ldwork = lldwork
    call compute_gaussian_latitudes_and_weights(nlat,dtheta,dwts,dummy_variable,ldwork,ier)
    do  i=1,nlat
        thetag(i) = dtheta(i)
    end do
    call name("compute_gaussian_latitudes_and_weights")
    call iout(ier," ier")
    call vecout(thetag,"thtg",nlat)
    !
    !     test all theta derivative subroutines
    !
    do icase=1,4
        !
        !     icase=1 test vtsec
        !     icase=2 test vtses
        !     icase=3 test vtsgc
        !     icase=4 test vtsgs
        !
        call name("****")
        call name("****")
        call iout(icase,"icas")
        !
        !
        !     set vector field v,w and compute theta derivatives in (vtsav,wtsav)
        !
        do k=1,nt
            do j=1,nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1,nlat
                    theta = real(i - 1, kind=wp)*dlat
                    if (icase ==3 .or. icase==4) theta = thetag(i)
                    cost = cos(theta)
                    sint = sin(theta)
                       !
                       !    set x,y,z and their theta derivatives at colatitude theta and longitude p
                       !
                    x = sint*cosp
                    dxdt = cost*cosp
                    y = sint*sinp
                    dydt = cost*sinp
                    z = cost
                    dzdt = -sint
                    !
                    !     set (v,w) field corresponding to stream function
                    !     S = exp(y)+exp(-z) and velocity potential function
                    !     P = exp(x)+exp(z)
                    !
                    ex = exp(x)
                    ey = exp(y)
                    ez = exp(z)
                    emz = exp(-z)
                    w(i,j,k) =-ex*sinp+emz*sint+ey*cost*sinp
                    v(i,j,k) =-ey*cosp-ez*sint+ex*cost*cosp
                    !
                    !     set theta derivatives differentiating w,v above
                    !
                    wtsav(i,j,k) = -ex*dxdt*sinp+emz*(-dzdt*sint+cost) &
                        +ey*sinp*(dydt*cost-sint)
                    vtsav(i,j,k) = -ey*dydt*cosp-ez*(dzdt*sint+cost) &
                        +ex*cosp*(dxdt*cost-sint)
                end do
            end do
        end do

        !     call a3out(wtsav,"wtsv",nlat,nlon,nt)
        !     call a3out(vtsav,"vtsv",nlat,nlon,nt)



        if (icase==1) then

            call name("**ec")

            call vhaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("vhai")
            call iout(ierror,"ierr")

            call vhaec(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,nlat, &
                nlat,wsave,lsave,work,lwork,ierror)
            call name("vha ")
            call iout(ierror,"ierr")

            !
            !     now compute theta derivatives of v,w
            !
            call vtseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)

            call name("vtsi")
            call iout(ierror,"ierr")

            call vtsec(nlat,nlon,ityp,nt,vt,wt,nlat,nlon,br,bi,cr,ci,nlat, &
                nlat,wsave,lsave,work,lwork,ierror)
            call name("vts ")
            call iout(ierror,"ierr")

        else if (icase==2) then

            call name("**es")

            call vhaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("vhai")
            call iout(ierror,"ierr")

            call vhaes(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,nlat, &
                nlat,wsave,lsave,work,lwork,ierror)
            call name("vha ")
            call iout(ierror,"ierr")

            call vtsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("vtsi")
            call iout(ierror,"ierr")

            call vtses(nlat,nlon,ityp,nt,vt,wt,nlat,nlon,br,bi,cr,ci,nlat, &
                nlat,wsave,lsave,work,lwork,ierror)
            call name("vts ")
            call iout(ierror,"ierr")

        else if (icase==3) then

            call name("**gc")

            call name("vhgi")
            call iout(nlat,"nlat")

            call vhagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("vhai")
            call iout(ierror,"ierr")

            call vhagc(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,nlat, &
                nlat,wsave,lsave,work,lwork,ierror)
            call name("vha ")
            call iout(ierror,"ierr")

            !
            !     now synthesize v,w from br,bi,cr,ci and compare with original
            !
            call vtsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("vtsi")
            call iout(ierror,"ierr")

            call vtsgc(nlat,nlon,ityp,nt,vt,wt,nlat,nlon,br,bi,cr,ci,nlat, &
                nlat,wsave,lsave,work,lwork,ierror)
            call name("vts ")
            call iout(ierror,"ierr")

        else if (icase==4) then

            call name("**gs")
            call vhagsi(nlat,nlon,wsave,lsave,work,lwork,ierror)
            call name("vhai")
            call iout(ierror,"ierr")

            call vhags(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,nlat, &
                nlat,wsave,lsave,work,lwork,ierror)
            call name("vha ")
            call iout(ierror,"ierr")

            call vtsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("vtsi")
            call iout(ierror,"ierr")

            call vtsgs(nlat,nlon,ityp,nt,vt,wt,nlat,nlon,br,bi,cr,ci,nlat, &
                nlat,wsave,lsave,work,lwork,ierror)
            call name("vts ")
            call iout(ierror,"ierr")

        end if

        !
        !     compute "error" in vt,wt
        !
        err2v = norm2(vt- vtsav)
        err2w = norm2(wt - wtsav)
        !
        !     set and print least squares error in v,w
        !
        call vout(err2v,"errv")
        call vout(err2w,"errw")
    !
    !     end of icase loop
    !
    end do


contains


    subroutine iout(ivar,nam)
        implicit none
        integer :: ivar
        character(len=*) nam
        write( stdout, 10) nam , ivar
10      format(1h a4, 3h = ,i8)
        return
    end subroutine iout



    subroutine vout(var,nam)
        implicit none
        real :: var
        character(len=*) nam
        write( stdout, 10) nam , var
10      format(1h a4,3h = ,e12.5)
        return
    end subroutine vout


    subroutine name(nam)
        implicit none
        character(len=*) nam
        write( stdout, 100) nam
100     format(1h a8)
        return
    end subroutine name



    subroutine vecout(vec,nam,vec_size)
        implicit none
        integer :: l
        integer :: vec_size
        real :: vec
        dimension vec(vec_size)
        character(len=*) nam
        write( stdout, 109) nam, (vec(l),l=1,vec_size)
109     format(1h a4,/(1h 8e11.4))
        return
    end subroutine vecout

end program tvts
