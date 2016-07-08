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
program tvtsgs

    use, intrinsic :: iso_fortran_env, only: &
        stdout => OUTPUT_UNIT

    use spherepack_library, only: &
        pi, shagsi, vhagsi, vhsgsi, vtsgsi, &
        gaqd, vhsgs, shags, shsgs, vtsgs

    implicit none
    real :: a
    real :: b
    real :: bi
    real :: br
    real :: c
    real :: ci
    real :: cr
    real :: da
    real :: db
    real :: dmax1
    real :: dmax2
    real :: dphi
    real :: dtr
    real :: hpi
    integer :: i
    integer :: idp
    integer :: ierror
    integer :: iq
    integer :: isym
    integer :: j
    integer :: jdp
    integer :: jq
    integer :: k
    integer :: ldwork
    integer :: lwork
    integer :: lwshi
    integer :: lwvha
    integer :: lwvhs
    integer :: lwvts
    integer :: mdab
    integer :: mp1
    integer :: ndab
    integer :: nlat
    integer :: nlon
    integer :: np1
    integer :: nt
    real :: phi
    real :: s
    real :: theta
    real :: tmp
    real :: u
    real :: v
    real :: vt
    real :: w
    real :: work
    real :: wshi
    real :: wt
    real :: wvha
    real :: wvhs
    real :: wvts
    real :: x
    real :: y
    real :: z
    real :: dummy_variable
    !
    !     this program checks vtsgs for computing the colatitudinal
    !     derivatives of the velocity ....
    !
    parameter (nlat=17,nlon=32)
    parameter (idp=nlat,jdp=nlon,mdab=nlat,ndab=nlat)
    !
    dimension u(idp,jdp),v(idp,jdp),w(idp,jdp), &
        x(idp,jdp),y(idp,jdp),z(idp,jdp), &
        vt(idp,jdp),wt(idp,jdp), &
        a(mdab,ndab,3),b(mdab,ndab,3), &
        da(mdab,ndab,3,3),db(mdab,ndab,3,3), &
        tmp(mdab,ndab), br(mdab,ndab),bi(mdab,ndab), &
        cr(mdab,ndab),ci(mdab,ndab), &
        c(idp,jdp,3,3),s(idp,jdp,3,3), &
        dthet(idp),dwts(idp)
    dimension wshi(608411),wvha(1098771),wvhs(1098771), &
        wvts(1098771),work(165765)
    real dthet,dwts,dwork(25741)

    write( *, '(/a/)') '     testvtsgs *** TEST RUN *** '

!    lwshi = 608411
!    lwvha = 1098771
!    lwvhs = 1098771
!    lwvts = 1098771
!    lwork = 165765
!    ldwork = 25741
!    !
!    hpi = pi/2.
!    dtr = pi/180
!    isym = 0
!    nt = 1
!    !
!    !     initialize spherepack routines
!    !
!    call shagsi(nlat,nlon,wshi,lwshi,work,lwork,dwork,ldwork,ierror)
!    if (ierror /= 0) write( stdout, 55) ierror
!55  format('testvtsgs:  error' i4 ' in shigs')
!    call vhagsi(nlat,nlon,wvha,lwvha,dwork,ldwork,ierror)
!    if (ierror /= 0) write( stdout, 57) ierror
!57  format('testvtsgs:  error' i4 ' in vhagsi')
!    call vhsgsi(nlat,nlon,wvhs,lwvhs,dwork,ldwork,ierror)
!    if (ierror /= 0) write( stdout, 58) ierror
!58  format(' testvtsgs: error' i4 ' in vhsgsi')
!    call vtsgsi(nlat,nlon,wvts,lwvts,work,lwork,dwork,ldwork,ierror)
!    if (ierror /= 0) write( stdout, 59) ierror
!59  format(' testvtsgs: error' i4 ' in vtsgsi')
!    !
!    !     compute gauss points and weights
!    !
!    call gaqd(nlat,dthet,dwts,dummy_variable,lwork,ierror)
!    !
!    !     zero vector harmonic coefficients
!    !
!    br = 0.0
!    bi = 0.0
!    cr = 0.0
!    ci = 0.0
!    !
!    !     initialize arrays with random numbers
!    !     old style non-portable commented out
!    !
!    !     call srand(227)
!    !     do np1=1,nlat-1
!    !     do mp1=1,np1
!    !     br(mp1,np1) = rand()
!    !     bi(mp1,np1) = rand()
!    !     cr(mp1,np1) = rand()
!    !     ci(mp1,np1) = rand()
!    !     end do
!    !     end do
!    !
!    !     set vector harmonic coefficients
!    !     (new style using standard Fortran90
!    !     intrinsics
!    !
!    call random_seed()
!
!    call random_number(tmp)
!    do np1=1,nlat-1
!        do mp1=1,np1
!            br(mp1,np1) = tmp(mp1,np1)
!        end do
!    end do
!    !
!    call random_number(tmp)
!    do np1=1,nlat-1
!        do mp1=1,np1
!            bi(mp1,np1) = tmp(mp1,np1)
!        end do
!    end do
!    !
!    call random_number(tmp)
!    do np1=1,nlat-1
!        do mp1=1,np1
!            cr(mp1,np1) = tmp(mp1,np1)
!        end do
!    end do
!    !
!    call random_number(tmp)
!    do np1=1,nlat-1
!        do mp1=1,np1
!            ci(mp1,np1) = tmp(mp1,np1)
!        end do
!    end do
!    !
!    call vhsgs(nlat,nlon,0,1,v,w,idp,jdp,br,bi, &
!        cr,ci,mdab,ndab,wvhs,lwvhs,work,lwork,ierror)
!    if (ierror /= 0) write( stdout, 79) ierror
!79  format(' testvtsgs: error' i4 ' in vhsgs at point 1')
!    !
!    !==> Initialize
!    u = 0.0
!    !
!    !     convert to cartesian coordinates
!    !
!    dphi = 2.*pi/nlon
!    do j=1,nlon
!        phi = (j-1)*dphi
!        do i=1,nlat
!            theta = dthet(i)
!            call sph2cart(theta,phi,u(i,j),v(i,j),w(i,j),x(i,j),y(i,j),z(i,j))
!        end do
!    end do
!    !
!    !     compute harmonic components of x component
!    !
!    call shags(nlat,nlon,0,1,x,idp,jdp,a(1,1,1),b(1,1,1), &
!        mdab,ndab,wshi,lwshi,work,lwork,ierror)
!    if (ierror /= 0) write( stdout, 16) ierror
!16  format(' testvtsgs: error' i4 ' in shags at point 2')
!    !
!    !     write harmonic coefficients for x
!    !
!    !      write( stdout, 20)
!    ! 20   format(//'  harmonic coefficients for x'//)
!    !      do mp1=1,nlat
!    !      do np1=mp1,nlat
!    !      write( stdout, 5) mp1,np1,a(mp1,np1,1),b(mp1,np1,1)
!    ! 5    format(2i5,1p2e15.6)
!    !      end do
!    !      end do
!    !
!    !     compute harmonic components of y component
!    !
!    call shags(nlat,nlon,0,1,y,idp,jdp,a(1,1,2),b(1,1,2), &
!        mdab,ndab,wshi,lwshi,work,lwork,ierror)
!    if (ierror /= 0) write( stdout, 17) ierror
!17  format(' testvtsgs: error' i4 ' in shags at point 3')
!    !
!    !     write harmonic coefficients for y
!    !
!    !      write( stdout, 21)
!    ! 21   format(//'  harmonic coefficients for y'//)
!    !      do mp1=1,nlat
!    !      do np1=mp1,nlat
!    !      write( stdout, 5) mp1,np1,a(mp1,np1,2),b(mp1,np1,2)
!    !      end do
!    !      end do
!    !
!    !     compute harmonic components of z component
!    !
!    call shags(nlat,nlon,0,1,z,idp,jdp,a(1,1,3),b(1,1,3), &
!        mdab,ndab,wshi,lwshi,work,lwork,ierror)
!    if (ierror /= 0) write( stdout, 18) ierror
!18  format(' testvtsgs: error' i4 ' in shags at point 4')
!    !
!    !     write harmonic coefficients for z
!    !
!    !      write( stdout, 22)
!    ! 22   format(//'  harmonic coefficients for z'//)
!    !      do mp1=1,nlat
!    !      do np1=mp1,nlat
!    !      write( stdout, 5) mp1,np1,a(mp1,np1,3),b(mp1,np1,3)
!    !      end do
!    !      end do
!    !
!    !     compute partials of x, y, and z wrt z
!    !
!    call dbdz(nlat,nlat,a(1,1,1),b(1,1,1),da(1,1,1,3),db(1,1,1,3))
!    call dbdz(nlat,nlat,a(1,1,2),b(1,1,2),da(1,1,2,3),db(1,1,2,3))
!    call dbdz(nlat,nlat,a(1,1,3),b(1,1,3),da(1,1,3,3),db(1,1,3,3))
!    !
!    !     compute partials of x, y, and z wrt y
!    !
!    call dbdy(nlat,nlat,a(1,1,1),b(1,1,1),da(1,1,1,2),db(1,1,1,2))
!    call dbdy(nlat,nlat,a(1,1,2),b(1,1,2),da(1,1,2,2),db(1,1,2,2))
!    call dbdy(nlat,nlat,a(1,1,3),b(1,1,3),da(1,1,3,2),db(1,1,3,2))
!    !
!    !     compute partials of x, y, and z wrt x
!    !
!    call dbdx(nlat,nlat,a(1,1,1),b(1,1,1),da(1,1,1,1),db(1,1,1,1))
!    call dbdx(nlat,nlat,a(1,1,2),b(1,1,2),da(1,1,2,1),db(1,1,2,1))
!    call dbdx(nlat,nlat,a(1,1,3),b(1,1,3),da(1,1,3,1),db(1,1,3,1))
!    !
!    !     transform cartesian jacobian to physical space
!    !
!    call shsgs(nlat,nlon,0,9,c,idp,jdp,da,db,idp,idp, &
!        wshi,lwshi,work,lwork,ierror)
!    if (ierror /= 0) write( stdout, 19) ierror
!19  format(' testvtsgs: error' i4 ' in shsgs at point 5')
!    !
!    !     convert to jacobian to spherical coordinates
!    !        (multiply cartesian jacobian by q)
!    !
!    do k=1,3
!        do j=1,nlon
!            phi = (j-1)*dphi
!            do i=1,nlat
!                theta = dthet(i)
!                call cart2sph(theta,phi,c(i,j,1,k),c(i,j,2,k),c(i,j,3,k), &
!                    s(i,j,1,k),s(i,j,2,k),s(i,j,3,k))
!            end do
!        end do
!    end do
!    !
!    !     form s = (q sq**T)**T
!    !
!    do k=1,3
!        do j=1,nlon
!            phi = (j-1)*dphi
!            do i=1,nlat
!                theta = dthet(i)
!                call cart2sph(theta,phi,s(i,j,k,1),s(i,j,k,2),s(i,j,k,3), &
!                    c(i,j,k,1),c(i,j,k,2),c(i,j,k,3))
!            end do
!        end do
!    end do
!
!    do j=1,nlon
!        do i=1,nlat
!            do iq=1,3
!                do jq=1,3
!                    s(i,j,iq,jq) = c(i,j,iq,jq)
!                end do
!            end do
!        end do
!    end do
!    !
!    !     check derivative program
!    !
!    call vtsgs(nlat,nlon,0,1,vt,wt,idp,jdp,br,bi,cr,ci, &
!        mdab,ndab,wvts,lwvts,work,lwork,ierror)
!    if (ierror /= 0) write( stdout, 4) ierror
!4   format(' testvtsgs: error' i4 ' in vtsgs during initialization')
!
!    dmax1 = maxval(abs(s(:,:,2,2)-vt))
!    dmax2 = maxval(abs(s(:,:,3,2)-wt))
!
!    write( stdout, 2) dmax1,dmax2
!2   format(' testvtsgs: error in vt '1pe15.6' error in wt '1pe15.6)

contains



    subroutine cart2sph(theta,phi,x,y,z,u,v,w)
        implicit none
        real :: cosp
        real :: cost
        real :: phi
        real :: sinp
        real :: sint
        real :: temp1
        real :: temp2
        real :: theta
        real :: u
        real :: v
        real :: w
        real :: x
        real :: y
        real :: z
        !
        !     this program computes the components of a vector
        !     field in spherical coordinates u, v, and w, from
        !     its components x, y, and z in cartesian coordinates
        !
        sint = sin(theta)
        cost = cos(theta)
        sinp = sin(phi)
        cosp = cos(phi)
        temp1 = cosp*x+sinp*y
        temp2 = cosp*y-sinp*x
        u = sint*temp1+cost*z
        v = cost*temp1-sint*z
        w = temp2

    end subroutine cart2sph



    subroutine sph2cart(theta,phi,u,v,w,x,y,z)
        implicit none
        real :: cosp
        real :: cost
        real :: phi
        real :: sinp
        real :: sint
        real :: temp1
        real :: temp2
        real :: theta
        real :: u
        real :: v
        real :: w
        real :: x
        real :: y
        real :: z
        !
        !     this program computes the components of a vector
        !     field in cartesian coordinates x, y, and z, from
        !     its components u, v, and w in spherical coordinates
        !
        sint = sin(theta)
        cost = cos(theta)
        sinp = sin(phi)
        cosp = cos(phi)
        temp1 = sint*u+cost*v
        temp2 = cost*u-sint*v
        x = cosp*temp1-sinp*w
        y = sinp*temp1+cosp*w
        z = temp2

    end subroutine sph2cart



    subroutine dbdx(l,mdim,a,b,dxa,dxb)
        implicit none
        real :: a
        real :: a1
        real :: a2
        real :: b
        real :: cn
        real :: dxa
        real :: dxb
        real :: fm
        real :: fn
        integer :: l
        integer :: lm1
        integer :: mdim
        integer :: mp1
        integer :: n
        integer :: np1
        !
        !     subroutine to compute the coefficients in the spherical
        !     harmonic representation of the derivative with respect to x
        !     of a scalar function. i.e. given the coefficients a and b
        !     in the spectral representation of a function, then dbdx
        !     computes the coefficients dxa and dxb in the spectral
        !     representation of the derivative of the function
        !     with respect to x.
        !
        !     the arrays a and dxa can be the same as well as the arrays
        !     b and dxb, i.e. the arrays a and b can be overwritten by
        !     dxa and dxb respectively.
        !
        !     dimension a(mdim,1),b(mdim,1),dxa(mdim,1),dxb(mdim,1)
        dimension a(mdim,*),b(mdim,*),dxa(mdim,*),dxb(mdim,*)
        dxa(1,1) = sqrt(6.)*a(2,2)
        dxb(1,1) = 0.
        lm1 = l-1
        do np1=2,lm1
            n = np1-1
            fn = real(n)
            cn = (fn+fn+3.)/(fn+fn+1.)
            dxa(1,np1) = sqrt(cn*(fn+2.)*(fn + 1.0))*a(2,np1+1)
            dxb(1,np1) = 0.
            do mp1=2,np1
                fm = real(mp1-1)
                a1 = .5*sqrt(cn*(fn+fm+2.)*(fn+fm+1.))
                a2 = .5*sqrt(cn*(fn-fm+2.)*(fn-fm+1.))
                dxa(mp1,np1) = a1*a(mp1+1,np1+1)-a2*a(mp1-1,np1+1)
                dxb(mp1,np1) = a1*b(mp1+1,np1+1)-a2*b(mp1-1,np1+1)
            end do
        end do
        do mp1=1,l
            dxa(mp1,l) = 0.
            dxb(mp1,l) = 0.
        end do

    end subroutine dbdx



    subroutine dbdy(l,mdim,a,b,dya,dyb)
        implicit none
        real :: a
        real :: a1
        real :: a2
        real :: b
        real :: cn
        real :: dya
        real :: dyb
        real :: fm
        real :: fn
        integer :: l
        integer :: lm1
        integer :: mdim
        integer :: mp1
        integer :: n
        integer :: np1
        !
        !     subroutine to compute the coefficients in the spherical
        !     harmonic representation of the derivative with respect to y
        !     of a scalar function. i.e. given the coefficients a and b
        !     in the spectral representation of a function, then dbdy
        !     computes the coefficients dya and dyb in the spectral
        !     representation of the derivative of the function
        !     with respect to y.
        !
        !     the arrays a and dya can be the same as well as the arrays
        !     b and dyb, i.e. the arrays a and b can be overwritten by
        !     dya and dyb respectively.
        !
        !     dimension a(mdim,1),b(mdim,1),dya(mdim,1),dyb(mdim,1)
        dimension a(mdim,*),b(mdim,*),dya(mdim,*),dyb(mdim,*)
        dya(1,1) = -sqrt(6.)*b(2,2)
        dyb(1,1) = 0.
        lm1 = l-1
        do np1=2,lm1
            n = np1-1
            fn = real(n)
            cn = (fn+fn+3.)/(fn+fn+1.)
            dya(1,np1) = -sqrt(cn*(fn+2.)*(fn + 1.0))*b(2,np1+1)
            dyb(1,np1) = 0.
            do mp1=2,np1
                fm = real(mp1-1)
                a1 = .5*sqrt(cn*(fn+fm+2.)*(fn+fm+1.))
                a2 = .5*sqrt(cn*(fn-fm+2.)*(fn-fm+1.))
                dya(mp1,np1) = -a1*b(mp1+1,np1+1)-a2*b(mp1-1,np1+1)
                dyb(mp1,np1) =  a1*a(mp1+1,np1+1)+a2*a(mp1-1,np1+1)
            end do
        end do

        do mp1=1,l
            dya(mp1,l) = 0.
            dyb(mp1,l) = 0.
        end do

    end subroutine dbdy



    subroutine dbdz(l,mdim,a,b,dza,dzb)
        implicit none
        real :: a
        real :: a1
        real :: b
        real :: cn
        real :: dza
        real :: dzb
        real :: fm
        real :: fn
        integer :: l
        integer :: lm1
        integer :: mdim
        integer :: mp1
        integer :: n
        integer :: np1
        !
        !     subroutine to compute the coefficients in the spherical
        !     harmonic representation of the derivative with respect to z
        !     of a scalar function. i.e. given the coefficients a and b
        !     in the spectral representation of a function, then dbdz
        !     computes the coefficients dza and dzb in the spectral
        !     representation of the derivative of the function
        !     with respect to z.
        !
        !     the arrays a and dza can be the same as well as the arrays
        !     b and dzb, i.e. the arrays a and b can be overwritten by
        !     dza and dzb respectively.
        !
        dimension a(mdim,*),b(mdim,*),dza(mdim,*),dzb(mdim,*)
        lm1 = l-1
        do np1=1,lm1
            n = np1-1
            fn = real(n)
            cn = (fn+fn+3.)/(fn+fn+1.)
            do mp1=1,np1
                fm = real(mp1-1)
                a1 = sqrt(cn*(fn-fm+1.)*(fn+fm+1.))
                dza(mp1,np1) = a1*a(mp1,np1+1)
                dzb(mp1,np1) = a1*b(mp1,np1+1)
            end do
        end do

        do mp1=1,l
            dza(mp1,l) = 0.
            dzb(mp1,l) = 0.
        end do

    end subroutine dbdz

end program tvtsgs
