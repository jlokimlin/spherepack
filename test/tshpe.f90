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
!    a program testing subroutines shpei and shpe in file shpe.f
!    also requires files  shaec.f shsec.f type_SpherepackAux.f and type_HFFTpack.f
!
!    shpe is the n**2 filter with complement, odd/even
!    factorization and zero truncation on an equally spaced
!    grid. The projection is defined in the JCP paper "Generalized
!    discrete spherical harmonic transforms" by 
!    Paul N. Swarztrauber and William F. Spotz
!    J. Comp. Phys., 159(2000) pp. 213-230.
!
!                     April 2002
!
program tshpe

    use, intrinsic :: ISO_Fortran_env, only: &
        sp => REAL32

    use spherepack_library

    implicit none

    real :: dmax1
    real :: g
    real :: ga
    real :: gb
    real :: gh
    real :: gw
    integer :: i
    integer :: idimg
    integer :: idp
    integer :: ierror
    integer :: iprint
    integer :: isym
    integer :: iwshp
    integer :: j
    integer :: jdimg
    integer :: k
    integer :: kdp
    integer :: liwshp
    integer :: lwork
    integer :: lwrk
    integer :: lwrk1
    integer :: lwsha
    integer :: lwshp
    integer :: lwshs
    integer :: mode
    integer :: mp1
    integer :: mtr
    integer :: mtrunc
    integer :: nlat
    integer :: nlon
    integer :: np1
    integer :: nt
    real :: sx
    real :: sy
    real :: thold
    real :: toe
    real :: tusl
    real :: wrk1
    real :: wrk2
    real :: wshaec
    real :: wshp
    real :: wshsec

    parameter (idp=32)
    parameter (kdp=idp+idp-2)
    parameter (lwshp=2*(idp+1)**2+kdp+20, &
        liwshp=4*(idp+1),lwrk=1.25*(idp+1)**2+7*idp+8)
    parameter (lwrk1=idp*kdp)
    parameter(lwork = 5*idp*(idp-1), &
        lwsha=idp*(idp+1)+3*(idp-2)*(idp-1)/2+kdp+15)
    real :: work(lwrk)
    dimension sx(idp,kdp),sy(idp,kdp), &
        wshp(lwshp),iwshp(liwshp),wrk1(lwrk1)
    dimension g(idp,kdp,2),ga(idp,idp,2),gb(idp,idp,2), &
        gh(idp,kdp,2),gw(idp,kdp,2), &
        wrk2(lwork),wshaec(lwsha),wshsec(lwsha)
    real (sp) :: t1(2)

    write( *, '(a)') ''
    write( *, '(a)') '     tshpe *** TEST RUN *** '
    write( *, '(a)') ''
    !
    iprint = 0
    nt = 1
    !
    do nlat=6,8
        do mtr=1,3
            mtrunc = nlat-mtr
            idimg = idp
            jdimg = kdp
            nlon = 2*(nlat-1)
            mtrunc = min(mtrunc,nlat-1,nlon/2)
            call shaeci(nlat,nlon,wshaec,lwsha,work,lwrk,ierror)
            if (ierror /= 0) write(6,70) ierror
70          format('   ierror0' ,i5)
            !
            lwshs = lwsha
            call shseci(nlat,nlon,wshsec,lwshs,work,lwrk,ierror)
            if (ierror /= 0) write(6,71) ierror
71          format('   ierror1' ,i5)
            !
            mode = 0
            if (iprint /= 0) write (6,5) mode,nlat,nlon
5           format(' mode =' ,i5,'  nlat =',i5,'  nlon =',i5)
            !
            do k=1,nt
                do i=1,nlat
                    do j=1,nlon
                        g(i,j,k) = cos(real(i*j))
                        gh(i,j,k) = g(i,j,k)
                    end do
                end do
            end do
            !
            call cpu_time(thold)
            
            call shaec(nlat,nlon,mode,nt,g,idimg,jdimg,ga,gb,idimg,idimg, &
                wshaec,lwsha,wrk2,lwork,ierror)
            if (ierror /= 0) write(6,72) ierror
72          format('   ierror2' ,i5)
            !
            if (mtrunc<nlat-1) then
                do np1=mtrunc+2,nlat
                    do mp1=1,np1
                        ga(mp1,np1,1) = 0.0
                        gb(mp1,np1,1) = 0.0
                    end do
                end do
            end if
            call shsec(nlat,nlon,mode,nt,gw,idimg,jdimg,ga,gb,idimg,idimg, &
                wshsec,lwshs,wrk2,lwork,ierror)
            call cpu_time(tusl)
            tusl = tusl-thold
            if (ierror /= 0) write(6,73) ierror
73          format('   ierror3' ,i5)
            !
            !     faster filter
            !
            isym = 0
            !      mtrunc = nlat-1
            call shpei(nlat,nlon,isym,mtrunc,wshp,lwshp,iwshp,liwshp, &
                work,lwrk,ierror)
            !
            if (ierror /= 0) write(*,429) ierror
429         format(' ierror =',i5,' at 429')

            do j=1,nlon
                do i=1,nlat
                    sx(i,j) = gh(i,j,1)
                end do
            end do
            !
            !      do i=1,nlat
            !       write(*,427) i,(sx(i,j),j=1,nlon)
            ! 427   format(i5,1p4e15.6/(5x,1p4e15.6))
            !      end do
            call cpu_time(thold)
            toe = t1(1)
            call shpe(nlat,nlon,isym,mtrunc,sx,sy,idp, &
                wshp,lwshp,iwshp,liwshp,wrk1,lwrk1,ierror)
            call cpu_time(thold)
            toe = t1(1)-toe
            if (ierror /= 0) write(*,428) ierror
428         format(' ierror =',i5,' at 428')
            if (iprint>0)write(*,431)
431         format(/' approx and exact solution'/)
            do j=1,nlon
                if (iprint>0)write(*,437) j,(sy(i,j),gw(i,j,1),i=1,nlat)
437             format(' j=',i5,1p4e15.6/(8x,1p4e15.6))
            end do
            dmax1 = 0.
            do j=1,nlon
                do i=1,nlat
                    dmax1 = max(dmax1,abs(sy(i,j)-gw(i,j,1)))
                end do
            end do
            write (6,134) nlat,mtrunc
134         format(/'case nlat =',i5,' and mtrunc =',i5/)
            write(*,135) dmax1,tusl,toe
135         format(/' error =',1pe15.6/ &
                ' tusl  =',1pe15.6/ &
                ' toe   =',1pe15.6)
        end do
    end do


end program tshpe

