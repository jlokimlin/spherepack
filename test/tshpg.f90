!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 1998 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                      SPHEREPACK                               *
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
!    a program testing subroutines shpgi and shpg in file shpg.f
!    also requires files  shagc.f, shsgc.f, gaqd.f 
!    sphcom.f and hrfft.f 
!
!    shpg is the n**2 projection with complement, odd/even
!    factorization and zero truncation on a Gaussian distributed
!    grid. The projection is defined in the JCP paper "Generalized
!    discrete spherical harmonic transforms" by 
!    Paul N. Swarztrauber and William F. Spotz
!    J. Comp. Phys., 159(2000) pp. 213-230.
!
!                     April 2002
!
program tshpg
    use spherepack_library
    implicit none
    real(wp) :: dmax1
    real(wp) :: g
    real(wp) :: ga
    real(wp) :: gb
    real(wp) :: gh
    integer(ip) :: i
    integer(ip) :: idimg
    integer(ip) :: idp
    integer(ip) :: ierror
    integer(ip) :: iprint
    integer(ip) :: isym
    integer(ip) :: iwshp
    integer(ip) :: j
    integer(ip) :: jdimg
    integer(ip) :: kdp
    integer(ip) :: liwshp
    integer(ip) :: lwork
    integer(ip) :: lwrk
    integer(ip) :: lwrk1
    integer(ip) :: lwsha
    integer(ip) :: lwshp
    integer(ip) :: lwshs
    integer(ip) :: mode
    integer(ip) :: mp1
    integer(ip) :: mtr
    integer(ip) :: mtrunc
    integer(ip) :: nlat
    integer(ip) :: nlon
    integer(ip) :: np1
    integer(ip) :: nt
    real(wp) :: sx
    real(wp) :: sy
    real(wp) :: thold
    real(wp) :: toe
    real(wp) :: tusl
    real(wp) :: wrk1
    real(wp) :: wrk2
    real(wp) :: wshagc
    real(wp) :: wshp
    real(wp) :: wshsgc
    parameter (idp=8)
    parameter (kdp=idp+idp-2)
    parameter (lwshp=2*(idp+1)**2+kdp+20, &
        liwshp=4*(idp+1),lwrk=1.25*(idp+1)**2+7*idp+8)
    parameter (lwrk1=idp*kdp)
    parameter(lwork = 4*idp*(idp-1), &
        lwsha=idp*(4*idp+1)+idp+idp+15)
    !     1 lwsha=idp*(idp+1)+3*(idp-2)*(idp-1)/2+kdp+15)

    real work(lwrk)
    dimension sx(idp,kdp),sy(idp,kdp), &
        wshp(lwshp),iwshp(liwshp),wrk1(lwrk1)
    dimension g(idp,kdp,2),ga(idp,idp,2),gb(idp,idp,2), &
        gh(idp,kdp,2), &
        wrk2(lwork),wshagc(lwsha),wshsgc(lwsha)
    real*4 t1(2)
    !
    iprint = 0
    nt = 1
    isym = 0
    mode = 0
    !
    do nlat=6,8
        do mtr=1,2
            nlon = 2*(nlat-1)
            mtrunc = nlat-mtr
            mtrunc = min(mtrunc,nlat-1,nlon/2)
            idimg = idp
            jdimg = kdp
            call shagci(nlat,nlon,wshagc,lwsha,work,lwrk,ierror)
            if(ierror /= 0) write(6,70) ierror
70          format('   ierror0' ,i5)
            !
            lwshs = lwsha
            call shsgci(nlat,nlon,wshsgc,lwshs,work,lwrk,ierror)
            if(ierror /= 0) write(6,71) ierror
71          format('   ierror1' ,i5)
            !
            !     initiate faster filter
            !
            call shpgi(nlat,nlon,isym,mtrunc,wshp,lwshp,iwshp,liwshp, &
                work,lwrk,ierror)
            if(ierror/=0) write(*,429) ierror
429         format(' ierror2 =',i5,' at 429')
            !
            if(iprint /= 0) write (6,5) mode,nlat,nlon
5           format(' mode =' ,i5,'  nlat =',i5,'  nlon =',i5)
            !
            !     initialize with pseudo random field
            !
            do i=1,nlat
                do j=1,nlon
                    sx(i,j) =  cos(real(i*j))
                    g(i,j,1) = sx(i,j)
                end do
            end do
            !
            thold = etime(t1)
            thold = t1(1)
            call shagc(nlat,nlon,mode,nt,g,idimg,jdimg,ga,gb,idimg,idimg, &
                wshagc,lwsha,wrk2,lwork,ierror)
            if(ierror /= 0) write(6,72) ierror
72          format('   ierror2' ,i5)
            !
            if(mtrunc<nlat-1) then
                do np1=mtrunc+2,nlat
                    do mp1=1,np1
                        ga(mp1,np1,1) = 0.
                        gb(mp1,np1,1) = 0.
                    end do
                end do
            end if
            call shsgc(nlat,nlon,mode,nt,gh,idimg,jdimg,ga,gb,idimg,idimg, &
                wshsgc,lwshs,wrk2,lwork,ierror)
            tusl = etime(t1)
            tusl = t1(1)-thold
            if(ierror /= 0) write(6,73) ierror
73          format('   ierror3' ,i5)
            !
            thold = etime(t1)
            toe = t1(1)
            call shpg(nlat,nlon,isym,mtrunc,sx,sy,idp, &
                wshp,lwshp,iwshp,liwshp,wrk1,lwrk1,ierror)
            thold = etime(t1)
            toe = t1(1)-toe
            if(ierror/=0) write(*,428) ierror
428         format(' ierror =',i5,' at 428')
            if(iprint>0)write(*,431)
431         format(/' approx and exact solution'/)
            do j=1,nlon
                if(iprint>0)write(*,437) j,(sy(i,j),gh(i,j,1),i=1,nlat)
437             format(' j=',i5,1p4e15.6/(8x,1p4e15.6))
            end do
            dmax1 = 0.
            do j=1,nlon
                do i=1,nlat
                    dmax1 = max(dmax1,abs(sy(i,j)-gh(i,j,1)))
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
end program tshpg

