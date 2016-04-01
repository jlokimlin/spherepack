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
!     a program for testing all scalar analysis and synthesis subroutines
!
program tsha
    !
    !     set dimensions with parameter statements
    !
    parameter(nnlat= 15,nnlon= 18, nnt = 3)
    !     parameter(nnlat=14,nnlon=20,nnt=3)
    parameter (lleng= 5*nnlat*nnlat*nnlon,llsav= 5*nnlat*nnlat*nnlon)
    parameter (lldwork = nnlat*(nnlat+4))
    real dwork(lldwork)
    dimension work(lleng),wsave(llsav)
    dimension a(nnlat,nnlat,nnt),b(nnlat,nnlat,nnt),s(nnlat,nnlon,nnt)
    dimension thetag(nnlat),dtheta(nnlat),dwts(nnlat)
    real dtheta, dwts
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
    isym = 0
    !
    !     set equally spaced colatitude and longitude increments
    !
    pi = acos( -1.0 )
    dphi = (pi+pi)/nlon
    dlat = pi/(nlat-1)
    !
    !     compute nlat gaussian points in thetag
    !
    ldwork = lldwork
    call gaqd(nlat,dtheta,dwts,dwork,ldwork,ier)
    do  i=1,nlat
        thetag(i) = dtheta(i)
    end do
    call name("gaqd")
    call iout(ier," ier")
    call vecout(thetag,"thtg",nlat)
    !
    !     test all analysis and synthesis subroutines
    !
    do icase=1,4
        !
        !     icase=1 test shaec,shsec
        !     icase=2 test shaes,shses
        !     icase=3 test shagc,shsgc
        !     icase=4 test shags,shsgs
        !
        call name("****")
        call name("****")
        call iout(icase,"icas")
        !
        !
        !     set scalar field as (x*y*z)**k) restricted to the sphere
        !
        do k=1,nt
            do j=1,nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1,nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta=thetag(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    xyzk = (sint*(sint*cost*sinp*cosp))**k
                    !           s(i,j,k) = exp(xyzk)
                    s(i,j,k) = xyzk
                end do
            end do
        !     call iout(k,"   k")
        !     call aout(s(1,1,k),"   s",nlat,nlon)
        end do

        do l=1,lsave
            wsave(l) = 0.0
        end do
        if (icase==1) then

            call name("**ec")
            call shaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)

            call name("shai")
            call iout(ierror,"ierr")

            call shaec(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave, &
                lsave,work,lwork,ierror)

            call name("sha ")
            call iout(ierror,"ierr")

            call shseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)

            call name("shsi")
            call iout(ierror,"ierr")

            call shsec(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave, &
                lsave,work,lwork,ierror)

            call name("shs ")
            call iout(ierror,"ierr")

        else if (icase==2) then

            call name("**es")
            call shaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)

            call name("shai")
            call iout(ierror,"ierr")

            call shaes(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave, &
                lsave,work,lwork,ierror)

            call name("sha ")
            call iout(ierror,"ierr")

            call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)

            call name("shsi")
            call iout(ierror,"ierr")

            call shses(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave, &
                lsave,work,lwork,ierror)

            call name("shs ")
            call iout(ierror,"ierr")

        else if (icase==3) then

            call name("**gc")

            call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)

            call name("shai")
            call iout(ierror,"ierr")

            call shagc(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave, &
                lsave,work,lwork,ierror)

            call name("sha ")
            call iout(ierror,"ierr")

            call shsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)

            call name("shsi")
            call iout(ierror,"ierr")

            call shsgc(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave, &
                lsave,work,lwork,ierror)

            call name("shs ")
            call iout(ierror,"ierr")

        else if (icase==4) then

            call name("**gs")

            call shagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)

            call name("shai")
            call iout(ierror,"ierr")

            call shags(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave, &
                lsave,work,lwork,ierror)

            call name("sha ")
            call iout(ierror,"ierr")

            call shsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call shsgs(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave, &
                lsave,work,lwork,ierror)

            call name("shs ")
            call iout(ierror,"ierr")
        end if
        !
        !     compute "error" in s
        !
        err2 = 0.0
        do k=1,nt
            do j=1,nlon
                phi = (j-1)*dphi
                sinp = sin(phi)
                cosp = cos(phi)
                do i=1,nlat
                    theta = (i-1)*dlat
                    if (icase>2) theta = thetag(i)
                    cost = cos(theta)
                    sint = sin(theta)
                    xyzk = (sint*(sint*cost*sinp*cosp))**k
                    !           err2 = err2+ (exp(xyzk)-s(i,j,k))**2
                    err2 = err2 + (xyzk-s(i,j,k))**2
                end do
            end do
        !     call iout(k,"   k")
        !     call aout(s(1,1,k),"   s",nlat,nlon)
        end do
        err2 = sqrt(err2/(nt*nlat*nlon))
        call vout(err2,"err2")
    end do
end program tsha
subroutine iout(ivar,nam)
    real nam
    write(6,10) nam , ivar
10  format(1h a4, 3h = ,i8)
    return
end subroutine iout
!
subroutine vout(var,nam)
    real nam
    write(6,10) nam , var
10  format(1h a4,3h = ,e12.5)
    return
end subroutine vout
!
subroutine name(nam)
    real nam
    write(6,100) nam
100 format(1h a8)
    return
end subroutine name
!
subroutine vecout(vec,nam,len)
    dimension vec(len)
    real nam
    write(6,109) nam, (vec(l),l=1,len)
109 format(1h a4,/(1h 8e11.4))
    return
end subroutine vecout
