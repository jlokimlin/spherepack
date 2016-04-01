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
!     3/97
!
!     a program for testing slap,islap (ec,es,gc,gs)
!
!     (1) set a scalar field s as poly in x,y,z restricted to sphere
!
!     (2) compute scalar laplacian in array sclp using slap(ec,es,gc,gs)
!
!     (3) compare (2) with analytic scalar laplacian in sclpe
!
!     (4) compute the inverse  of (2) and compare with (1)
!
program tslap
    parameter(nnlat=15,nnlon= 22, nnt = 3)
    parameter (mmdab = (nnlon+2)/2)

    parameter (lleng=15*nnlat*nnlat*nnlon,llsav=5*nnlat*nnlat*nnlon)
    dimension work(lleng),wsave(llsav)
    parameter (lldwork = nnlat*(nnlat+4))
    real dwork(lldwork)
    dimension a(mmdab,nnlat,nnt),b(mmdab,nnlat,nnt),s(nnlat,nnlon,nnt)
    dimension sclp(nnlat,nnlon,nnt)
    dimension sclpe(nnlat,nnlon,nnt)
    dimension ptrb(nnt),xlm(nnt)
    dimension thetag(nnlat),dtheta(nnlat),dwts(nnlat)
    real dtheta, dwts
    !
    !     set dimension variables
    !
    nlat = nnlat
    nlon = nnlon
    mdab = mmdab
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
    call gaqd(nlat,dtheta,dwts,work,lwork,ier)
    do  i=1,nlat
        thetag(i) = dtheta(i)
    end do
    call name("gaqd")
    call iout(ier," ier")
    call vecout(thetag,"thtg",nlat)
    !
    !     set helmholtz constant zero for laplacian inversion
    !
    do k=1,nt
        xlm(k) = 0.0
    end do
    !
    !     test all analysis and synthesis subroutines
    !
    do icase=1,4
        call name("****")
        call name("****")
        call iout(icase,"icas")
        !
        !
        !     set scalar field as poly in x,y,z restricted to the sphere
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
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    dxdt = cost*cosp
                    d2xdt2 = -x
                    dxdp = -cost*sinp
                    d2xdp2 = -x
                    dydt = cost*sinp
                    d2ydt2 = -y
                    dydp = cost*cosp
                    d2ydp2 = -y
                    dzdt = -sint
                    d2zdt2 = -z
                    dzdp = 0.
                    d2zdp2 = 0.
                    if (k==1) then
                        s(i,j,k) = x+y
                        sclpe(i,j,k) = -2.*(x+y)
                    else if (k==2) then
                        s(i,j,k) = x+z
                        sclpe(i,j,k) = -2.*(x+z)
                    else if (k==3) then
                        s(i,j,k) = y+z
                        sclpe(i,j,k) = -2.*(y+z)
                    end if
                end do
            end do
        end do
        !     do k=1,nt
        !     call iout(k,"   k")
        !     call aout(s(1,1,k),"   s",nlat,nlon)
        !     call aout(sclpe(1,1,k),"sclp",nlat,nlon)
        !     end do

        !     call aout(s,"   s",nlat,nlon)


        if (icase==1) then

            call name("**ec")

            call shaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("shai")
            call iout(ierror,"ierr")

            call shaec(nlat,nlon,isym,nt,s,nlat,nlon,a,b,mdab,nlat,wsave, &
                lsave,work,lwork,ierror)
            call name("sha ")
            call iout(ierror,"ierr")

            call shseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call slapec(nlat,nlon,isym,nt,sclp,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ierror)
            call name("slap")
            call iout(ierror,"ierr")

        else if (icase==2) then

            call name("**es")

            call shaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("shai")
            call iout(ierror,"ierr")

            call shaes(nlat,nlon,isym,nt,s,nlat,nlon,a,b,mdab,nlat,wsave, &
                lsave,work,lwork,ierror)
            call name("sha ")
            call iout(ierror,"ierr")

            call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call slapes(nlat,nlon,isym,nt,sclp,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ierror)
            call name("slap")
            call iout(ierror,"ierr")

        else if (icase==3) then

            call name("**gc")

            call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("shai")
            call iout(ierror,"ierr")

            call shagc(nlat,nlon,isym,nt,s,nlat,nlon,a,b,mdab,nlat,wsave, &
                lsave,work,lwork,ierror)
            call name("sha ")
            call iout(ierror,"ierr")

            call shsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call slapgc(nlat,nlon,isym,nt,sclp,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ierror)
            call name("slap")
            call iout(ierror,"ierr")

        else if (icase==4) then

            call name("**gs")

            call shagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("shai")
            call iout(ierror,"ierr")

            call shags(nlat,nlon,isym,nt,s,nlat,nlon,a,b,mdab,nlat,wsave, &
                lsave,work,lwork,ierror)
            call name("sha ")
            call iout(ierror,"ierr")

            call shsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call slapgs(nlat,nlon,isym,nt,sclp,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ierror)
            call name("slap")
            call iout(ierror,"ierr")
        end if
        !
        !     compute "error" in sclp
        !
        err2 = 0.0
        do k=1,nt
            do j=1,nlon
                do i=1,nlat
                    err2 = err2 + (sclpe(i,j,k)-sclp(i,j,k))**2
                end do
            end do
        !     call iout(k,"   k")
        !     call aout(sclp(1,1,k),"sclp",nlat,nlon)
        end do
        err2 = sqrt(err2/(nt*nlat*nlon))
        call vout(err2,"err2")
        !
        !     invert sclp
        !
        if (icase==1) then

            call shaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call shaec(nlat,nlon,isym,nt,sclp,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ierror)

            call shseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call islapec(nlat,nlon,isym,nt,xlm,s,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ptrb,ierror)
            call name("isla")
            call iout(ierror,"ierr")
            call vecout(ptrb,"ptrb",nt)

        else if (icase==2) then

            call shaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call shaes(nlat,nlon,isym,nt,sclp,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ierror)

            call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call islapes(nlat,nlon,isym,nt,xlm,s,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ptrb,ierror)
            call name("isla")
            call iout(ierror,"ierr")
            call vecout(ptrb,"ptrb",nt)

        else if (icase==3) then

            call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call shagc(nlat,nlon,isym,nt,sclp,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ierror)

            call shsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call islapgc(nlat,nlon,isym,nt,xlm,s,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ptrb,ierror)
            call name("isla")
            call iout(ierror,"ierr")
            call vecout(ptrb,"ptrb",nt)

        else if (icase==4) then

            call name("**gs")

            call shagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call shags(nlat,nlon,isym,nt,sclp,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ierror)

            call shsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
            call name("shsi")
            call iout(ierror,"ierr")

            call islapgs(nlat,nlon,isym,nt,xlm,s,nlat,nlon,a,b,mdab,nlat, &
                wsave,lsave,work,lwork,ptrb,ierror)
            call name("isla")
            call iout(ierror,"ierr")
            call vecout(ptrb,"ptrb",nt)

        end if

        !     call aout(s,"   s",nlat,nlon)


        !
        !     compare s with original
        !
        err2s = 0.0
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
                    x = sint*cosp
                    y = sint*sinp
                    z = cost
                    if (k==1) then
                        se = x+y
                    else if (k==2) then
                        se = x+z
                    else if (k==3) then
                        se = y+z
                    end if
                    err2s = err2s+(s(i,j,k) - se)**2
                end do
            end do
        end do
        err2s = sqrt(err2s/(nlat*nlon*nt))
        call vout(err2s,"errs")
    !
    !     end of icase loop
    !
    end do
end program tslap
!
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
