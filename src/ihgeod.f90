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
module module_ihgeod

    use spherepack_precision, only: &
        wp, & ! working precision
        ip, & ! integer precision
        PI

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: ihgeod
    public :: sph2cart
    public :: cart2sph

contains

    subroutine ihgeod(m, idp, jdp, x, y, z)
        !
        !     m         is the number of points on the edge of a
        !               single geodesic triangle
        !
        !     x, y, z     the coordinates of the geodesic points on
        !               the sphere are x(i, j, k), y(i, j, k), z(i, j, k)
        !               where i=1, ..., m+m-1; j=1, ..., m; and k=1, ..., 5.
        !               the indices are defined on the unfolded
        !               icosahedron as follows for the case m=3
        !
        !                north pole
        !
        !                 (5, 1)          0      l
        !        i     (4, 1) (5, 2)              a    (repeated for
        !           (3, 1) (4, 2) (5, 3)  theta1   t    k=2, 3, 4, 5 in
        !        (2, 1) (3, 2) (4, 3)              i        -->
        !     (1, 1) (2, 2) (3, 3)        theta2   t    the longitudinal
        !        (1, 2) (2, 3)                    u    direction)
        !           (1, 3)                pi     d
        !      j                                e
        !         south pole
        !
        !                total number of points is 10*(m-1)**2+2
        !                total number of triangles is 20*(m-1)**2
        !                total number of edges is 30*(m-1)**2
        !
        !-------------------------------------------------------------
        ! Dictionary: calling arguments
        !-------------------------------------------------------------
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: idp
        integer (ip), intent (in)     :: jdp
        real (wp),    intent (in out) :: x(idp, jdp, 5)
        real (wp),    intent (in out) :: y(idp, jdp, 5)
        real (wp),    intent (in out) :: z(idp, jdp, 5)
        !-------------------------------------------------------------
        ! Dictionary: local variables
        !-------------------------------------------------------------
        real (wp) :: beta
        real (wp) :: dphi
        real (wp) :: dxi
        real (wp) :: dxj
        real (wp) :: dyi
        real (wp) :: dyj
        real (wp) :: dzi
        real (wp) :: dzj
        real (wp) :: hdphi
        integer (ip) :: i
        integer (ip) :: j
        integer (ip) :: k
        real (wp) :: phi
        real (wp) :: rad
        real (wp) :: tdphi
        real (wp) :: theta
        real (wp) :: theta1
        real (wp) :: theta2
        real (wp) :: x1
        real (wp) :: x2
        real (wp) :: x3
        real (wp) :: x4
        real (wp) :: x5
        real (wp) :: x6
        real (wp) :: xs
        real (wp) :: y1
        real (wp) :: y2
        real (wp) :: y3
        real (wp) :: y4
        real (wp) :: y5
        real (wp) :: y6
        real (wp) :: ys
        real (wp) :: z1
        real (wp) :: z2
        real (wp) :: z3
        real (wp) :: z4
        real (wp) :: z5
        real (wp) :: z6
        real (wp) :: zs
        !-------------------------------------------------------------

        dphi = 0.4_wp*pi
        beta = cos(dphi)
        theta1 = acos(beta/(1.0_wp-beta))
        theta2 = pi-theta1
        hdphi = dphi/2
        tdphi = 3.0_wp*hdphi

        do k=1, 5
            phi = (k-1)*dphi
            call sph2cart(1.0_wp, theta2, phi, x1, y1, z1)
            call sph2cart(1.0_wp, pi, phi+hdphi, x2, y2, z2)
            call sph2cart(1.0_wp, theta2, phi+dphi, x3, y3, z3)
            dxi = (x2-x1)/(m-1)
            dyi = (y2-y1)/(m-1)
            dzi = (z2-z1)/(m-1)
            dxj = (x3-x2)/(m-1)
            dyj = (y3-y2)/(m-1)
            dzj = (z3-z2)/(m-1)

            do i=1, m
                xs = x1 + (i-1)*dxi
                ys = y1 + (i-1)*dyi
                zs = z1 + (i-1)*dzi
                do j=1, i
                    x(j, i, k) = xs + (j-1)*dxj
                    y(j, i, k) = ys + (j-1)*dyj
                    z(j, i, k) = zs + (j-1)*dzj
                end do
            end do

            call sph2cart(1.0_wp, theta1, phi+hdphi, x4, y4, z4)

            dxi = (x3-x4)/(m-1)
            dyi = (y3-y4)/(m-1)
            dzi = (z3-z4)/(m-1)
            dxj = (x4-x1)/(m-1)
            dyj = (y4-y1)/(m-1)
            dzj = (z4-z1)/(m-1)
            do j=1, m
                xs = x1 + (j-1)*dxj
                ys = y1 + (j-1)*dyj
                zs = z1 + (j-1)*dzj
                do i=1, j
                    x(j, i, k) = xs + (i-1)*dxi
                    y(j, i, k) = ys + (i-1)*dyi
                    z(j, i, k) = zs + (i-1)*dzi
                end do
            end do

            call sph2cart(1.0_wp, theta1, phi+tdphi, x5, y5, z5)

            dxj = (x5-x3)/(m-1)
            dyj = (y5-y3)/(m-1)
            dzj = (z5-z3)/(m-1)
            do i=1, m
                xs = x4 + (i-1)*dxi
                ys = y4 + (i-1)*dyi
                zs = z4 + (i-1)*dzi
                do j=1, i
                    x(j+m-1, i, k) = xs + (j-1)*dxj
                    y(j+m-1, i, k) = ys + (j-1)*dyj
                    z(j+m-1, i, k) = zs + (j-1)*dzj
                end do
            end do

            call sph2cart(1.0_wp, 0.0_wp, phi+dphi, x6, y6, z6)

            dxi = (x5-x6)/(m-1)
            dyi = (y5-y6)/(m-1)
            dzi = (z5-z6)/(m-1)
            dxj = (x6-x4)/(m-1)
            dyj = (y6-y4)/(m-1)
            dzj = (z6-z4)/(m-1)
            do j=1, m
                xs = x4 + (j-1)*dxj
                ys = y4 + (j-1)*dyj
                zs = z4 + (j-1)*dzj
                do i=1, j
                    x(j+m-1, i, k) = xs + (i-1)*dxi
                    y(j+m-1, i, k) = ys + (i-1)*dyi
                    z(j+m-1, i, k) = zs + (i-1)*dzi
                end do
            end do
        end do

        do k=1, 5
            do j=1, m+m-1
                do i=1, m
                    call cart2sph(x(j, i, k), y(j, i, k), z(j, i, k), rad, theta, phi)
                    call sph2cart(1.0_wp, theta, phi, x(j, i, k), y(j, i, k), z(j, i, k))
                end do
            end do
        end do

    end subroutine ihgeod


    pure subroutine cart2sph(x, y, z, r, theta, phi)
        !-------------------------------------------------------------
        ! Dictionary: calling arguments
        !-------------------------------------------------------------
        real (wp), intent (in)   :: x
        real (wp), intent (in)   :: y
        real (wp), intent (in)   :: z
        real (wp), intent (out)  :: r
        real (wp), intent (out)  :: theta
        real (wp), intent (out)  :: phi
        !-------------------------------------------------------------
        real (wp) :: radial
        !-------------------------------------------------------------

        radial = hypot(x,y)**2

        if (radial == 0.0_wp) then
            phi = 0.0_wp
            theta = 0.0_wp

            if (z < 0.0_wp) theta = PI

        else
            r = sqrt(radial+z**2)
            radial = sqrt(radial)
            phi = atan2(y, x)
            theta = atan2(radial, z)
        end if

    end subroutine cart2sph


    pure subroutine sph2cart(r, theta, phi, x, y, z)
        !-------------------------------------------------------------
        ! Dictionary: calling arguments
        !-------------------------------------------------------------
        real (wp), intent (in)  :: r
        real (wp), intent (in)  :: theta
        real (wp), intent (in)  :: phi
        real (wp), intent (out) :: x
        real (wp), intent (out) :: y
        real (wp), intent (out) :: z
        !-------------------------------------------------------------

        x = r*sin(theta)*cos(phi)
        y = r*sin(theta)*sin(phi)
        z = r*cos(theta)

    end subroutine sph2cart

end module module_ihgeod
