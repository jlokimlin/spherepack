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
subroutine vsurf(xeye, yeye, zeye, ntri, x1, y1, z1, x2, y2, z2, &
    x3, y3, z3, itype, work, iwork)
    implicit none
    integer :: itype
    integer :: ntri
    real :: work
    real :: x1
    real :: x2
    real :: x3
    real :: xeye
    real :: y1
    real :: y2
    real :: y3
    real :: yeye
    real :: z1
    real :: z2
    real :: z3
    real :: zeye
    !
    !    subroutine vsurf is like subroutine hidel except the triangles
    !    are categorized. vsurf is also like solid except triangles rather
    !    than lines are covered.
    !
    !     written by paul n. swarztrauber, national center for atmospheric
    !     research, p.o. box 3000, boulder, colorado, 80307
    !
    !    this program plots visible lines for the surface defined
    !    by the input 3-d triangles with corners at (x1, y1, z1), (x2, y2, z2)
    !    and (x3, y3, z3). the sides of these these triangles may or
    !    may not be plotted depending on itype. if itype is 1 then the
    !    side between points (x1, y1, z1) and (x2, y2, z2) is plotted if it
    !    is visible. if itype is 2 then the side between (x2, y2, z2)
    !    and (x3, y3, z3) is plotted. if itype is 3 then the visible portion
    !    of the side between (x3, y3, z3) and (x1, y1, z1) is plotted.
    !    any combination is possible by specifying itype to be one
    !    of the following values: 0, 1, 2, 3, 12, 13, 23, 123.
    !
    !    the length of real    array  work must be at least 19*ntri
    !
    !    the length of integer array iwork must be at least 19*ntri
    !
    !
    !    the vertices of the triangles are renumbered by vsurf so that
    !    their projections are orientated counterclockwise. the user need
    !    only be aware that the vertices may be renumbered by vsurf.
    !
    dimension x1(ntri), y1(ntri), z1(ntri), x2(ntri), y2(ntri), z2(ntri), &
        x3(ntri), y3(ntri), z3(ntri), itype(ntri), work(19*ntri)
    integer iwork(19*ntri)
    !
    call vsurf1(xeye, yeye, zeye, ntri, x1, y1, z1, x2, y2, z2, x3, y3, z3, &
        itype, work, work(ntri+1), work(2*ntri+1), work(3*ntri+1), &
        work(4*ntri+1), work(5*ntri+1), work(6*ntri+1), work(7*ntri+1), &
        work(8*ntri+1), work(9*ntri+1), work(10*ntri+1), work(11*ntri+1), &
        work(12*ntri+1), work(13*ntri+1), iwork(14*ntri+1), iwork(6*ntri+1), &
        iwork(15*ntri+1), iwork(17*ntri+1))

end subroutine vsurf
