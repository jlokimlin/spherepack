SPHEREPACK Version 3.2
----------------------
A Package for Modeling Geophysical Processes

12/06/2011

Author:  John  Adams and Paul Swarztrauber

CAVEAT
------
SPHEREPACK 3.2 is written in Fortran 77 and 90 but it may not fully comply 
with either Fortran standard. We have not run the code through a rigorous
standards checker, other than the various compilers we have used in testing
the code. Users whose applications require strict adherence to the standard 
must provide this assurance themselves.

Documentation
-------------
Documentation for this software package is provided in file SPHEREPACK3.2.html
and companion image file SPHERE.gif.  We intend the document for browsing on 
your local computer.

Recent bugfixes
---------------
None. This library is regarded as legacy.

Compiling the Library and Test Programs
---------------------------------------
Our Makefile requires the following: a Unix or Linux or Mac operating system, gmake, 
and a Fortran 90 compiler.  If successful, this Makefile builds a static library 
libspherepack.a and several binary executables.

Examine file make.inc to see if your OS and compiler are represented.  If they 
are not, you should modify file make.inc and the Makefile in each directory 
(main, src and test) so that they are represented.

If you desire double precision floating-point arithmetic, you should modify Makefile 
and make.inc with compiler options for promoting single to double precision.  Your 
compiler's user guide will have information on the required compiler options.
