SIU with pelagic and discontinuous layers

% ------------------------------------------------------------------------------
%   MILAMIN: MATLAB-based FEM solver for large problems
%
%   Version 1.0
%
%   Copyright (C) 2007, Marcin Dabrowski, Marcin Krotkiewski, Daniel W. Schmid
%                       University of Oslo, Physics of Geological Processes
%
%   http://milamin.org
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; version 2 of the License.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------
   
   MILAMIN Readme

    Below we provide instructions how to obtain and install Triangle mesh generator
    and SuiteSparse package that are required by MILAMIN. It is recommended to 
    install an optional reordering package METIS and change BLAS used by MATLAB to
    GotoBLAS in order to achieve best performance. Finally, we provide a list of 
	MILAMIN files.
	
   
------------------------------------------
1. Triangle

1.1 Simple binary output support for triangle

    The triangle_bin.diff patch adds basic binary output support for Shewchuk's 2D mesh
	generator, Triangle. The original code can be downloaded from:

	http://www.cs.cmu.edu/~quake/triangle.html

    In order to apply the patch on Linux/Unix systems, download triangle v1.6 
    from the above website and unpack it. Copy the triangle_bin.diff to the same
    directory, 'cd' to that directory and execute the below command:

        patch -p1 < triangle_bin.diff
	
	Patch is available for windows systems from: 
	
	http://gnuwin32.sourceforge.net/packages/patch.htm


1.2 Usage

    By default, triangle outputs mesh information in ASCII. To enable binary
    output supply '-b' flag on the command line. When this is done, triangle
    will output nodes, elements, edges and neighbor information in binary format.
    The output file names do not change. Also, the output format for all of these
    files is the same as defined for original triangle. All data written is
    of type REAL (defaults to double) defined in triangle.c. 

    Note that when using binary output, meshes may not be directly movable between
    different computer architectures due to different endianness.


------------------------------------------
2. SuiteSparse

   MILAMIN requires several function that are part of SuiteSparse package developed 
   by Tim Davis. This package can be downloaded from

      http://www.cise.ufl.edu/research/sparse/SuiteSparse/

   The packages required are CHOLMOD, AMD, and CSparse. For compiling instructions 
   see README files distributed with SuiteSparse. It is recommended to compile SuiteSparse
   package with separete reordering package METIS that can be obtained here

      http://glaros.dtc.umn.edu/gkhome/metis/metis/download


------------------------------------------
3. GotoBLAS
   
   It is recommended to change BLAS used by MATLAB to GotoBLAS that can be downloaded 
   from

      http://www.tacc.utexas.edu/resources/software/

   First, compile GotoBLAS into a dynamically linked library (see GotoBLAS installation
   instructions). Next, before running MATLAB set the following environment variable

      BLAS_VERSION=<full path to GotoBLAS library>/<full name of .so library>

   e.g. on linux you could write:

      $ export BLAS_VERSION=/home/milamin/GotoBLAS/libgoto_opteron-r1.19.so

   To verify that the desired BLAS version is indeed used, also set LAPACK_VERBOSITY=1.
   During startup, MATLAB will print the names of libraries used.
