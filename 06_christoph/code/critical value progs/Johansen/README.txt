This directory contains data and programs for the paper

  James G. MacKinnon, Alfred A. Haug, and Leo Michelis, "Numerical 
  distribution functions of likelihood ratio tests for cointegration,"
  Journal of Applied Econometrics, 14, 1999, 563-577.

All programs and data are copyright by the authors.  They may be used
freely for non-commercial purposes, provided proper attribution is
given; please cite the paper in all publications.  However, no part of
them may be incorporated in any form into any non-free publication or
computer program without the express, written consent of the authors.

The data files have been zipped to save space and bandwidth. They are in
the zip file tabs-all.zip. The files in this zip file are in DOS format,
with CR/LF pairs at the ends of lines. Unix users should use the -a option
of unzip when unzipping this file, or they should use some other method to
strip off the unneeded CRs.

The zip file contains 121 files: probs.tab, joh-1.tab through joh-12.tab,
and lrc-xx-y.tab for xx=01 to 12 and y=0 to 8. On a DOS, Windows, or OS/2
system, it is recommended that these be put in the directory \URCDIST. On
a Unix system, it is recommended that they be put in /usr/local/urcdist.

For those who are interested only in the case k=0 (no exogenous I(1) 
variables), the file tabs-k0.zip contains a subset of the files in
tabs-all.zip. The programs will work correctly for the case k=0 if only
these files are present. The file tabs-extra.zip contains the remaining
files that are in tabs-all.zip and not in tabs-k0.zip. 

The program lrcdist.f is written in Fortran 77.  It should be compilable
on any Unix system or on a PC with a 32-bit Fortran compiler.  On a
typical Unix system, the command

      f77 -O lrcdist.f -o lrcdist

will compile the program and create an executable file called lrcdist; 
the compiler may be called g77, fort77, or something else instead of f77.
Users of non-Unix system will have to change the default location for the
.tab files before compiling.

For the benefit of PC users who do not have access to a Fortran compiler, an
executable file, lrcdist.exe (which is also zipped in lrcdist.zip), is
provided. It requires a modern computer running Windows. Precisely which
versions of Windows it works on, I do not know, since I never run Windows
myself. It was compiled using CVF 6.6 on a machine running Windows 2000. The
default location for the .tab files is \urcdist. NOTE: This file contains a
minor error, which will affect reported values to a small extent.

The lrcdist program searches for the .tab files in the current directory
and in the default location.  If it cannot find these files, it will
stop.

For those who wish to incorporate these routines into their own programs,
two separate files called johrouts.f and lrcrouts.f are provided. These
are also included in the file lrc-source.zip.  The former computes P
values and critical values for cases 0, 1, 1*, 2, and 2*. The latter
computes them for cases I, II, III, IV, and V. Both sets of routines are
incorporated into the lrcdist program. Your program must call either
johval or lrcval with the correct arguments; see the file itself for
documentation.  The routine will read the necessary .tab files and call
other routines as needed. 

List of files:

	f77l3.eer	needed by lrcdist.exe if there is an error
	johrouts.f	routines for cases 0, 1, 1*, 2, and 2*.
	lrc-source.zip	contains all .f files
	lrcdist.exe	executable for Windows systems
	lrcdist.zip	zipped executable for Windows systems
	lrcdist.f	source file for program for all cases
	lrcdist.for     source file for Windows version
	lrcrouts.f	routines for cases I, II, III, IV, and V.
	tabs-all.zip	contains all .tab file in DOS format
        tabs-k0.zip     contains only .tab files for k=0 only
        tabs-extra.zip	contains files in tabs-all that are not in tabs-k0.

Because I have no access to a Fortran compiler for Windows, the Windows
executable has not been updated to correct an error that was found on
2003-5-05. This error may cause very small discrepancies in P values
and critical values. All source files have been updated. Users are
therefore urged to use the source, not the Windows executable.

Sept. 14, 2004
