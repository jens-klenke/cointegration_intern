James G. MacKinnon, "Numerical Distribution Functions for Unit Root and
Cointegration Tests", Journal of Applied Econometrics, Vol. 11, No. 6, 
1996, pp. 601-618.

The computer programs and data files for this paper are copyright by the
author.  They may be used freely for non-commercial purposes, provided
proper attribution is given; please cite the paper in all publications. 
However, no part of them may be incorporated in any form into any non-free
publication or computer program without the express, written consent of
the author. 

The data files have been zipped to save space and bandwidth.  There are
two zip files, one called tabs-dos.zip and one called tabs-unix.zip.
The only difference between them is that the files in tabs-dos.zip are
in DOS format (with CR/LF pairs at the ends of lines), while the files
in tabs-unix are in Unix format (with LF at the ends of lines).

Each zip file contains thirteen data files:  probs.tab, and urc-1.tab
through urc-12.tab.  If you will be using a DOS, Windows, or OS/2
system, you should unzip the file tabs-dos.zip on such a system, using
pkunzip or some other program.  Do not unzip this file on a Unix system.
If you will be using a Unix system, it should not matter where or how
you unzip tabs-unix.zip.

For reference, the size of probs.tab is 4420 bytes in the DOS version
and 4199 bytes in the Unix version.

The program urcdist.f is written in Fortran 77.  It should be compilable
on any Unix system or on a PC with a 32-bit Fortran compiler.  On a
typical Unix system, the command

      f77 -O urcdist.f -o urcdist

will compile the program and create an executable file called urcdist.
This program should also compile correctly using f2c + gcc. The default
location for the .tab files is /usr/local/urcdist. Users of non-Unix
systems will have to change this, and users of Unix systems may wish to
change it.

For the benefit of PC users who do not have access to a Fortran compiler, an
executable file, urcdist.exe (which is zipped in urcdist.zip), is also
provided.  It requires a modern computer running Windows. Precisely which
versions of Windows it works on, I do not know, since I never run Windows
myself. It was compiled using CVF 6.6 on a machine running Windows 2000. The
default location for the .tab files is \urcdist. 

The urcdist program searches for the .tab files in the current directory
and in the default location.  If it cannot find these files, it will
stop.

For those who wish to incorporate these routines into their own program,
a separate file called urcrouts.f is provided.  Simply call urcval with
the correct arguments; see the file itself for documentation.  It will
read the necessary .tab files and call other routines as needed.

The file jgm-data.dat contains the data used in the empirical example of
Section 7. The first 43 lines of the file contain the data, as six columns
of numbers, and the remainder of the file contains documentation. Because 
this file is very small, it has not been zipped.

Warning: An error was discovered in all programs on 2003-5-5. This error
causes small discrepancies in P values and critical values in some cases.
It has been corrected in the source files, but not yet in the Windows
executable.
