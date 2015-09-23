The aim of this program is to simulate the evolution of an open quantum system. The structure of the system is described by a mathematical graph.
External particle currents and dephasing are being considered.
In order to successfully compile the code, you should have installed
the libraries BLAS, LAPACK, LAPACK95 and EXPOKIT on your machine. Some
parallelizations have also been done using the OpenMP-Standard.

The modules have been documented by using Doxygen. Running 'doxygen' in your
shell generates the documentation as LaTex and HTML-files.



How to install:

1.Check the Makefile and specify at the top where your libraries are.
2.Run 'make' in your shell.


'make clean' removes previous results and object-files from the compilation
which are no longer needed.

'make realclean' also removes the compiled executables.


NOTE:
Running 'make res' generates an executable file which runs a simplified
version of the code for a variety of different dephasing-strengths. The
range and sampling points for the dephasing-strength can be set when running
the script. Use the 'res' program if only the resistance over dephasing
strength is of interest, as it is much faster than the whole code.
So far, only the simple sampling works.
