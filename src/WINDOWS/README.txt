INSTRUCTIONS FOR COMPILING LIGGGHTS WITH VISUAL STUDIO 2010/2012/2013
(Ultimate, Professional or Express Versions)

#######################################################################################
IMPORTANT: If you used Git on Windows to get this code, ensure that Unix-style line 
endings are used throughout the project. "Git for Windows" will ask during installation
how it should handle newlines, we recommend Option 2. "Checkout as-is, commit Unix-style
line endings". The simplest way to check if everything is alright is opening
Make.sh in the source folder using a text editor such as Notepad++. Verify that lines
end with LF and not CR LF.
#######################################################################################

#######################################################################################
IMPORTANT: Prior to opening the Visual Studio project, some files must be generated and
updated. This is why the project file only contains lammps.cpp and lammps.h in the
beginning.

LIGGGHTS build routine uses GNU tools before compilation to generate headers.
The easiest way to generate these files is to install Cygwin, a utility which allows 
many Unix utilities to run on Windows. Using these ported unix utilties one can
trigger the file generation.

1. Download the Cygwin installer (https://www.cygwin.com/)
2. Install Cygwin: Beside the core installation, also install Python
2. Open a Cygwin shell and go to source folder (parent folder of this one)

   $ cd /cygdrive/c/your-windows-path-to-your-liggghts-folder/

3. Run the following commands in that folder

   $ sh Make.sh style
   $ sh Make.sh models

4. Verify if files were generated:

   $ ls style_*
   should output a list of style file headers

   $ cat style_contact_model.h
   should output a long list of GRAN_MODEL(....) lines

5. Finally one must update the Visual Studio project in the WINDOWS folder. To do this,
   run the following Python script inside of the WINDOWS folder

   $ cd WINDOWS/
   $ python update_project.py LIGGGHTS.vcxproj

   This will update all headers and implementation files from the LIGGGHTS source
   directory and insert them into the Visual Studio project.
#######################################################################################

To compile LIGGGHTS open the LIGGGHTS_VS2013 Solution

The LIGGGHTS project has configurations to compile either with MPI 
support or with MPI stubs. *

To compile WITH MPI:

1.  Install MS-MPI by downloading the HPC Pack MS-MPI Redistributable Package
    from http://www.microsoft.com/en-us/download/details.aspx?id=41634
    Validate corresponding include and lib directories in the project properties 
    of LIGGGHTS: LIGGGHTS/Properties/Configuration Properties/VC++ Directories **
	Here is a tutorial on MS-MPI: http://www.cs.ucla.edu/~zhu/tutorial/Using_MS-MPI.pdf

2.  Compile LIGGGHTS using Debug or Release configurations from the
    provided projects (use x64 for 64bit binary)

To compile WITHOUT MPI, but instead using MPI STUBS
   
1.  Compile STUBS.vcproj 

2.  Compile LIGGGHTS using Debug_STUBS or Release_STUBS configurations
from the provided project (use x64 for 64bit binary)


* For Visual Studio versions prior to 2013  the Platform Toolset setting has to
be adjusted for each project. This setting can be changed under:
Properties/Configuration Properties/General/General/Platform Toolset