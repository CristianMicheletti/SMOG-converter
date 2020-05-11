# SMOG converter
converter for SMOG input scripts and data files from GROMACS to LAMMPS
Version 1.0 - May 9, 2020 

    Please cite in publications as:
    A. Suma, L. Coronel, G. Bussi and C. Micheletti
    "Directional translocation resistance of Zika xrRNA".
_____________________________________________________

AUTHORS:
--------
A. Suma, L. Coronel and C. Micheletti

CONTENT DESCRIPTION:
--------------------

The package contains:
    
    - the conversion script
    - the source code of LAMMPS extensions required to run SMOG simulations.


Conversion script
-----------------
The code of the conversion tool can be freely used, changed and redistributed for Academic use as long as this file and the header of the code are included.


The script converts SMOG input and data files from GROMACS format to LAMMPS format.
Notice that SMOG input files must be created with the "cut off" contact matrix option.

The output script/data files are for LAMMPS simulations with the style "real", where physical quantities are expressed as:

    distances=Angstroms
    masses=g/mol
    temperature=K
    energies=Kcal/mol
    forces=Kcal/(mol Angstroms)
    time=fs

By default the integration time step is set equal to 2fs, the damp coefficient is equal to 2ps and the masses, which are expressed in g/mol or amu units, are all set to the uniform value of 16.0.
This occurs even if the "heterogeneous masses" SMOG option is chosen. 
Appropriate heterogeneous masses can be directly set in the data file.


Sample usage of the conversion script

The shell command:
    bash GROMACS_to_LAMMPS_conversion.x test
will process files *test.gro* and *test.top* and will convert them to the output files *input_script.lammps* and *data_file.lammps*.



Custom LAMMPS extensions
------------------------

The provided custom Lammps extension can be compiled along with the other source files of LAMMPS.

_____________________________________________________
END
