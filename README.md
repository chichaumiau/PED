
#PED: Point Electron Density

===============Install======================================================
1. Download chimera from https://www.cgl.ucsf.edu/chimera/ and install
Find the binary execution of chimera in the 'bin' folder.

2. Replace the line 'chimera_bin=' in 'PED' file with your chimera dir


===============Usage======================================================
./PED pdb_file omap_file
example:	./PED 1ehz.pdb 1ehz.omap

input files: 
1. PDB file #standard PDB format, only 'ATOM' and 'HETATM' lines are considered
2. Electron density map file, omap, available at http://eds.bmc.uu.se/eds/

output files:
1. '.dens' file, list of atoms with point electron density
example:
G	A	1	OP3	99.85	1.16
G	A	1	P	100.19	1.88
G	A	1	OP1	100.19	0.95
G	A	1	OP2	99.21	1.18
G	A	1	O5'	99.82	1.54
G	A	1	C5'	98.63	0.92
G	A	1	C4'	97.84	1.03
G	A	1	O4'	97.10	0.93
G	A	1	C3'	98.07	1.08
G	A	1	O3'	99.39	1.05
G	A	1	C2'	96.96	1.53
G	A	1	O2'	96.77	0.99
Residue=Ch===#Res======Atom===B factor=PointElecDensity

2. '.ent' file, PDB format file, the SegId column is replaced with point electron density
example:
ATOM      1  N   GLN A   1      26.544  15.247  15.544  1.00 40.54      0.48 N  
ATOM      2  CA  GLN A   1      27.006  15.744  16.843  1.00 23.59      2.12 C  
ATOM      3  C   GLN A   1      25.837  16.192  17.710  1.00 19.71      2.69 C  
ATOM      4  O   GLN A   1      25.930  17.211  18.394  1.00 20.35      3.41 O  
ATOM      5  CB  GLN A   1      27.962  16.905  16.627  1.00 27.73      1.56 C  
ATOM      6  CG  GLN A   1      27.432  17.948  15.648  1.00 48.19      1.02 C  
ATOM      7  CD  GLN A   1      27.620  19.361  16.180  1.00 62.99      0.81 C  
ATOM      8  OE1 GLN A   1      28.461  19.615  17.051  1.00 52.44      0.80 O  
ATOM      9  NE2 GLN A   1      26.824  20.288  15.653  1.00 76.44      0.29 N
======================================================================-segId==



Please contact Chichau by email chichaumiau@gmail.com if you have any problem.
