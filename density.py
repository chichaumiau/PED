#!/usr/bin/python

#===========================================================
#Copyright(c)2015, IBMC, CNRS
#All rights reserved.
#NAME:		density.py
#ABSTRACT:	
#DATE:		Thu Nov 12 22:19:33 2015
#Usage:		
#VERSION: 	0.01
#AUTHOR: 	Miao Zhichao
#CONTACT: 	chichaumiau AT gmail DOT com
#NOTICE: This is free software and the source code is freely
#available. You are free to redistribute or modify under the
#conditions that (1) this notice is not removed or modified
#in any way and (2) any modified versions of the program are
#also available for free.
#		** Absolutely no Warranty **
#===========================================================

import sys,os
import chimera

import ringerParms
from parms import UserParameters
import densityMap
import pdbResidues

##########################################################################
#read information for structure file
def mapFile(parameters, mapfilename, pdb):
	try: #catch if map is too big
		density = densityMap.open_map(parameters, mapfilename, pdb[0])
	except MemoryError:
		parameters.outfile.write("\nERROR: Memory allocation exceeded " 
				"while reading in %s map." %mapfilename)
		chimera.openModels.close(pdb)
		sys.exit()
	return density

def density(fpdb,fmap):
	parameters = UserParameters(None, None)
	pdb = chimera.openModels.open(fpdb)
	map1 = mapFile(parameters, fmap, pdb)
	out=''
	dd=[]
	for res in pdb[0].residues:
		for atm in res.atoms:
			d=map1.interpolated_values([atm.coord()])[0]
			out+='%s\t%s\t%s\t%s\t%.2f\t%.2f\n'%(res.type, res.id.chainId, res.id.position,atm.name,atm.bfactor,d)
			dd.append(d)
	open(fpdb.replace('.pdb','.dens'),'w').write(out)
	ll=[]
	for line in file(fpdb):
		if line.startswith('ATOM') or line.startswith('HETATM'):
			ll.append(line)
	if len(dd) != len(ll):
		sys.stderr.write("ERROR:	Length different: %d - %d\n"%(len(dd),len(ll)))
		exit(-1)
	out=''
	for i,j in zip(ll,dd):
		out+='%s%6.2f%s'%(i[:70],j,i[76:])
	open(fpdb.replace('.pdb','.ent'),'w').write(out)
	parameters=None
	pdb=None
	map1=None

if __name__ == '__main__':
	density(sys.argv[1],sys.argv[2])
	
	
