"""
##########################################################################
#  CLASS:  Sigma
#     contains sigma values
##########################################################################
"""

class Sigma:

	def __init__(self):
		self.sigma1 = []
		self.sigma2 = []
		self.sigma3 = []
		self.sigma4 = []
		self.peak1 = []
		self.peak2 = []
		self.peak3 = []
		self.peak4 = []

"""
##########################################################################
#  CLASS:  Atom
#     contains information about atom

##########################################################################
"""

class Atom:

	def __init__(self, cAtom):
	      	self.coord = cAtom.coord()
		self.occupancy = cAtom.occupancy
		self.element = cAtom.element
		self.type = cAtom.idatmType
		self.name = cAtom.name
		
"""
##########################################################################
#  CLASS:  Residue
#     contains information about each residue
##########################################################################
"""

class Residue:

	def __init__(self, cRes, pdbCode):
		self.name = cRes.type
		self.pdbCode = pdbCode
		self.fullName = name(cRes)
		self.chiList = []
		self.peakList = []
		self.noiseList = []
		self.altLocs = []
		self.chi1 = False
		self.ring = False
		self.chi2 = False
		self.chi3 = False
		self.chi4 = False
		self.methyl = False
		self.other = True
		
	##########################################################################	
	# set up rotamer class
	def parse_residue(self, movable, pdbRes, parms):

	    # count how many chi angles are in residue
	    confs = count_chi(movable, pdbRes, parms)
		
	    #if residue is not in the parameter set, return
	    #an empty list
	    if confs == None:
    		return
		
	    #categorize residue
	    self.sampleType(pdbRes.type)

	    #create list of chi atoms 
	    self.createChiList(pdbRes, confs, movable)
	    
	    #identify backbone atoms
	    if parms.writeConformers:
	    	self.BBoneCoords[0] = pdbRes.atomsMap["N"][0].coord()
		self.BBoneCoords[1] = pdbRes.atomsMap["CA"][0].coord()
		self.BBoneCoords[2] = pdbRes.atomsMap["C"][0].coord()

	##########################################################################
	# removes individual residues with multiple conformations
	def removeMulti(self, parms):

		removed = []
		if len(self.chiList) > 1:
			#check to see if same residue occurs in multiple tables
			if parms.skipMultiConf:
				removed.append('%s' %(self.fullName))
				self.chiList = []
			else:
				self.chiList = self.chiList[0]
			removed.append('%s' %(self.fullName))
		elif len(self.chiList) > 0:
			self.chiList = self.chiList[0]

		return removed
		
	##########################################################################
	# categorize each residue by type of sampling
	def sampleType(self, name):

		if name in ['SER', 'GLN', 'ASN', 'GLU', 'ASP', 'ARG', \
				'LYS', 'MET', 'CYS', 'LEU']:
			self.other = False
			self.chi1 = True
		if name in ['PHE', 'TYR', 'TRP', 'HIS']:
			self.other = False
			self.ring = True
		if name in ['GLN', 'GLU', 'ARG', 'LYS', 'MET', 'ILE']:
			self.other = False
			self.chi2 = True
		if name in ['LYS', 'ARG', 'MET']:
			self.other = False
			self.chi3 = True
		if name in ['LYS', 'ARG']:
			self.other = False
			self.chi4 = True
		if name in ['ILE', 'LEU', 'THR']:
			self.other = False
			self.methyl = True
			self.chi3 = True
			self.chi4 = True
		if name == 'VAL':
			self.other = False
			self.methyl = True
			self.chi2 = True
			self.chi3 = True
		if name == 'ALA':
			self.other = False
			self.methyl = True
			self.chi1 = True
			
	##########################################################################
	# create list of chi angles
	def createChiList(self, pdbRes, confs, movable):
	
	    #create empty list to hold chi atoms
	    for i in range(confs):
        	self.chiList.append([])

	    #collect atoms to generate chi angles
	    for name in movable[pdbRes.type]:
        	tmp_chi = []
		for i in range(confs):
		    tmp_chi.append([])
		for ialtConf, altConf in enumerate(tmp_chi):
		    for atom in name:
		    	try:
				#temporarily reset the index for atoms that have 
				#single positions in a multi-conformer residue
				if pdbRes.atomsMap[atom][0].element.name == 'H':
					continue
				if ialtConf >= len(pdbRes.atomsMap[atom]):
				    tmp_index = -1
				    tmp_atom = Atom(pdbRes.atomsMap[atom][tmp_index])
				    tmp_chi[ialtConf].append(tmp_atom)
				else:
				    tmp_atom = Atom(pdbRes.atomsMap[atom][ialtConf])
				    tmp_chi[ialtConf].append(tmp_atom)
			except KeyError: 
				continue
		    self.chiList[ialtConf].append(tmp_chi[ialtConf])

	    #if residue has more than one conformation
	    #sort conformations by occupancy 
	    if len(self.chiList) > 1:
	    	self.chiList.sort(lambda a,b: 
			  	    cmp(b[-1][-1].occupancy, a[-1][-1].occupancy))
	##########################################################################
	# create a new chi list from chi1 peak
	def createNewChiList(self, parms):
	
		#append atoms from current chi2 angle
		newChiList = []
		tmpNewChiList = []	

		if len(self.chiList[1]) > 3:
			for a in self.chiList[1][:-1]:
				tmpNewChiList.append(a)
		else:
			for a in self.chiList[1]:
				tmpNewChiList.append(a)
		newChiList.append(tmpNewChiList)

		#create new coordinates for third atom from 
		#secondary peak of chi1
		from chimera import molEdit
		try:
			values = parms.AMBERtype['%s:%s:%s' 
						%(self.chiList[0][1].type, 
						self.chiList[0][2].type, 
						self.chiList[0][3].type)] 
		except KeyError:
			values = [1.53, 111.10]

		coords = molEdit.findPt(self.chiList[0][2].coord,
							self.chiList[0][1].coord,
							self.chiList[0][0].coord,
							values[0], values[1],
							self.peakList[0].peak1[1][1])
		newChiList[0][2].coord = coords

		return newChiList
	
############################################################################
############################################################################
############################################################################
# count the number of chi angles in each residue type

def count_chi(movable, res, parms):
    
    try:
        confs = 0
        for name in movable[res.type]:
            #determine if residue has multiple modeled residues
	    for atom in name[:-1]:
	        if len(res.atomsMap[atom]) > confs:
	            confs = len(res.atomsMap[atom])
    except KeyError:
	#if residue does not exist in full_chi.lib parameter file,
	#if water or glycine, just ignore
	if res.type == 'HOH' or res.type == 'GLY':
	    return None
	else:
            # otherwise, print warning message and continue
	    parms.outfile.write("      WARNING:  No chi angles for ")
	    parms.outfile.write("residue: %s %s %s\n" %(res.type, 
	                                              res.id.chainId, 
						      res.id.position))
            return None

    return confs

##########################################################################
#generate labels for verbose list of peaks
def name(r):
	
	if r.id.chainId != '':
		name = '%s_%s_%s' %(r.type, r.id.chainId, r.id.position)
	else:
		name = '%s_%s' %(r.type, r.id.position)
		
	return name

##########################################################################
# check if built chi1 is in same rotamer well as primary peak
def eval_built_conformer(res):
	
	from chimera import dihedral
	builtChi1 = dihedral(res.chiList[0][3].coord, 
				res.chiList[0][2].coord,
				res.chiList[0][1].coord, 
				res.chiList[0][0].coord)
	
	try:
		badFlag = in_well(res.peakList[0].peak1[0][1], builtChi1)
	except IndexError:
		badFlag = True
	
	badPeakFlag = False
	try:
		if res.peakList[0].peak1[0][0] < 1.0:
			badPeakFlag = True
	except IndexError:
		badFlag = True
	
	return badFlag, badPeakFlag, builtChi1
	
##########################################################################
# compare peak and chi angle
def in_well(primary, peak):
	
	flag = True
	boundaries = []
	#collect bracketed rotamer values
	pad = 10
	boundaries.append(primary + pad)
	boundaries.append(primary - pad)
	
	#shift into proper degree range
	for high_low in boundaries:
		if high_low > 360:
			high_low -= 360
		if high_low < 0:
			high_low += 360
	if peak < 0:
		peak += 360
	if peak > 360:
		peak -= 360
		
	#if in rotamer well
	if peak <= boundaries[0] and peak >= boundaries[1]:
		flag = False
	#deal with edge effects
	if boundaries[0] < boundaries[1]:	
		if peak <= boundaries[0]:
			flag = False
		if peak >= boundaries[1]:
			flag = False
	return flag

##########################################################################
# write out summary from removing multiple conformers

def write_removeMulti(parms, removed):
	
	if len(removed) > 0:
		parms.outfile.write('\n      Removed multi-conformer residues:\n')
		start = 0
		end = 5
		last = False
		for y in range((len(removed)/5)+1):
			if end > len(removed):
				end = len(removed)
				last = True
			comb = ", ".join(removed[start:end])
			if not last:
				comb += ','
			parms.outfile.write('        %s\n' %comb)
			start += 5
			end += 5

##########################################################################
# write out summary of chi1 conformers not built at primary ringer peak

def write_badBuiltConformers(parms, badBuiltRes, badPeakRes):
	
	if len(badBuiltRes) > 0: 
		parms.outfile.write('\n      WARNING:  Following chi1 atoms do not match\n' \
		 			'      Ringer peaks.  Recommend double checking\n' \
					'      backbone position!\n')
		parms.outfile.write('        Residue     Built Chi 1    Ringer Peak Chi\n')
		for r in badBuiltRes:
			try:
				peak = r[0].peakList[0].peak1[0][1]
			except IndexError:
				continue
			if r[1] < 0:
				r[1] += 360
			parms.outfile.write('        %-10s%10.2f%16s\n' \
							%(r[0].fullName, r[1], peak))
			
	else:
		parms.outfile.write('\n      All built sides chains match chi1 Ringer peaks!\n')

	if len(badPeakRes) > 0:
		parms.outfile.write('\n      WARNING:  Following chi1 atoms have electron density\n' \
		 			'      values lower than 1.0 sigma.  Recommend double checking\n' \
					'      backbone position!\n')
		parms.outfile.write('        Residue     Sigma Value\n')
		for r in badPeakRes:
			try:
				peak = '%.2f' %r.peakList[0].peak1[0][0]
			except IndexError:
				peak = 'NO PEAK'
			parms.outfile.write('        %-10s%10s\n' 
							%(r.fullName, peak))

	else:
		parms.outfile.write('\n      Sigma values for all chi1 angles are > 1.0 sigma!\n')
		
		
