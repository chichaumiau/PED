import sys, os

"""
##########################################################################
#  CLASS:  UserParameters
# 1)  opens main output file
# 2)  collects, parses, and validates user parameters
##########################################################################
"""
class UserParameters:

	def __init__(self, infileName, outfileName):
		self.mapfileName = None
		self.noisemapfileName = None
		self.pdbfileName = None
		self.pdbCode = None
		self.writePlot = False
		self.writePeakList = False
		self.skipMultiConf = False
		self.mapType = 'sigma'
		self.atomSample = 'dynamic'
		self.chiSample = 10
		self.lowerCutoff = '0.3'
		self.upperCutoff = False
		self.upperCutoffValue = '0.8'
		self.buildChi2Chi1 = False
		self.writePeakCoords = False
		self.writeConformers = False
		self.outfileName = outfileName
		self.outfile = None
		self.verbOutfileName = 'ringer'
		self.infileName = infileName
	
	########################################################################	
	# opens main output file
	def open_general_output_file(self):
		
		#open output file
		if self.outfileName == None:
			self.outfileName = 'ringer.out'
		try:
			self.outfile = open('%s' %self.outfileName, 'w')
		except:
			print "ERROR:  Cannot open output files for writing!"
			sys.exit()

		#print header
		self.outfile.write("\nWelcome to Ringer!\n\n")

	########################################################################	
	# reads parameters and parses them
	def parm_reader(self):

		#open output file for writing
		self.open_general_output_file()
		
		#open user parameter file and collect data
		parameters = self.read_input_file()

		#print out header for user parameters section
		self.outfile.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		self.outfile.write("~~~~~~~~~~~~~~~~~~~~~~\n")
		self.outfile.write("Reading in User Parameters...\n\n")
		self.outfile.write("   PARAMETERS\n\n")
		
		#set flag for verbose output
		verbose = False
		
		#check for existence of user parameters
		valid, parm = self.validate(parameters, 'map_name')		
		if valid: self.mapfileName = parm
		else:
			print "ERROR:  Invalid or missing value for map_name (parameter",
			print "for density map file name).  Please add variable to run Ringer."
			sys.exit()
		
		valid, parm = self.validate(parameters, 'pdb_name')
		if valid: self.pdbfileName = parm
		else:
			print "ERROR:  Invalid or missing value for pdb_name (parameter",
			print "for pdb file name).  Please add variable to run Ringer."
			sys.exit()
		tmpName = self.pdbfileName.strip().split('/')
		self.pdbCode = tmpName[-1][:-4]
		
		valid, parm = self.validate(parameters, 'noise_map_name')
		if valid: self.noisemapfileName = parm
		
		#NOTE:  test to convert user parameter to bool
		#USAGE:  bool(myVar == a > 0 and b or c)
		# where myVar will be set to b if a is greater than zero (or True), 
		# and to c otherwise		
		valid, parm = self.validate(parameters, 'map_type', ['sigma', 'volume'])
		if valid:  self.mapType = parm
		else: self.outfile.write('      %-30s%-30s\n'
					  %('map_type', 'sigma'))
		if self.mapType == 'sigma':
			if self.noisemapfileName:
				print "ERROR:  Noise map cannot be used with sigma scaling!!!"

		valid, parm = self.validate(parameters, 'skip_multi_conf', ['on', 'off'])
		if valid: self.skipMultiConf = bool(parm == 'on' > 0 and True or False)
		else: self.outfile.write('      %-30s%-30s\n'
					  %('skip_multi_conf', 'off'))

                valid, parm = self.validate(parameters, 
						'atom_sample_type', ['constant',
						'dynamic', 'experimental'])	
		if valid: self.atomSample = parm
		else: self.outfile.write('      %-30s%-30s\n'
					  %('atom_sample_type', 'dynamic'))
		if self.atomSample == 'dynamic':
			self.read_dynamic_bond_parm()

		valid, parm = self.validate(parameters, 'chi_sample_degree')
		if valid: self.chiSample = parm
		else: self.outfile.write('      %-30s%-30s\n'
					  %('chi_sample_degree', '10'))
		try:
			self.chiSample = int(self.chiSample)
		except ValueError:
			print "ERROR:  chi_sample_degree must be an integer"
			sys.exit()
		
		valid, parm = self.validate(parameters, 'lower_sigma_cutoff')
		if valid: self.lowerCutoff = parm
		else: self.outfile.write('      %-30s%-30s\n'
					  %('lower_sigma_cutoff', '0.3')) 
		try:
			self.lowerCutoff = float(self.lowerCutoff)
		except ValueError:
			print "ERROR: Lower sigma cutoff must be a number"
			sys.exit()
		
		valid, parm = self.validate(parameters, 'upper_sigma_cutoff', ['on', 'off'])
		if valid: self.upperCutoff = bool(parm == 'on' > 0 and True or False)
		else: self.outfile.write('      %-30s%-30s\n'
					  %('upper_sigma_cutoff', 'off'))
		
		if self.upperCutoff:
			valid, parm = self.validate(parameters, 'upper_sigma_cutoff_value')
			if valid: self.upperCutoffValue = float(parm)
			else: self.outfile.write('      %-30s%-30s\n'
						  %('upper_sigma_cutoff_value', '0.8')) 
			try:
				self.upperCutoffValue = float(self.upperCutoffValue)
			except ValueError:
				print "ERROR: Upper sigma cutoff value must be a number"
				sys.exit()		
					  
		valid, parm = self.validate(parameters, 'write_plot', ['on', 'off'])
		if valid: self.writePlot = bool(parm == 'on' > 0 and True or False)
		else: self.outfile.write('      %-30s%-30s\n'
					  %('write_plot', 'off'))
		if self.writePlot:
			verbose = True	  
			
		valid, parm = self.validate(parameters, 'write_chi2chi1', ['on', 'off'])
		if valid: self.buildChi2Chi1 = bool(parm == 'on' > 0 and True or False)
		else:  self.outfile.write('      %-30s%-30s\n'
					  %('write_chi2chi1', 'off'))
		if self.buildChi2Chi1:
			verbose = True
					  			  
		valid, parm = self.validate(parameters,'write_peak_list', ['on', 'off'])
		if valid: self.writePeakList = bool(parm == 'on' > 0 and True or False)
		else: self.outfile.write('      %-30s%-30s\n'
					  %('write_peak_list', 'off'))
		if self.writePeakList:
			verbose = True
			  		  
		if verbose:
			valid, parm = self.validate(parameters, 'verbose_outfile_prefix')
			if valid: self.verbOutfileName = parm
			else: self.outfile.write('      %-30s%-30s\n'
					  %('verbose_outfile_prefix', 'ringer'))
		
		self.outfile.write('\n')

		#print out parameters that were not used
		if len(parameters):
			self.outfile.write("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			self.outfile.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
			self.outfile.write("   INVALID/UNUSED PARAMETERS:\n\n")
			for key,value in parameters.items():
				self.outfile.write("      %-30s%-30s\n" 
						   %(key, value))
			self.outfile.write('\n')
	##########################################################################
	# determine if parameter is valid
	def validate(self, parameters, key, valid_opt=''):

		value = None
		valid_parm = False
		#if parameter does not exist, use default
		if parameters.has_key(key):
		
			#check if parameter is valid
			if valid_opt != '':
				for vo in valid_opt:
					if parameters[key] == vo:
						valid_parm = True
			else:
				valid_parm = True
			#if parameter is valid
			if valid_parm:
				
				value = parameters[key]
				#print parameter to logfile
				self.outfile.write("      %-30s%-30s\n" 
						    %(key, value))
						    
				#remove parameter from total list of parameters
				del parameters[key]

		return valid_parm, value
				
	########################################################################
	# close output files
	def close_general_output_file(self):

		#clean up file
		self.outfile.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		self.outfile.write("~~~~~~~~~~~~~~~~~~~~~~\n")
		self.outfile.write("Calculation complete!\n")
		self.outfile.close()
		
	
	##########################################################################
	# read in user input file
	def read_input_file(self):

		tmp_parms = {}
		try:
		   f1 = open(self.infileName, 'r')
		except:
		   print "ERROR:  Cannot open input parameter file!"
		   sys.exit()

		for next in f1:
		   #skip comments
		   if next[0] == '#':
		       continue
		   tmp = next.strip().split()
		   #skip newlines
		   if len(tmp) < 1:
		       continue
		   #process correct parameters
		   if len(tmp) == 2:
		       tmp_parms[tmp[0]] = tmp[1]	
		   #skip improperly formated parameters and print warning
		   else:
		       self.outfile.write("   WARNING:  Incorrect format ")
		       self.outfile.write("for parameters %s.\n" %next.strip())
		       self.outfile.write("             Parameter will not be ")
		       self.outfile.write("used.\n")
		       continue

		return tmp_parms

	##########################################################################
	# read in atom type, bond, and angle parameters for dynamic bond length
	def read_dynamic_bond_parm(self):

		self.AMBERtype = {}
		path = os.getenv('RINGER')
		f1 = open('%s/parameters/bonds.parm' %path, 'r')
		for next in f1:
			if next[0] == '#':
				continue
			aType = next.strip().split()
			if len(aType) > 2:
				self.AMBERtype[aType[0]] = [float(aType[1]), float(aType[2])]
		f1.close()


