import sys
#########################################################################
##MODULE FUNCTIONS########################################################
###########################################################################
# Make a copy of the active volume with values shifted so that the mean is 0
# and scaled so the standard deviation is 1.0.
# author: Tom Goddard

def make_normalized_map(data):
    
	from numpy import mean, std
	m = data.full_matrix()
     	mn = (m - m.mean()) / m.std() 
		
	from VolumeData import Array_Grid_Data
	dn = Array_Grid_Data(mn, data.data.origin, data.data.step, data.data.cell_angles,
    	                	data.data.rotation, name = data.data.name + ' normalized')

	from VolumeViewer import volume_from_grid_data
	vn = volume_from_grid_data(dn, show_data = False)
	vn.initialize_thresholds()
        vn.copy_settings_from(data)
        vn.data.symmetries = data.data.symmetries
	
        return vn 

############################################################################
# Use a crystallographic unit cell map to create a new map that covers a
# PDB model plus some padding.  Written by Tom Goddard.
#
def bounding_xray_map_symmetry(atoms, pad, volume):

     # Get atom positions.
    from _multiscale import get_atom_coordinates
    xyz = get_atom_coordinates(atoms, transformed = True)

    # Transform atom coordinates to volume ijk indices.
    from Matrix import multiply_matrices, xform_matrix
    tf = multiply_matrices(volume.data.xyz_to_ijk_transform,
    			   xform_matrix(volume.model_transform().inverse()))
    from _contour import affine_transform_vertices
    affine_transform_vertices(xyz, tf)
    ijk = xyz

    # Find integer bounds.
    from math import floor, ceil
    ijk_min = [int(floor(i-pad)) for i in ijk.min(axis=0)]
    ijk_max = [int(ceil(i+pad)) for i in ijk.max(axis=0)]

    #gather information about unit cell and symmetry
    ijk_cell_size = getattr(volume.data, 'unit_cell_size', volume.data.size)
    syms = volume.data.symmetries
    step = volume.data.step
    from VolumeFilter import cover
    cg = cover.map_covering_box(volume, ijk_min, ijk_max, ijk_cell_size, syms, step)
    
    from VolumeViewer import volume_from_grid_data
    v = volume_from_grid_data(cg, show_data = False)
    v.copy_settings_from(volume, copy_region = False)
    
    return v	
	
##########################################################################
#opens map file and normalized data
def open_map(parms, mapFileName, structure):

    from chimera import openModels
    mapCURR = openModels.open(mapFileName)[0]

    if parms.mapType == 'sigma':
        #normalize map to mean of 0 and stdev of 1
        scaledDen = make_normalized_map(mapCURR)
    elif parms.mapType == 'volume':
        scaledDen = mapCURR
    else:
        print 'ERROR:  Unrecognized map type!'
        sys.exit()
    
    #extend map to 5.0 A beyond pdb coordinates
    return bounding_xray_map_symmetry(structure.atoms, 5.0, scaledDen)
