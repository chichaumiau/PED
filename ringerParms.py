###########################################################################
#identify chi angles that should be rotated for each residue

def read_res_templates(cdir):
    f1 = open('%s/parameters/full_chi.lib' %cdir, 'r')
    movable = {}
    for next in f1:
	tmp = next.strip().split()
	if len(tmp) > 1:
	    if tmp[0] == '#':
	        continue
	    if movable.has_key(tmp[0]):
	        movable[tmp[0]].append(tmp[2:6])
	    else:
                movable[tmp[0]] = [tmp[2:6]]
    return movable
