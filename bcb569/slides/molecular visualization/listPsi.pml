python
import pymol
for i in range(3,7):
	cmd.select("myr","chain a and resi %d" % (i))
	x = cmd.get_dihedral("myr and n. N","myr and n. C","myr and n. CA","myr and n. O")
	print "residue %d psi = %f" % (i,x)
python end

