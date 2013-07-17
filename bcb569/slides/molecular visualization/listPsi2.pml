python
import pymol

f = open('Psi.dat', 'w')

for i in range(1,99):
	cmd.select("myr","chain a and resi %d" % (i))
	x = cmd.get_dihedral("myr and n. N","myr and n. C","myr and n. CA","myr and n. O")
	f.write("residue %d psi = %f\n" % (i,x))
f.close()

python end

