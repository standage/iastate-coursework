all:	
	gcc -Wall -O3 -g -lm -o hw01 hw01.c Array.c AminoAcid.c Atom.c Geometry.c Vector.c
	./hw01 2GB1.pdb 2GB1-new.pdb

clean:	
	rm -rf hw01 hw01.dSYM

mem:	
	gcc -Wall -O3 -g -lm -o hw01 hw01.c Array.c AminoAcid.c Atom.c Geometry.c Vector.c
	valgrind ./hw01 2GB1.pdb 2GB1-new.pdb
