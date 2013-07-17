BIN=hw00

hw:	test4 test8
	
test:	test2
	
clean:	
	rm hw00
	rm hw00.o

test2:	$(BIN)
	mpirun -np 2 hw00

test4:	$(BIN)
	@echo ""
	@echo "! [Run test p=4] mpirun -np 4 hw00"
	@-mpirun -np 4 hw00

test8:	$(BIN)
	@echo ""
	@echo "! [Run test p=8] mpirun -np 8 hw00"
	@-mpirun -np 8 hw00

$(BIN):	hw00.f90
	mpif90 -g -C -CB -o hw00 hw00.f90
