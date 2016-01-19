F90 = mpifort
INC = /usr/include/

OBJ = typedef.o mpi_transpose.o rbmat.o ffts.o dnsdata.o channel.o
#flags =  -fall-intrinsics -ffree-line-length-none -std=f2008 -I$(INC) -cpp  -g -fcheck=all -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -lfftw3
flags =   -fall-intrinsics -ffree-line-length-none -std=f2008 -I$(INC) -cpp -lfftw3  -O3

channel: $(OBJ)
	$(F90) $(flags) -o $@ $(OBJ)
%.o : %.f90
	$(F90) $(flags) -c $<
clean: 
	rm *.mod *.o
