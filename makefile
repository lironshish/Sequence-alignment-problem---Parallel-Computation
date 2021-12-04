build:
	mpicxx -fopenmp -c main.c -o main.o
	nvcc -I./inc -c cuda.cu -o cuda.o
	mpicxx -fopenmp -o exec  main.o cuda.o /usr/local/cuda/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./exec

run:
	mpiexec -np 4 ./exec < input.txt
