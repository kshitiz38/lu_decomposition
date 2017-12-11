all: convertin convertout computeLU matgen

convertin: convertin.c
	gcc -o convertin convertin.c

convertout: convertout.c
	gcc -o convertout  convertout.c

computeLU: computeLU.c
	mpicc -o computeLU computeLU.c

matgen: matgen.c
	gcc -o matgen matgen.c

clean:
	rm matgen convertin convertout computeLU
