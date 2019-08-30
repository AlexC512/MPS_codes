LIBS = -lm -llapack -lnetcdf

all: ft

ft.o: ft.c 
	cc -c ft.c

ft: ft.o
	cc -oft $^ $(LIBS)
