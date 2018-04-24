#C makefile
all: serial

serial: serial_diffusion.o
	$(CC) serial_diffusion.o -o serial
  
serial_diffusion.o: serial_diffusoin.c
	$(CC) $(CFLAGS) gcc serial_diffusion.c -c
  
clean:
	rm *.o
	rm serial
