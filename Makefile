CC=g++
DEPS = input_new.h
OBJ = cocluster+.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< 

cocluster+: $(OBJ)
	$(CC) -o $@ $^ 

clean:
	rm -f cocluster+.o 

