CC = nvcc
LD = nvcc
CFLAGS = 
LIBS = -lcudart
INCLUDE = -Iinclude
OBJS  = $(patsubst src/%.cc,lib/%.o,$(wildcard src/*.cc))
CUOBJS = $(patsubst src/%.cu,lib/%.o,$(wildcard src/*.cu))
EXECS = $(patsubst exe/%.cc,bin/%,$(wildcard exe/*.cc))
EXEOBJS  = $(patsubst exe/%.cc,lib/%.o,$(wildcard exe/*.cc))
CUEXECS = $(patsubst exe/%.cu,bin/%,$(wildcard exe/*.cu))
CUEXEOBJS  = $(patsubst exe/%.cu,lib/%.o,$(wildcard exe/*.cu))


all: $(OBJS) $(CUOBJS) $(EXEOBJS) $(EXECS) $(CUEXEOBJS) $(CUEXECS)

lib/%.o : src/%.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : src/%.cu
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : exe/%.cc
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

lib/%.o : exe/%.cu
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@


bin/% : $(OBJS) $(CUOBJS) lib/%.o
	$(LD) $(LIBS) $(OBJS) $(CUOBJS) lib/$*.o -o bin/$*





clean:
	rm -f $(EXECS) $(CUEXECS) lib/*.o

