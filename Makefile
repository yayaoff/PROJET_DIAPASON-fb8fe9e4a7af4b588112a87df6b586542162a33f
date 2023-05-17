GMSH_DIR := ../gmsh-sdk
CC := gcc
LIB := $(GMSH_DIR)/lib/libgmsh.so

PROG := project.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c src/RCM.c src/opti.c
PROG_PY := src/opti_for_py.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c src/RCM.c
PROG_HALF := project_half.c src_half/elasticity.c src_half/lu.c src_half/matrix.c src_half/design.c src_half/eigen.c src_half/RCM.c src_half/opti.c
OBJS := project.o elasticity.o lu.o matrix.o eigen.o design.o RCM.o opti.o
OBJS_HALF := project_half.o elasticity.o lu.o matrix.o eigen.o design.o RCM.o opti.o
OBJS_PY := opti_for_py.o elasticity.o lu.o matrix.o eigen.o design.o RCM.o
LDFLAGS := -Wl,-rpath,$(GMSH_DIR)/lib -lm

.PHONY: clean

fec:
	make clean
	$(CC) -c $(PROG)
	$(CC) -o project $(OBJS) $(LIB) $(LDFLAGS)
	rm -f *.o

opti_py:
	make clean
	$(CC) -c $(PROG_PY)
	$(CC) -o opti_for_py $(OBJS_PY) $(LIB) $(LDFLAGS)
	rm -f *.o

half:

	$(CC) -c $(PROG_HALF)
	$(CC) -o project $(OBJS_HALF) $(LIB) $(LDFLAGS)

clean:
	rm -f project
	rm -f *.o

# GMSH_DIR := gmsh-4.11.1-Windows64-sdk
# CC := gcc
# LIB := $(GMSH_DIR)/lib/gmsh-4.11.dll

# PROG := project.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c src/RCM.c src/opti.c
# PROG_HALF := project_half.c src_half/elasticity.c src_half/lu.c src_half/matrix.c src_half/design.c src_half/eigen.c src_half/RCM.c src_half/opti.c
# OBJS := project.o elasticity.o lu.o matrix.o eigen.o design.o RCM.o opti.o
# OBJS_HALF := project_half.o elasticity.o lu.o matrix.o eigen.o design.o RCM.o opti.o
# LDFLAGS := -Wl,-rpath,$(GMSH_DIR)/lib -lm


# .PHONY: clean

# fec:
# 	$(CC) -c $(PROG)
# 	$(CC) -o project $(OBJS) $(LIB) $(LDFLAGS)

# half:

# 	$(CC) -c $(PROG_HALF)
# 	$(CC) -o project $(OBJS_HALF) $(LIB) $(LDFLAGS)
# clean:
# 	del /F project.exe
# 	del /F *.o