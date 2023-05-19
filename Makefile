GMSH_DIR := ../gmsh-sdk
CC := gcc
LIB := $(GMSH_DIR)/lib/libgmsh.so
LDFLAGS := -Wl,-rpath,$(GMSH_DIR)/lib -lm

PROG := project.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c src/RCM.c src/optimize.c src/freq.c
OBJS := project.o elasticity.o lu.o matrix.o design.o eigen.o RCM.o optimize.o freq.o

PROG_PY := src/run_optimize.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c src/RCM.c src/optimize.c src/freq.c
OBJS_PY := run_optimize.o elasticity.o lu.o matrix.o design.o eigen.o RCM.o optimize.o freq.o


.PHONY: clean

fec:
	make clean
	$(CC) -c $(PROG)
	$(CC) -o project $(OBJS) $(LIB) $(LDFLAGS)
	rm -f *.o

opti_py:
	make clean
	$(CC) -c $(PROG_PY)
	$(CC) -o opti_py $(OBJS_PY) $(LIB) $(LDFLAGS)
	rm -f *.o

clean:
	rm -f opti_py
	rm -f project
	rm -f *.o

# GMSH_DIR := ../gmsh-4.11.1-Windows64-sdk
# CC := gcc
# LIB := $(GMSH_DIR)/lib/gmsh-4.11.dll
# LDFLAGS := -Wl,-rpath,$(GMSH_DIR)/lib -lm

# PROG := project.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c src/RCM.c src/optimize.c src/freq.c
# OBJS := project.o elasticity.o lu.o matrix.o design.o eigen.o RCM.o optimize.o freq.o

# SIM := simulation.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c src/RCM.c src/optimize.c src/freq.c
# OBJS_SIM := simulation.o elasticity.o lu.o matrix.o design.o eigen.o RCM.o optimize.o freq.o

# PROG_PY := src/opti_for_py.c src/elasticity.c src/lu.c src/matrix.c src/design.c src/eigen.c src/RCM.c
# OBJS_PY := opti_for_py.o elasticity.o lu.o matrix.o eigen.o design.o RCM.o


# .PHONY: clean

# fec:
# 	$(CC) -c $(PROG)
# 	$(CC) -o project $(OBJS) $(LIB) $(LDFLAGS)
# 	del /F *.o
# sim:
# 	$(CC) -c $(SIM)
# 	$(CC) -o project $(OBJS_SIM) $(LIB) $(LDFLAGS)
# 	del /F *.o

# opti_py:
# 	$(CC) -c $(PROG_PY)
# 	$(CC) -o opti_for_py $(OBJS_PY) $(LIB) $(LDFLAGS)
# 	del /F *.o

# clean:
# 	del /F project
# 	del /F *.o


