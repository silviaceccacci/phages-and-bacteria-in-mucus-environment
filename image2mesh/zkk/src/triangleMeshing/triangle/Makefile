  #
  #
  #----------------------------------------------
  #----------------------------------------------
  #  TRIANGLE compilation
  #----------------------------------------------
  #----------------------------------------------
  #	
CC = cc
CSWITCHES = -O -I/usr/X11R6/include -L/usr/X11R6/lib

triangleDelaunay: triangle.c
	$(CC) $(CSWITCHES) -o triangle.x triangle.c -lm
	@echo '---------------------------->>> END OF TRIANGLE DELAUNAY COMPILATION'
new:
	@rm -rf *.o *.mod core*
	@make

clean:
	@rm -f core core.* *~ *.o *.mod
#------------------------------------------------------------------------------
# All is not set as default make option to avoid compiling tetgen each time
all: triangleDelaunay
tri: triangleDelaunay
