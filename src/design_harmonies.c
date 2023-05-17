#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "../../gmsh-sdk/include/gmshc.h"
#include "../gmsh-4.11.1-Windows64-sdk/include/gmshc.h"

void designTuningFork(double r1, double r2, double e, double l, double meshSizeFactor, char * filename) {
  /**
   * r1 = inner radius (half-distance between prongs)
   * r2 = outer radius (half-width of fork)
   * e  = length of handle
   * l  = length of prongs
   * meshSizeFactor = meshSize / width of prongs
   * if `filename` is not NULL, save to file
  */
  
  int ierr;

  gmshClear(&ierr);

  double h = r2 - r1; // width of prongs
  double h2=(r2-r1)/4;
  double meshSize = h * meshSizeFactor;

  // Add points
  double x = 0;
  double y = 0;
  double z = 0;
  gmshModelOccAddPoint(x,y,z,meshSize,1,&ierr);
  x += h;
  gmshModelOccAddPoint(x,y,z,meshSize,2,&ierr);
  y += e;
  gmshModelOccAddPoint(x,y,z,meshSize,3,&ierr);
  x += r2;
  y += r2;
  gmshModelOccAddPoint(x,y,z,meshSize,4,&ierr);
  y += l;
  gmshModelOccAddPoint(x,y,z,meshSize,5,&ierr);
  x += r2/2;
  y+=r2/2;
  gmshModelOccAddPoint(x,y,z,meshSize,6,&ierr);
  y += l/2;
  gmshModelOccAddPoint(x,y,z,meshSize,7,&ierr);
  x -= h2;
  gmshModelOccAddPoint(x,y,z,meshSize,8,&ierr);
  y -= l/2;
  gmshModelOccAddPoint(x,y,z,meshSize,9,&ierr);
  x -= r1/2;
  y -= r1/2;
  gmshModelOccAddPoint(x,y,z,meshSize,10,&ierr);
  x-= 6*h2;
  gmshModelOccAddPoint(x,y,z,meshSize,11,&ierr);
  x -= r1/2;
  y+=r1/2;
  gmshModelOccAddPoint(x,y,z,meshSize,12,&ierr);
  y += l/2;
  gmshModelOccAddPoint(x,y,z,meshSize,13,&ierr);
  x -= h2;
  gmshModelOccAddPoint(x,y,z,meshSize,14,&ierr);
  y -= l/2;
  gmshModelOccAddPoint(x,y,z,meshSize,15,&ierr);
  x += r2/2;
  y-=r2/2;
  gmshModelOccAddPoint(x,y,z,meshSize,16,&ierr);
  y-=l;
  gmshModelOccAddPoint(x,y,z,meshSize,17,&ierr);
  x-=r1;
  y-=r1;
  gmshModelOccAddPoint(x,y,z,meshSize,18,&ierr);
  x-=h;
  gmshModelOccAddPoint(x,y,z,meshSize,19,&ierr);
  x-=r1;
  y+=r1;
  gmshModelOccAddPoint(x,y,z,meshSize,20,&ierr);
  y+=l;
  gmshModelOccAddPoint(x,y,z,meshSize,21,&ierr);
  x+=r2/2;
  y+=r2/2;
  gmshModelOccAddPoint(x,y,z,meshSize,22,&ierr);
  y+=l/2;
  gmshModelOccAddPoint(x,y,z,meshSize,23,&ierr);
  x-=h2;
  gmshModelOccAddPoint(x,y,z,meshSize,24,&ierr);
  y-=l/2;
  gmshModelOccAddPoint(x,y,z,meshSize,25,&ierr);
  x-=r1/2;
  y-=r1/2;
  gmshModelOccAddPoint(x,y,z,meshSize,26,&ierr);
  x-=6*h2;
  gmshModelOccAddPoint(x,y,z,meshSize,27,&ierr);
  x-=r1/2;
  y+=r1/2;
  gmshModelOccAddPoint(x,y,z,meshSize,28,&ierr);
  y+=l/2;
  gmshModelOccAddPoint(x,y,z,meshSize,29,&ierr);
  x-=h2;
  gmshModelOccAddPoint(x,y,z,meshSize,30,&ierr);
  y-=l/2;
  gmshModelOccAddPoint(x,y,z,meshSize,31,&ierr);
  x+=r2/2;
  y-=r2/2;
  gmshModelOccAddPoint(x,y,z,meshSize,32,&ierr);
  y-=l;
  gmshModelOccAddPoint(x,y,z,meshSize,33,&ierr);
  x+=r2;
  y-=r2;
  gmshModelOccAddPoint(x,y,z,meshSize,34,&ierr);
  y+=r2;
  gmshModelOccAddPoint(x,y,z,meshSize,35,&ierr);
  x+=h;
  gmshModelOccAddPoint(x,y,z,meshSize,36,&ierr);
  x+=r1;
  y+=l+r2/2;
  gmshModelOccAddPoint(x,y,z,meshSize,37,&ierr);
  x-=2*r1+h;
  gmshModelOccAddPoint(x,y,z,meshSize,38,&ierr);
  x-=h;
  gmshModelOccAddPoint(x,y,z,meshSize,39,&ierr);
  x-=h2;
  gmshModelOccAddPoint(x,y,z,meshSize,40,&ierr);
  x+=2*h2+h;
  gmshModelOccAddPoint(x,y,z,meshSize,41,&ierr);
  x+=2*r1+h-2*h2;
  gmshModelOccAddPoint(x,y,z,meshSize,42,&ierr);
  x+=h2+h;
  gmshModelOccAddPoint(x,y,z,meshSize,43,&ierr);
  x+=h2;
  gmshModelOccAddPoint(x,y,z,meshSize,44,&ierr);




  
  // Add curves
  gmshModelOccAddLine(1,2,1,&ierr);
  gmshModelOccAddLine(2,3,2,&ierr);
  gmshModelOccAddCircleArc(3,36,4,3,&ierr);
  gmshModelOccAddLine(4,5,4,&ierr);
  gmshModelOccAddCircleArc(5,43,6,5,&ierr);
  gmshModelOccAddLine(6,7,6,&ierr);
  gmshModelOccAddLine(7,8,7,&ierr);
  gmshModelOccAddLine(8,9,8,&ierr);
  gmshModelOccAddCircleArc(9,44,10,9,&ierr);
  gmshModelOccAddLine(10,11,10,&ierr);
  gmshModelOccAddCircleArc(11,42,12,11,&ierr);
  gmshModelOccAddLine(12,13,12,&ierr);
  gmshModelOccAddLine(13,14,13,&ierr);
  gmshModelOccAddLine(14,15,14,&ierr);
  gmshModelOccAddCircleArc(15,37,16,15,&ierr);
  gmshModelOccAddLine(16,17,16,&ierr);
  gmshModelOccAddCircleArc(17,36,18,17,&ierr);
  gmshModelOccAddLine(18,19,18,&ierr);
  gmshModelOccAddCircleArc(19,35,20,19,&ierr);
  gmshModelOccAddLine(20,21,20,&ierr);
  gmshModelOccAddCircleArc(21,38,22,21,&ierr);
  gmshModelOccAddLine(22,23,22,&ierr);
  gmshModelOccAddLine(23,24,23,&ierr);
  gmshModelOccAddLine(24,25,24,&ierr);
  gmshModelOccAddCircleArc(25,41,26,25,&ierr);
  gmshModelOccAddLine(26,27,26,&ierr);
  gmshModelOccAddCircleArc(27,40,28,27,&ierr);
  gmshModelOccAddLine(28,29,28,&ierr);
  gmshModelOccAddLine(29,30,29,&ierr);
  gmshModelOccAddLine(30,31,30,&ierr);
  gmshModelOccAddCircleArc(31,39,32,31,&ierr);
  gmshModelOccAddLine(32,33,32,&ierr);
  gmshModelOccAddCircleArc(33,35,34,33,&ierr);
  gmshModelOccAddLine(34,1,34,&ierr);

  // Add wire (closed curve)
  int curveTags[34] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};
  gmshModelOccAddWire(curveTags, 34, 1, 1, &ierr);

  // Add surface
  int wireTags[1] = {1};
  gmshModelOccAddPlaneSurface(wireTags, 1, 100, &ierr);

  // Sync
  gmshModelOccSynchronize(&ierr);

  // Create physical group for surface
  int surfaceTags[1] = {100};
  gmshModelAddPhysicalGroup(2, surfaceTags, 1, -1, "bulk", &ierr);

  // Create physical group for clamped curves
  int clampedCurveTags[3] = {1, 2, 34};
  gmshModelAddPhysicalGroup(1, clampedCurveTags, 3, -1, "clamped", &ierr);

  gmshModelMeshGenerate(2, &ierr);
  // gmshFltkInitialize(&ierr);
  // gmshFltkRun(&ierr);

  if(filename != NULL) gmshWrite(filename, &ierr);
}