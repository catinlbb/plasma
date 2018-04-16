#include "phg2D.h"
#include "phg.h"

int main(int argc, char* argv[])
{

  GRID *g;

  char *fn = "../mesh_data/plasma2D.mesh";

  phgInit(&argc, &argv);

  printf("@@@Dim=%d\n", Dim);
  g = phgNewGrid(-1);

  printf("@@@@@@@@@ NVert=%d  NEdge=%d\n", NVert, NEdge);
  
  if (!phgImport2D(g, fn, FALSE))
    phgError(1, "can't read file \"%s\".\n", fn);
  
  phgPrintf("Final mesh written to \"%s\".\n",
	    phgExportVTK(g, "equil2D.vtk", NULL));



  
  phgFreeGrid(&g);
  phgFinalize();
  return(0);
}
