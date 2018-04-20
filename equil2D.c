#include "phg2D.h"
#include "phg.h"

int main(int argc, char* argv[])
{

  GRID *g;

  char *fn = "../mesh_data/plasma2D.mesh";

  RegisterParallelFunction(phgExportVTK2D_);
  phgInit(&argc, &argv);






  
  printf("@@@Dim=%d\n", Dim);
  g = phgNewGrid(-1);

  
  if (!phgImport2D(g, fn, FALSE))
    phgError(1, "can't read file \"%s\".\n", fn);
  
  
  printf("@@@@@@@@@ gele=%d  gvert=%d\n", g->nelem, g->nvert);



#if 0
  phgPrintf("Final mesh written to \"%s\".\n",
	    phgExportVTK2D(g, "equil2D.vtk", NULL));

#endif

  
  phgFreeGrid(&g);
  phgFinalize();
  return(0);
}
