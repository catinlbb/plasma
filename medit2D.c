#include "phg.h"
#include "phg2D.h"

#include <string.h>
#include <stdlib.h>


static int get_token_ret = 0;
static BOOLEAN
get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (TRUE) {
	if (fscanf(fp, "%s", token) != 1)
	    return get_token_ret = FALSE;
	if (token[0] != '#')
	    break;
	/* skip to newline */
	do
	    if ((c = fgetc(fp)) == EOF)
		return get_token_ret = FALSE;
	while (c != '\n');
    }
    if ((p = strchr(token, '#')) != NULL)
	*p = '\0';

    phgInfo(4, "  got token \"%s\"\n", token);

    return get_token_ret = TRUE;
}



static GRID *g_;

static void
process_triangle(GRID *g, INT nlist, INT (*list)[3], INT *vtypes)
/* assign boundary types (e->bound_type[]) using the list of triangles */
{
    /* Note: each face is stored as (elem, face, type) */
    INT i;
    INT key[3], (*found)[3];
    ELEMENT *e;
    int j;

    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	key[0] = i;
	/* element type was saved in mark */
	for (j = 0; j < NFace; j++) {
	    e->bound_type[j] = UNDEFINED;
	    key[1] = j;
	    /* adjust bdry type according to types of triangles */
	    if (nlist == 0)
		continue;
	    found = bsearch(key, list, nlist, sizeof(*list), comp_face2);
	    if (found != NULL) {
		e->bound_type[j] = _phg_bcmap(g_, (*found)[2]);
		/* printf(" ! set e:[%d,%d], %d:%d\n", e->index, j,  */
		/*        (*found)[3], _phg_bcmap(g_, (*found)[3])); */
	    }
	}
    }
}




int
phgImportMedit2D(GRID *g, const char *filename, BOOLEAN parallel)
{
  FILE *fp;
  char token[1024]; 
  int i, n;
  INT *vtypes = NULL;			/* vertex types */
  INT nlist2 = 0, (*list2)[3] = NULL;	/* list of edges */

  INT nlist4 = 0, (*list4)[5] = NULL;	/* list of quadrilaterals */
  INT nlistt = 0;			/* number of triangles */

#define READ_NUMBER						\
  if (!get_token(fp, token)) strcpy(token, "End");		\
  if (isalpha((int)(token[0]))) {				\
    phgWarning("fewer entries (%d) than expected.\n", i);	\
    break;							\
  }

  FunctionEntry;

  (void)(parallel);	/* avoid gcc warning */

  if (phgRank >= g->nprocs)
    Return 0;

  if ((fp = phgOpenInputFile_(filename)) == NULL) {
    phgWarning("can't open mesh file \"%s\"!\n", filename);
    phgFreeGrid(&g);
    Return -1;
  }

  g_ = g;

  phgInfo(1, "loading file \"%s\".\n", filename);

  token[0] = '\0';
  if (!get_token(fp, token) || strcasecmp(token, "MeshVersionFormatted") ||
      !get_token(fp, token) || strcmp(token, "1") ||
      !get_token(fp, token) || strcasecmp(token, "Dimension") ||
      !get_token(fp, token) || strcmp(token, "2")) {
    n = __LINE__;
  error:
    phgError(1, "(%s:%d) invalid input file: ret=%s, token=%s\n",
	     __FILE__, n, get_token_ret ? "TRUE" : "FALSE", token);
  }



#undef ERROR
#define ERROR	{n = __LINE__; goto error;}

  while (TRUE) {
    if (!get_token(fp, token))
      break;
  next_token:
    if (!strcasecmp(token, "End")) {
      break;
    }
    else if (!strcasecmp(token, "Vertices")) {
      BOOLEAN flag = FALSE;
      if (!get_token(fp, token))
	ERROR
	  n = atoi(token);
      phgInfo(2, "number of vertices: %d\n", n);
      g->verts = phgNewVertices(n);
      vtypes = phgAlloc(n * sizeof(*vtypes)); 
      for (i = 0; i < n; i++) {
	READ_NUMBER
	  g->verts[i][0] = atof(token);
	if (!get_token(fp, token))
	  ERROR
	    g->verts[i][1] = atof(token);
	if (!get_token(fp, token))
	  ERROR
	    /* the type of the vertex */
	    if (!get_token(fp, token))
	      ERROR
		if ((vtypes[i] = atoi(token)) != 0)
		  flag = TRUE;
      }
      g->nvert_global = g->nvert = i;
      if (!flag) {
	phgFree(vtypes);
	vtypes = NULL;
      }
      if (i < n)
	goto next_token;
    }
    else if (!strcasecmp(token, "Edges")) {
      if (!get_token(fp, token))
	ERROR
	  n = atoi(token);
      phgInfo(2, "number of edges: %d \n", n);
      list2 = phgRealloc_(list2, (n + nlist2) * sizeof(*list2),
			  nlist4 * sizeof(*list2));
      for (i = 0; i < n; i++) {
	READ_NUMBER
	  list2[nlist2][0] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR
	    list2[nlist2][1] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR
	    list2[nlist2][2] = atoi(token);	/* boundary type */
	/* sort the two vertices */
	qsort(list2[nlist2++], 2, sizeof(INT), phgCompINT);

      }
      if (nlist2 == 0) {
	phgFree(list2);
	list2 = NULL;
      }
      if (i < n)
	goto next_token;
    }
    else if (!strcasecmp(token, "Triangles")) {
      if (!get_token(fp, token))
	ERROR
      n = atoi(token);
      phgInfo(2, "number of triangles: %d\n", n);
      g->roots = phgReallocElements(g->roots, nlistt, n + nlistt);
      for (i = 0; i < n; i++) {
	READ_NUMBER
        g->roots[nlistt].verts[0] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR
	g->roots[nlistt].verts[1] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR
	g->roots[nlistt].verts[2] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR

	/* region_mark*/    
	g->roots[nlistt].region_mark =
	    g->roots[nlistt].mark = atoi(token);
	nlistt++;
      }
      if (i < n)
	goto next_token;
    }
    else if (!strcasecmp(token, "Quadrilaterals")) {
#if 0 /*!!!! Warning: Not implemented yet for 2D elements*/
      if (!get_token(fp, token))
	ERROR
	  n = atoi(token);
      phgInfo(2, "number of quadrilaterals: %d\n", n);
      list4 = phgRealloc_(list4, (n + nlist4) * sizeof(*list4),
			  nlist4 * sizeof(*list4));
      for (i = 0; i < n; i++) {
	READ_NUMBER
	  list4[nlist4][0] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR
	    list4[nlist4][1] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR
	    list4[nlist4][2] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR
	    list4[nlist4][3] = atoi(token) - 1;
	if (!get_token(fp, token))
	  ERROR
	    list4[nlist4][4] = atoi(token);	/* boundary type */
	/* sort the four vertices */
	qsort(list4[nlist4++], 4, sizeof(INT), phgCompINT);
      }
      if (nlist4 == 0) {
	phgFree(list4);
	list4 = NULL;
      }
      if (i < n)
	goto next_token;
#endif      
    }
    else if (!strcasecmp(token, "Corners")) {
      if (!get_token(fp, token))
	ERROR
	  n = atoi(token);
      for (i = 0; i < n; i++) {
	READ_NUMBER
      }
      if (i < n)
	goto next_token;
    }
    else if (!strcasecmp(token, "RequiredVertices")) {
      if (!get_token(fp, token))
	ERROR
	  n = atoi(token);
      for (i = 0; i < n; i++) {
	READ_NUMBER
      }
      if (i < n)
	goto next_token;
    }
    else if (!strcasecmp(token, "Ridges")) {
      if (!get_token(fp, token))
	ERROR
	  n = atoi(token);
      for (i = 0; i < n; i++) {
	READ_NUMBER
	  }
      if (i < n)
	goto next_token;
    }
    else if (!strcasecmp(token, "RequiredEdges")) {
      if (!get_token(fp, token))
	ERROR
	  n = atoi(token);
      for (i = 0; i < n; i++) {
	READ_NUMBER
	  }
      if (i < n)
	goto next_token;
    }
    else {
#if 0
      int lineno = 1;
      size_t pos, current_pos = ftell(fp);
      fseek(fp, 0, SEEK_SET);
      for (pos = 0; pos < current_pos; pos++)
	if (fgetc(fp) == '\n')
	  lineno++;
      phgCloseInputFile_(fp);
      phgError(1, "unkown token \"%s\" near line %d.\n", token, lineno);
#else
      /* tetgen 1.4.0 seems to generate inconsistent (smaller) number of
       * triangles, as a workaround extra numbers are silently ignored */
      static BOOLEAN warned = FALSE;
      if (!warned) {
	int lineno = 1;
	size_t pos, current_pos = ftell(fp);
	warned = TRUE;
	fseek(fp, 0, SEEK_SET);
	for (pos = 0; pos < current_pos; pos++)
	  if (fgetc(fp) == '\n')
	    lineno++;
	phgWarning("extra token \"%s\" near line %d ignored.\n",
		   token, lineno);
      }
      continue;
#endif
    }
  }

  phgCloseInputFile_(fp);

  g->nelem_global = g->ntree = g->nleaf = g->nelem = g->nroot = nlistt;

  if (g->nroot > 0)
    process_triangle(g, nlist2, list2, vtypes);

  if (list2 != NULL) {
    phgFree(list2);
    list2 = NULL;
    nlist2 = 0;
  }

  if (list4 != NULL) {
    phgFree(list4);
    list4 = NULL;
    nlist4 = 0;
  }

  if (vtypes != NULL) {
    phgFree(vtypes);
    vtypes = NULL;
  }

  Return 0;
  return 0;	/* for avoiding MIPS cc warning */
}
