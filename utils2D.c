#include "phg.h"
#include "phg2D.h"
#include <math.h>
#include <string.h>
#include <limits.h>	/* PATH_MAX */



static int btype_index = -1;

static const char *btype_names[] = {
    "dirichlet", "neumann", 
    "user0", "user1", "user2", "user3", "user4", 
    "user5", "user6", "user7", "user8", "user9",
    "undefined", NULL
};
static BTYPE btype_values[] = {
    DIRICHLET, NEUMANN, 
    BDRY_USER0, BDRY_USER1, BDRY_USER2, BDRY_USER3, BDRY_USER4,
    BDRY_USER5, BDRY_USER6, BDRY_USER7, BDRY_USER8, BDRY_USER9, 
    UNDEFINED
};

/* struct mapping input bc in the interval [min,max] to PHG type b */
typedef struct {
    int		min, max, lineno;
    BTYPE	btype;
} bcmap_t;

static bcmap_t *bcmap = NULL;
static int bcmap_n = 0, bcmap_alloc = 0;
static char *bcmap_file = NULL;


static int
map_comp(const void *p0, const void *p1)
/* compares two maps in the order: min, max, btype */
{
    int i;

    if ((i = bcmap[*(const int *)p0].min - bcmap[*(const int *)p1].min) != 0)
	return i;
    if ((i = bcmap[*(const int *)p0].max - bcmap[*(const int *)p1].max) != 0)
	return i;
    return bcmap[*(const int *)p0].btype - bcmap[*(const int *)p1].btype;
}


static BOOLEAN
load_bcmap_file(void)
/* Note: each line in the bc_map_file is of the form:
 * 	val[:val] PHG bc
 * where 'val' can be '*' or omitted, meaning +-infty. For examples:
 * 	-1	DIRICHLET
 * 	0:*	BDRY_USER1
 * The keywords in the second column are case insensitive.
 *
 * Characters between '#' and EOL are ignored (treated as comments) */
{
    FILE *f;
    char *p, *q, *r, line[1024];
    const char *s;
    BOOLEAN star;
    int i, j, k, lineno;
    int *ordering;
    bcmap_t *bcmap1;

    if ((f = fopen(bcmap_file, "r")) == NULL)
	return FALSE;

    lineno = 0;
    while (TRUE) {
	if (fgets(line, sizeof(line), f) == NULL)
	    break;
	lineno++;
	/* delete comment part */
	if ((p = strchr(line, '#')) != NULL)
	    *p = '\0';
	/* delete trailing spaces */
	for (p = line + strlen(line) - 1;
	     p >= line && isspace(*(BYTE *)p);
	     p--);
	p[1] = '\0';
	/* skip leading spaces */
	for (p = line; isspace(*(BYTE *)p); p++);
	/* ignore empty line */
	if (*p == '\0')
	    continue;
	/* find next token */
	for (r = p; *r != '\0' && !isspace(*(BYTE *)r); r++);
	/* skip spaces */
	for (q = r; isspace(*(BYTE *)q); q++);
	if (*q == '\0')
	    phgError(1, "%s: syntax error in %s, line %d.\n", __func__, 
			bcmap_file, lineno);
	*r = '\0';
	/* now p points to the first token, q points to the second token */
	if (bcmap_n >= bcmap_alloc) {
	    bcmap = phgRealloc_(bcmap, (bcmap_alloc + 16) * sizeof(*bcmap),
				bcmap_alloc * sizeof(*bcmap));
	    bcmap_alloc += 16;
	}
	bcmap[bcmap_n].lineno = lineno;
	/* get btype */
	if (!strncasecmp(q, "BDRY_", 5))
	    q += 5;			/* skip optional "BDRY_" prefix */
	for (i = 0; (s = btype_names[i]) != NULL && strcasecmp(s, q); i++);
	if (s == NULL)
	    phgError(1, "%s: syntax error in %s, line %d.\n", __func__, 
			bcmap_file, lineno);
	bcmap[bcmap_n].btype = btype_values[i];
	/* get min value */
	if (*p == '*' || *p == ':') {
	    star = TRUE;
	    bcmap[bcmap_n].min = INT_MIN;
	    r = (*p == '*' ? p + 1 : p);
	}
	else {
	    star = FALSE;
	    bcmap[bcmap_n].min = strtol(p, &r, 10);
	}
	/* get max value */
	if (*r == ':') {
	    if (*(++r) == '*' || *r == '\0') {
		bcmap[bcmap_n].max = INT_MAX;
		if (*r == '*')
		    r++;
	    }
	    else {
		bcmap[bcmap_n].max = strtol(r, &r, 10);
	    }
	}
	else {
	    bcmap[bcmap_n].max = (star ? INT_MAX : bcmap[bcmap_n].min);
	}
	if (*r != '\0')
	    phgError(1, "%s: syntax error in %s, line %d.\n", __func__, 
			bcmap_file, lineno);
	bcmap_n++;
    }

    fclose(f);

    if (bcmap_n == 0)
	return TRUE;

    /* sort and merge intervals in bcmap[] */
    ordering = phgAlloc(bcmap_n * sizeof(*ordering));
    for (i = 0; i < bcmap_n; i++)
	ordering[i] = i;
    qsort(ordering, bcmap_n, sizeof(*ordering), map_comp);
    bcmap1 = phgAlloc(bcmap_n * sizeof(*bcmap1));
    i = k = 0;
    while (i < bcmap_n) {
	bcmap1[k] = bcmap[ordering[i]];
	j = i + 1;
	while (j < bcmap_n && bcmap[ordering[j]].min <= bcmap1[k].max) {
	    /* intervals ordering[i] and ordering[j] overlap */
	    if (bcmap[ordering[i]].btype != bcmap[ordering[j]].btype)
		phgError(1, "%s: inconsistent mapping, line %d <=> line %d\n",
			 __func__, bcmap[ordering[i]].lineno,
				   bcmap[ordering[j]].lineno);
	    if (bcmap[ordering[j]].max > bcmap1[k].max)
		bcmap1[k].max = bcmap[ordering[j]].max;
	    j++;
	}
	i = j;
	k++;
    }
    phgFree(ordering);
    phgFree(bcmap);
    bcmap = bcmap1;
    bcmap_n = k;

    return TRUE;
}



BOOLEAN
phgImport2D(GRID *g, const char *filename, BOOLEAN distr)
{
  FILE *fp;
  char line[128], magic[128];
  BOOLEAN ret = FALSE;
  INT i;
  int flag;
  static BOOLEAN initialized = FALSE;
  static BOOLEAN auto_bcmap = TRUE;

  if(g == NULL)
    phgError(1, "invalid grid pointer.\n");
  
    if (g->nleaf_global > 0) {
	phgWarning("g already contains a valid grid, new grid not imported.\n");
	return ret;
    }

  
    if (phgRank > 0)
	goto hamilton;

    if (phgVerbosity > 0)
	phgPrintf("Default boundary type: %s.\n",
			BTypeName(btype_values[btype_index]));

    if ((fp = phgOpenInputFile_(filename)) == NULL) {
	phgError(1, "cannot read file \"%s\".\n", filename);
	goto hamilton;
    }


    if (fgets(line, sizeof(line), fp) == NULL ||
	sscanf(line, "%s", magic) != 1) {
	phgCloseInputFile_(fp);
	goto hamilton;
    }
    phgCloseInputFile_(fp);

    FreeAtExit(bcmap);
    if (bcmap_file != NULL) {
	phgFree(bcmap);
	bcmap_n = bcmap_alloc = 0;
	if (!load_bcmap_file())
	    phgError(1, "cannot open BC map file \"%s\".\n", bcmap_file);
	phgFree(bcmap_file);
	bcmap_file = NULL;
    }
    else if (auto_bcmap) {
	/* construct bcmap filename by replacing file extension with .bcmap */
	char *p, *q;
	bcmap_file = phgAlloc(PATH_MAX);
	strncpy(bcmap_file, filename, PATH_MAX);
	bcmap_file[PATH_MAX - 1] = '\0';
	p = strrchr(bcmap_file, '.');
	q = strrchr(bcmap_file, '/');
	if (p == NULL || ((q = strrchr(bcmap_file, '/')) != NULL && p < q))
	    p = bcmap_file + strlen(bcmap_file);
	if (PATH_MAX - (p - bcmap_file) > 0)
	    strncpy(p, ".bcmap", PATH_MAX - (p - bcmap_file));
	bcmap_file[PATH_MAX - 1] = '\0';
	phgFree(bcmap);
	bcmap_n = bcmap_alloc = 0;
	if (load_bcmap_file())
	    phgPrintf("Load bcmap file \"%s\" ...\n", bcmap_file);
	phgFree(bcmap_file);
	bcmap_file = NULL;
    }

    
    /* Note:
     *	flag < 0:	input error,
     *	flag == 0: 	the elements don't have refinement types,
     *	flag > 0:  	the elements already have refinement types */
    /*!!! Warning: Only Medit format is implemented right now in 2D*/
    
    if (!strncmp(magic, "DIM:", 4) || !strncmp(magic, "DIM_OF_WORLD:", 13)) {
	flag = phgImportALBERT(g, filename, FALSE);
    }
    else if (!strcmp(magic, "MeshVersionFormatted")) {
	flag = phgImportMedit2D(g, filename, FALSE);
    }
    else if (!strcmp(magic, "CONTROL")) {
        flag = phgImportGambit(g, filename, FALSE);
    }
    else {
	flag = -1;
	phgWarning("can't determine type of input file, mesh not imported.\n");
    }

    ret = (flag >= 0);

    _phg_bcmap(NULL, 0);	/* reset input_list[] */

    if (ret) {
	g->filename = strdup(filename);

	if (phgRank == 0)
	    phgInfo(1, "%"dFMT" elements loaded.\n", g->nleaf);\

	/* compute BoundingBox and volume of the mesh */
	g->volume = 0.;
	if (g->nvert == 0) {
	    for (i = 0; i < Dim; i++) {
		g->bbox[0][i] = g->bbox[1][i] = 0.;
	    }
	}
	else {
	    int j;
	    FLOAT f;

	    for (j = 0; j < Dim; j++) {
		g->bbox[0][j] = g->bbox[1][j] = g->verts[0][j];
	    }
	    for (i = 1; i < g->nvert; i++) {
		COORD *c = g->verts + i;
		for (j = 0; j < Dim; j++) {
		if ((f = (*c)[j]) < g->bbox[0][j]) g->bbox[0][j] = f;
		if (f > g->bbox[1][j]) g->bbox[1][j] = f;
		}
	    }

	    for (i = 0; i < g->nroot; i++) {
		_phg_compute_elem_data(g, g->roots + i, &f, NULL);
		g->volume += f;
	    }
	}

	phgPeriodInit(g);
	phgUpdateNeighbours(g);
	phgRemoveDegeneratedElements(g);

	phgUpdateEdges(g);
#if Dim == 3
	phgUpdateFaces(g);
#endif

	phgCheckConformity(g);

#if ALLOW_CURVED_BOUNDARY
	phgInitBoundaryFunctions(g);
#endif

	phgRefineInit(g, !flag);

hamilton:
#if USE_MPI
	if (phgNProcs > 1) {
	    if (phgRank == 0) {
		/* compute initial order and reorder the root elements.
		 * if fails, then avoid to use RTK repartitioner
                 * serial version */
		flag = phgGridInitOrder(g, FALSE);
		if (flag != 0) {
        	    phgWarning("Error in computing initial  grid order,\n");
        	    phgWarning("Grid order not created.\n");
		}
	    }
	    if (!phgMasterSlave)
		MPI_Bcast(&flag, 1, MPI_INT, 0, phgComm);
	    if (flag != 0) {
		if (!strcmp(phgOptionsGetKeyword("-partitioner"), "rtk")) {
#if USE_PARMETIS
		    if (!phgMasterSlave) {
			if (phgRank == 0)
        		    phgWarning("switching to ParMETIS partitioner.\n");
			phgOptionsSetKeyword("-partitioner", "metis");
		    }
#elif USE_ZOLTAN
		    if (!phgMasterSlave) {
			if (phgRank == 0)
        		    phgWarning("switching to Zoltan partitioner.\n");
			phgOptionsSetKeyword("-partitioner", "zoltan");
		    }
#else
		    if (phgRank == 0)
			phgWarning("poor partitioning results expected.\n");
#endif
		}
            }
	}
	if (phgRank > 0)
	    goto end;
#endif	/* USE_MPI */

	/* phgUpdateBoundaryTypes might have been called in RefineInit
	 * if pre-refinement is enabled */
	if ((g->flags & (VERT_FLAG | EDGE_FLAG | FACE_FLAG | ELEM_FLAG))
	    && g->types_vert == NULL)
	    phgUpdateBoundaryTypes(g);
    }

#if USE_MPI
end:
    if (!phgMasterSlave && phgNProcs > 1) {
	MPI_Bcast(&ret, sizeof(ret), MPI_BYTE, 0, phgComm);
	if (distr) /* force grid distribution */
	    phgRedistributeGrid(g);
    }
#endif

    if ((g->flags & GEOM_FLAG))
	phgGeomInit(g);

    if (ret)
	g->serial_no++;

    return ret;
}

#if ALLOW_CURVED_BOUNDARY

void
phgInitBoundaryFunctions(GRID *g)
{
    EXPR **funcs;
    BYTE *flags;
    ELEMENT *e;
    INT i;
    int j;

    if ((funcs = g->bdry_funcs) == NULL)
	return;

    flags = phgCalloc(g->nvert, sizeof(*flags));

    /* mark boundary vertices (bit 0) */
    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	for (j = 0; j < NFace; j++) {
	    if (e->bound_type[j] & BDRY_MASK) {
		int k;
		for (k = 0; k < NFace; k++)
		    if (k != j)
			flags[e->verts[k]] = 1;
	    }
	}
    }

    while (*funcs != NULL) {
	if (phgVerbosity > 1) {
	    int k = (INT)(funcs - g->bdry_funcs);
	    for (j = 0; j < 4; j++)
		phgPrintf("Boundary function %d = %p, \"%s\"\n", j + k,
				funcs[k], phgDump3DFunction(funcs[k]));
	}
	for (i = 0; i < g->nvert; i++) {
	    COORD *c;
	    if (!flags[i])
		continue;
	    c = g->verts + i;
	    if (fabs(phgEvaluate3DFunction(*funcs, (*c)[0], (*c)[1], (*c)[2],
			0.0)) < 1e-4)
		flags[i] |= 2;
	    else
		flags[i] &= 1;
	}
	for (i = 0; i < g->nroot; i++) {
	    e = g->roots + i;
	    for (j = 0; j < NEdge; j++) {
		if (e->bound_func[j] != -1)
		    continue;
		if ((flags[e->verts[GetEdgeVertex(j, 0)]] & 2) != 0 &&
		    (flags[e->verts[GetEdgeVertex(j, 1)]] & 2) != 0)
		    e->bound_func[j] = (funcs - g->bdry_funcs) / 4;
	    }
	}
	funcs += 4;
    }

    phgFree(flags);
}
FLOAT
phgGetBoundaryError(GRID *g, ELEMENT *e)
/* returns the maximum square of the distances between the middle points of
 * the edges and their projections onto curved boundary */
{
    EXPR **f;
    COORD *v0, *v1;
    FLOAT d, dmax, x, y, z, dx, dy, dz;
    int edge;

    dmax = 0.0;
    for (edge = 0; edge < NEdge; edge++) {
	if (e->bound_func[edge] == -1)
	    continue;
    
	f = g->bdry_funcs + e->bound_func[edge] * 4;
	v0 = g->verts + e->verts[GetEdgeVertex(edge, 0)];
	v1 = g->verts + e->verts[GetEdgeVertex(edge, 1)];
	x = ((*v0)[0] + (*v1)[0]) * .5;
	y = ((*v0)[1] + (*v1)[1]) * .5;
	z = ((*v0)[2] + (*v1)[2]) * .5;
	dx = phgEvaluate3DFunction(f[1], x, y, z, 0.0) - x;
	dy = phgEvaluate3DFunction(f[2], x, y, z, 0.0) - y;
	dz = phgEvaluate3DFunction(f[3], x, y, z, 0.0) - z;
	/* FIXME: sum or max? */
	d = dx * dx + dy * dy + dz * dz;
	if (dmax < d)
	    dmax = d;
    }

    return dmax;
}

#else	/* ALLOW_CURVED_BOUNDARY */

FLOAT
phgGetBoundaryError(GRID *g, ELEMENT *e)
{
    Unused(g);
    Unused(e);
    return 0.0;
}

#endif	/* ALLOW_CURVED_BOUNDARY */

