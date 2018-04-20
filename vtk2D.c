#include <stdarg.h>
#include <limits.h>	/* PATH_MAX */

#include "phg2D.h"
#include "phg.h"


static BOOLEAN vtk_nonlinear_elements = TRUE;

enum {VTK_SINGLE = 0, VTK_DOUBLE = 1};
static const char *vtk_precision_list[] = {"single", "double", NULL};
static int vtk_precision = VTK_SINGLE;


/* buffer for broadcasting filename */
static const char *filename;
/* buffer for broadcasting list of DOF ids */
static unsigned short *dof_ids = NULL;	/* last element for filename length */
static unsigned short dof_ids_len = 0;

static unsigned char *buffer = NULL;
static size_t buffer_pos, buffer_size;



void
phgExportVTK2D_(GRID *g)
/* Output grid in VTK unstructured grid format. */
{
  int k, a, dim;
    size_t size, vtk_float_size;
    INT i;
    INT nleaf, nvert, nedge;
    ELEMENT *s;
    DOF *dof;
    BOOLEAN has_cell_data = FALSE, has_point_data = FALSE;
    FLOAT *values;
    static char fn[PATH_MAX];
    unsigned char c, output_order;

    ParallelFunction(phgExportVTK_, TRUE);

    if (g == NULL || g->comm == MPI_COMM_NULL)
	return;

    vtk_float_size = (vtk_precision == VTK_SINGLE ?
				sizeof(float) : sizeof(double));
    
    printf("@@@ %d %d %d    %d  %d\n",nleaf,nvert,nedge, g->rank, g->flags );
    nleaf = g->nleaf;
    nvert = g->nvert;
    nedge = g->nedge;

    assert((g->flags & VERT_FLAG));
    assert(g->rank != 0 || sizeof(a) == 4);	/* to be fixed later */



#if XJ

#if USE_MPI
    if (phgMasterSlave) {
	if (g->rank == 0) {
	    memcpy(fn, filename, a = strlen(filename) + 1);
	    dof_ids[dof_ids_len - 1] = a ;
	}
	filename = fn;
	MPI_Bcast(&dof_ids_len, sizeof(dof_ids_len), MPI_BYTE, 0, g->comm);
	MPI_Bcast(dof_ids, sizeof(*dof_ids)*dof_ids_len, MPI_BYTE, 0, g->comm);
	a = dof_ids[dof_ids_len - 1];
	MPI_Bcast((void *)filename, a, MPI_BYTE, 0, g->comm);
    }
    else
#else
    {
	strcpy(fn, filename);
	filename = fn;
    }
#endif

    filename = phgIOAddExtension(filename, "vtk");
    if (!phgIOOpen(g, (char *)filename)) {
	phgError(0, "cannot create VTK file \"%s\".\n", filename);
	return;
    }
    if (g->rank == 0)
	phgInfo(1, "creating VTK file \"%s\".\n", filename);

    /* determine output order */
    if (vtk_nonlinear_elements == FALSE) {
	output_order = 1;
    }
    else {
	output_order = 0;
	for (a = 0; a < dof_ids[0]; a++) {
	    dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	    if (DofIsHP(dof)) {
		if (output_order < dof->hp->max_order)
		    output_order = dof->hp->max_order;
	    }
	    else {
		if (output_order < dof->type->order)
		    output_order = dof->type->order;
	    }
	}
    }

    if (output_order > 1 && g->period != NULL) {
	static BOOLEAN warned = FALSE;
	if (warned == FALSE && g->rank == 0)
	    phgWarning("nonlinear element for periodic boundaries unimplemented yet, using linear element instead.\n");
	warned = TRUE;
	output_order = 1;
    }

    if (output_order <= 1) {
	/* use linear tetrahedra */
	output_order = 1;
	buffer_size = nvert * vtk_float_size * 3;
	if ((size = nleaf * 5 * sizeof(a)) > buffer_size)
	    buffer_size = size;
    }
    else {
	/* use quadratic tetrahedra */
	output_order = 2;
	buffer_size = (nvert > nedge ? nvert : nedge) * vtk_float_size * 3;
	if ((size = nleaf * 11 * sizeof(a)) > buffer_size)
	    buffer_size = size;
    }

    /* compute size of DOF data */
    for (a = 0; a < dof_ids[0]; a++) {
	dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
        dim = DofDim(dof);
	dim = dim > 4 ? 4 : dim;
	if (dof->type == DOF_P0 || dof->type == DOF_DG0) {
	    size = nleaf * dim * vtk_float_size;
	    has_cell_data = TRUE;
	}
	else {
	    if (output_order <= 1)
		size = nvert * dim * vtk_float_size;
	    else
		size = (nvert + nedge) * dim * vtk_float_size;
	    has_point_data = TRUE;
	    if (DofIsHP(dof)) {
		has_cell_data = TRUE;	/* hp orders are saved as cell data */
		if (buffer_size < size)
		    buffer_size = size;
		size = nleaf * sizeof(unsigned char);
	    }
	}
	if (buffer_size < size)
	    buffer_size = size;
    }
    buffer = phgAlloc(buffer_size);

    /* File header and vertices */
    phgIOPrintf(g, "# vtk DataFile Version 2.0\n"
		"Tetrahedral grid created by PHG %s.\n"
		"BINARY\nDATASET UNSTRUCTURED_GRID\n"
		"POINTS %"dFMT" %s\n", PHG_VERSION,
		g->nvert_global + (output_order > 1 ? g->nedge_global : 0),
		vtk_precision == VTK_DOUBLE ? "double" : "float");
    buffer_pos = 0;
    for (i = 0; i < nvert; i++) {
	for (a = 0; a < 3; a++) {
	    VTK_PUT_FLOAT(g->verts[i][a])
	}
    }
#if USE_MPI
    phgIOWrite(g, buffer, vtk_float_size * 3, nvert, g->nvert_global,
		g->types_vert, g->L2Gmap_vert, TRUE);
#else
    phgIOWrite(g, buffer, vtk_float_size * 3, nvert, g->nvert_global,
		NULL, NULL, FALSE);
#endif

    if (output_order > 1) {
	/* save list of nodes at edge centers */
	VEF_MAP *vef = phgDofSetupVEFMap(g, NULL, EDGE_FLAG);
	assert(vef != NULL);
	buffer_pos = 0;
	for (i = 0; i < nedge; i++) {
	    COORD *p0, *p1;
	    s = vef->Emap[i];
	    if (s == NULL) {
		/* unreferenced edge */
		VTK_PUT_FLOAT(0.)
		VTK_PUT_FLOAT(0.)
		VTK_PUT_FLOAT(0.)
		continue;
	    }
	    k = vef->Eind[i];
	    p0 = g->verts + s->verts[GetEdgeVertex(k, 0)];
	    p1 = g->verts + s->verts[GetEdgeVertex(k, 1)];
	    for (a = 0; a < 3; a++) {
		VTK_PUT_FLOAT(((*p0)[a] + (*p1)[a]) * 0.5)
	    }
	}
#if USE_MPI
	phgIOWrite(g, buffer, vtk_float_size * 3, nedge, g->nedge_global,
		   g->types_edge, g->L2Gmap_edge, TRUE);
#else
	phgIOWrite(g, buffer, vtk_float_size * 3, nedge, g->nedge_global,
		   NULL, NULL, FALSE);
#endif
	phgDofFreeVEFMap(&vef);
    }
    
    /* tetrahedra */
    phgIOPrintf(g, "\nCELLS %"dFMT" %"dFMT"\n", g->nelem_global,
		g->nelem_global * (output_order > 1 ? 11 : 5));
    buffer_pos = 0;
    ForAllElements(g, s) {
	/* Note: nodes in a VTK_QUADRATIC_TETRA are defined as:
	 * 	v0, v1, v2, v3, e01, e12, e02, e03, e13, e23 */
	static int edge_order[] = {0, 3, 1, 2, 4, 5};
	a = (output_order > 1 ? 10 : 4);
	BigEndianAppend(&a, sizeof(a), 1);
	/* vertices */
	for (i = 0; i < 4; i++) {
	    a = GlobalVertex(g, s->verts[i]);
	    BigEndianAppend(&a, sizeof(a), 1);
	}
	if (output_order <= 1)
	    continue;
	/* edge centers */
	for (i = 0; i < 6; i++) {
	    a = GlobalEdge(g, s->edges[edge_order[i]]) + g->nvert_global;
	    BigEndianAppend(&a, sizeof(a), 1);
	}
    }
    phgIOWrite(g, buffer, sizeof(a) * (output_order > 1 ? 11 : 5), nleaf,
		g->nelem_global, NULL, NULL, FALSE); 

    phgIOPrintf(g, "\nCELL_TYPES %"dFMT"\n", g->nelem_global);
    buffer_pos = 0;
    a = (output_order > 1 ? /* VTK_QUADRATIC_TETRA = */24 : /*VTK_TETRA = */10);
    for (i = 0; i < nleaf; i++)
	BigEndianAppend(&a, sizeof(a), 1);
    phgIOWrite(g, buffer, sizeof(a), nleaf, g->nelem_global, NULL, NULL, FALSE);

    if (has_point_data) {	/* point data */
	phgIOPrintf(g, "\nPOINT_DATA %"dFMT"\n",
		g->nvert_global + (output_order > 1 ? g->nedge_global : 0));
	for (a = 0; a < dof_ids[0]; a++) {
	    DOF *dof0;
	    dof = dof0 = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	    if (dof->type == DOF_P0 || dof->type == DOF_DG0)
		continue;
	    if (dof->type != DOF_Pn[output_order])
		dof = phgDofCopy(dof0, NULL, DOF_Pn[output_order], dof0->name);
	    if (g->rank == 0)
		phgInfo(2, "writing DOF %d: \"%s\"\n", a, dof->name);
	    dim = dof->dim;
	    if (dim > 4) {
		if (g->rank == 0)
		    phgWarning("only saving the first 4 components of \"%s\"\n",
				dof->name);
		dim = 4;
	    }
	    if (USE_VECTORS && dim == Dim)
	        phgIOPrintf(g, "\nVECTORS %s %s\n", parse_name(g, a),
			    vtk_precision == VTK_DOUBLE ? "double" : "float");
	    else
	        phgIOPrintf(g, "\nSCALARS %s %s %d\nLOOKUP_TABLE default\n",
			    parse_name(g, a),
			    vtk_precision == VTK_DOUBLE ? "double" : "float",
			    dim);
	    /* collect vertex data */
	    buffer_pos = 0;
	    values = DofVertexData(dof, 0);
	    for (i = 0; i < dim * nvert; i++) {
		VTK_PUT_FLOAT(*(values++))
	    }
#if USE_MPI
	    phgIOWrite(g, buffer, vtk_float_size * dim, nvert,
			g->nvert_global, g->types_vert, g->L2Gmap_vert, TRUE); 
#else
	    phgIOWrite(g, buffer, vtk_float_size * dim, nvert,
			g->nvert_global, NULL, NULL, FALSE); 
#endif
	    if (output_order > 1) {
		/* collect edge data */
		buffer_pos = 0;
		values = DofEdgeData(dof, 0);
		for (i = 0; i < dim * nedge; i++) {
		    VTK_PUT_FLOAT(*(values++))
		}
#if USE_MPI
		phgIOWrite(g, buffer, vtk_float_size * dim, nedge,
			g->nedge_global, g->types_edge, g->L2Gmap_edge, TRUE); 
#else
		phgIOWrite(g, buffer, vtk_float_size * dim, nedge,
			g->nedge_global, NULL, NULL, FALSE); 
#endif
	    }
	    if (dof != dof0)
		phgDofFree(&dof);
	}
    }

    phgIOPrintf(g, "\nCELL_DATA %"dFMT"\n", g->nelem_global);
    if (has_cell_data) {
	for (a = 0; a < dof_ids[0]; a++) {
	    dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);

	    if (DofIsHP(dof)) {
		/* variable order DOF, save orders as cell data */
		phgIOPrintf(g, "\nSCALARS %s_orders unsigned_char 1\n"
			   "LOOKUP_TABLE default\n", parse_name(g, a));
		buffer_pos = 0;
		ForAllElements(g, s) {
		    c = dof->hp->elem_order[s->index];
		    BigEndianAppend(&c, sizeof(c), 1);
		}
		phgIOWrite(g, buffer, sizeof(c), nleaf, g->nelem_global,
				NULL, NULL, FALSE); 
		continue;
	    }

	    if (dof->type != DOF_P0 && dof->type != DOF_DG0)
		continue;
	    if (g->rank == 0)
		phgInfo(2, "writing DOF %d: \"%s\"\n", a, dof->name);
	    dim = dof->dim;
	    if (dim > 4) {
		if (g->rank == 0)
		    phgWarning("only saving the first 4 components of \"%s\"\n",
				dof->name);
		dim = 4;
	    }
	    if (USE_VECTORS && dim == Dim)
	 	phgIOPrintf(g, "\nVECTORS %s %s\n", parse_name(g, a),
			    vtk_precision == VTK_DOUBLE ? "double" : "float");
	    else
	 	phgIOPrintf(g, "\nSCALARS %s %s %d\nLOOKUP_TABLE default\n",
			    parse_name(g, a),
			    vtk_precision == VTK_DOUBLE ? "double" : "float",
			    dim);
	    buffer_pos = 0;
	    ForAllElements(g, s) {
		values = DofElementData(dof, s->index);
		for (k = 0; k < dim; k++) {
		    VTK_PUT_FLOAT(*(values++))
		}
	    }
	    phgIOWrite(g, buffer, vtk_float_size * dim, nleaf,
			g->nelem_global, NULL, NULL, FALSE); 
	}
    }

    /* save regional mark as cell data */
    phgIOPrintf(g, "\nSCALARS region_mark int 1\nLOOKUP_TABLE default\n");
    buffer_pos = 0;
    ForAllElements(g, s) {
	a = s->region_mark;
	BigEndianAppend(&a, sizeof(a), 1);
    }
    phgIOWrite(g, buffer, sizeof(a), nleaf, g->nelem_global, NULL, NULL, FALSE);

    /* save submesh indices as cell data */
    phgIOPrintf(g, "\nSCALARS submesh_no int 1\nLOOKUP_TABLE default\n");
    buffer_pos = 0;
    ForAllElements(g, s) {
	a = g->rank;
	BigEndianAppend(&a, sizeof(a), 1);
    }
    phgIOWrite(g, buffer, sizeof(a), nleaf, g->nelem_global, NULL, NULL, FALSE);

    phgIOClose();

    phgFree(buffer);
    buffer = NULL;
#endif/*XJ*/
    return;
}


const char *
phgExportVTKn2D(GRID *g, const char *fn, int ndof, DOF **dofs)
/* Outputs mesh in VTK format, returns the output filename.
 * The variable part of the arguments is a NULL terminated list of DOFs
 * to save to the VTK file */
{
    int i, j, k;
    DOF *dof;
    BTYPE *types_vert = NULL;

    FreeAtExit(dof_ids);

    assert(fn != NULL && Dim == 2);
    if (g == NULL)
	return NULL;
    
    /* output fn */
    filename = fn;
    /* initialize dof_ids */
    dof_ids_len = ndof + 2;
    phgFree(dof_ids);
    dof_ids = phgAlloc(sizeof(*dof_ids) * dof_ids_len);

    /* list of DOFs */
    for (i = 0, k = 0; i < ndof; i++) {
	if ((dof = dofs[i])->type == DOF_ANALYTIC) {
	    if (g->rank == 0) {
		phgWarning("DOF type 'DOF_ANALYTIC' is not supported by "
			   "phgExportVTK(). \n");
		phgWarning("DOF '%s' will not be saved.\n", dof->name);
	    }
	    continue;
	}
	if (dof->g != g)
	    phgError(1, "%s: DOF \"%s\" doesn't belong to the grid.\n",
			__func__, dof->name);
	/* test for missing NULL terminator */
	if (!phgDofIsValid(g, dof)) {
	    if (g->rank == 0)
		phgWarning("phgExportVTK: missing NULL terminator"
			   " in argument list.\n");
	    break;
	}
	if ((j = phgDofPointer2Index(dof)) < 0)
	    break;
	dof_ids[++k] = j;
    }
    dof_ids[0] = k;

    if (g->period != NULL) {
	/* reconstruct types_vert[] without periodicity */
	BYTE	flags;
	PERIOD	*period;
	INT	nvert_owned, *owner_index_vert;
	int	*owner_rank_vert;

	period = g->period;
	g->period = NULL;
	types_vert = g->types_vert;
	g->types_vert = NULL;
	owner_index_vert = g->owner_index_vert;
	g->owner_index_vert = NULL;
	owner_rank_vert = g->owner_rank_vert;
	g->owner_rank_vert = NULL;
	flags = g->flags;
	g->flags = VERT_FLAG;
	nvert_owned = g->nvert_owned;

	phgUpdateBoundaryTypes(g);

	g->period = period;
	phgFree(g->owner_index_vert);
	g->owner_index_vert = owner_index_vert;
	phgFree(g->owner_rank_vert);
	g->owner_rank_vert = owner_rank_vert;
	g->flags = flags;
	g->nvert_owned = nvert_owned;
    }

    phgExportVTK2D_(g);

    if (g->period != NULL) {
	phgFree(g->types_vert);
	g->types_vert = types_vert;
    }
    
    
    return filename;
}



const char *
phgExportVTK2D(GRID *g, const char *fn, DOF *dof, ...)
{
    int ndof;
    DOF **dofs;
    va_list ap;

    assert(g != NULL || phgInitialized == FALSE);

    if (g == NULL) {
	/* register options for VTK export */
	phgOptionsRegisterTitle("\nVTK export options:", "\n", "vtk");
	phgOptionsRegisterKeyword("-vtk_precision", "VTK data precision",
				  vtk_precision_list, &vtk_precision);
	phgOptionsRegisterNoArg("-vtk_nonlinear_elements",
		"Use VTK's non linear elements for high order elements",
		&vtk_nonlinear_elements);

	return NULL;
    }

    dofs = phgAlloc(256 * sizeof(*dofs));

    va_start(ap, dof);
    for (ndof = 0; ndof < 256; ndof++) {
	if (dof == NULL)
	    break;
	dofs[ndof] = dof;
	dof = va_arg(ap, DOF *);
    }
    va_end(ap);

    fn = phgExportVTKn2D(g, fn, ndof, dofs);
    phgFree(dofs);

    return fn;
}
