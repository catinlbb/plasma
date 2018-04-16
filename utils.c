/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: utils.c,v 1.386 2017/03/01 08:21:01 zlb Exp $ */

#define _GNU_SOURCE	/* for fcloseall() */

#include "phg.h"
#include "phg/io.h"

#if USE_MPI
#include "phg/partition-utils.h"
#endif	/* USE_MPI */
#if USE_ZOLTAN
#include "phg/partition-zoltan.h"
#endif	/* USE_ZOLTAN */

#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>	/* usleep */
#include <limits.h>	/* PATH_MAX */
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include <sys/time.h>
#ifdef __WIN__
#elif defined(USETIMES)
# include <sys/times.h>
# include <unistd.h>
#else
# include <sys/resource.h>
#endif
#define TRUNK_SIZE 8

#if HAVE_FEENABLEEXCEPT
# include <fenv.h>
#endif	/* HAVE_FEENABLEEXCEPT */
#include <signal.h>

FLOAT _phg_round_tmp;	/* scratch variable used by the 'Round_' macro */

/* list of valid grids */
static size_t g_list_count = 0;
static GRID **g_list = NULL;

int phgMaxThreads = 1;
int phgThreadId = 0;
INT phgVerbosity = 0;
const char *phgCurrentFunctionName_ = "PHG";
INT _phg_i = -1;	/* scratch INT variable used by ForAllElements */

/* function used by ForAllElements macro */
ELEMENT *
_phg_step_element(GRID *g)
{
    ELEMENT *e = NULL;

    while (_phg_i < g->nelem && (e = g->elems[_phg_i]) == NULL)
	_phg_i++;

    if (e == NULL)
	_phg_i = -1;

    return e;
}

static GRID *g_traverse;	/* used by phgTraverseElements, other functions
				   shouldn't use it */
static GRID *g_;		/* for traversal callback functions */
BOOLEAN phgInitialized = FALSE;
BOOLEAN phgMasterSlave = FALSE;

char *phgLogFilename = NULL;
char *phgOutFilename = NULL;
static FILE *log_file = NULL;
static FILE *out_file = NULL;

#define LogFile	(log_file == NULL ? stderr : log_file)
#define OutFile	(out_file == NULL ? stdout : out_file)

static void
init_log_file(void)
{
    char s[PATH_MAX];
    int rank = 0, nprocs = 1;

#if USE_MPI
    MPI_Initialized(&rank);
    if (!rank)
	return;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif	/* USE_MPI */

    if (phgLogFilename == NULL || log_file != NULL)
	return;

    if (strcmp(phgLogFilename, "stdout") == 0) {
	log_file = stdout;
	return;
    }
    else if (strcmp(phgLogFilename, "stderr") == 0) {
	log_file = stderr;
	return;
    }

    if (nprocs > 1)
	sprintf(s, "%s.%04d", phgLogFilename, rank);
    else
	strcpy(s, phgLogFilename);
    
    log_file = fopen(s, "w+t");
    if (log_file == NULL) {
	phgFree(phgLogFilename);
	phgLogFilename = NULL;
    }
}

static void
close_log_file(void)
{
    int rank = 0, nprocs = 1;

#if USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif	/* USE_MPI */

    if (log_file != NULL && log_file != stdout && log_file != stderr) {
	fclose(log_file);
	if (rank == 0) {
	    if (nprocs > 1)
		phgPrintf("Log files \"%s.####\" created.\n", phgLogFilename);
	    else
		phgPrintf("Log file \"%s\" created.\n", phgLogFilename);
	}
    }
    log_file = NULL;

    fflush(OutFile);
    if (out_file != NULL && out_file != stdout && out_file != stderr) {
	fclose(out_file);
    }
    out_file = NULL;
}

int
phgCompFLOAT(const void *p0, const void *p1)
/* compares two 'FLOAT's (used in qsort and bsearch) */
{
    FLOAT d = *((FLOAT *)p0) - *((FLOAT *)p1);

    return (d < 0.0) ? -1 : ((d == 0.0) ? 0 : 1);
}

int
phgCompINT(const void *p0, const void *p1)
/* compares two 'INT's (used in qsort and bsearch) */
{
    INT i = *((INT *)p0) - *((INT *)p1);

    return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
}

int
phgCompint(const void *p0, const void *p1)
/* compares two 'int's (used in qsort and bsearch) */
{
    return *((int *)p0) - *((int *)p1);
}

INT
phgBinarySearch(INT n, const void *table, const INT *ordering, const void *key,
		size_t size, int (*comp)(const void *, const void *))
/* performs binary search in table[ordering[]] which is a sorted table of
 * size n, returns i in [0, n) such that (assume table[ordering[-1]] == -infty)
 * 	(key == table[ordering[i]] ||
 * 	 (key < table[ordering[i]] && key > table[ordering[i - 1]])),
 * or n if key > table[ordering[n - 1]]
 *
 * Note:
 * 1. comp() is always called with the first argument pointing to table[]
 *    and the second argument pointing to key.
 *
 * 2. ordering may be a NULL pointer, in which case ordering[i] == i is assumed.
 */
{
    INT a, b, c;
    int i;

    if (n <= 0)
	return n;

    if (ordering != NULL) {
	if (comp((const char *)table + size * ordering[a = 0], key) >= 0)
	    return 0;
	if (n == 1)
	    return 1;
	if ((i = comp((const char *)table + size * ordering[b = n - 1],
			key)) < 0)
	    return n;
	if (i == 0)
	    return n - 1;
	while (b > a + 1) {
	    if ((i = comp((const char *)table +
				size * ordering[c = (a + b) / 2], key)) == 0)
		return c;
	    else if (i > 0)
		b = c;
	    else
		a = c;
	}
    }
    else {
	if (comp(table, key) >= 0)
	    return 0;
	a = 0;
	if (n == 1)
	    return 1;
	if ((i = comp((const char *)table + size * (b = n - 1), key)) < 0)
	    return n;
	if (i == 0)
	    return n - 1;
	while (b > a + 1) {
	    if ((i = comp((const char *)table + size * (c = (a + b) / 2),
				key)) == 0)
		return c;
	    else if (i > 0)
		b = c;
	    else
		a = c;
	}
    }

    return b;
}

#ifdef phgBinarySearchINT
# undef phgBinarySearchINT
#endif	/* defined(phgBinarySearchINT) */
INT
phgBinarySearchINT(INT size, const INT *list, const INT *ordering, INT index)
/* performs binary search in list[ordering[]] which is a sorted table of
 * size size, returns i in the range [0, size) such that
 * 	(index == list[ordering[i]] ||
 * 	 (index < list[ordering[i]] && index > list[ordering[i - 1]])),
 * (assume table[ordering[-1]] == -infty), or n if key > table[ordering[n - 1]]
 *
 * This function is approximately equivalent to:
 * 	for (i = 0; i < size && list[i] < index; i++); return i;
 *
 * Note: This function may be defined as a macro using phgBinarySearch above,
 *	 but since it's a frequently used one we leave it here as a separate
 *	 function for better performance (see utils.h). */
{
    INT a, b, c, i;

    if (size <= 0)
	return size;

    if (ordering != NULL) {
	if (index <= list[ordering[a = 0]])
	    return 0;
	if (index >= (i = list[ordering[b = size - 1]]))
	    return index == i ? size - 1 : size;
	while (b > a + 1) {
	    if (index == (i = list[ordering[c = (a + b) / 2]]))
		return c;
	    else if (index < i)
		b = c;
	    else
		a = c;
	}
    }
    else {
	if (index <= list[a = 0])
	    return 0;
	if (index >= (i = list[b = size - 1]))
	    return index == i ? size - 1 : size;
	while (b > a + 1) {
	    if (index == (i = list[c = (a + b) / 2]))
		return c;
	    else if (index < i)
		b = c;
	    else
		a = c;
	}
    }

    return b;
}

#if HAVE_FEENABLEEXCEPT
static void
fpe_handler(int signum, siginfo_t *siginfo, void *dummy)
{
    switch (siginfo->si_code) {
	case FPE_FLTDIV:
	    phgError(1, "floating point exception: divide by zero.\n");
	    break;
	case FPE_FLTOVF:
	    phgError(1, "floating point exception: overflow.\n");
	    break;
	case FPE_FLTUND:
	    phgError(1, "floating point exception: underflow.\n");
	    break;
	case FPE_FLTRES:
	    phgError(1, "floating point exception: inexact result.\n");
	    break;
	case FPE_FLTINV:
	    phgError(1, "floating point exception: invalid operation.\n");
	    break;
	default:
	    phgError(1, "program received signal SIGFPE, abort.\n");
    }
}
#endif	/* HAVE_FEENABLEEXCEPT */

void
phgTrapSignals_(BOOLEAN init)
/* Note: should be called once with 'init == TRUE' once before phgInit,
 * and once with 'init == FALSE' to set FPE trap after initialization
 * of all external packages (like PETSc) */
{
#ifdef __WIN__
    Unused(init);
#else	/* defined(__WIN__) */
#if HAVE_FEENABLEEXCEPT
    int flags = 0;
#endif	/* HAVE_FEENABLEEXCEPT */
    static const char *fpetrap_keywords[] = {"none", "fatal", "all", NULL};
    static int fpetrap = 1;
    static BOOLEAN reset_signals = TRUE;
    struct sigaction sa;

    if (init) {
	/* Note: the function is called once with 'init == TRUE' for
	 * registering options, then once with 'init == FALSE' after
	 * initialization of external solvers, to make sure the sigactions
	 * are called after external packages */
	phgOptionsRegisterNoArg("reset_signals", "Reset signal handlers",
		&reset_signals);
#if HAVE_FEENABLEEXCEPT
	phgOptionsRegisterKeyword("fpetrap", "Trap floating point exceptions",
		fpetrap_keywords, &fpetrap);
#else	/* HAVE_FEENABLEEXCEPT */
	phgOptionsRegisterKeyword("fpetrap",
		"This option is ignored on this system",
		fpetrap_keywords, &fpetrap);
#endif	/* HAVE_FEENABLEEXCEPT */
	return;
    }

    if (reset_signals) {
	/* reset some signal handlers to default */
	sa.sa_handler = SIG_DFL;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = 0;
	if (sigaction(SIGPIPE, &sa, NULL))
	    if (phgVerbosity > 0)
		perror("sigaction");
	if (sigaction(SIGINT, &sa, NULL))
	    if (phgVerbosity > 0)
		perror("sigaction");
	if (sigaction(SIGQUIT, &sa, NULL))
	    if (phgVerbosity > 0)
		perror("sigaction");
	if (sigaction(SIGILL, &sa, NULL))
	    if (phgVerbosity > 0)
		perror("sigaction");
	if (sigaction(SIGSEGV, &sa, NULL))
	    if (phgVerbosity > 0)
		perror("sigaction");
	if (sigaction(SIGTERM, &sa, NULL))
	    if (phgVerbosity > 0)
		perror("sigaction");
    }

#if HAVE_FEENABLEEXCEPT
    switch (fpetrap) {
	case 2:	/* all */
	    flags = FE_UNDERFLOW /*| FE_INEXACT*/;
	    /* fall through the next case */
	case 1:	/* fatal */
	    flags |= FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW;
	    /* Clear current exceptions queue.
	     * Note: the initialization code of the Anasazi package from
	     * Trilinos-7.0.5 generates 'invalid operation' traps on i386 */
	    feclearexcept(flags);
	    if (feenableexcept(flags) == -1) {
		if (phgVerbosity > 0)
		    perror("");
		exit(1);
	    }
	    break;
	case 0:	/* none */
	    return;
    }

    /* install SIGFPE handler */
    sa.sa_sigaction = fpe_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_SIGINFO;
    if (sigaction(SIGFPE, &sa, NULL))
	if (phgVerbosity > 0)
	    perror("sigaction");
#endif	/* HAVE_FEENABLEEXCEPT */
#endif	/* defined(__WIN__) */
}

BOOLEAN _phg_pause = FALSE;

void
phgInit(int *argc, char ***argv)
{
    if (phgInitialized)
	return;

    /* test validity of basic datatypes */
    assert(sizeof(int) >= 4);
    assert(sizeof(unsigned int) >= 4);
    assert(sizeof(short) >= 2);
    assert(sizeof(unsigned short) >= 2);

#if USE_OMP
    /* set default # threads to 1 if OMP_NUM_THREADS is undefined or empty */
    {
	const char *t = getenv("OMP_NUM_THREADS");
	if (t != NULL)
	    while (isspace((int)(*t))) t++;
	if (t == NULL || *t == '\0') {
	    putenv("OMP_NUM_THREADS=1");
	    omp_set_num_threads(1);
	}
	omp_set_dynamic(0);
	omp_set_nested(0);
	phgMaxThreads = omp_get_max_threads();
#pragma omp parallel
	phgThreadId = omp_get_thread_num();
    }
#endif	/* USE_OMP */

    /* Initialize: register internal options */
    phgOptionsRegisterInit_();

    /* Register partitioner options */
    phgRedistributeGrid(NULL);
    /* Register marking options */
    phgMarkElements(MARK_DEFAULT, NULL, 0., NULL, 0., 0,
		    MARK_DEFAULT, NULL, 0., 0, 0.);
    /* Register solver options */
    phgSolverRegisterOptions();

    /* Initialize pre-defined DOF_TYPEs */
    phgDofTypesInit();

#if USE_MPI
    phgInitMPI(argc, argv);
#else	/* USE_MPI */
    phgPrintf("Parallel Hierarchical Grid (version " PHG_VERSION ").\n");
    phgOptionsParseCmdline(argc, argv);

    /* Initialize modules after parsing all options.
     * Note: also add the calls below in phgInitMPI() in phg-mpi.c */
    phgSolverInitialize(argc, argv);	/* init phgSolver */
#if USE_ZOLTAN
    phgPartitionZoltanInitialize(*argc, *argv);
#endif	/* USE_ZOLTAN */

    phgOptionsHelp();

    if (phgVerbosity >= 2) {
	phgPrintf("sizeof(INT)        = %d\n", sizeof(INT));
	phgPrintf("sizeof(FLOAT)      = %d\n", sizeof(FLOAT));
	phgPrintf("sizeof(ELEMENT)    = %d\n", sizeof(ELEMENT));
	phgPrintf("sizeof(GRID)       = %d\n", sizeof(GRID));
    }
    fflush(OutFile);

    phgTrapSignals_(FALSE);
    phgPerfInit();
#endif	/* USE_MPI */

    if (phgVerbosity > 0)
	phgPrintf("phgMaxThreads = %d\n", phgMaxThreads);

    phgInitialized = TRUE;
    if (_phg_pause)
	phgPause(0);
    if (phgMasterSlave)
	atexit(phgFinalize);
}

#if USE_MPI
void
phgSetVerbosity_(GRID *g)
{
    ParallelFunction(phgSetVerbosity_, TRUE);

    if (phgNProcs <= 0)
	return;

    MPI_Bcast(&phgVerbosity, 1, MPI_INT, 0, phgComm);

    return;
}
#endif

int
phgSetVerbosity(int verbosity)
{
    int old = phgVerbosity;
    if (verbosity >= 0) phgVerbosity = verbosity;
#if USE_MPI
    phgSetVerbosity_(NULL);
#endif
    return old;
}

void
phgFinalize(void)
{
    if (!phgInitialized)
	return;

    /* Free all grids (to clean up memory blocks) */
    while (g_list != NULL) {
	GRID *g = *g_list;
	phgFreeGrid(&g);
    }

    phgPerfFinalize();
    phgSolverFinalize();
    phgQuadReset();

#if USE_MPI
    if (!phgMasterSlave)
	phgCreateComm_(-1, NULL);	/* tell idle processes to exit */
#endif	/* USE_MPI */

    close_log_file();

    phgOptionsReset();

    phgInitialized = FALSE;

#if USE_MPI
    phgFinalizeMPI();
#endif
 
    phgFreeAtExit(NULL, 0);

    if (phgRank > 0)
	exit(0);
}

void
phgPause(int seconds)
{
    char s[128];

    if (seconds <= 0) {
	int i;
	pid_t pid = getpid();
	FLOAT speed = phgGetProcessorSpeed(NULL);
#if USE_MPI
	char name[MPI_MAX_PROCESSOR_NAME + sizeof(pid_t) + sizeof(speed)],
	     (*names)[MPI_MAX_PROCESSOR_NAME + sizeof(pid_t) + sizeof(speed)];
	int len;
	MPI_Get_processor_name(name, &len);
	len = MPI_MAX_PROCESSOR_NAME;
	memcpy(name + len, &pid, sizeof(pid));
	len += sizeof(pid);
	memcpy(name + len, &speed, sizeof(speed));
	len += sizeof(speed);
	names = phgAlloc(phgNProcs * len);
	MPI_Gather(name, len, MPI_BYTE, names, len, MPI_BYTE, 0, phgComm);
#endif	/* USE_MPI */
	if (phgRank == 0) {
            fprintf(stderr, "Program paused.\n");
	    for (i = 0; i < phgNProcs; i++) {
		char *p;
#if USE_MPI
		len = MPI_MAX_PROCESSOR_NAME;
		p = names[i];
		memcpy(&pid, p + len, sizeof(pid));
		len += sizeof(pid);
		memcpy(&speed, p + len, sizeof(speed));
#else	/* USE_MPI */
		p = "localhost";
#endif	/* USE_MPI */
		fprintf(stderr, "  proc %d (%s), pid = %d, speed = %lg\n",
				i, p, (int)pid, (double)speed);
	    }
	}
#if USE_MPI
	phgFree(names);
#endif	/* USE_MPI */
    }
    if (phgRank == 0) {
	if (seconds <= 0) {
            fflush(stderr);
            fprintf(stderr, "<<< Press enter to continue >>>\n");
            fflush(stderr);
            if (fgets(s, sizeof(s), stdin) == NULL)
		phgError(1, "unexpected error.\n");
            fflush(stderr);
	}
	else {
	    sleep(seconds);
	}
    }
#if USE_MPI
    MPI_Barrier(phgComm);
#endif
}

void
phgInitElements(ELEMENT *s, INT count)
{
    INT i;
    int j;

    if (count == 0)
	return;
    memset(s, 0, count * sizeof(*s));

    for (i = 0; i < count; i++) {
	for (j = 0; j < NEdge; j++) {
	    s[i].edges[j] = -1;
#if ALLOW_CURVED_BOUNDARY
	    s[i].bound_func[j] = -1;
#endif
	}
	for (j = 0; j < NFace; j++) {
	    s[i].bound_type[j] = UNDEFINED;
#if (Dim == 3)
	    s[i].faces[j] = -1;
#endif
	}
	s[i].index = -1;
	s[i].region_mark = 0;
    }
}

ELEMENT *
phgNewElements(INT count)
{
    ELEMENT *s = phgAlloc(count * sizeof(ELEMENT));

    phgInitElements(s, count);
    return s;
}

ELEMENT *
phgReallocElements(ELEMENT *s, INT oldcount, INT newcount)
{
    if (newcount == oldcount)
	return s;
    s = phgRealloc_(s, newcount * sizeof(ELEMENT), oldcount * sizeof(ELEMENT));
    if (newcount <= oldcount)
	return s;
    phgInitElements(s + oldcount, newcount - oldcount);

    return s;
}

COORD *
phgNewVertices(INT count)
{
    COORD *v;

    v = phgAlloc(count * sizeof(COORD));
    return v;
}

static int flags_tmp;
static GRID *g_new;

void
phgNewGrid_(GRID *g)
{
    ParallelFunction(phgNewGrid_, TRUE);

#if USE_MPI
    if (phgMasterSlave && phgNProcs > 0) {
	MPI_Bcast(&flags_tmp, sizeof(flags_tmp), MPI_BYTE, 0, phgComm);
    }
#endif

    /* for Master/Slave mode in which only one grid is allowed */
    g_new = g;
    if (g != NULL)
	return;

    g_new = g = phgCalloc(1, sizeof(GRID));
    g->lif = 1.0;
    g->nelem = g->nleaf = g->nroot = g->nvert = g->nedge = g->ntree = 0;
#if Dim == 3
    g->nface = 0;
    g->nface_global = 0;
#endif
    g->verts = NULL;
    g->roots = NULL;
    g->elems = NULL;
    g->dof   = NULL;
    g->flags = flags_tmp;
    g->flags |= VERT_FLAG;	/* always needed */
    g->flags |= ELEM_FLAG;	/* always needed */
#if USE_ZOLTAN
    g->flags |= FACE_FLAG;	/* needed by Zoltan */
#endif	/* USE_ZOLTAN */

    assert(sizeof(int) >= sizeof(BTYPE));
    g->bc_n = g->bc_alloc = 0;
    g->bc_list = NULL;
    g->bc_rmap = NULL;

#if ALLOW_CURVED_BOUNDARY
    g->bdry_funcs = NULL;
#endif	/* ALLOW_CURVED_BOUNDARY */
    g->types_vert = NULL;
    g->types_edge = NULL;
    g->types_face = NULL;
    g->types_elem = NULL;
    g->alien = NULL;
    g->alien_map = NULL;
    g->period = NULL;
    g->nelem_global = g->nvert_global = g->nedge_global = 0;
#if USE_MPI
    g->L2Gmap_vert = NULL;
    g->L2Gmap_edge = NULL;
    g->L2Gmap_face = NULL;
    g->L2Gmap_elem = NULL;
#if NO_NEW_COMMUNICATOR
    g->comm = MPI_COMM_SELF;
#else	/* NO_NEW_COMMUNICATOR */
    MPI_Comm_dup(MPI_COMM_SELF, &g->comm);
#endif	/* NO_NEW_COMMUNICATOR */

    g->neighbours.counts = g->neighbours.displs = NULL;
    g->neighbours.list = NULL;
    g->neighbours.count = g->neighbours.allocated = 0;
#else	/* USE_MPI */
    g->comm = phgComm;
#endif	/* USE_MPI */
    g->rank = 0;
    g->nprocs = 1;
    g->serial_no = 0;

    g->last_partitioner = -1;

    /* add g to g_list (for phgGridIsValid) */
    g_list = phgRealloc_(g_list, (g_list_count + 1) * sizeof(*g_list),
				 g_list_count * sizeof(*g_list));
    g_list[g_list_count++] = g;

    if (_phg_g_sys == NULL)
	_phg_g_sys = g;

    return;
}

GRID *
phgNewGrid(int flags)
{
    flags_tmp = flags;
    phgNewGrid_(NULL);

    return g_new;
}

int
phgGridId(const GRID *g)
/* returns id of the given grid, -1 if grid is invalid */
{
    int i;

    for (i = 0; i < g_list_count; i++)
	if (g_list[i] == g)
	    return i;

    return -1;
}

static INT index_ = 0;

static BOOLEAN
elems_callback CB_ARGS(e)
{
    ELEMENT *e0, *e1;
    int i, j, v[NVert];
    INT V[NVert];

    e0 = e->children[0];
    e1 = e->children[1];

    assert(e->index >= 0 && e->index < g_->nelem);

    if (e0 == NULL && e1 == NULL) {
	g_->elems[e->index] = e;
    }
    else {
	/* point e->parent to parent */
	if (e0 != NULL)
	    e0->parent = e;
	if (e1 != NULL)
	    e1->parent = e;
    }
    if (e->generation == 0)
	e->parent = NULL;

    if (e->ordering != 0)
	return TRUE;

    /* compute ordering of vertices with increasing global indices */
    for (i = 0; i < NVert; i++) {
	V[i] = GlobalVertexP(g_, e->verts[i]);
	v[i] = 0;
    }
    for (i = 0; i < NVert; i++) {
	for (j = i+1; j < NVert; j++) {
            if (V[i] > V[j])
		v[i]++;
            else
		v[j]++;
        }
    }
    for (i = 0; i < NVert; i++) 
	e->ordering |= v[i] << (i * 3);

    return TRUE;
}

static BOOLEAN
types_callback CB_ARGS(e)
{
    int i, type;

    if ((g_->flags & ELEM_FLAG))
	g_->types_elem[e->index] |= OWNER;

    for (i = 0; i < NFace; i++) {
	type = e->bound_type[i];
#if 0
	/* use DIRICHLET for UNDEFINED boundary type */
	if (type == UNDEFINED)
	    type = DIRICHLET;
#endif
	type |= OWNER;

	if ((g_->flags & VERT_FLAG)) {
	    switch (i) {
		case 0:
		    g_->types_vert[e->verts[1]] |= type;
		    g_->types_vert[e->verts[2]] |= type;
		    g_->types_vert[e->verts[3]] |= type;
		    break;
		case 1:
		    g_->types_vert[e->verts[0]] |= type;
		    g_->types_vert[e->verts[2]] |= type;
		    g_->types_vert[e->verts[3]] |= type;
		    break;
		case 2:
		    g_->types_vert[e->verts[0]] |= type;
		    g_->types_vert[e->verts[1]] |= type;
		    g_->types_vert[e->verts[3]] |= type;
		    break;
		case 3:
		    g_->types_vert[e->verts[0]] |= type;
		    g_->types_vert[e->verts[1]] |= type;
		    g_->types_vert[e->verts[2]] |= type;
		    break;
	    }
	}

	if ((g_->flags & EDGE_FLAG)) {
	    switch (i) {
		case 0:
		    g_->types_edge[e->edges[3]] |= type;
		    g_->types_edge[e->edges[4]] |= type;
		    g_->types_edge[e->edges[5]] |= type;
		    break;
		case 1:
		    g_->types_edge[e->edges[1]] |= type;
		    g_->types_edge[e->edges[2]] |= type;
		    g_->types_edge[e->edges[5]] |= type;
		    break;
		case 2:
		    g_->types_edge[e->edges[0]] |= type;
		    g_->types_edge[e->edges[2]] |= type;
		    g_->types_edge[e->edges[4]] |= type;
		    break;
		case 3:
		    g_->types_edge[e->edges[0]] |= type;
		    g_->types_edge[e->edges[1]] |= type;
		    g_->types_edge[e->edges[3]] |= type;
		    break;
	    }
	}

	if ((g_->flags & FACE_FLAG)) {
	    g_->types_face[e->faces[i]] |= type;
	}

	if ((g_->flags & ELEM_FLAG)) {
	    g_->types_elem[e->index] |= (type & REMOTE);
	}
    }

    return TRUE;
}

#if USE_MPI
typedef struct {
    INT		index;	/* global index */ 
    INT		lindex;	/* local index */ 
    int		rank;	/* process owning the entry */
    GTYPE	gtype;	/* whether a vertex, edge, face or element */
    BTYPE	btype;	/* boundary type flags */
} B_LIST;

static B_LIST *blist;

static int
blist_comp_indices(const void *p0, const void *p1)
{
    int ii;
    INT i;
    B_LIST *b0 = blist + *((const INT *)p0), *b1 = blist + *((const INT *)p1);

    if ((ii = b0->gtype - b1->gtype) != 0)
	return ii;

    return (i = b0->index - b1->index) > 0 ? 1 : (i < 0 ? -1 : 0);
}
#endif

static void set_default_bdry_type(GRID *g);

static int
comp_pvertex(const void *p0, const void *p1)
{
    return g_->period->L2Gmap_vert[*(const INT *)p0] - 
	   g_->period->L2Gmap_vert[*(const INT *)p1];
}

void
phgUpdateBoundaryTypes(GRID *g)
/* Sets boundary types (g->types_xxxx), updates g->elems, g->owner_xxxx_xxxx
 * and the peer_face member of g->neighbours. */
{
    ELEMENT *e;
    INT i, j, nd, nm, n0;
#if USE_MPI
    BTYPE b;
    int ii, jj, n, rank, *scnts, *sdsps, *rcnts, *rdsps;
    INT nsend, nrecv, *rind;		/* sorted indices of brecv[] */
    MPI_Datatype type;
    B_LIST *bsend, *brecv, *bl = NULL;
#endif
    double time0 = phgGetTime(NULL);

    /*assert((g->flags & ELEM_FLAG))*/; /* elems_callback requires e->index */

    phgFree(g->elems);
    g->elems = NULL;

    if (g == NULL || g->comm == MPI_COMM_NULL)
	return;

    set_default_bdry_type(g);

    if ((g->flags & VERT_FLAG)) {
	phgFree(g->types_vert);
	g->types_vert = phgCalloc(g->nvert, sizeof(*g->types_vert));
    }

    if ((g->flags & EDGE_FLAG)) {
	phgFree(g->types_edge);
	g->types_edge = phgCalloc(g->nedge, sizeof(*(g->types_edge)));
    }

    if ((g->flags & FACE_FLAG)) {
	phgFree(g->types_face);
	g->types_face = phgCalloc(g->nface, sizeof(*(g->types_face)));
    }

    if ((g->flags & ELEM_FLAG)) {
	phgFree(g->types_elem);
	g->types_elem = phgCalloc(g->nelem, sizeof(*(g->types_elem)));
    }

    g_ = g;
    g->elems = phgCalloc(g->nelem, sizeof(*g->elems));
    index_ = 0;
    /* construct g->elems, update e->parent and e->odering */
    phgTraverseAllElements(g, elems_callback);

    if (!(g->flags & (VERT_FLAG | EDGE_FLAG | FACE_FLAG | ELEM_FLAG)))
	return;

    if ((g->flags & VERT_FLAG))
	g->nvert_owned = g->nvert;
    if ((g->flags & EDGE_FLAG))
	g->nedge_owned = g->nedge;
    if ((g->flags & FACE_FLAG))
	g->nface_owned = g->nface;
    if ((g->flags & ELEM_FLAG))
	g->nelem_owned = g->nelem;

    ForAllElements(g, e)
	types_callback(e);

    if (g->period != NULL) {
	/* construct g->period->ordering of period->L2Gmap_vert */
	phgFree(g->period->ordering);
	g->period->ordering = phgAlloc(g->nvert * sizeof(*g->period->ordering));
	for (i = 0; i < g->nvert; i++)
	    g->period->ordering[i] = i;
	qsort(g->period->ordering, g->nvert,
			sizeof(*g->period->ordering), comp_pvertex);
	/* set REMOTE bit on periodic vertices */
	if ((g->flags & VERT_FLAG)) {
	    for (i = 1; i < g->nvert; i++) {
		if (g->period->L2Gmap_vert[n0 = g->period->ordering[i]] !=
		    g->period->L2Gmap_vert[nm = g->period->ordering[i - 1]])
		    continue;
		if (g->types_vert[n0] != UNREFERENCED)
		    g->types_vert[n0] |= REMOTE;
		if (g->types_vert[nm] != UNREFERENCED)
		    g->types_vert[nm] |= REMOTE;
	    }
	}
    }

    if (g->nprocs == 1) {
	if (g->period == NULL || !(g->flags & VERT_FLAG))
	    return;

	/* update OWNER bit in types_vert[] and set up owner_xxxxx_vert[] */
	phgFree(g->owner_rank_vert);
	g->owner_rank_vert = phgAlloc(g->nvert * sizeof(*g->owner_rank_vert));
	phgFree(g->owner_index_vert);
	g->owner_index_vert = phgAlloc(g->nvert * sizeof(*g->owner_index_vert));
	g->nvert_owned = 0;
	nd = n0 = -1;
	for (i = 0; i < g->nvert; i++) {
	    g->owner_rank_vert[nm = g->period->ordering[i]] = g->rank;
	    if (nd != (j = g->period->L2Gmap_vert[nm])) {
		nd = j;
		n0 = nm;
		g->nvert_owned++;
	    }
	    else {
		g->types_vert[nm] &= ~OWNER;
	    }
	    g->owner_index_vert[nm] = n0;
	}
	return;
    }

#if USE_MPI
    if ((g->flags & VERT_FLAG)) {
	phgFree(g->owner_rank_vert);
	g->owner_rank_vert = phgAlloc(g->nvert * sizeof(*g->owner_rank_vert));
	phgFree(g->owner_index_vert);
	g->owner_index_vert = phgAlloc(g->nvert * sizeof(*g->owner_index_vert));
	for (i = 0; i < g->nvert; i++) {
	    g->owner_rank_vert[i] = -1;
	    g->owner_index_vert[i] = -1;
	}
    }

    if ((g->flags & EDGE_FLAG)) {
	phgFree(g->owner_rank_edge);
	g->owner_rank_edge = phgAlloc(g->nedge * sizeof(*g->owner_rank_edge));
	phgFree(g->owner_index_edge);
	g->owner_index_edge = phgAlloc(g->nedge * sizeof(*g->owner_index_edge));
	for (i = 0; i < g->nedge; i++) {
	    g->owner_rank_edge[i] = -1;
	    g->owner_index_edge[i] = -1;
	}
    }

    /* set up the OWNER bit in types_xxxx[] and other bits such that
     * they are consistent across all processes, update owner_xxxx_xxxx
     * arrays and the peer_face member in g->neighbours */

    /* Algorithm for updating the types_xxxx array:
     *    a vertex or edge is blockwisely (or cyclically?) assigned to a
     *    process using its global index. All vertices and edges with REMOTE
     *    flag set are sent to assigned processes, processed there, and then
     *    sent back with updated data. */

#define IsRemote(t)    ((b = (t)) & REMOTE)

    /* count entries which have REMOTE bit set */
    scnts = phgCalloc(4 * g->nprocs, sizeof(*scnts));
    sdsps = scnts + g->nprocs;
    rcnts = sdsps + g->nprocs;
    rdsps = rcnts + g->nprocs;
    if ((g->flags & VERT_FLAG)) {
	if (g->period != NULL) {
	    nd = g->period->nvert_global / g->nprocs;
	    nm = g->period->nvert_global % g->nprocs;
	}
	else {
	    nd = g->nvert_global / g->nprocs;
	    nm = g->nvert_global % g->nprocs;
	}
	n0 = nm * (nd + 1);
	for (i = 0; i < g->nvert; i++) {
	    if (IsRemote(g->types_vert[i])) {
		j = GlobalVertexP(g, i);
		rank = (j < n0 ? j / (nd + 1) : nm + (j - n0) / nd);
		scnts[rank]++;
	    }
	    else if (b == UNREFERENCED) {
		g->nvert_owned--;
	    }
	}
    }

    if ((g->flags & EDGE_FLAG)) {
	nd = g->nedge_global / g->nprocs;
	nm = g->nedge_global % g->nprocs;
	n0 = nm * (nd + 1);
	for (i = 0; i < g->nedge; i++) {
	    if (IsRemote(g->types_edge[i])) {
		j = GlobalEdge(g, i);
		rank = (j < n0 ? j / (nd + 1) : nm + (j - n0) / nd);
		scnts[rank]++;
	    }
	    else if (b == UNREFERENCED) {
		g->nedge_owned--;
	    }
	}
    }

    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, g->comm);

    nsend = nrecv = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	sdsps[rank] = nsend;
	rdsps[rank] = nrecv;
	assert(nsend <= INT_MAX - scnts[rank]);
	nsend += scnts[rank];
	assert(nrecv <= INT_MAX - rcnts[rank]);
	nrecv += rcnts[rank];
    }
    phgInfo(4, "nsend = %"dFMT", nrecv = %"dFMT", ntotal = %"dFMT"\n",
		nsend, nrecv, g->nvert + g->nedge + g->nface + g->nelem);

    /* collect entries which have REMOTE bit */
    bsend = phgAlloc((nsend + nrecv) * sizeof(*bsend));
    brecv = bsend + nsend;

#if DEBUG
    for (ii = 0; ii < nsend; ii++)
	bsend[ii].gtype = MIXED;
#endif	/* DEBUG */

    if ((g->flags & VERT_FLAG)) {
	if (g->period != NULL) {
	    nd = g->period->nvert_global / g->nprocs;
	    nm = g->period->nvert_global % g->nprocs;
	}
	else {
	    nd = g->nvert_global / g->nprocs;
	    nm = g->nvert_global % g->nprocs;
	}
	n0 = nm * (nd + 1);
	for (i = 0; i < g->nvert; i++) {
	    if (IsRemote(g->types_vert[i])) {
		j = GlobalVertexP(g, i);
		rank = (j < n0 ? j / (nd + 1) : nm + (j - n0) / nd);
		bl = bsend + (sdsps[rank]++);
		bl->index  = j;
		bl->lindex = i;
		bl->rank   = g->rank;
		bl->gtype  = VERTEX;
		bl->btype  = b & ~(/*UNDEFINED |*/ INTERIOR | OWNER);
		phgInfo(4, "%s vertex %"dFMT", sent to proc %d\n", 
				b==UNREFERENCED ? "unref" : "shared", j, rank);
	    }
	}
    }

    if ((g->flags & EDGE_FLAG)) {
	nd = g->nedge_global / g->nprocs;
	nm = g->nedge_global % g->nprocs;
	n0 = nm * (nd + 1);
	for (i = 0; i < g->nedge; i++) {
	    if (IsRemote(g->types_edge[i])) {
		j = GlobalEdge(g, i);
		rank = (j < n0 ? j / (nd + 1) : nm + (j - n0) / nd);
		bl = bsend + (sdsps[rank]++);
		bl->index  = j;
		bl->lindex = i;
		bl->rank   = g->rank;
		bl->gtype  = EDGE;
		bl->btype  = b & ~(/*UNDEFINED |*/ INTERIOR | OWNER);
	    }
	}
    }

    assert(sdsps[g->nprocs - 1] == nsend);

    /* restore sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
	sdsps[rank] -= scnts[rank];

    /* send lists */
    MPI_Type_contiguous(sizeof(*bsend), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(bsend, scnts, sdsps, type, brecv, rcnts, rdsps, type,
		  g->comm);

    /* process received lists */
    rind = phgAlloc(nrecv * sizeof(*rind));
    for (ii = 0; ii < nrecv; ii++)
	rind[ii] = ii;
    blist = brecv;
    qsort(rind, nrecv, sizeof(*rind), blist_comp_indices);
    for (ii = 0; ii < nrecv; ) {
	for (n = ii + 1;
	     n < nrecv && blist_comp_indices(rind + ii, rind + n) == 0;
	     n++);
	/* the entries in the range [ii, n) correspond to the same entity */
	rank = -1;		/* the owner process rank */
	n0 = -1;		/* the local index in the owner process */
#if SIZEOF_PHG_INT >= 4
	i = 99999999;		/* hash value */
#else	/* SIZEOF_PHG_INT >= 4 */
	i = 9999;		/* hash value */
#endif	/* SIZEOF_PHG_INT >= 4 */
	b = 0;			/* BTYPE */
	for (jj = ii; jj < n; jj++) {
	    bl = brecv + rind[jj];
	    b |= bl->btype;
	    j = (bl->index % g->nprocs + bl->rank) % g->nprocs;
	    phgInfo(4, "brecv[%"dFMT"]: gtype=%d, rank=%d, index=%"dFMT
			", type=0x%02x\n",
			rind[jj], bl->gtype, bl->rank, bl->index, bl->btype);
	    assert(bl->btype != UNREFERENCED);
	    if (i > j) {
		i = j;
		rank = bl->rank;
		n0 = bl->lindex;
	    }
	}
	if (rank == -1)
	    phgError(1, "gtype=%d, rank=%d, index=%"dFMT"\n",
				bl->gtype, bl->rank, bl->index);
	/* Note: at this stage data in the range [ii, jj) are no longer
	 * needed, they are replaced with the owner index, owner rank and
	 * updated BTYPE value which will be sent back to incoming process */
	for (jj = ii; jj < n; jj++) {
	    bl = brecv + rind[jj];
	    bl->rank = rank;	/* replace with owner's rank */
	    bl->index = n0;	/* replace with owner's local index */
	    bl->btype = b;	/* replace with updated BTYPE */
	}
	ii = n;
    }
    phgFree(rind);

    /* send back lists */
    MPI_Alltoallv(brecv, rcnts, rdsps, type, bsend, scnts, sdsps, type,
		  g->comm);
    MPI_Type_free(&type);

    /* process received list */
    for (i = 0; i < nsend; i++) {
	bl = bsend + i;
	b = bl->btype;
	if (bl->rank == g->rank && bl->index == bl->lindex)
	    b |= OWNER;
	else
	    b &= ~OWNER;
	switch (bl->gtype) {
	    case VERTEX:
		g->owner_rank_vert[bl->lindex] = bl->rank;
		g->owner_index_vert[bl->lindex] = bl->index;
		g->types_vert[bl->lindex] = b;
		if (!(b & OWNER))
		    g->nvert_owned--;
		break;
	    case EDGE:
		g->owner_rank_edge[bl->lindex] = bl->rank;
		g->owner_index_edge[bl->lindex] = bl->index;
		g->types_edge[bl->lindex] = b;
		if (!(b & OWNER))
		    g->nedge_owned--;
		break;
	    default:
		phgError(1, "unexpected.\n");
	}
    }

    phgFree(bsend);
    phgFree(scnts);

    if ((g->flags & FACE_FLAG)) {
	RNEIGHBOUR *rn;
	struct {
	    ELEMENT	*e;
	    INT		face_no;
	    BYTE	vertex;
	} *sbuf, *rbuf;
	MPI_Datatype type;

	/* update the peer_face member of g->neighbours */
	sbuf = phgAlloc(2 * g->neighbours.count * sizeof(*sbuf));
	rbuf = sbuf + g->neighbours.count;
	for (i = 0; i < g->neighbours.count; i++) {
	    rn = g->neighbours.list + i;
	    sbuf[i].e	    = rn->remote;
	    sbuf[i].vertex  = rn->op_vertex;
	    sbuf[i].face_no = rn->local->faces[rn->vertex];
	}
	MPI_Type_contiguous(sizeof(*sbuf), MPI_BYTE, &type);
	MPI_Type_commit(&type);
	MPI_Alltoallv(sbuf, g->neighbours.counts, g->neighbours.displs, type,
		      rbuf, g->neighbours.counts, g->neighbours.displs, type,
		      g->comm);
	MPI_Type_free(&type);
	for (i = 0; i < g->neighbours.count; i++) {
	    rn = GetRNeighbour(g, rbuf[i].e, rbuf[i].vertex);
	    rn->peer_face = rbuf[i].face_no;
	}
	phgFree(sbuf);

	/* update OWNER flag for faces */
	for (i = 0; i < g->neighbours.count; i++) {
	    rn = g->neighbours.list + i;
	    n0 = GlobalFace(g, j = rn->local->faces[rn->vertex]);
	    assert(g->types_face[j] != UNREFERENCED);
	    if ((g->rank  + (n0 % g->nprocs)) % g->nprocs >
		(rn->rank + (n0 % g->nprocs)) % g->nprocs) {
		g->types_face[j] &= ~OWNER;
		g->nface_owned--;
	    }
	}

	for (i = 0; i < g->nface; i++)
	    if (g->types_face[i] == UNREFERENCED)
		g->nface_owned--;

	g->nelem_owned = g->nleaf;
    }

#if 0 * DEBUG
#warning comment out me!
    /* check g->nxxx_owned members */
    if ((g->flags & VERT_FLAG)) {
	n0 = 0;
	for (i = 0; i < g->nvert; i++) {
	    if ((g->types_vert[i] & OWNER))
		n0++;
	    if (g->types_vert[i] == UNREFERENCED)
		assert(g->owner_rank_vert[i] == -1 &&
			g->owner_index_vert[i] == -1);
	}
	assert(n0 == g->nvert_owned);
    }

    if ((g->flags & EDGE_FLAG)) {
	n0 = 0;
	for (i = 0; i < g->nedge; i++) {
	    if ((g->types_edge[i] & OWNER))
		n0++;
	    if (g->types_edge[i] == UNREFERENCED)
		assert(g->owner_rank_edge[i] == -1 &&
			g->owner_index_edge[i] == -1);
	}
	assert(n0 == g->nedge_owned);
    }

    if ((g->flags & FACE_FLAG)) {
	n0 = 0;
	for (i = 0; i < g->nface; i++) {
	    if ((g->types_face[i] & OWNER))
		n0++;
	    if (g->types_face[i] == UNREFERENCED)
		assert(g->owner_rank_face == NULL ||
			(g->owner_rank_face[i] == -1 &&
			 g->owner_index_face[i] == -1));
	}
	assert(n0 == g->nface_owned);
    }

    if ((g->flags & ELEM_FLAG)) {
	n0 = 0;
	for (i = 0; i < g->nelem; i++) {
	    if ((g->types_elem[i] & OWNER))
		n0++;
	    if (g->types_elem[i] == UNREFERENCED)
		assert(g->owner_rank_elem == NULL ||
			(g->owner_rank_elem[i] == -1 &&
			 g->owner_index_elem[i] == -1));
	}
	assert(n0 == g->nelem_owned);
    }
#endif	/* DEBUG */

#if 0 * DEBUG
    if ((g->flags & EDGE_FLAG) && g->nedge_global > 0) {
	n = 0;
	for (i = 0; i < g->nedge; i++) {
	    if (g->types_edge[i] & OWNER)
		n++;
	}
	ii = g->nedge_global;
	MPI_Allreduce(&n, &jj, 1, MPI_INT, MPI_SUM, g->comm); 
	if (ii != jj) {
	    phgInfo(-1, "total edges: %d, edges owned: %"dFMT"\n", ii, n);
#if 0
	    phgDumpGrid(g);
	    for (i = 0; i < g->nedge; i++) {
		phgInfo(-1, "types_edge[%"dFMT"] = %d (owner = %d)\n",
			GlobalEdge(g, i), g->types_edge[i],
			(g->types_edge[i] & OWNER) != 0);
	    }
#endif
	    phgError(1, "%s:%d, wrong # of edges owned: %d\n",
			__FILE__, __LINE__, jj);
	}
    }

    if ((g->flags & FACE_FLAG) && g->nface_global > 0) {
	n = 0;
	for (i = 0; i < g->nface; i++) {
	    if (g->types_face[i] & OWNER)
		n++;
	}
	ii = g->nface_global;
	MPI_Allreduce(&n, &jj, 1, MPI_INT, MPI_SUM, g->comm); 
	if (ii != jj) {
	    phgInfo(-1, "total faces: %d, faces owned: %"dFMT"\n", ii, n);
#if 0
	    phgDumpGrid(g);
	    for (i = 0; i < g->nface; i++) {
		phgInfo(-1, "types_face[%"dFMT"] = %d (owner = %d)\n",
			GlobalFace(g, i), g->types_face[i],
			(g->types_face[i] & OWNER) != 0);
	    }
#endif
	    phgError(1, "%s:%d, wrong # of faces owned: %d\n",
			__FILE__, __LINE__, jj);
	}
    }
#endif	/* DEBUG */

#endif	/* USE_MPI */

    if (phgVerbosity > 0)
	phgPrintf("phgUpdateBoundaryTypes: wall time = %lg\n", 
			phgGetTime(NULL) - time0);
}

void *
phgAlloc(size_t size)
{
    void *ptr;

    if (size != 0) {
#if OMP_ALLOC_HACK && USE_OMP
# pragma omp critical (alloc)
#endif	/* OMP_ALLOC_HACK && USE_OMP */
	ptr = malloc(size);
    }
    else {
	ptr = NULL;
    }
    
    if (ptr == NULL && size != 0)
	phgError(1, "failed to malloc %llu bytes memory.\n", (long long)size);

    return ptr;
}

void *
phgCalloc(size_t nmemb, size_t size)
{
    void *ptr;

    if (nmemb != 0 && size != 0) {
#if OMP_ALLOC_HACK && USE_OMP
# pragma omp critical (alloc)
#endif	/* OMP_ALLOC_HACK && USE_OMP */
	ptr = calloc(nmemb, size);
    }
    else {
	ptr = NULL;
    }

    if (ptr == NULL && nmemb != 0 && size != 0)
	phgError(1, "%s:%d, failed to calloc %llu*%llu bytes memory.\n",
			__FILE__, __LINE__, (long long)nmemb, (long long)size);

    return ptr;
}

void *
phgRealloc0(void *ptr, size_t size)
{
    void *p;

#if OMP_ALLOC_HACK && USE_OMP
# pragma omp critical (alloc)
#endif	/* OMP_ALLOC_HACK && USE_OMP */
    p = realloc(ptr, size);
    
    if (p == NULL && size != 0)
	phgError(1, "failed to reallocate %llu bytes memory at %p.\n",
			(long long)size, ptr);

    return p;
}

void *
phgRealloc__(void *ptr, size_t size, size_t oldsize, const char *file, int line)
/* Note: this function is a workaround on DeepComp7000 when using MVAPICH2
 * 	 (if realloc() fails, then malloc() + memcpy() is used).
 *
 * 	 'oldsize' needs not to be the actual allocated size of the buffer,
 * 	 it should be >= size_of_data_to_copy and <= allocated_size */
{
    void *p;

#if OMP_ALLOC_HACK && USE_OMP
# pragma omp critical (alloc)
#endif	/* OMP_ALLOC_HACK && USE_OMP */
    p = realloc(ptr, size);
    
    if (p == NULL && size != 0) {
	static BOOLEAN warned = FALSE;
	if (phgRank == 0 && !warned) {
	    phgWarning("(%s:%d) realloc() failed, try malloc().\n", file, line);
	    warned = TRUE;
	}
#if OMP_ALLOC_HACK && USE_OMP
# pragma omp critical (alloc)
#endif	/* OMP_ALLOC_HACK && USE_OMP */
	p = malloc(size);
	if (p == NULL)
	    phgError(1, "(%s:%d) malloc() failed (size=%llu), abort.\n",
			file, line, (long long)size);
	if (oldsize > 0 && size > 0 && ptr != NULL)
	    memcpy(p, ptr, oldsize < size ? oldsize : size);
	if (ptr != NULL) {
#if OMP_ALLOC_HACK && USE_OMP
# pragma omp critical (alloc)
#endif	/* OMP_ALLOC_HACK && USE_OMP */
	    free(ptr);
	}
    }

    return p;
}

void
phgFree(void *ptr)
{
    if (ptr != NULL) {
#if 1 && USE_OMP
#pragma omp critical (alloc)
#endif	/* ?? && USE_OMP */
	free(ptr);
    }
}

int
phgMemcmp(const void *v0, const void *v1, size_t n)
/* replaces buggy memcmp() on many systems */
{
    const BYTE *s0 = v0, *s1 = v1;
    int i = 0;

    for (; n--; )
	if ((i = *(s0++) - *(s1++)) != 0)
	    break;

    return i;
}

void
phgFreeAtExit(void **ptr, int level)
{
    static void ***ptr_table = NULL;
    static size_t table_n = 0, allocated = 0;
    static int level_max = 0, *level_table = NULL;

    if (ptr == NULL) {	/* free all registered pointers */
	int i;

	CheckThread
	for (level = 0; level <= level_max; level++) {
	    for (i = table_n - 1; i >= 0; i--) {
		if (level_table[i] != level)
		    continue;
		phgFree(*(ptr_table[i]));
		*(ptr_table[i]) = NULL;	/* avoid double frees */
	    }
	}
	phgFree(ptr_table);
	ptr_table = NULL;
	phgFree(level_table);
	level_table = NULL;
	table_n = allocated = 0;
	level_max = 0;

	return;
    }

#if USE_OMP
#pragma omp critical (FreeAtExit)
{
#endif	/* USE_OMP */
    if (table_n >= allocated) {
	ptr_table = (void***)phgRealloc_(ptr_table,
				sizeof(*ptr_table) * (allocated + 128),
				sizeof(*ptr_table) * allocated);
	level_table = phgRealloc_(level_table,
				sizeof(*level_table) * (allocated + 128),
				sizeof(*level_table) * allocated);
	allocated += 128;
    }
    level_table[table_n] = level;
    ptr_table[table_n++] = ptr; 
#if USE_OMP
}
#endif	/* USE_OMP */
    if (level_max < level)
	level_max = level;

    assert(table_n <= 1000);	/* for trapping leaks */

    return;
}

double
phgGetTime(double tarray[])
/* tarray[0]=User time, tarray[1]=Sys time, tarray[2]=Elapsed time,
   all expressed in seconds.

   DAWN1000 nf77
	double F77_FUNC(dclock,DCLOCK)();
	tarray[0]=tarray[1]=tarray[2]=F77_FUNC(dclock,DCLOCK)(); */
{
    struct timeval tv;
  
    gettimeofday(&tv, (struct timezone *)0);

    if (tarray == NULL)
	return tv.tv_sec + (double)tv.tv_usec * 1e-6;

    tarray[2] = tv.tv_sec + (double)tv.tv_usec * 1e-6;

#if defined(__WIN__) || defined(__SUNWAY__)
    /* TODO: figure out how to get user/sys time under Windows */
    tarray[0] = tarray[1] = tarray[2];
#elif defined(USETIMES) /*==================================================*/
    struct tms ts;
    static double ClockTick = 0.0;

    if (ClockTick == 0.0)
	    ClockTick = (double)sysconf(_SC_CLK_TCK);
    if ( times(&ts) == -1 ) {
	    tarray[0]=tarray[1] = 0.0;
    } else {
	    tarray[0] = (double)ts.tms_utime / ClockTick;
	    tarray[1] = 0.0;
    }
#else		 /*==================================================*/
    /* For Solaris */
    /*#include </usr/ucbinclude/sys/rusage.h> */
    struct rusage RU;
    getrusage(RUSAGE_SELF, &RU);
    tarray[0] = RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec * 1e-6;
    tarray[1] = RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec * 1e-6;
#endif		 /*==================================================*/

    return tarray[2];
}

void
phgWarning(const char *fmt, ...)
{
    va_list ap;
    char s[1024];
    const char *n;

    init_log_file();

    n = phgCurrentFunctionName_;
    if (strlen(n) >= 3 && !phgMemcmp(n, "phg", 3))
	n += 3;

    if (phgNProcs > 1)
	sprintf(s, "*%d %s: *** WARNING: %s", phgRank, n, fmt);
    else
	sprintf(s, "* %s: *** WARNING: %s", n, fmt);
    va_start(ap, fmt);
    vfprintf(LogFile, s, ap);
    va_end(ap);
    fflush(LogFile);
}

void
phgAbort(int code)
{
    close_log_file();

    /* leave time for processes to clean up */
    if (phgNProcs > 1)
	usleep(500000);

#if HAVE_FCLOSEALL
    fcloseall();
#endif	/* HAVE_FCLOSEALL */

#if USE_MPI
    if (phgInitialized && phgNProcs > 1) {
	MPI_Abort(phgComm, code);
    }
#endif

    exit(code);
}

void 
phgError(int code, const char *fmt, ...)
{
    va_list ap;
#if USE_MPI
    char s[4096 + MPI_MAX_PROCESSOR_NAME], pn[MPI_MAX_PROCESSOR_NAME + 3];
    int len;
#else
    char s[4096], pn[1];
#endif
    const char *n;

    init_log_file();

    n = phgCurrentFunctionName_;
    if (strlen(n) >= 3 && !phgMemcmp(n, "phg", 3))
	n += 3;

#if USE_MPI
    pn[0] = ' ';
    pn[1] = '(';
    MPI_Get_processor_name(pn + 2, &len);
    pn[2 + len] = ')';
    pn[2 + len + 1] = '\0';
#else
    pn[0] = '\0';
#endif
    if (phgNProcs > 1)
	sprintf(s, "\n*%d%s %s: *** ERROR: %s\n", phgRank, pn, n, fmt);
    else
	sprintf(s, "\n*%s %s: *** ERROR: %s\n", pn, n, fmt);
    va_start(ap, fmt);
    vfprintf(LogFile, s, ap);
    va_end(ap);
    fflush(LogFile);

    if (code == 0)
	return;
    phgAbort(code);
}

void
phgInfo(int verbose_level, const char *fmt, ...)
/* writes formatted message to log_file when log_file != NULL,
 * if verbose_level < 0 then writes to stdout if log_file == NULL */
{
    va_list ap;
    char s[1024];
    const char *n;

    init_log_file();

#if USE_MPI
    /*if (phgRank != 0)
	return;*/
#endif
    if ((verbose_level >= 0 && log_file == NULL) ||
	verbose_level > phgVerbosity)
	return;

    n = phgCurrentFunctionName_;
    if (strlen(n) >= 3 && !phgMemcmp(n, "phg", 3))
	n += 3;

    if (phgNProcs > 1)
	sprintf(s, "*%d %s: %s", phgRank, n, fmt);
    else
	sprintf(s, "* %s: %s", n, fmt);
    va_start(ap, fmt);
    vfprintf(LogFile, s, ap);
    va_end(ap);
    fflush(LogFile);
}

#undef Log

int
phgPrintf(const char *fmt, ...)
{
    va_list ap;
    int ret = 0;

    if (phgRank == 0) {
	if (out_file == NULL && phgOutFilename != NULL) {
	    if (strcmp(phgOutFilename, "stdout") == 0) {
		out_file = stdout;
	    }
	    else if (strcmp(phgOutFilename, "stderr") == 0) {
		out_file = stderr;
	    }
	    else {
		out_file = fopen(phgOutFilename, "w+t");
		if (out_file == NULL) {
		    phgWarning("cannot open output file \"%s\".",
							phgOutFilename);
		    phgFree(phgOutFilename);
		    phgOutFilename = NULL;
		}
	    }
	}
	va_start(ap, fmt);
	ret = vfprintf(OutFile, fmt, ap);
	va_end(ap);
	fflush(OutFile);
    }

#if USE_MPI
    /*if (g->nprocs > 0)
	MPI_Bcast(&ret, sizeof(ret), MPI_BYTE, 0, g->comm);*/
#endif

    return ret;
}


int
phgPrintf2(GRID *g, const char *fmt, ...)
/* Same as phgPrint, but use rank 0 of grid to print. */
{
    va_list ap;
    int ret = 0;

    if (g->rank == 0) {
	if (out_file == NULL && phgOutFilename != NULL) {
	    if (strcmp(phgOutFilename, "stdout") == 0) {
		out_file = stdout;
	    }
	    else if (strcmp(phgOutFilename, "stderr") == 0) {
		out_file = stderr;
	    }
	    else {
		out_file = fopen(phgOutFilename, "w+t");
		if (out_file == NULL) {
		    phgWarning("cannot open output file \"%s\".",
							phgOutFilename);
		    phgFree(phgOutFilename);
		    phgOutFilename = NULL;
		}
	    }
	}
	va_start(ap, fmt);
	ret = vfprintf(OutFile, fmt, ap);
	va_end(ap);
	fflush(OutFile);
    }

#if USE_MPI
    /*if (g->nprocs > 0)
	MPI_Bcast(&ret, sizeof(ret), MPI_BYTE, 0, g->comm);*/
#endif

    return ret;
}


FILE *
phgFOpen(GRID *g, const char *path, const char *mode)
{
#if USE_MPI
    FILE *fp;

    if (g->rank == 0)
	fp = fopen(path, mode);
# if 1
    else
	fp = (void *)-1;
# else
    if (g->nprocs > 1)
	MPI_Bcast(&fp, sizeof(fp), MPI_BYTE, 0, g->comm);
# endif
    return fp;
#else
    Unused(g);
    return fopen(path, mode);
#endif
}

int
phgFPrintf(GRID *g, FILE *fp, char *fmt, ...)
{
    va_list ap;
    int ret;

#if USE_MPI
    if (g->rank == 0) {
	va_start(ap, fmt);
	ret = vfprintf(fp, fmt, ap);
	va_end(ap);
	fflush(fp);
    }
# if 1
    else
	ret = 1;
# else
    if (g->nprocs > 1)
	MPI_Bcast(&ret, sizeof(ret), MPI_BYTE, 0, g->comm);
# endif
#else
    Unused(g);
    va_start(ap, fmt);
    ret = vfprintf(fp, fmt, ap);
    va_end(ap);
    fflush(fp);
#endif

    return ret;
}

int
phgFClose(GRID *g, FILE *fp)
{
    int ret;

#if USE_MPI
    if (g->rank == 0)
	ret = fclose(fp);
# if 1
    else
	ret = 0;
# else
    if (g->nprocs > 1)
	MPI_Bcast(&ret, sizeof(ret), MPI_BYTE, 0, g->comm);
# endif
#else
    Unused(g);
    ret = fclose(fp);
#endif

    return ret;
}

BYTE
phgOppositeFace(GRID *g, ELEMENT *e, BYTE v, ELEMENT *e_op)
/* returns the vertex e_op opposite to the vertex v of e.
   This function does not rely on the neighbours links */
{
    INT i, u0 = 0, u1 = 0, u2 = 0;

#if DEBUG
    if (e_op == NULL)
	phgError(1, "unexpected error in %s (possibly caused by inconsistent "
		 "boundary types in the input mesh file).\n", __func__);
#endif	/* DEBUG */

    switch (v) {
	case 0: u0 = e->verts[1]; u1 = e->verts[2]; u2 = e->verts[3]; break;
	case 1: u0 = e->verts[0]; u1 = e->verts[2]; u2 = e->verts[3]; break;
	case 2: u0 = e->verts[0]; u1 = e->verts[1]; u2 = e->verts[3]; break;
	case 3: u0 = e->verts[0]; u1 = e->verts[1]; u2 = e->verts[2]; break;
    }

    if (g->period == NULL) {
	if ((i = e_op->verts[0]) != u0 && i != u1 && i != u2)
	    return 0;
	if ((i = e_op->verts[1]) != u0 && i != u1 && i != u2)
	    return 1;
	if ((i = e_op->verts[2]) != u0 && i != u1 && i != u2)
	    return 2;
	if ((i = e_op->verts[3]) != u0 && i != u1 && i != u2)
	    return 3;
    }
    else {
#define VERTEX(e, i)	g->period->L2Gmap_vert[e->verts[i]]
	u0 = g->period->L2Gmap_vert[u0];
	u1 = g->period->L2Gmap_vert[u1];
	u2 = g->period->L2Gmap_vert[u2];
	if ((i = VERTEX(e_op, 0)) != u0 && i != u1 && i != u2)
	    return 0;
	if ((i = VERTEX(e_op, 1)) != u0 && i != u1 && i != u2)
	    return 1;
	if ((i = VERTEX(e_op, 2)) != u0 && i != u1 && i != u2)
	    return 2;
	if ((i = VERTEX(e_op, 3)) != u0 && i != u1 && i != u2)
	    return 3;
#undef VERTEX
    }

    return 4;	/* make gcc happy */
}

void
phgFreeElement(ELEMENT **eptr)
{
    ELEMENT *e;

    if (eptr == NULL || (e = *eptr) == NULL)
	return;
#if 0
    if (e->userdata != NULL)
	phgFree(e->userdata);
#endif
    phgFree(e);
    *eptr = NULL;
}

static void
free_grid_traverse(ELEMENT *s)
{
    if (s->children[0] != NULL)
	free_grid_traverse(s->children[0]);
    if (s->children[1] != NULL)
	free_grid_traverse(s->children[1]);
    phgFreeElement(&s);
}

void
phgFreeGrid_(GRID *g)
{
    INT i, n;

    /* FIXME: meshes with comm = MPI_COMM_SELF won't be freed on non root
     *	      processes, changing the condition to 'TRUE' fixes it, but
     *	      causes phg_tcl to hang when switching displaying of submeshes
     *	      on non root process (no problem when phgMasterSlave == FALSE). */
    ParallelFunction(phgFreeGrid_, g != NULL && g->nprocs > 1);

    if (g == NULL)
	return;

    for (i = 0; i < g->nroot; i++) {
	ELEMENT *s = g->roots + i;
	if (s->children[0] != NULL)
	    free_grid_traverse(s->children[0]);
	if (s->children[1] != NULL)
	    free_grid_traverse(s->children[1]);
#if 0
	if (s->userdata != NULL)
	    phgFree(s->userdata);
#endif
    }
    phgFree(g->filename);
    phgFree(g->roots);
    phgFree(g->verts);
    phgFree(g->elems);
#if ALLOW_CURVED_BOUNDARY
    if (g->bdry_funcs != NULL) {
	EXPR **p = g->bdry_funcs;
	while (*p != NULL)
	    phgFree3DFunction(*(p++));
	phgFree(g->bdry_funcs);
    }
#endif
    phgFree(g->types_vert);
    phgFree(g->types_edge);
    phgFree(g->types_face);
    phgFree(g->types_elem);
    phgFree(g->alien_map);
    phgPeriodFree(g);
    phgFree(g->owner_index_vert);
    phgFree(g->owner_rank_vert);
#if USE_MPI
# if !NO_NEW_COMMUNICATOR
    if (g->comm != MPI_COMM_NULL)
	MPI_Comm_free(&g->comm);
# endif	/* !NO_NEW_COMMUNICATOR */
    phgFree(g->L2Gmap_vert);
    phgFree(g->L2Gmap_edge);
# if Dim == 3
    phgFree(g->L2Gmap_face);
# endif	/* Dim == 3 */
    phgFree(g->L2Gmap_elem);

    phgFree(g->owner_index_edge);
    phgFree(g->owner_rank_edge);
    phgFree(g->owner_index_face);
    phgFree(g->owner_rank_face);
    phgFree(g->owner_index_elem);
    phgFree(g->owner_rank_elem);

    phgFree(g->neighbours.counts);
    phgFree(g->neighbours.displs);
    phgFree(g->neighbours.list);
#endif	/* USE_MPI */

    /* avoid warning on unfreed DOF */
    if (g->geom != NULL)
	phgDofFree(&g->geom);

    /* Free DOFs attached to the grid */
    if (g->dof != NULL) {
	DOF **p;
	for (p = g->dof; *p != NULL; p++);
	i = (INT)(p - g->dof);
	while (i > 0) {
	    DOF *d = g->dof[--i];
#if DEBUG
	    /* For helping finding out memory leaks due to unfreed DOFs */
	    phgPrintf("*** WARNING: unfreed DOF \"%s\" (created at %s:%d)\n",
			d->name, d->srcfile, d->srcline);
#endif
	    d->g = NULL; /* make phgGridIsValid return FALSE in phgDofFree */
	    phgDofFree(&d);
	}
	phgFree(g->dof);
    }

    /* Free HP_TYPEs attached to the grid */
    if (g->hp != NULL) {
	HP_TYPE **php;
	for (php = g->hp; *php != NULL; php++) {
#if DEBUG
	    /* Warn unfreed HP_TYPEs */
	    phgPrintf("*** WARNING: unfreed HP_TYPE (created at %s:%d)\n",
			(*php)->srcfile, (*php)->srcline);
#endif
	    (*php)->refcount = 0;	/* force freeing *php */
	    phgHPFree(php);
	}
	phgFree(g->hp);
    }

    phgFree(g->bc_list);
    phgFree(g->bc_rmap);

    phgFree(g);

    /* remove g from g_list */
    for (i = 0; i < g_list_count; i++)
	if (g_list[i] == g) break;

    if (i >= g_list_count)
	phgError(1, "phgFreeGrid: invalid grid!\n");

    n = g_list_count - i - 1;
    if (n > 0) {
	memmove(g_list + i, g_list + i + 1, n * sizeof(*g_list));
    }
    if (--g_list_count > 0) {
	g_list = phgRealloc_(g_list, g_list_count * sizeof(*g_list),
				     (g_list_count + 1) * sizeof(*g_list));
    }
    else {
	phgFree(g_list);
	g_list = NULL;
    }

    /* for Master/Slave mode in which only one grid is allowed */
    if (_phg_g_sys == g)
	_phg_g_sys = NULL;

    return;
}

void
phgFreeGrid(GRID **gptr)
{
    GRID *g = (gptr == NULL ? NULL : *gptr);

    if (g == NULL)
	return;

    if (gptr != NULL)
	*gptr = NULL;

    phgFreeGrid_(g);
}

/* flag controlling whether to call callback function for non-leaf elements */
static BOOLEAN traverse_all_flag = FALSE;

/* Depth of the element in the tree.
   This variable is made global and may be used by the callback functions */
int phgTraverseDepth = 0;
INT phgTraverseIndex = 0;

static BOOLEAN
traverse_element(ELEMENT *e, BOOLEAN(*cb) CB_ARGS(e))
{
    BOOLEAN is_leaf = TRUE;

    if (e->children[0] != NULL) {
        is_leaf = FALSE;
	phgTraverseDepth++;
	if (!traverse_element(e->children[0], cb))
	    return FALSE;
	phgTraverseDepth--;
    }
    if (e->children[1] != NULL) {
        is_leaf = FALSE;
	phgTraverseDepth++;
	if (!traverse_element(e->children[1], cb))
	    return FALSE;
	phgTraverseDepth--;
    }
    if (!is_leaf && !traverse_all_flag) {
	return TRUE;
    }
    else {
	BOOLEAN ret = cb CB_ARGS0(g_traverse, e);
	Unused(g_traverse);
	phgTraverseIndex++;
	return ret;
    }
}

void
phgTraverseElements(GRID *g, BOOLEAN(*cb) CB_ARGS(e))
/* function traversing (traverse_all_flag ? leaf : all) elements.
 * Stop traversing if callback function returns FALSE */
{
    INT i;

    if (g->nroot == 0)
	return;

    g_traverse = g;
    phgTraverseIndex = 0;
    for (i = 0; i < g->nroot; i++) {
	phgTraverseDepth = 0;
	if (!traverse_element(g->roots + i, cb)) break;
    }
}

void
phgTraverseAllElements(GRID *g, BOOLEAN(*cb) CB_ARGS(e))
/* function traversing all (including non leaf) elements */
{
    traverse_all_flag = TRUE;
    phgTraverseElements(g, cb);
    traverse_all_flag = FALSE;
}

#define EXTENSIVE_TEST 0

#if DEBUG && EXTENSIVE_TEST
static BOOLEAN
singular_solve(int n, int m, double a[n][m], double b[n])
/* solves the singular system AX == B (A is a singular NxM dense matrix),
 * returns FALSE if no solution */
{
    int i, j, k;
    double d, e;

    assert(m <= n);

    /* normalize the coefficients */
    for (i = 0; i < n; i++) {
	d = fabs(a[i][0]);
	for (j = 1; j < m; j++) {
	    if (d < (e = fabs(a[i][j])))
		d = e;
	}
	if (d == 0.)
	    continue;
	d = 1. / d;
	for (j = 1; j < m; j++)
	    a[i][j] *= d;
	b[i] *= d;
    }

    /* Gauss elimination */
    for (i = 0; i < m; i++) {
	k = i;
	d = fabs(a[i][i]);
	for (j = i + 1; j < n; j++) {
	    if (d < (e = fabs(a[j][i]))) {
		d = e;
		k = j;
	    }
	}
	if (d < 1e-12)
	    return FALSE;
	if (k != i) {
	    /* exchange row i with row k */
	    d = b[i];
	    b[i] = b[k];
	    b[k] = d;
	    for (j = 0; j < m; j++) {
		d = a[i][j];
		a[i][j] = a[k][j];
		a[k][j] = d;
	    }
	}
	/* elimination */
	d = 1. / a[i][i];
	for (j = 0; j < m; j++)
	    a[i][j] *= d;
	b[i] *= d;
	for (j = i + 1; j < n; j++) {
	    for (k = i + 1; k < m; k++)
		a[j][k] -= a[j][i] * a[i][k];
	    b[j] -= a[j][i] * b[i];
	}
    }

    /* check consistency of the system */
    for (; i < n; i++) {
	if (fabs(b[i]) > 1e-12)
	    return FALSE;
    }

    /* back substitute */
    for (i = m - 2; i >= 0; i--)
	for (j = i + 1; j < m; j++)
	    b[i] -= a[i][j] * b[j];

    return TRUE;
}

static BOOLEAN
check_intersect(COORD *v0, COORD *v1, COORD *u0, COORD *u1)
/* returns TRUE if the segments v0-v1 and u0-u1 intersect */
{
    int i;
    double a[Dim][2], b[Dim];

    for (i = 0; i < Dim; i++) {
	b[i] = (*u0)[i] - (*v0)[i];
	a[i][0] = (*v1)[i] - (*v0)[i];
	a[i][1] = -((*u1)[i] - (*u0)[i]);
    }

    if (!singular_solve(Dim, 2, a, b))
	return FALSE;

    if (b[0] < -1e-12 || b[0] > 1 + 1e-12 || b[1] < -1e-12 || b[1] > 1 + 1e-12)
	return FALSE;

    return TRUE;
}

static void
check2(ELEMENT *s, int face, int fa, int fb, int fc, ELEMENT *s1)
/* check the case in which 's1' has two vertices on the give face of 's'
 * (vertices 'fa' and 'fb' are shared by 's' and 's1').
 *
 * This check handles a special case in which two edges of two faces
 * sharing two common vertices intersect. This typically happens when
 * a hexahedron is adjacent to a tetrahedron and is tetrahedralized
 * inconsistently */
{
    int j;
    INT v0, v1, v2, v3;
    COORD *coor;

    if (IS_BDRY(s->bound_type[face]))
	return;

    v0 = s->verts[GetFaceVertex(face, fa)],
    v1 = s->verts[GetFaceVertex(face, fb)],
    v2 = s->verts[GetFaceVertex(face, fc)];
    coor = g_->verts;

    /* make sure edges of s1 does not cross the face */
    for (j = 0; j < NVert; j++) {
	if ((v3 = s1->verts[j]) == v0 || v3 == v1)
	    continue;
	assert(v3 != v2);
	/* (v3-v0) should not intersect with (v2-v1) and
	   (v3-v1) should not intersect with (v2-v0) */
	if (check_intersect(coor + v3, coor + v0, coor + v2, coor + v1) ||
	    check_intersect(coor + v3, coor + v1, coor + v2, coor + v0))
	    phgError(1, "nonconforming mesh, abort.\n");
    }
}

typedef struct {
    ELEMENT *e;
    BYTE v;
} ELIST;

static int
comp_elist(const void *p0, const void *p1)
{
    const ELIST *el0 = p0, *el1 = p1;
    return (el0->e > el1->e ? 1 : (el0->e == el1->e ? 0 : -1));
}

#endif	/* DEBUG && EXTENSIVE_TEST */

static BOOLEAN conformity; /* used to check conformity during traversal */

static struct {
    ELEMENT **elements;
    SHORT   nelem, nalloc;
} *vertices, *pv, *pv0[Dim];
static INT *L2Lmap;

static BOOLEAN
find_neighbours_callback0 CB_ARGS(s)
{
    INT k;
    int j;

    /* clear all neighbour links */
    memset(s->neighbours, 0, sizeof(s->neighbours));

    for (j = 0; j < NVert; j++) {
	k = s->verts[j];
	/* map periodic vertex */
	if (L2Lmap != NULL)
	    k = L2Lmap[k];
	pv = vertices + k;
	if (pv->nelem >= pv->nalloc) {
	    pv->elements = phgRealloc_(pv->elements,
				(pv->nalloc + 32) * sizeof(s),
				pv->nalloc * sizeof(s));
	    pv->nalloc += 32;
	}
	pv->elements[pv->nelem++] = s;
    }

    if (s->index == -1)
	s->index = phgTraverseIndex;

    return TRUE;
}

static BOOLEAN
find_neighbours_callback1 CB_ARGS(s)
{
#if 0
    /* new code based on sorting the list of neighbours of local vertices.
     * FIXME: the code is much slower and breaks mesh coarsening algorithm */
    ELEMENT *s1;
    struct {
	INT	eno;
	BYTE	v;
    } *elist = NULL;
    int i, j, n, m;
    INT eno;

    /* count number of elements in the lists of the 3 vertices */
    n = 0;
    for (j = 0; j < NVert; j++) {
	s->neighbours[j] = NULL;
	n += vertices[s->verts[j]].nelem;
    }
    n -= NVert;
    if (n <= 0)
	goto cont;
    elist = phgAlloc(n * sizeof(*elist));
    /* collect element indices */
    n = 0;
    for (j = 0; j < NVert; j++) {
	for (i = 0; i < vertices[s->verts[j]].nelem; i++) {
	    if ((s1 = vertices[s->verts[j]].elements[i]) == s)
		continue;
	    elist[n].eno = s1 - g_->roots;
	    elist[n++].v = j;
	}
    }
    qsort(elist, n, sizeof(*elist), phgCompINT);

    eno = -1;	/* current element */
    m = 0;	/* count of duplicate elements */
    for (j = 0; j <= n; j++) {
	if (j < n && eno == elist[j].eno) {
	    m++;
	    continue;
	}
	assert(m <= 3);
	if (m == 3) {
	    BYTE v0 = elist[j - 3].v, v1 = elist[j - 2].v, v2 = elist[j - 1].v;
	    /* find the face number */
	    if (v0 != 0 && v1 != 0 && v2 != 0)
		i = 0;
	    else if (v0 != 1 && v1 != 1 && v2 != 1)
		i = 1;
	    else if (v0 != 2 && v1 != 2 && v2 != 2)
		i = 2;
	    else
		i = 3;
	    s->bound_type[i] &= ~(REMOTE | UNDEFINED);
	    s->bound_type[i] |= INTERIOR;
	    s->neighbours[i] = g_->roots + eno;
	}
	else if (m == 2) {
	    /* TODO: further conformity check */
	}
	if (j >= n)
	    break;
	eno = elist[j].eno;
	m = 1;
    }

cont:
    for (i = 0; i < NVert; i++) {
	if (s->neighbours[i] != NULL)
	    continue;
	if (s->bound_type[i] & INTERIOR) {
	    if (g_->nprocs == 1) {
		conformity = FALSE;
	    }
	    else {
		s->bound_type[i] |= REMOTE;
		s->neighbours[i] = NULL;
	    }
	}
    }

    phgFree(elist);
#else	/* #if 0|1 */
    ELEMENT *s1 = NULL;
    INT i;
    int face, j, k, l;

    for (face = 0; face < NFace; face++) {
	k = 0;
	for (j = 0; j < NVert; j++) {
	    if (j == face)
		continue;
	    i = s->verts[j];
	    /* map periodic vertex */
	    if (L2Lmap != NULL)
		i = L2Lmap[i];
	    pv0[k++] = vertices + i;
	}
	/* element appearing in all the Dim lists is the neighbour */
	for (k = 0; k < pv0[0]->nelem; k++) {
	    s1 = pv0[0]->elements[k];
	    if (s1 == s)
		continue;
	    for (l = 0; l < pv0[1]->nelem; l++)
		if (s1 == pv0[1]->elements[l])
		    break;
	    if (l >= pv0[1]->nelem)
		continue;
#if Dim == 3
	    for (l = 0; l < pv0[2]->nelem; l++)
		if (s1 == pv0[2]->elements[l])
		    break;
	    if (l >= pv0[2]->nelem)
		continue;
#endif	/* Dim == 3 */
	    break;
	}
	if (k < pv0[0]->nelem) {
	    s->bound_type[face] &= ~(REMOTE | UNDEFINED);
	    s->bound_type[face] |= INTERIOR;
	    s->neighbours[face] = s1;
	}
	else {
	    if (s->bound_type[face] & INTERIOR) {
		if (g_->nprocs == 1) {
		    /* Note: can only detect non conformity after coarsening */
		    conformity = FALSE;
		}
		else {
		    s->bound_type[face] |= REMOTE;
		    s->neighbours[face] = NULL;
		}
	    }
	    /* already cleared by find_neighbours_callback0() */
	    /* s->neighbours[face] = NULL;*/
	}
#if DEBUG && EXTENSIVE_TEST	/* Slow! don't enable unless for debugging */
	/* Check intersecting edges */
	{
	    static BOOLEAN warned = FALSE;
	    ELIST *elist = NULL;
	    int n, m;

	    if (!warned) {
		warned = TRUE;
		phgWarning("slow code! don't enable except for debugging.\n");
	    }

	    /* count number of elements in the lists of the 3 vertices */
	    n = pv0[0]->nelem + pv0[1]->nelem + pv0[2]->nelem - 3;
	    if (n <= 0)
		continue;
	    elist = phgAlloc(n * sizeof(*elist));
	    /* collect element indices */
	    n = 0;
	    for (j = 0; j < 3; j++) {
		for (m = 0; m < pv0[j]->nelem; m++) {
		    if ((s1 = pv0[j]->elements[m]) == s)
			continue;
		    elist[n].e = s1;
		    elist[n++].v = j;
		}
	    }
	    qsort(elist, n, sizeof(*elist), comp_elist);
	
	    s1 = NULL;	/* current element */
	    m = 0;	/* count of duplicate elements */
	    for (j = 0; j <= n; j++) {
		if (j < n && s1 == elist[j].e) {
		    m++;
		    continue;
		}
		assert(m <= 3);
		if (m == 2) {
		    int v0 = elist[j - 1].v, v1 = elist[j - 2].v, v2; 
		    if (((v2 = 0) != v0 && v2 != v1) ||
			((v2 = 1) != v0 && v2 != v1) ||
			((v2 = 2) != v0 && v2 != v1))
			check2(s, face, v0, v1, v2, s1);
		}
		if (j >= n)
		    break;
		s1 = elist[j].e;
		m = 1;
	    }
	}
#endif /* DEBUG && EXTENSIVE_TEST */
    }
#endif	/* #if 0|1 */

    return TRUE;
}

void
phgUpdateNeighbours(GRID *g)
/* fills 'neighbours' array of all simplices. */
{
    INT i, j, k, n;

    FunctionEntry;

    if (g == NULL)
	return;

    g_ = g;

    if (g->period == NULL) {
	L2Lmap = NULL;
    }
    else if (g->nprocs == 1) {
	L2Lmap = g->period->L2Gmap_vert;
    }
    else {
	/* build L2Lmap */
	L2Lmap = phgAlloc(g->nvert * sizeof(*L2Lmap));
	for (i = 0; i < g->nvert; i++)
	    L2Lmap[i] = i;
	qsort(L2Lmap, g->nvert, sizeof(*L2Lmap), comp_pvertex);
	n = -1;
	k = -1;
	for (i = 0; i < g->nvert; i++) {
	    if ((j = g->period->L2Gmap_vert[L2Lmap[i]]) != k) {
		k = j;
		n++;
	    }
	    L2Lmap[i] = n;
	}
	n++;
    }

#if 0
    if (g->nprocs == 1 && g->nroot == g->nleaf) {
	/* check for duplicate tetrahedra in the initial mesh */
	INT j;
	qsort(g->roots, g->nroot, sizeof(*g->roots), elem_comp);
	j = -1;
	for (i = 0; i < g->nroot; i++) {
	    if (j < 0 || elem_comp(g->roots + j, g->roots + i) != 0)
		j++;
	    if (j != i)
		g->roots[j] = g->roots[i];
	}
	if (j < g->nroot)
	    phgWarning("duplicate tetrahedra found and removed.\n");
    }
#endif

    vertices = phgCalloc(g->nvert, sizeof(*vertices));

    /* For each vertex, a list of elements to which it belongs is built */
    phgTraverseElements(g, find_neighbours_callback0);

    /* then search neighbours in the small lists of elements */
    conformity = TRUE;
    phgTraverseElements(g, find_neighbours_callback1);

    /* reset lists of elements */
    for (i = 0; i < g->nvert; i++) {
	if (vertices[i].elements != NULL) {
	    phgFree(vertices[i].elements);
	    /*vertices[i].elements = NULL;
	    vertices[i].nelem = vertices[i].nalloc = 0;*/
	}
    }
    phgFree(vertices);

    if (L2Lmap != NULL && L2Lmap != g->period->L2Gmap_vert)
	phgFree(L2Lmap);

    PrintTime(1);
    Return;
}

static BOOLEAN
check_conformity_callback CB_ARGS(e)
{
    INT i, j, k;
    ELEMENT *e1;

    for (i = 0; i < NVert; i++) {
	/* check boundary types */
#if 0
	if (e->bound_type[i] == 0) {
	    phgInfo(-1, "undefined boundray type for face %"dFMT"\n", i);
	    phgDumpElement(g_, e);
	    return conformity = FALSE;
	}
#endif
	/* check indices of vertices */
	if (e->verts[i] < 0 || e->verts[i] >= g_->nvert) {
	    phgInfo(-1, "invalid vertex index %"dFMT"\n", e->verts[i]);
	    phgDumpElement(g_, e);
	    return conformity = FALSE;
	}
	/* check links to neighbours */
	if (!HasLocalNeighbour(e, i))
	    continue;
	e1 = e->neighbours[i];
	if (e1 == NULL) {
	    phgInfo(-1, "interior face (%"dFMT") without neighbour:\n", i);
	    phgDumpElement(g_, e);
	    if (i > 1)
		phgDumpPatch(g_, e, 0, 1);
	    else
		phgDumpPatch(g_, e, 2, 3);
	    return conformity = FALSE;
	}
	if (!IsLeaf(e1)) {
	    phgInfo(-1, "non-root neighbour:\n");
	    phgDumpElement(g_, e);
	    phgDumpElement(g_, e1);
	    if (e1->children[0] != NULL)
		phgDumpElement(g_, e1->children[0]);
	    if (e1->children[1] != NULL)
		phgDumpElement(g_, e1->children[1]);
	    return conformity = FALSE;
	}
	/* check that the vertices of the common face (edge in 2D) match */
	for (j = 0; j < NVert; j++) {
	    INT v;
	    if (j == i)
		continue;
#define VERTEX(e, i)	\
    (g_->period == NULL ? e->verts[i] : g_->period->L2Gmap_vert[e->verts[i]])
	    v = VERTEX(e, j);
	    for (k = 0; k < NFace; k++) {
		if (v == VERTEX(e1, k))
		    break;
	    }
#undef VERTEX
	    if (k >= NFace) {
		phgInfo(-1, "non conforming face:\n");
		phgDumpElement(g_, e);
		phgDumpElement(g_, e1);
		return conformity = FALSE;
	    }
	}
#if Dim == 3
	if (HasLocalNeighbour(e, i) && (g_->flags & FACE_FLAG)) {
	    j = phgOppositeFace(g_, e, i, e1 = e->neighbours[i]);
	    if (e->faces[i] != e1->faces[j]) {
		phgInfo(-1, "non conforming face indices: %"dFMT" <-> %"
			    dFMT"\n", i, j);
		phgDumpElement(g_, e);
		phgDumpElement(g_, e1);
		return conformity = FALSE;
	    }
	}
#endif
    }

    /* TODO: check edge indices */

    return TRUE;
}

static BOOLEAN
check_conformity_nonleaf_callback CB_ARGS(e)
{
    int fmap[NVert], emap[NEdge];
    ELEMENT *e0, *e1;

    /* TODO: check vertices */

    if ((e0 = e->children[0]) != NULL) {
	if ((g_->flags & ELEM_FLAG) && e->index != e0->index &&
	    GlobalVertexP(g_, e->verts[0]) < GlobalVertexP(g_, e->verts[1])) {
	    phgDumpElement(g_, e);
	    phgDumpElement(g_, e0);
	    phgInfo(-1, "%s:%d, index of child 0 != index of parent\n",
			__FILE__, __LINE__);
	    return conformity = FALSE;
	}

	if ((g_->flags & (EDGE_FLAG | FACE_FLAG)))
	    phgMapP2C(e, fmap, emap, 0);

	if ((g_->flags & EDGE_FLAG) &&
	    (e->edges[1] != e0->edges[emap[1]] ||
	     e->edges[2] != e0->edges[emap[2]] ||
	     e->edges[5] != e0->edges[emap[5]] ||
	     (GlobalVertexP(g_, e->verts[0]) < GlobalVertexP(g_, e->verts[1]) &&
	      e->edges[0] != e0->edges[emap[0]]))) {
	    phgDumpElement(g_, e);
	    phgDumpElement(g_, e0);
	    phgInfo(-1, "edge indices of parent and child 0 mismatch\n");
	    return conformity = FALSE;
	}

	if ((g_->flags & FACE_FLAG) &&
	    (e->faces[1] != e0->faces[3] ||
	     (GlobalVertexP(g_, e->verts[0]) < GlobalVertexP(g_, e->verts[1]) &&
	      (e->faces[2] != e0->faces[fmap[2]] ||
	       e->faces[3] != e0->faces[fmap[3]])))) {
	    phgDumpElement(g_, e);
	    phgDumpElement(g_, e0);
	    phgInfo(-1, "face indices of parent and child 0 mismatch\n");
	    return conformity = FALSE;
	}
    }

    if ((e1 = e->children[1]) != NULL) {
	if ((g_->flags & ELEM_FLAG) && e->index != e1->index &&
	    GlobalVertexP(g_, e->verts[0]) > GlobalVertexP(g_, e->verts[1])) {
	    phgDumpElement(g_, e);
	    phgDumpElement(g_, e1);
	    phgInfo(-1, "%s:%d, index of child 1 != index of parent\n",
			__FILE__, __LINE__);
	    return conformity = FALSE;
	}

	if ((g_->flags & (EDGE_FLAG | FACE_FLAG)))
	    phgMapP2C(e, fmap, emap, 1);

	if ((g_->flags & EDGE_FLAG) &&
	    (e->edges[3] != e1->edges[emap[3]] ||
	     e->edges[4] != e1->edges[emap[4]] ||
	     e->edges[5] != e1->edges[emap[5]] ||
	     (GlobalVertexP(g_, e->verts[0]) > GlobalVertexP(g_, e->verts[1]) &&
	      e->edges[0] != e1->edges[emap[0]]))) {
	    phgDumpElement(g_, e);
	    phgDumpElement(g_, e1);
	    phgInfo(-1, "edge indices of parent and child 1 mismatch\n");
	    return conformity = FALSE;
	}

	if ((g_->flags & FACE_FLAG) &&
	    (e->faces[0] != e1->faces[3] ||
	     (GlobalVertexP(g_, e->verts[0]) > GlobalVertexP(g_, e->verts[1]) &&
	      (e->faces[2] != e1->faces[fmap[2]] ||
	       e->faces[3] != e1->faces[fmap[3]])))) {
	    phgDumpElement(g_, e);
	    phgDumpElement(g_, e1);
	    phgInfo(-1, "face indices of parent and child 1 mismatch\n");
	    return conformity = FALSE;
	}
    }

    return TRUE;
}

#if USE_MPI
static int
get_edge_no(int v0, int v1)
/* returns local number of an edge, given its two local vertices */
{
    int v;
 
    assert(v0 != v1);
    if (v0 > v1) {
	v = v0;
	v0 = v1;
	v1 = v;
    }

    switch (v0 * 10 + v1) {
	case 01:
	    return 0;
	case 02:
	    return 1;
	case 03:
	    return 2;
	case 12:
	    return 3;
	case 13:
	    return 4;
	case 23:
	    return 5;
    }

    return -1;
}

static INT (*dup_list)[NVert] = NULL;

static BOOLEAN
count_duplicate_elements(ELEMENT *e)
{
    if (IsLeaf(e)) {
	e->flag = 0;
	return FALSE;
    }

    if (e->children[0] != NULL)
	e->flag = count_duplicate_elements(e->children[0]);
    if (e->children[1] != NULL)
	e->flag |= count_duplicate_elements(e->children[1]);

    if ((e->children[0] == NULL || e->children[1] == NULL))
	e->flag = TRUE;

    if (e->flag)
	index_++;

    return e->flag;
}

static BOOLEAN
collect_duplicate_elements CB_ARGS(e)
{
    if (e->flag) {
	dup_list[index_][0] = GlobalVertexP(g_, e->verts[0]);
	dup_list[index_][1] = GlobalVertexP(g_, e->verts[1]);
	dup_list[index_][2] = GlobalVertexP(g_, e->verts[2]);
	dup_list[index_][3] = GlobalVertexP(g_, e->verts[3]);
	index_++;
    }
    return TRUE;
}

static int
duplicate_comp(const void *p0, const void *p1)
{
    INT i;
    const INT *v0 = p0, *v1 = p1;
    if ((i = *(v0++) - *(v1++)) != 0)
	return i > 0 ? 1 : -1;
    if ((i = *(v0++) - *(v1++)) != 0)
	return i > 0 ? 1 : -1;
    if ((i = *(v0++) - *(v1++)) != 0)
	return i > 0 ? 1 : -1;
    return (i = *v0 - *v1) > 0 ? 1 : (i < 0 ? -1 : 0);
}
#endif

#if DEBUG && USE_MPI
static void
check_L2Gmap(INT *map, INT n, INT nglobal, const char *file, int line)
{
    INT *buf;
    INT i;

    if (map == NULL)
	return;

    for (i = 0; i < n; i++) {
	INT j = map[i];
	if (j < 0 || j >= nglobal)
	    phgError(1, "%s:%d, invalid entry: map[%"dFMT"]=%"dFMT", "
			"nglobal=%"dFMT"\n", file, line, i, j, nglobal);
    }

    if (n <= 1)
	return;
    buf = phgAlloc(n * sizeof(*map));
    memcpy(buf, map, n * sizeof(*map));
    qsort(buf, n, sizeof(*map), phgCompINT);

    for (i = 0; i < n - 1; i++) {
	if (buf[i] == buf[i + 1])
	    phgError(1, "%s:%d, duplicate value: %"dFMT"\n",
			file, line, buf[i]);
    }
    free(buf);
}

# define CheckL2Gmap(map, n, ng) check_L2Gmap(map, n, ng, __FILE__, __LINE__)
#else	/* DEBUG && USE_MPI */
# define CheckL2Gmap(map, n, ng)
#endif	/* DEBUG && USE_MPI */

void
phgCheckConformity(GRID *g)
{
    ParallelFunction(phgCheckConformity, g != NULL && g->nprocs > 1);

    if (g == NULL)
	return;

#if USE_MPI
    CheckL2Gmap(g->L2Gmap_vert, g->nvert, g->nvert_global);
    if ((g->flags & EDGE_FLAG))
	CheckL2Gmap(g->L2Gmap_edge, g->nedge, g->nedge_global);
    if ((g->flags & FACE_FLAG))
	CheckL2Gmap(g->L2Gmap_face, g->nface, g->nface_global);
    if ((g->flags & ELEM_FLAG))
	CheckL2Gmap(g->L2Gmap_elem, g->nelem, g->nelem_global);
#endif

    g_ = g;
    conformity = TRUE;
    phgTraverseElements(g, check_conformity_callback);

    if (conformity)
	phgTraverseAllElements(g, check_conformity_nonleaf_callback);

#if USE_MPI
    if (g->nprocs > 1) {
	int c0 = conformity, c;
	MPI_Allreduce(&c0, &c, 1, MPI_INT, MPI_SUM, g->comm);
	conformity = (c == g->nprocs);
    }
#endif

    if (!conformity) {
	phgError(1, "local non conformity.\n");
    }

#if USE_MPI
    if (conformity && g->nprocs > 1) {
	struct {
	    RNEIGHBOUR rn;	/* the RNEIGHBOUR data */
	    INT verts[NFace - 1]; /* global indices of the 3 shared vertices */
	    INT edges[NEdge];
#if Dim == 3
	    INT face_index;
#endif
	} *local, *remote, *p;
	RNEIGHBOUR *rlist = phgAlloc(g->neighbours.count * sizeof(*rlist));
	MPI_Datatype type;
	INT *verts;
	int i, j, k;
	int *counts = g->neighbours.counts, *displs = g->neighbours.displs;

	local = phgAlloc(g->neighbours.count * sizeof(*local));
	remote = phgAlloc(g->neighbours.count * sizeof(*remote));

	for (i = 0; i < g->neighbours.count; i++) {
	    int vertex;
	    ELEMENT *e;

	    local[i].rn = g->neighbours.list[i];
	    vertex = local[i].rn.vertex;
	    e = local[i].rn.local;
	    verts = e->verts;
	    for (k = 0, j = 0; j < NFace; j++) {
		if (j == vertex) continue;
		local[i].verts[k++] = GlobalVertexP(g, verts[j]);
	    }
	    local[i].rn.depth = 0;
	    for (j = 0; j < NEdge; j++) {
		local[i].edges[j] = GlobalEdge(g, e->edges[j]);
	    }
#if Dim == 3
	    local[i].face_index = GlobalFace(g, e->faces[vertex]);
#endif
	}

	MPI_Type_contiguous(sizeof(*local), MPI_BYTE, &type);
	MPI_Type_commit(&type);
	MPI_Alltoallv(local,  counts, displs, type,
		      remote, counts, displs, type, g->comm);
	MPI_Type_free(&type);

	for (i = 0; i < g->neighbours.count && conformity; i++) {
	    RNEIGHBOUR *rn = &remote[i].rn;
	    RNEIGHBOUR *rn0 = GetRNeighbour(g, rn->remote, rn->op_vertex);
	    p = remote + i;
	    if (!IsLeaf(rn->remote)) {
		phgWarning("element linked to non-leaf element.\n");
		phgWarning("remote data: local=%p, remote=%p, vertex=%d, "
			   "op=%d, rface=%d %d %d\n", rn->remote, rn->local,
			   rn->vertex, rn->op_vertex, rn->rface[0],
			   rn->rface[1], rn->rface[2]);
		phgWarning("local  data: local=%p, remote=%p, vertex=%d, "
			   "op=%d, rface=%d %d %d\n", rn0->remote, rn0->local,
			   rn0->vertex, rn0->op_vertex, rn0->rface[0],
			   rn0->rface[1], rn0->rface[2]);
		break;
	    }
	    if (local[i].rn.depth != 0) {
		phgWarning("different entries linked to same element.\n");
		break;
	    }
	    local[i].rn.depth = 1;
	    if (rn0->vertex != rn->op_vertex || rn0->op_vertex != rn->vertex ||
		rn0->local != rn->remote || rn0->remote != rn->local ||
		rn->rank != g->rank) {
		conformity = FALSE;
		phgWarning("global non conformity: linked data don't match.\n");
		break;
	    }
	    /* check whether the vertices match */
	    verts = rn0->local->verts;
	    for (j = 0; j < (NFace - 1) && conformity; j++) {
		INT v = GlobalVertexP(g, verts[rn->rface[j]]);
		if (v == p->verts[j])
		    continue;
		conformity = FALSE;
		phgWarning("global non conformity: vertices don't match.\n");
		phgWarning("local=(%p %d %d%d%d), remote=(%p %d %d%d%d)\n",
			rn0->local, rn0->vertex, rn0->rface[0], rn0->rface[1],
			rn0->rface[2], rn0->remote, rn0->op_vertex,
			rn->rface[0], rn->rface[1], rn->rface[2]);
		phgWarning("remote vertices=%"dFMT" %"dFMT" %"dFMT"\n",
			p->verts[0], p->verts[1], p->verts[2]);
		v = verts[rn->rface[j]];
		phgWarning("unmatched vertices: (%lf %lf %lf)\n",
			(double)g->verts[v][0], (double)g->verts[v][1],
			(double)g->verts[v][2]);
		phgDumpElement(g, rn0->local);
		break;
	    }
	    if ((g->flags & EDGE_FLAG)) {
	      /* check whether the edges match */
	      int a = 0, b = 1;
	      for (j = 0; j < NVert - 1 && conformity; j++) {
		if (j == rn0->vertex) continue;
		for (k = j + 1; k < NVert; k++) {
		    int u, v;
		    u = rn0->rface[a];
		    v = rn0->rface[b];
		    if (k == rn0->vertex) continue;
		    if (GlobalEdge(g, rn0->local->edges[get_edge_no(j, k)]) !=
			remote[i].edges[get_edge_no(u, v)]) {
			conformity = FALSE;
			phgWarning("non conformity: edges don't match.\n");
			phgWarning("unmatched edge: %d-%d %d <==> %d-%d %d\n",
				j, k, get_edge_no(j, k),
				u, v, get_edge_no(u, v));
			phgWarning("remote edges: %d %d %d %d %d %d\n",
				remote[i].edges[0], remote[i].edges[1],
				remote[i].edges[2], remote[i].edges[3],
				remote[i].edges[4], remote[i].edges[5]);
			phgWarning("remote vertices=%"dFMT" %"dFMT" %"dFMT"\n",
				p->verts[0], p->verts[1], p->verts[2]);
			phgDumpElement(g, rn0->local);
			break;
		    }
		    if (++b >= 3) {
			a++;
			b = a + 1;
		    }
		}
	      }
	    }
#if Dim == 3
	    if ((g->flags & FACE_FLAG)) { /* check whether the faces match */
		ELEMENT *e = remote[i].rn.remote;
		int vertex = remote[i].rn.op_vertex;
		if (GlobalFace(g, e->faces[vertex]) != remote[i].face_index) {
		    phgWarning("non conformity: face indices don't match.\n");
		    phgWarning("local(%p,%d)=%"dFMT", remote(%p,%d)=%"dFMT"\n",
			e, vertex, GlobalFace(g, e->faces[vertex]),
			remote[i].rn.local, remote[i].rn.vertex,
			remote[i].face_index);
		    break;
		}
	    }
#endif
	}

	for (i = 0; i < g->neighbours.count && conformity; i++) {
	    if (local[i].rn.depth == 1) continue;
	    phgWarning("global non conformity: unmatched link(s)!\n");
	    break;
	}

	phgFree(rlist);
	phgFree(local);
	phgFree(remote);
	
	i = conformity;
	MPI_Allreduce(&i, &j, 1, MPI_INT, MPI_MIN, g->comm);
	conformity = j;
    }

    /* check vertex indices on duplicate non-leaf elements */
    if (FALSE && conformity && g->nprocs > 1) {
	INT k, l, *counts, nmax;
	INT (*buffer)[NVert], (*p0)[NVert], (*p1)[NVert];
	MPI_Status status;
	MPI_Datatype type;
	int i, j, m, src, dst;

	index_ = 0;
	for (k = 0; k < g->nroot; k++)
	    count_duplicate_elements(g->roots + k);
	dup_list = phgAlloc(index_ * sizeof(*dup_list));
	counts = phgAlloc(g->nprocs * sizeof(*counts));
	MPI_Allgather(&index_, sizeof(index_), MPI_BYTE,
			counts, sizeof(index_), MPI_BYTE, g->comm);
	nmax = 0;
	for (i = 0; i < g->nprocs; i++)
	    if (nmax < counts[i])
		nmax = counts[i];
	buffer = phgAlloc(2 * nmax * sizeof(*dup_list));
	index_ = 0;
	phgTraverseAllElements(g_, collect_duplicate_elements);
	qsort(dup_list, index_, sizeof(*dup_list), duplicate_comp);
	/* circulate and compare the lists */
	j = g->rank;
	dst = (g->rank == 0) ? g->nprocs - 1 : g->rank - 1;
	src = (g->rank == g->nprocs - 1) ? 0 : g->rank + 1;
	MPI_Type_contiguous(sizeof(*buffer), MPI_BYTE, &type);
	MPI_Type_commit(&type);
	p0 = dup_list;
	p1 = buffer;
	for (i = 0; i < g->nprocs - 1; i++) {
	    p1 = buffer + (nmax - (p1 - buffer));
	    if ((l = j + 1) >= g->nprocs)
		l = 0;
	    MPI_Sendrecv(p0, counts[j], type, dst, 111,
			 p1, counts[l], type, src, 111, g->comm, &status);
	    j = l;
	    p0 = p1;
	    /* compare dup_list with p1 */
	    k = 0;
	    l = 0;
	    while (TRUE) {
		while (k < index_ && dup_list[k][0] < 0)
		    k++;
		if (k >= index_)
		    break;
		while (l < counts[j]) {
		    if ((m = duplicate_comp(p0 + l, dup_list + k)) >= 0)
			break;
		    l++;
		}
		if (l >= counts[j])
		    break;
		while (m > 0 && k < index_) {
		    if (dup_list[k][0] >=0 &&
			(m = duplicate_comp(p0 + l, dup_list + k)) <= 0)
			break;
		    k++;
		}
		if (k >= index_)
		    break;
		/* a matching pair of tetrahedra found */
		dup_list[k][0] = dup_list[k][1] = dup_list[k][2]
							= dup_list[k][3] = -1;
		k++;
		l++;
	    }
	}
	MPI_Type_free(&type);
	phgFree(buffer);
	phgFree(counts);
	/* check whether we have unmatched tetrahedra */
	for (k = 0; k < index_; k++) {
	    if (dup_list[k][0] >= 0)
		phgError(1, "%s:%d, tetrahedra: "
			"%"dFMT" %"dFMT" %"dFMT" %"dFMT"\n",
			__FILE__, __LINE__, dup_list[k][0],
			dup_list[k][1], dup_list[k][2], dup_list[k][3]);
	}
	phgFree(dup_list);
    }
#endif	/* USE_MPI */

    if (conformity) {
	if (IS_RANK0) phgInfo(1, "the grid is conforming.\n");
    }
    else {
#if 0
	phgDumpGrid(g);
#endif
	phgError(1, "%s:%d: abort.", __FILE__, __LINE__);
    }

    return;
}

#if 0
static int
elem_comp(const void *p0, const void *p1)
{
    return phgMemcmp((const ELEMENT *)p0, (const ELEMENT *)p1, sizeof(ELEMENT));
}
#endif

/* Note: enabling REMOVE_DEGENERATED_ELEMENTS causes PHG to crash
 * on some bad input meshes (e.g., cask1.mesh, cask.mesh, etc.)
 *
 * It is simpler to move the code to before updating neighbour links,
 * but the present code allows transferring boundary types of deleted
 * elements. */
#define REMOVE_DEGENERATED_ELEMENTS 1

void
phgRemoveDegeneratedElements(GRID *g)
/* checks and removes degenerated elements and unused vertices.
 *
 * Note: this function should be called after updating the neighbour links,
 *	 necessary for transfering boundary types of deleted elements */
{
    INT *vmap;
    INT i, j;
    ELEMENT *e;
#if REMOVE_DEGENERATED_ELEMENTS
    INT k, u, v, w;
    ELEMENT *e0, *e1;
    int a, b, c;
    BOOLEAN flag;
#endif	/* REMOVE_DEGENERATED_ELEMENTS */

    FunctionEntry;

#if USE_MPI
    assert(g->nprocs <= 1);
    assert(g->L2Gmap_vert == NULL);
#endif
    assert(g->serial_no == 0);

    if (g->nvert <= 0)
	Return;

    vmap = phgAlloc(g->nvert * sizeof(*vmap));

    /* check and remove tetrahedra with zero length edges.
     * TODO: check for other types of degenerated tetrahedra */

#if REMOVE_DEGENERATED_ELEMENTS
    for (j = 0; j < g->nvert; j++)
	vmap[j] = j;
    flag = FALSE;
    k = 0;
    for (j = 0; j < g->nroot; j++) {
	e = g->roots + j;
	for (i = 0; i < NEdge; i++) {
	    u = e->verts[a = GetEdgeVertex(i, 0)];
	    v = e->verts[b = GetEdgeVertex(i, 1)];
	    if (g->verts[u][0] == g->verts[v][0] &&
		g->verts[u][1] == g->verts[v][1] &&
		g->verts[u][2] == g->verts[v][2])
		break;
	}
	if (i >= NEdge) {
	    /* no zero length edge in this element */
	    e->index = k++;
	    continue;
	}
	e->index = -1;	/* mark e as deleted */
	/* mark the vertex with larger value as deleted */
	w = (u > v) ? v : u;
	while (vmap[w] != w)
	    w = vmap[w];
	vmap[u] = vmap[v] = w;
	if (flag)
	    continue;
	/* update neighbour links and boundary types of the neighbours */
	e0 = e->neighbours[a];
	e1 = e->neighbours[b];
	if (e0 != NULL) {
	    for (c = 0; c < NFace; c++)
		if (e0->neighbours[c] == e)
		    break;
	    if (c >= NVert) {
		flag = TRUE;
		continue;
	    }
	    e0->neighbours[c] = e1;
	    e0->bound_type[c] = e->bound_type[b];
	}
	if (e1 != NULL) {
	    for (c = 0; c < NFace; c++)
		if (e1->neighbours[c] == e)
		    break;
	    if (c >= NVert) {
		flag = TRUE;
		continue;
	    }
	    e1->neighbours[c] = e0;
	    e1->bound_type[c] = e->bound_type[a];
	}
    }

    /* update vmap[] */
    if (k < g->nroot) {
	for (j = 0; j < g->nvert; j++) {
	    v = j;
	    while (v != vmap[v])
		v = vmap[v];
	    vmap[j] = v;
	}
    }

    /* remove tetrahedra with duplicate vertices */
    k = 0;
    for (j = 0; j < g->nroot; j++) {
	INT verts[NVert];
	e = g->roots + j;
	if (e->index == -1)
	    continue;
	for (i = 0; i < NVert; i++)
	    verts[i] = e->verts[i] = vmap[e->verts[i]];
	qsort(verts, NVert, sizeof(*verts), phgCompINT);
	if (verts[0] < 0 || verts[NVert - 1] >= g->nvert) {
	    flag = TRUE;
	    e->index = -1;
	    continue;
	}
	u = verts[0];
	for (i = 1; i < NVert; i++) {
	    v = verts[i];
	    if (u == v)
		break;
	    u = v;
	}
	if (i < NVert) {
	    flag = TRUE;
	    e->index = -1;
	    continue;
	}
	e->index = k++;
    }

    if (k < g->nroot) {
	/* some elements are to be deleted, update neighbours and vertices */
	i = g->nroot - k;
	phgInfo(2, "deleting %"dFMT" degenerated element%s\n",
			i, i > 1 ? "s" : "");
	for (j = 0; j < g->nroot; j++) {
	    e = g->roots + j;
	    if (e->index == -1)
		continue;
	    for (i = 0; i < NVert; i++) {
		if ((e0 = e->neighbours[i]) != NULL) {
		    if (!flag) {
			assert(e0->index >= 0);
			e->neighbours[i] = g->roots + (e0->index);
		    }
		    else {
			e->neighbours[i] = NULL;
		    }
		}
	    }
	}
	for (j = k = 0; j < g->nroot; j++) {
	    if (g->roots[j].index == -1)
		continue;
	    if (j != k)
		g->roots[k] = g->roots[j];
	    k++;
	}
	g->ntree = g->nleaf = g->nelem = g->nroot = g->nelem_global = k;
	if (flag) {
	    phgWarning("bdry type of deleted elements may not be "
			"correctly transferred.\n");
	    phgUpdateNeighbours(g);
	}
    }
#endif	/* REMOVE_DEGENERATED_ELEMENTS */

    /* remove unused vertices, they may either come from the imported file
     * or from the deleted elements */

    memset(vmap, 0, g->nvert * sizeof(*vmap));
    for (j = 0; j < g->nroot; j++) {
	e = g->roots + j;
	for (i = 0; i < NVert; i++)
	    vmap[e->verts[i]] = 1;
    }

    /* compute new index in vmap and update g->verts */
    j = 0;
    for (i = 0; i < g->nvert; i++) {
	if (vmap[i] == 0) {
	    /* vertex unused */
	    vmap[i] = -1;
	    continue;
	}
	vmap[i] = j;
	if (j != i)
	    memcpy(g->verts + j, g->verts + i, sizeof(COORD));
	j++;
    }

    if (j == g->nvert) {
	/* no unused vertices */
	phgFree(vmap);
	Return;
    }

    i = g->nvert - j;
    phgInfo(2, "deleting %"dFMT" unused vert%s\n", i, i > 1 ? "ices" : "ex");
    g->nvert = g->nvert_global = j;
    if (j == 0)
	g->verts = NULL;
    else {
	g->verts = phgRealloc_(g->verts, j * sizeof(*g->verts),
					 j * sizeof(*g->verts));
	for (j = 0; j < g->nroot; j++) {
	    e = g->roots + j;
	    for (i = 0; i < NVert; i++)
		e->verts[i] = vmap[e->verts[i]];
	}
    }

    phgFree(vmap);
    Return;
}

static ELEMENT *
find_edge_neighbour(GRID *g, ELEMENT *e, INT u0, INT u1, INT *u)
/* Returns the neighbour not containing the face {u0,u1,u}, but containing
 * the edge u0-u1, where u0, u1, *u are vertex indices. On return
 * *u is set to the element index of the third vertex of the common face. */
{
    int i, j;
    INT *v, verts[Dim + 1];

    if (g->period != NULL) {
	v = verts;
	for (i = 0; i < Dim + 1; i++)
	    v[i] = g->period->L2Gmap_vert[e->verts[i]];
    }
    else {
	v = e->verts;
    }

    /* one of the vertices must equal to "*u" */
    if (v[i = 0] != *u && v[++i] != *u && v[++i] != *u && v[++i] != *u) {
	phgError(-1, "%s:%d, unexpected error, u0 u1 = %d %d, *u = %"dFMT", "
		     "v = %"dFMT" %"dFMT" %"dFMT" %"dFMT" (%"dFMT" %"dFMT" %"dFMT" %"dFMT")\n", __FILE__, __LINE__,
		     u0, u1, *u, v[0], v[1], v[2], v[3],
		     e->verts[0], e->verts[1], e->verts[2], e->verts[3]);
    }

    for (j = 0; j == i || v[j] == u0 || v[j] == u1; j++);
#if DEBUG
    if (j >= Dim + 1) {
	phgError(-1, "%s:%d, unexpected error, u0 u1 = %d %d, *u = %"dFMT", "
		     "v = %"dFMT" %"dFMT" %"dFMT" %"dFMT" (%"dFMT" %"dFMT" %"dFMT" %"dFMT")\n", __FILE__, __LINE__,
		     u0, u1, *u, v[0], v[1], v[2], v[3],
		     e->verts[0], e->verts[1], e->verts[2], e->verts[3]);
    }
#endif	/* DEBUG */
    *u = v[j];

    /* it is assumed that e->neighbours[i] == NULL for all non-interior and
     * REMOTE neighbours */

    return IsLeaf(e) ? GetNeighbour(e, i) : e->neighbours[i];
}

int
phgGetPatch(GRID *g, ELEMENT *e, BYTE v0, BYTE v1, ELEMENT ***p)
/* constructs the list of elements constituting the patch around
 * the edge v0-v1 (local indices) of e. returns number of elements
 * in the patch and list of elements through 'p' */
/* TODO: patch composed of several non contiguous parts */
{
    INT u0 = e->verts[v0], u1 = e->verts[v1], u;
    /* gcc is unhappy if v2, v3 are not initialized :( */
    BYTE v2 = 0, v3 = 0;
    int count = 0, size;
    ELEMENT *e1;

    if (g->period != NULL) {
	u0 = g->period->L2Gmap_vert[u0];
	u1 = g->period->L2Gmap_vert[u1];
    }

    *p = phgAlloc((size = 6) * sizeof(**p));
    (*p)[count++] = e;

    /* find out the two opposite vertices, v2 and v3, of the edge */
    switch (v0 * 4 + v1) {
	case 0 * 4 + 1: v2 = 2; v3 = 3; break;
	case 0 * 4 + 2: v2 = 1; v3 = 3; break;
	case 0 * 4 + 3: v2 = 1; v3 = 2; break;
	case 1 * 4 + 2: v2 = 0; v3 = 3; break;
	case 1 * 4 + 3: v2 = 0; v3 = 2; break;
	case 2 * 4 + 3: v2 = 0; v3 = 1; break;
    }
    
    /* Note: refine_patch depends on the direction, don't change! */

    /* search in first direction */
    e1 = IsLeaf(e) ? GetNeighbour(e, v2) : e->neighbours[v2];
    u = e->verts[v3];
    if (g->period != NULL)
	u = g->period->L2Gmap_vert[u];
    while (e1 != e && e1 != NULL) {
	if (count >= size) {
	    if (count > 100)
		phgError(1, "%s:%d, unexpected error.\n", __FILE__, __LINE__);
	    *p = phgRealloc_(*p, (size + 6) * sizeof(**p), size * sizeof(**p));
	    size += 6;
	}
	(*p)[count++] = e1;
	e1 = find_edge_neighbour(g, e1, u0, u1, &u);
    }
    if (e1 == e) {
	phgDebug((4, "circular patch, size = %"dFMT"\n", count));
	return count;
    }
    /* search in second direction */
#if 0
{
    phgInt count0 = count;	/* save position */
#endif
    e1 = IsLeaf(e) ? GetNeighbour(e, v3) : e->neighbours[v3];
    u = e->verts[v2];
    if (g->period != NULL)
	u = g->period->L2Gmap_vert[u];
    /* Note: e1 may come back to e when dealing with non-leaf elements
     * on a distributed mesh */
    while (e1 != e && e1 != NULL) {
	if (count >= size) {
	    if (count > 100)
		phgError(1, "%s:%d, unexpected error.\n", __FILE__, __LINE__);
	    *p = phgRealloc_(*p, (size + 6) * sizeof(**p), size * sizeof(**p));
	    size += 6;
	}
	(*p)[count++] = e1;
	e1 = find_edge_neighbour(g, e1, u0, u1, &u);
#if DEBUG
	if (e1 == e) {
	    int i;
	    phgInfo(-1, "unexpected case: e1 == e,\n");
	    for (i = 0; i < count; i++)
		phgDumpElement(g, (*p)[i]);
	    phgError(1, "abort.\n");
	}
#endif
    }

#if 0
    /* We have the list of elements around the patch in the following order:
     *		*p[0] -> *p[1] -> ... -> *p[count0-1]
     *		*p[0] <- *p[count0] <- ... <- *p[count-1]
     *
     * We first commute (*p[0], *p[count0-1]), (*p[1], *p[count0-2]), etc.,
     * the we commute (*p[0], *p[count-1]), (*p[1], *p[count-2]), etc.,
     * such that the elements in the list are chained in the right direction */
    for (u0 = 0, u1 = count0 - 1; u0 < u1; u0++, u1--) {
	e = (*p)[u0];
	(*p)[u0] = (*p)[u1];
	(*p)[u1] = e;
    }
    for (u0 = 0, u1 = count - 1; u0 < u1; u0++, u1--) {
	e = (*p)[u0];
	(*p)[u0] = (*p)[u1];
	(*p)[u1] = e;
    }
}
#endif

    phgDebug((4, "open patch, size = %"dFMT"\n", count));

    return count;    
}

#if USE_MPI
/* functions for managing global edge numbers */

static INT nedge_save, nedge_global_save;

static struct {
    ELEMENT *e;		/* first element on the patch */
    BYTE v0, v1;	/* vertices of the edge */
    BYTE shared;	/* whether the edge is shared */
} *edges_local = NULL;

static INT edges_local_count = 0, edges_local_allocated = 0;

static struct {
    INT index;		/* index to edges_local */
    INT v0, v1;		/* global indices of the two vertices */
} *edges_shared = NULL, *edges_gathered = NULL;
static INT edges_shared_count = 0, edges_shared_allocated = 0;

static void
add_local_edge(ELEMENT *e, int v0, int v1, BOOLEAN shared)
{
    if (edges_local_count >= edges_local_allocated) {
	edges_local = phgRealloc_(edges_local,
		sizeof(*edges_local) * (edges_local_allocated + 4096),
		sizeof(*edges_local) * edges_local_allocated);
	edges_local_allocated += 4096;
    }
    edges_local[edges_local_count].e = e;
    edges_local[edges_local_count].v0 = v0;
    edges_local[edges_local_count].v1 = v1;
    edges_local[edges_local_count++].shared = shared;
}

static int
comp_edge0(const void *p0, const void *p1)
{
    INT i = *((const INT *)p0), j = *((const INT *)p1);
    INT v0 = edges_local[i].e->verts[edges_local[i].v0];
    INT v1 = edges_local[i].e->verts[edges_local[i].v1];
    INT u0 = edges_local[j].e->verts[edges_local[j].v0];
    INT u1 = edges_local[j].e->verts[edges_local[j].v1];
    INT c;

    if (v0 > v1) {
	c = v0; v0 = v1; v1 = c;
    }

    if (u0 > u1) {
	c = u0; u0 = u1; u1 = c;
    }

    if ((c = v0 - u0) != 0)
	return (c > 0) ? 1 : -1;
    return ((c = v1 - u1) > 0) ? 1 : (c < 0) ? -1 : 0;
}

static int
comp_edge(const void *p0, const void *p1)
{
    int i;

    if ((i = comp_edge0(p0, p1)) != 0)
	return i;
    return edges_local[*((const INT *)p1)].shared -
	   edges_local[*((const INT *)p0)].shared;
}

static void
add_shared_edge(INT index, INT v0, INT v1)
{
    if (edges_shared_count >= edges_shared_allocated) {
	edges_shared = phgRealloc_(edges_shared,
		sizeof(*edges_shared) * (edges_shared_allocated + 4096),
		sizeof(*edges_shared) * edges_shared_allocated);
	edges_shared_allocated += 4096;
    }

    if (v0 < v1) {
	edges_shared[edges_shared_count].v0 = v0;
	edges_shared[edges_shared_count].v1 = v1;
    }
    else {
	edges_shared[edges_shared_count].v0 = v1;
	edges_shared[edges_shared_count].v1 = v0;
    }
    edges_shared[edges_shared_count++].index = index;
}

static int
comp_edge_index(const void *p0, const void *p1)
{
    INT i = *((const INT *)p0), j = *((const INT *)p1);
    INT v0 = edges_gathered[i].v0, v1 = edges_gathered[i].v1;
    INT u0 = edges_gathered[j].v0, u1 = edges_gathered[j].v1;
    INT c;

    if ((c = v0 - u0) != 0)
	return (c > 0) ? 1 : -1;
    return ((c = v1 - u1) > 0) ? 1 : ((c < 0) ? -1 : 0);
}

static void
update_edge_no(INT index, INT no)
/* update an edge number */
{
    int j, k, n;
    ELEMENT **p, *e;
    INT v0, v1;

    e = edges_local[index].e;
    if (g_->period != NULL) {
	v0 = g_->period->L2Gmap_vert[e->verts[edges_local[index].v0]];
	v1 = g_->period->L2Gmap_vert[e->verts[edges_local[index].v1]];
    }
    else {
	v0 = e->verts[edges_local[index].v0];
	v1 = e->verts[edges_local[index].v1];
    }
    n = phgGetPatch(g_, e, edges_local[index].v0, edges_local[index].v1, &p);
    for (j = 0; j < n; j++) {
	e = p[j];
	for (k = 0; k < NEdge; k++) {
	    INT u0, u1;
	    if (g_->period != NULL) {
		u0 = g_->period->L2Gmap_vert[e->verts[GetEdgeVertex(k, 0)]];
		u1 = g_->period->L2Gmap_vert[e->verts[GetEdgeVertex(k, 1)]];
	    }
	    else {
		u0 = e->verts[GetEdgeVertex(k, 0)],
		u1 = e->verts[GetEdgeVertex(k, 1)];
	    }
	    if ((v0 == u0 && v1 == u1) || (v0 == u1 && v1 == u0))
		break;
	}
	assert(k < NEdge);
	e->edges[k] = no;
	phgDebug((4, "update_edge_no: e=%p, edge=%d, new=%"dFMT"\n", e, k, no));
    }
    phgFree(p);
}

static void
sync_edges(void)
/* update saved edge numbers. */
{
    int flag, *counts, *displs;
    INT tmp[2], *index, *nbuffer, i, j, local_count, total_local, total_shared;
    INT edge_no_local = 0, edge_no_global = 0;
    MPI_Datatype type;
    void *p;

    /* sort indices of edges_local */
    index = NULL;
    if (edges_local_count > 0) {
	index = phgAlloc(edges_local_count * sizeof(*index));
	for (i = 0; i < edges_local_count; i++)
	    index[i] = i;
	qsort(index, edges_local_count, sizeof(*index), comp_edge);
    }

    /* count # of private edges and build list of shared edges */
    i = local_count = 0;
    while (i < edges_local_count) {
	if (edges_local[j = index[i]].shared) {
 	    add_shared_edge(j,
		GlobalVertexP(g_, edges_local[j].e->verts[edges_local[j].v0]),
		GlobalVertexP(g_, edges_local[j].e->verts[edges_local[j].v1]));
	    i++;
	    continue;
	}
	j = i;
	while (++i < edges_local_count
		&& comp_edge0(index + j, index + i) == 0);
	local_count++;
    }

    phgDebug((3, "edges: local=%d (%d), shared=%d, total=%d\n", local_count,
			edges_local_count,
			edges_shared_count, local_count + edges_shared_count));

    nbuffer = phgAlloc(2 * g_->nprocs * sizeof(*nbuffer));
    counts = phgAlloc(2 * g_->nprocs * sizeof(*counts));;
    displs = counts + g_->nprocs;

    p = g_->L2Gmap_edge;
    g_->L2Gmap_edge = phgRealloc_(g_->L2Gmap_edge, (nedge_save + local_count +
			edges_shared_count) * sizeof(*(g_->L2Gmap_edge)),
		(p == NULL ? 0 : nedge_save) * sizeof(*(g_->L2Gmap_edge)));
    if (p == NULL) {
	for (i = 0; i < nedge_save; i++)
	    g_->L2Gmap_edge[i] = i;
    }

    tmp[0] = local_count;
    tmp[1] = edges_shared_count;
    MPI_Allgather(tmp, 2, PHG_MPI_INT, nbuffer, 2, PHG_MPI_INT, g_->comm);

    edge_no_local = nedge_save;
    total_local = total_shared = 0;
    for (i = 0; i < g_->nprocs; i++) {
	if (i == g_->rank) edge_no_global = total_local + nedge_global_save;
	total_local += nbuffer[2 * i];
	assert(nbuffer[2 * i + 1] <= INT_MAX);
	total_shared += (counts[i] = nbuffer[2 * i + 1]);
    }

    j = 0;
    for (i = 0; i < g_->nprocs; i++) {
	displs[i] = j;
	assert(j <= INT_MAX - counts[i]);
	j += counts[i];
    }

    /* update private edges */
    phgDebug((3, "private edges range: %"dFMT"-%"dFMT"\n", edge_no_global,
		edge_no_global + local_count - 1));
    i = 0;
    while (i < edges_local_count) {
	if (edges_local[index[i]].shared) {
	    i++;
	    continue;
	}
	j = i;
	while (++j < edges_local_count
		&& comp_edge0(index + i, index + j) == 0);
	while (i < j)
	    update_edge_no(index[i++], edge_no_local);
	g_->L2Gmap_edge[edge_no_local++] = edge_no_global++;
    }

    if (index != NULL)
	phgFree(index);

    edge_no_global = nedge_global_save + total_local;

    if (total_shared == 0)
	goto end;

    /* update shared edges */
    edges_gathered = phgAlloc(total_shared * sizeof(*edges_gathered));
    index = phgAlloc(total_shared * sizeof(*index));
    MPI_Type_contiguous(sizeof(*edges_gathered), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Allgatherv(edges_shared, edges_shared_count, type,
		   edges_gathered, counts, displs, type, g_->comm);
    MPI_Type_free(&type);
    for (i = 0; i < total_shared; i++)
	index[i] = i;
    qsort(index, total_shared, sizeof(*index), comp_edge_index);

    phgDebug((3, "update shared edges\n"));
    i = 0;
    flag = 0;
    while (i < total_shared) {
	j = index[i];
	if (j >= displs[g_->rank] && j < displs[g_->rank] + counts[g_->rank]) {
	    flag = 1;
	    update_edge_no(edges_gathered[j].index, edge_no_local);
	}
	j = i++;
	if (i >= total_shared || comp_edge_index(index + i, index + j) > 0) {
	    if (flag) {
		g_->L2Gmap_edge[edge_no_local++] = edge_no_global;
		flag = 0;
	    }
	    edge_no_global++;
	}
    }

    phgFree(index);
    phgFree(edges_gathered);

end:
    g_->nedge = edge_no_local;
    g_->nedge_global = edge_no_global;
    phgDebug((3, "shared edges range: %"dFMT"-%"dFMT"\n",
		nedge_global_save + total_local, edge_no_global - 1));

    phgFree(nbuffer);
    phgFree(counts);

    edges_local_count = edges_shared_count = 0;
    if (edges_local_allocated > 0) {
	phgFree(edges_local);
	edges_local = NULL;
	edges_local_allocated = 0;
    }
    if (edges_shared_allocated > 0) {
	phgFree(edges_shared);
	edges_shared = NULL;
	edges_shared_allocated = 0;
    }
}

#endif	/* USE_MPI */

static BOOLEAN
update_edges_callback CB_ARGS(e)
{
    int i, j, k, count;
    ELEMENT **p, *e1;
#if USE_MPI
    BOOLEAN flag;	/* indicate whether a shared edge */
#endif	/* USE_MPI */

    /* loop on all edges of the element */
    for (i = 0; i < NEdge; i++) {
	int v0 = GetEdgeVertex(i, 0), v1 = GetEdgeVertex(i, 1);
	INT g0, g1;	/* global indices of the two vertices */
	if (e->edges[i] != -1)
	    continue;
	/* new edge, assign it to all elements sharing it */
	g0 = e->verts[v0];
	g1 = e->verts[v1];
	if (g_->period != NULL) {
	    g0 = g_->period->L2Gmap_vert[g0];
	    g1 = g_->period->L2Gmap_vert[g1];
	}
	/* loop on elements in the patch */
	count = phgGetPatch(g_, e, v0, v1, &p);
#if USE_MPI
	flag = FALSE;
#endif	/* USE_MPI */
	for (j = 0; j < count; j++) {
#if USE_MPI
	    int u0 = 0, u1 = 0;
#endif	/* USE_MPI */
	    /* loop on edges of p[j] to find the one matching current edge */
	    e1 = p[j];
	    for (k = 0; k < NEdge; k++) {
		INT h0, h1;
		if (e1->edges[k] != -1)
		    continue;
		h0 = e1->verts[GetEdgeVertex(k, 0)];
		h1 = e1->verts[GetEdgeVertex(k, 1)];
		if (g_->period != NULL) {
		    h0 = g_->period->L2Gmap_vert[h0];
		    h1 = g_->period->L2Gmap_vert[h1];
		}
		if ((h0 == g0 && h1 == g1) || (h0 == g1 && h1 == g0))
		    break;
	    }
	    if (k >= NEdge) {
#if 0
/* save all elements on the patch for debugging */
FILE *f = fopen("debug.mesh", "w+t");
ELEMENT *s;
fprintf(f, "MeshVersionFormatted 1\nDimension\n3\n\n");
fprintf(f, "Vertices\n%"dFMT"\n", g_->nvert);
for (i = 0; i < g_->nvert; i++)
    fprintf(f, "%lf %lf %lf 0\n",
		g_->verts[i][0], g_->verts[i][1], g_->verts[i][2]);
count = 0;
for (s = g_->roots; s < g_->roots + g_->nroot; s++)
    for (j = 0; j < NVert - 1; j++)
	for (k = j + 1; k < NVert; k++)
	    if ((s->verts[j] == g0 && s->verts[k] == g1) ||
		(s->verts[j] == g1 && s->verts[k] == g0))
		count++;
fprintf(f, "\nTetrahedra\n%"dFMT"\n", count);
for (s = g_->roots; s < g_->roots + g_->nroot; s++)
    for (j = 0; j < NVert - 1; j++)
	for (k = j + 1; k < NVert; k++)
	    if ((s->verts[j] == g0 && s->verts[k] == g1) ||
		(s->verts[j] == g1 && s->verts[k] == g0))
		fprintf(f, "%"dFMT" %"dFMT" %"dFMT" %"dFMT" 0\n",
				s->verts[0] + 1, s->verts[1] + 1,
				s->verts[2] + 1, s->verts[3] + 1);
fprintf(f, "\nEnd\n");
fclose(f);
#endif
		phgInfo(-1, "===== element 'e'\n");
		phgDumpElement(g_, e);
		phgInfo(-1, "===== element 'e1'\n");
		phgDumpElement(g_, e1);
		phgInfo(-1, "===== the edge: %"dFMT" %"dFMT"\n", g0, g1);
		phgError(1, "%s:%d, unexpected error.\n", __FILE__, __LINE__);
	    }
	    else {
		e1->edges[k] = g_->nedge;
	    }
#if USE_MPI
	    if (g_->nprocs > 1) {
		/* determine whether it's an interface edge */
		switch (k) {
		    case 0:	/* 0-1 */
			u0 = 2; u1 = 3;
			break;
		    case 1:	/* 0-2 */
			u0 = 1; u1 = 3;
			break;
		    case 2:	/* 0-3 */
			u0 = 1; u1 = 2;
			break;
		    case 3:	/* 1-2 */
			u0 = 0; u1 = 3;
			break;
		    case 4:	/* 1-3 */
			u0 = 0; u1 = 2;
			break;
		    case 5:	/* 2-3 */
			u0 = 0; u1 = 1;
			break;
		}
		if ((e1->bound_type[u0] | e1->bound_type[u1]) & REMOTE)
		    flag = TRUE;
	    }
#endif
	}
	phgFree(p);
	g_->nedge++;
#if USE_MPI
	if (g_->nprocs > 1) {	/* no discontinuous patches if undistributed */
	    add_local_edge(e, v0, v1, flag);
	}
#endif
    }

    return TRUE;
}

void
phgUpdateEdges(GRID *g)
{
    FunctionEntry;

    if (g == NULL)
	return;

#if USE_MPI
    if (g->nprocs > 1) {
	nedge_global_save = g->nedge_global;
	nedge_save = g->nedge;
    }
#endif

    /* traverse all edges */
    phgDebug((3, "collecting edges\n"));
    g_ = g;
    phgTraverseElements(g, update_edges_callback);

#if USE_MPI
    phgDebug((3, "number of edges collected: %"dFMT"\n", edges_local_count));
    if (g->nprocs > 1)
	sync_edges();
    else
#endif
    g->nedge_global = g->nedge;

    PrintTime(1);
    phgInfo(2, "number of edges = %"dFMT"\n", g->nedge);

    Return;
}

#if Dim == 3
static INT faces_local, faces_shared, gindex;

static BOOLEAN
update_faces_callback0 CB_ARGS(e)
{
    int i;

    /* loop on all faces of the element */
    for (i = 0; i < NFace; i++) {
	if (e->faces[i] != -1)
	    continue;
	if (HasLocalNeighbour(e, i)) {
	    ELEMENT *e1 = e->neighbours[i];
	    if (e->verts[i] < e1->verts[phgOppositeFace(g_, e, i, e1)])
		faces_local++;
	}
	else if (!HasRemoteNeighbour(e, i)) {
	    faces_local++;
	}
	else {
	    faces_shared++;
	}
    }

    return TRUE;
}

static BOOLEAN
update_faces_callback1 CB_ARGS(e)
{
    int i, j;
    ELEMENT *e1;

    /* loop on all faces of the element */
    for (i = 0; i < NFace; i++) {
	if (e->faces[i] != -1 || (e->bound_type[i] & REMOTE))
	    continue;
#if USE_MPI
	if (g_->L2Gmap_face != NULL)
	    g_->L2Gmap_face[g_->nface] = gindex;
#endif
	if (HasLocalNeighbour(e, i)) {
	    j = phgOppositeFace(g_, e, i, e1 = e->neighbours[i]);
	    if (e->verts[i] < e1->verts[j]) {
		e->faces[i] = e1->faces[j] = g_->nface++;
		gindex++;
	    }
	}
	else {
	    e->faces[i] = g_->nface++;
	    gindex++;
	}
    }

    return TRUE;
}

void
phgUpdateFaces(GRID *g)
{
    FunctionEntry;

    if (g == NULL)
	Return;

    g_ = g;

    /* count number of unassigned faces (those having index -1) */
    faces_local = faces_shared = 0;
    phgTraverseElements(g, update_faces_callback0);
    phgDebug((3, "faces: local=%"dFMT", shared=%"dFMT"\n",
		faces_local, faces_shared));
    gindex = g->nface_global;
# if USE_MPI
    if (g->nprocs > 1) {
	int ii, *counts, *displs;
	INT *nbuffer, i, j, total_local, local0 = 0, total_shared;
	struct {
	    ELEMENT	*local, *remote;
	    INT		gindex;		/* global vertex index */
	    int		rank;
	    BYTE	vertex, op_vertex;
	} *local, *remote;
	MPI_Datatype type;

	nbuffer = phgAlloc(g->nprocs * sizeof(*nbuffer));
	counts = phgAlloc(2 * g->nprocs * sizeof(*counts));
	displs = counts + g->nprocs;
	local = (void *)g->L2Gmap_face;
	g->L2Gmap_face = phgRealloc_(g->L2Gmap_face, sizeof(*(g->L2Gmap_face)) *
				(g->nface + faces_local + faces_shared),
		(local == NULL ? 0 : g->nface) * sizeof(*(g->L2Gmap_face)));
	if (local == NULL) {
	    for (i = 0; i < g->nface; i++)
		g_->L2Gmap_face[i] = i;
	}
	/* gather counts of local faces in nbuffer */
	MPI_Allgather(&faces_local, 1, PHG_MPI_INT, nbuffer, 1, PHG_MPI_INT,
			g->comm);
	total_local = 0;
	for (ii = 0; ii < g->nprocs; ii++) {
	    if (ii == g->rank)
		local0 = total_local;
	    total_local += nbuffer[ii];
	    counts[ii] = 0;
	}
	local = phgAlloc(sizeof(*local) * faces_shared);
	remote = phgAlloc(sizeof(*remote) * faces_shared);
	/* copy unassigned shared faces and count the counts for each process */
	j = 0;
	for (i = 0; i < g->neighbours.count; i++) {
	    RNEIGHBOUR *rn = g->neighbours.list + i;
	    if ((rn->local)->faces[rn->vertex] >= 0)
		continue;
	    local[j].local = rn->local;
	    local[j].remote = rn->remote;
	    counts[local[j].rank = rn->rank]++;
	    local[j].gindex = GlobalVertexP(g, (rn->local)->verts[rn->vertex]);
	    local[j].vertex = rn->vertex;
	    local[j].op_vertex = rn->op_vertex;
	    j++;
	}
	j = 0;
	for (ii = 0; ii < g->nprocs; ii++) {
	    displs[ii] = j;
	    assert(j <= INT_MAX - counts[ii]);
	    j += counts[ii];
	}
	MPI_Type_contiguous(sizeof(*local), MPI_BYTE, &type);
	MPI_Type_commit(&type);
	MPI_Alltoallv(local,  counts, displs, type,
		      remote, counts, displs, type, g->comm);
	/* get number of shared faces which should be assigned by the process */
	j = 0;
	for (i = 0; i < faces_shared; i++) {
	    if (GlobalVertexP(g, remote[i].remote->verts[remote[i].op_vertex])
		< remote[i].gindex)
		j++;
	}
	MPI_Allgather(&j, 1, PHG_MPI_INT, nbuffer, 1, PHG_MPI_INT, g->comm);
	total_shared = 0;
	for (ii = 0; ii < g->nprocs; ii++) {
	    if (ii == g->rank)
		gindex += total_local + total_shared;
	    total_shared += nbuffer[ii];
	}
	/* assign indices to shared faces */
	for (i = 0; i < faces_shared; i++) {
	    if (GlobalVertexP(g, remote[i].remote->verts[remote[i].op_vertex])
		< remote[i].gindex) {
		if (g->L2Gmap_face != NULL)
		    g->L2Gmap_face[g->nface] = gindex;
		remote[i].gindex = gindex;
		remote[i].remote->faces[remote[i].op_vertex] = g->nface++;
		gindex++;
	    }
	    else {
		remote[i].gindex = -1;
	    }
	}
	MPI_Alltoallv(remote, counts, displs, type,
		      local,  counts, displs, type, g->comm);
	MPI_Type_free(&type);
	for (i = 0; i < faces_shared; i++) {
	    if (local[i].gindex == -1) continue;
	    if (g->L2Gmap_face != NULL)
		g->L2Gmap_face[g->nface] = local[i].gindex;
	    local[i].local->faces[local[i].vertex] = g->nface++;
	}
	gindex = g->nface_global + local0;
	g->nface_global += total_local + total_shared;
	phgFree(remote);
	phgFree(local);
	phgFree(counts);
	phgFree(nbuffer);
    }
# endif

    /* assign local indices */
    phgTraverseElements(g, update_faces_callback1);

    if (g->nprocs <= 1)
	g->nface_global = g->nface;

    PrintTime(1);
    Return;
}
#endif	/* Dim == 3 */

void
phgUpdateElementIndices(GRID *g, INT nelem_old, INT nelem_global_old)
{
#if USE_MPI
    INT i, index0, n;
#endif

    FunctionEntry;

#if USE_MPI
    if (g == NULL || g->nprocs <= 1)
	Return;

    g_ = g;

    n = g->nelem - nelem_old;
    MPI_Scan(&n, &index0, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    index0 += nelem_global_old - n;

    if (g->L2Gmap_elem == NULL) {
	g->L2Gmap_elem = phgAlloc(g->nelem * sizeof(*g->L2Gmap_elem));
	for (i = 0; i < nelem_old; i++)
	    g->L2Gmap_elem[i] = i;
    }
    else {
	g->L2Gmap_elem = phgRealloc_(g->L2Gmap_elem,
				g->nelem * sizeof(*g->L2Gmap_elem),
				nelem_old * sizeof(*g->L2Gmap_elem));
    }

    /* assign global indices */
    for (i = 0; i < n; i++)
	g->L2Gmap_elem[nelem_old + i] = index0 + i;

    PrintTime(1);
#else	/* USE_MPI */
    (void)g;
    (void)nelem_old;
    (void)nelem_global_old;
#endif	/* USE_MPI */

    Return;
}

void
phgDumpElement(GRID *g, ELEMENT *e)
{
    int i;
    char s[128];
#ifdef VL_
# undef VL_
#endif
#define VL_(i)	(e->verts[i])
#ifdef EL_
# undef EL_
#endif
#define EL_(i)	(e->edges[i])
#ifdef FL_
# undef FL_
#endif
#define FL_(i)	(e->faces[i])
#if USE_MPI
    RNEIGHBOUR *rn;
# ifdef VG_
#  undef VG_
# endif
# define VG_(i)	(e->verts[i] < 0 ? -1 : GlobalVertexP(g, e->verts[i]))
# ifdef EG_
#  undef EG_
# endif
# define EG_(i)	(e->edges[i] < 0 ? -1 : GlobalEdge(g, e->edges[i]))
# ifdef FG_
#  undef FG_
# endif
# define FG_(i)	(e->faces[i] < 0 ? -1 : GlobalFace(g, e->faces[i]))
#endif

    if (e == NULL) {
	phgWarning("phgDumpElement: e == NULL\n", e);
	return;
    }

#if Dim == 2
    phgInfo(-1, "E(%"dFMT"): %p, %s, mark=%d, bdry=%x:%x:%x, " "gen=%d\n",
		(g->flags & ELEM_FLAG) && e->index >= 0 ?
				GlobalElement(g, e->index) : (INT)(-1),
		e, GTypeName(e->type), e->mark,
		e->bound_type[0], e->bound_type[1], e->bound_type[2],
		e->generation);
#if USE_MPI
    phgInfo(-1, "  children=%p %p, parent=%p, flag=%d\n",
		e->children[0], e->children[1], e->parent, e->flag);
    phgInfo(-1, "  verts=%"dFMT" %"dFMT" %"dFMT" (%"dFMT" %"dFMT" %"dFMT")\n",
		VL_(0), VL_(1), VL_(2), VG_(0), VG_(1), VG_(2));
    phgInfo(-1, "  edges=%"dFMT" %"dFMT" %"dFMT" (%"dFMT" %"dFMT" %"dFMT")\n",
		EL_(0), EL_(1), EL_(2), EG_(0), EG_(1), EG_(2));
#else
    phgInfo(-1, "  children=%p %p, parent=%p\n",
		e->children[0], e->children[1]);
    phgInfo(-1, "  verts=%"dFMT" %"dFMT" %"dFMT"\n", VL_(0), VL_(1), VL_(2));
    phgInfo(-1, "  edges=%"dFMT" %"dFMT" %"dFMT"\n",
			e->edges[0], e->edges[1], e->edges[2]);
#endif
#else	/* Dim */
    phgInfo(-1, "E(%"dFMT"): %p, %s, mark=%d, bdry=%x:%x:%x:%x, " "gen=%d\n",
		(g->flags & ELEM_FLAG) && e->index >= 0 ?
				GlobalElement(g, e->index) : -1,
		e, GTypeName(e->type), e->mark,
		e->bound_type[0], e->bound_type[1], e->bound_type[2],
		e->bound_type[3], e->generation);
#if USE_MPI
    phgInfo(-1, "  children=%p %p, parent=%p, flag=%d\n",
		e->children[0], e->children[1], e->parent, e->flag);
    phgInfo(-1, "  verts=%"dFMT" %"dFMT" %"dFMT" %"dFMT" (%"dFMT" %"dFMT" %"dFMT" %"dFMT")\n",
		VL_(0), VL_(1), VL_(2), VL_(3), VG_(0), VG_(1), VG_(2), VG_(3));
    phgInfo(-1, "  edges=%"dFMT" %"dFMT" %"dFMT" %"dFMT" %"dFMT" %"dFMT" (%"dFMT" %"dFMT" %"dFMT" %"dFMT" %"dFMT" %"dFMT")\n",
		EL_(0), EL_(1), EL_(2), EL_(3), EL_(4), EL_(5),
		EG_(0), EG_(1), EG_(2), EG_(3), EG_(4), EG_(5));
    phgInfo(-1, "  faces=%"dFMT" %"dFMT" %"dFMT" %"dFMT" (%"dFMT" %"dFMT" %"dFMT" %"dFMT")\n",
		FL_(0), FL_(1), FL_(2), FL_(3),
		FG_(0), FG_(1), FG_(2), FG_(3));
#else
    phgInfo(-1, "  children=%p %p, parent=%p\n",
		e->children[0], e->children[1]);
    phgInfo(-1, "  verts=%"dFMT" %"dFMT" %"dFMT" %"dFMT"\n", VL_(0), VL_(1), VL_(2), VL_(3));
    phgInfo(-1, "  edges=%"dFMT" %"dFMT" %"dFMT" %"dFMT" %"dFMT" %"dFMT"\n", e->edges[0], e->edges[1],
		e->edges[2], e->edges[3], e->edges[4], e->edges[5]);
    phgInfo(-1, "  faces=%"dFMT" %"dFMT" %"dFMT" %"dFMT"\n", e->faces[0], e->faces[1],
		e->faces[2], e->faces[3]);
#endif
#endif	/* Dim */
#undef V_
    s[0] = '\0';
    for (i = 0; i < NFace; i++) {
	if (HasLocalNeighbour(e, i))
	    sprintf(s+strlen(s), "%p ", e->neighbours[i]);
	else if ((e->bound_type[i] & DIRICHLET))
	    sprintf(s+strlen(s), "D ");
	else if ((e->bound_type[i] & NEUMANN))
	    sprintf(s+strlen(s), "N ");
	else if ((e->bound_type[i] & UNDEFINED))
	    sprintf(s+strlen(s), "U ");
#if USE_MPI
	else if ((e->bound_type[i] & REMOTE) && IsLeaf(e)) {
	    if (g != NULL && (rn = GetRNeighbour(g, e, i)) != NULL)
		sprintf(s+strlen(s), "%p@p%d ", rn->remote, rn->rank);
	    else
		sprintf(s+strlen(s), "unset ");
	}
#endif
	else
	    sprintf(s+strlen(s), "undef ");
    }
    phgInfo(-1, "  neighbours=%s\n", s);

    if (phgVerbosity < 2)
	return;

#if USE_MPI
    if (IsLeaf(e))
      for (i = 0; i < NFace; i++) {
	if (!(e->bound_type[i] & REMOTE) || g == NULL ||
		(rn = GetRNeighbour(g,e,i)) == NULL)
	    continue;
	phgInfo(-1, "  RNeighbour%d: vertex=%d, op_vertex=%d, face=%d %d %d\n",
			i, rn->vertex, rn->op_vertex, 
			rn->rface[0], rn->rface[1], rn->rface[2]);
    }
#endif
}

void
phgDumpPatch(GRID *g, ELEMENT *e, int v0, int v1)
/* dumps out a patch */
{
    int j, n;
    ELEMENT **p;

    phgInfo(-1, "========= the patch:\n");
    n = phgGetPatch(g, e, v0, v1, &p);
    for (j = 0; j < n; j++)
	phgDumpElement(g, p[j]);
    phgFree(p);
    phgInfo(-1, "=========\n");
}

static BOOLEAN
dump_callback CB_ARGS(e)
{
    phgDumpElement(g_, e);

    return TRUE;
}

static void
dump_head(GRID *g)
{
    phgInfo(-1, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
    phgInfo(-1, "nroot = %"dFMT"\n", g->nroot);
    phgInfo(-1, "nleaf = %"dFMT"\n", g->nleaf);
    phgInfo(-1, "ntree = %"dFMT"\n", g->ntree);
    phgInfo(-1, "nvert = %"dFMT"\n", g->nvert);
    phgInfo(-1, "nedge = %"dFMT"\n", g->nedge);
#if Dim == 3
    phgInfo(-1, "nface = %"dFMT"\n", g->nface);
#endif
    phgInfo(-1, "nelem = %"dFMT"\n", g->nelem);
    /*if (g->nprocs > 1)*/ {
	if (g->period != NULL)
	    phgInfo(-1, "nvert_global = %"dFMT" (%"dFMT")\n", g->nvert_global,
						    g->period->nvert_global);
	else
	    phgInfo(-1, "nvert_global = %"dFMT"\n", g->nvert_global);
	phgInfo(-1, "nedge_global = %"dFMT"\n", g->nedge_global);
#if Dim == 3
	phgInfo(-1, "nface_global = %"dFMT"\n", g->nface_global);
#endif
	phgInfo(-1, "nelem_global = %"dFMT"\n", g->nelem_global);

	phgInfo(-1, "nvert_owned = %"dFMT"\n", g->nvert_owned);
	phgInfo(-1, "nedge_owned = %"dFMT"\n", g->nedge_owned);
#if Dim == 3
	phgInfo(-1, "nface_owned = %"dFMT"\n", g->nface_owned);
#endif
	phgInfo(-1, "nelem_owned = %"dFMT"\n", g->nelem_owned);
    }
#if USE_MPI
    phgInfo(-1, "neighbours.count = %"dFMT"\n", g->neighbours.count);
#endif
}

static void
dump_tail(GRID *g)
{
    Unused(g);
    phgInfo(-1, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
}

void
phgDumpGridInfo(GRID *g)
{
    ParallelFunction(phgDumpGridInfo, g != NULL && g->nprocs > 1);
    if (g == NULL)
	return;
    dump_head(g_ = g);
    dump_tail(g);
}

void
phgDumpGrid(GRID *g)
{
    ParallelFunction(phgDumpGrid, g != NULL && g->nprocs > 1);
    if (g == NULL)
	return;
    dump_head(g_ = g);
    phgTraverseElements(g, dump_callback);
    dump_tail(g);
}

void
phgDumpTree(GRID *g)
{
    ParallelFunction(phgDumpTree, g != NULL && g->nprocs > 1);
    if (g == NULL)
	return;
    dump_head(g_ = g);
    phgTraverseAllElements(g, dump_callback);
    dump_tail(g);
}

void
phgDumpBranch(GRID *g, ELEMENT *e)
{
    phgDumpElement(g, e);
    if (e->children[0] != NULL)
	phgDumpBranch(g, e->children[0]);
    if (e->children[1] != NULL)
	phgDumpBranch(g, e->children[1]);
}

#if defined(GZIP_PROG) || defined(BZIP2_PROG)
#define MAX_FILES 16
static struct {
    FILE *fp;
    BOOLEAN is_pipe;
} open_files[MAX_FILES];
#endif

void
phgCloseInputFile_(FILE *fp)
{
#if defined(GZIP_PROG) || defined(BZIP2_PROG)
    int i;
    BOOLEAN is_pipe = FALSE;
    for (i = 0; i < MAX_FILES; i++) {
	if (open_files[i].fp == fp) {
	    is_pipe = open_files[i].is_pipe;
	    open_files[i].fp = NULL;
	    break;
	}
    }
    (is_pipe ? pclose : fclose)(fp);
#else
    if (fp != NULL) fclose(fp);
#endif
}

FILE *
phgOpenInputFile_(const char *filename)
{
    FILE *fp;
#if defined(GZIP_PROG) || defined(BZIP2_PROG)
    char cmd[PATH_MAX + 16];
    int len = strlen(filename);
    BOOLEAN is_pipe = FALSE;
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	int i;
	initialized = TRUE;
	for (i = 0; i < MAX_FILES; i++)
	    open_files[i].fp = NULL;
    }
#endif

#ifdef GZIP_PROG
    if (len > 3 && !phgMemcmp(filename + len - 3, ".gz", 3)) {
	if ((fp = fopen(filename, "r")) != NULL) {
	    fclose(fp);
	    phgInfo(1, "gzipped file.\n");
	    sprintf(cmd, GZIP_PROG " -dc %s 2>/dev/null", filename);
	    fp = popen(cmd, "r");
	    is_pipe = TRUE;
	}
    } else
#endif
#ifdef BZIP2_PROG
    if (len > 4 && !phgMemcmp(filename + len - 4, ".bz2", 4)) {
	if ((fp = fopen(filename, "r")) != NULL) {
	    fclose(fp);
	    phgInfo(1, "bzip2 file.\n");
	    sprintf(cmd, BZIP2_PROG " -dc %s 2>/dev/null", filename);
	    fp = popen(cmd, "r");
	    is_pipe = TRUE;
	}
    } else
#endif
    {
	fp = fopen(filename, "r");
    }
    if (fp == NULL)
	phgWarning("can't open file \"%s\"!\n", filename);
#if defined(GZIP_PROG) || defined(BZIP2_PROG)
    else {
	int i;
	for (i = 0; i < MAX_FILES; i++)
	    if (open_files[i].fp == NULL) break;
	if (i >= MAX_FILES) {
	    phgWarning("too many open files.\n");
	    (is_pipe ? pclose : fclose)(fp);
	    fp = NULL;
	}
	else {
	    open_files[i].fp = fp;
	    open_files[i].is_pipe = is_pipe;
	}
    }
#endif

    return fp;
}

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

BTYPE
phgImportSetDefaultBdryType(BTYPE type)
{
    BTYPE t, oldtype;
    int i;

    oldtype = (btype_index < 0 ? UNDEFINED : btype_values[btype_index]);
    btype_index = -1;
    t = type & BDRY_MASK;
    for (i = 0; i < sizeof(btype_values) / sizeof(btype_values[0]); i++) {
	if (t == btype_values[i]) {
	    btype_index = i;
	    break;
	}
    }

    if (btype_index == -1)
	phgError(1, "%s: invalid boundary type %d.\n", __func__, type);

    return oldtype;
}

static BOOLEAN
default_bdry_callback CB_ARGS(e)
{
    int i;

    for (i = 0; i < NFace; i++) {
	if (e->bound_type[i] == UNDEFINED) {
	    e->bound_type[i] &= ~BDRY_MASK;
	    e->bound_type[i] |= btype_values[btype_index];
	}
    }

    return TRUE;
}

static void
set_default_bdry_type(GRID *g)
/* sets UNDEFINED boundary type to default type */
{
    if (btype_values[btype_index] == UNDEFINED)
	return;

    /* Note: can't use ForAllElements because g->elems is not up to date */
    phgTraverseElements(g, default_bdry_callback);
}

static int
map_comp1(const void *p0, const void *p1)
/* compares two maps, if the intervals overlap than returns 0, otherwise
 * returns p0->min - p1->min */
{
    if (((const bcmap_t *)p0)->min > ((const bcmap_t *)p1)->max)
	return 1;
    if (((const bcmap_t *)p0)->max < ((const bcmap_t *)p1)->min)
	return -1;
    return 0;
}

static int
default_bcmap(int bc)
{
    if (bcmap_n > 0) {
	INT i;
#if 0
	for (i = 0; i < bcmap_n; i++)
	    if (bc >= bcmap[i].min && bc <= bcmap[i].max)
		return bcmap[i].btype;
	return -1;
#else	/* 0 */
	bcmap_t key;
	key.min = key.max = bc;
	i = phgBinarySearch(bcmap_n, bcmap, NULL, &key, sizeof(key), map_comp1);
	if (i >= bcmap_n || map_comp1(bcmap + i, &key) != 0)
	    return -1;
	return bcmap[i].btype;
#endif	/* 0 */
    }

    /* fall to built-in default bcmap */
    switch (bc) {
	case 1:
	    return DIRICHLET;
	case 2:
	    return NEUMANN;
	case 3:
	    return BDRY_USER1;
	case 4:
	    return BDRY_USER2;
	case 0:
	    return UNDEFINED;
	default:
	    return -1;		/* invalid/undefined bctype */
    }
}

static BDRY_MAP_FUNC user_bcmap = NULL;

void
phgImportSetBdryMapFunc(BDRY_MAP_FUNC func)
/* Defines an user-function for mapping boundary types for Medit/Gambit.
 * Note: should be called before phgImport. */
{
    user_bcmap = func;
}

int
_phg_bcmap(GRID *g, int input_bc)
/* This function is used by phgImport{Medit/Gambit} to map a bc type from
 * the input file. This function also collects all types used in the input
 * file in input_list[], and creates/updates g->bc_list[] */
{
    static int *input_list = NULL, input_alloc = 0, input_n = 0;
    int i, j, bc = UNDEFINED;
    int *found;

    if (g == NULL) {
	/* reset input_list[] */
	phgFree(input_list);
	input_list = NULL;
	input_alloc = input_n = 0;
	return 0;
    }

    assert(sizeof(BTYPE) < sizeof(int));
    if (input_n == 0)
	found = NULL;
    else
	found = bsearch(&input_bc, input_list, input_n, sizeof(*input_list),
			phgCompint);

    bc = (*(user_bcmap == NULL ? default_bcmap : user_bcmap))(input_bc);

    if (found == NULL) {
	if (input_n >= input_alloc) {
	    input_list = phgRealloc_(input_list,
				(input_alloc + 16) * sizeof(*input_list),
				input_n * sizeof(*input_list));
	    input_alloc += 16;
	}
	for (i = input_n; i > 0; i--) {
	    if (input_list[i - 1] < input_bc)
		break;
	    input_list[i] = input_list[i - 1];
	}
	input_list[i] = input_bc;
	input_n++;

	/* adjust bc and add it to g->bc_list[] and g->bc_rmap[] */
	if (bc < 0 || (bc >> BDRY_SHIFT) >= BDRY_UB) {
	    bc = UNDEFINED;
	    phgWarning("bdry type %d mapped to UNDEFINED\n", input_bc);
	}
	else {
	    phgInfo(1, "bdry type %d mapped to %d (%s)\n",
				input_bc, bc, BTypeName(bc));
	}

	j = bc >> BDRY_SHIFT;
	if (g->bc_n == 0)
	    found = NULL;
	else
	    found = bsearch(&j, g->bc_list, g->bc_n, sizeof(*g->bc_list),
				phgCompint);
	if (found == NULL) {
	    /* add a new entry in g->bc_list[] and g->bc_rmap[] */
	    if (g->bc_n >= g->bc_alloc) {
		g->bc_list = phgRealloc_(g->bc_list,
				(g->bc_alloc + 16) * sizeof(*g->bc_list),
				g->bc_n * sizeof(*g->bc_list));
		g->bc_rmap = phgRealloc_(g->bc_rmap,
				(g->bc_alloc + 16) * sizeof(*g->bc_rmap),
				g->bc_n * sizeof(*g->bc_rmap));
		g->bc_alloc += 16;
	    }
	    for (i = g->bc_n; i > 0; i--) {
		if (g->bc_list[i-1] < j)
		    break;
		g->bc_list[i] = g->bc_list[i-1];
		g->bc_rmap[i] = g->bc_rmap[i-1];
	    }
	    g->bc_list[i] = j;
	    g->bc_rmap[i] = input_bc;
	    g->bc_n++;
	}
    }
    else if (bc < 0 || (bc >> BDRY_SHIFT) >= BDRY_UB) {
	bc = UNDEFINED;
    }

    return bc;
}

BOOLEAN
phgImport(GRID *g, const char *filename, BOOLEAN distr)
{
    FILE *fp;
    char line[128], magic[128];
    BOOLEAN ret = FALSE;
    INT i;
    int flag;
    static BOOLEAN initialized = FALSE;
    static BOOLEAN auto_bcmap = TRUE;

    if (!initialized) {
	initialized = TRUE;
	phgImportSetDefaultBdryType(UNDEFINED);
	phgOptionsRegisterKeyword("-default_bdry_type", "Default boundary type",
		btype_names, &btype_index);
	phgOptionsRegisterFilename("-bcmap_file", "Load bcmap from the "
				   "specified file. Note: this option is only "
				   "effective for Medit and Gambit formats",
				   &bcmap_file);
	phgOptionsRegisterNoArg("-auto_bcmap", "When loading an Medit mesh, "
				"load bcmap from the file with the same "
				"basename and extension .bcmap if the latter "
				"exists. Note: this option is only effective "
				"for Medit and Gambit formats", &auto_bcmap);
	return TRUE;
    }

    if (g == NULL)
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
    if (!strncmp(magic, "DIM:", 4) || !strncmp(magic, "DIM_OF_WORLD:", 13)) {
	flag = phgImportALBERT(g, filename, FALSE);
    }
    else if (!strcmp(magic, "MeshVersionFormatted")) {
	flag = phgImportMedit(g, filename, FALSE);
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

#if USE_MPI

static int
comp_rneighbour(const void *p1, const void *p2)
/* Note: the function is made global since it's also used by refine.c */
{
    int i;
    const RNEIGHBOUR *ptr1 = p1, *ptr2 = p2;

#if 0	/* no longer necessary */
    if (ptr1->local->children[0] != NULL || ptr2->local->children[0] != NULL) {
	/* move non-leaf element(s) to the end of list */
	if (ptr1->local->children[0] == NULL)
	    return -1;
	else if (ptr2->local->children[0] == NULL)
	    return 1;
	else
	    return 0;
    }
#endif

    if ((i = ptr1->rank - ptr2->rank) != 0)
	return i;

#if 0
    /* The next two lines are commented out for deterministic results.
     * 
     * Note: if in the future they turn out to be necessary, then a workaround
     * is compare the vertices of ptr1->local and ptr2->local */
    if (ptr1->local->index != ptr2->local->index)
	return (ptr1->local->index < ptr2->local->index) ? -1 : 1;
#endif

    return ptr1->vertex - ptr2->vertex;
}

INT
phgCountRNeighbours(int nprocs, RNEIGHBOUR *list, INT count, int *counts,
			int *displs)
/* counts number of elements for different processes, and removes entries on
   non leaf elements */
{
    int i;
    INT j;
    RNEIGHBOUR *p0, *p1;
    ELEMENT *e;

    if (counts == NULL)
	return count;

    /* First, remove entries for non leaf elements */
    p0 = list;
    while (p0 < list + count) {
	e = p0->local;
	if (!IsLeaf(e))
	    break;
	p0++;
    }
    p1 = p0;
    while (p1 < list + count) {
	do {
	    p1->local->neighbours[p1->vertex] = (void *)-1;
	} while (++p1 < list + count && ((e = p1->local)->children[0] != NULL
		|| e->children[1] != NULL));
	do {
	    *(p0++) = *(p1++);
	} while (p1 < list + count && (e = p1->local)->children[0] == NULL
		&& e->children[1] == NULL);
    }
    count = p0 - list;

    /* Then sort the list */
    if (count > 0) qsort(list, count, sizeof(RNEIGHBOUR), comp_rneighbour);

    /* Finally, count number of entries for each process */
    p0 = list;
    j = 0;
    for (i = 0; i < nprocs; i++) {
	while (p0 - list < count && p0->rank < i /*&& 
		p0->local->children[0] == NULL*/) p0++;
	p1 = p0;
	while (p1 - list < count && p1->rank == i /*&&
		p1->local->children[0] == NULL*/) p1++;
	counts[i] = p1 - p0;
	displs[i] = j;
	assert(j <= INT_MAX - (p1 - p0));
	j += p1 - p0;
	p0 = p1;
    }

    return j;
}

#endif	/* USE_MPI */

GRID *
phgDupGrid(GRID *g, BOOLEAN flatten)
{
    GRID *g1;
    ELEMENT *e, *e0 = NULL, *e1;
    INT i, j;
    INT *map_vert, *map_edge, *map_face, *map_elem;
    int ii, *scnts, *sdsps, *rcnts, *rdsps;

    if (g->period != NULL)
	phgError(1, "%s: unimplemented for periodic boundaries.\n", __func__);

    if (!flatten)
	phgError(1, "%s: unimplemented for flatten == FALSE.\n", __func__);

    g1 = phgNewGrid(g->flags);
    memcpy(g1->bbox, g->bbox, sizeof(g->bbox));
    if (g->bc_n > 0) {
	g1->bc_n = g1->bc_alloc = g->bc_n;
	g1->bc_rmap = phgAlloc(g->bc_n * sizeof(*g1->bc_rmap));
	memcpy(g1->bc_rmap, g->bc_rmap, g->bc_n * sizeof(*g->bc_rmap));
	g1->bc_list = phgAlloc(g->bc_n * sizeof(*g1->bc_list));
	memcpy(g1->bc_list, g->bc_list, g->bc_n * sizeof(*g->bc_list));
    }
    g1->volume = g->volume;

#if ALLOW_CURVED_BOUNDARY
    if (g->bdry_funcs != NULL) {
	EXPR **p = g->bdry_funcs;
	for (ii = 0; *p != NULL; ii++, p++);
	g1->bdry_funcs = phgAlloc((ii + 1) * sizeof(*g1->bdry_funcs));
	for (ii = 0, p = g->bdry_funcs; *p != NULL; ii++, p++)
	    g1->bdry_funcs[ii] = phgDup3DFunction(*p);
	g1->bdry_funcs[ii] = NULL;
    }
#endif	/* ALLOW_CURVED_BOUNDARY */

    g1->nvert = g->nvert;
    g1->nedge = g->nedge;
    g1->nface = g->nface;
    g1->nelem = g->nelem;

    map_vert = map_edge = map_face = map_elem = NULL;
    scnts = sdsps = rcnts = rdsps = NULL;
    if (g->nprocs > 1) {
	scnts = phgAlloc(4 * g->nprocs * sizeof(*scnts));
	sdsps = scnts + g->nprocs;
	rcnts = sdsps + g->nprocs;
	rdsps = rcnts + g->nprocs;
	memset(scnts, 0, g->nprocs * sizeof(*scnts));

	if ((g1->flags & VERT_FLAG)) {
	    map_vert = phgAlloc(g->nvert * sizeof(*map_vert));
	    for (i = j = 0; i < g->nvert; i++) {
		map_vert[i] = j;
		if (g->types_vert[i] != UNREFERENCED)
		    j++;
	    }
	    g1->nvert = j;
	}

	if ((g1->flags & EDGE_FLAG)) {
	    map_edge = phgAlloc(g->nedge * sizeof(*map_edge));
	    for (i = j = 0; i < g->nedge; i++) {
		map_edge[i] = j;
		if (g->types_edge[i] != UNREFERENCED)
		    j++;
	    }
	    g1->nedge = j;
	}

	if ((g1->flags & FACE_FLAG)) {
	    map_face = phgAlloc(g->nface * sizeof(*map_face));
	    for (i = j = 0; i < g->nface; i++) {
		map_face[i] = j;
		if (g->types_face[i] != UNREFERENCED)
		    j++;
	    }
	    g1->nface = j;
	}

	if ((g1->flags & ELEM_FLAG)) {
	    map_elem = phgAlloc(g->nelem * sizeof(*map_elem));
	    for (i = j = 0; i < g->nelem; i++) {
		map_elem[i] = j;
		if (g->types_elem[i] != UNREFERENCED)
		    j++;
	    }
	    g1->nelem = j;
	}
    }

#if USE_MPI
# if !NO_NEW_COMMUNICATOR
    MPI_Comm_free(&g1->comm);
    MPI_Comm_dup(g->comm, &g1->comm);
#else
    g1->comm = g->comm;
# endif /* !NO_NEW_COMMUNICATOR */
#endif	/* USE_MPI */

    g1->nprocs = g->nprocs;
    g1->rank = g->rank;

    g1->lif = g->lif;
    g1->filename = (g->filename != NULL ? strdup(g->filename) : NULL);
    

    g1->nvert_global = g->nvert_global;
    g1->nedge_global = g->nedge_global;
    g1->nface_global = g->nface_global;
    g1->nelem_global = g->nelem_global;

    g1->nvert_owned = g->nvert_owned;
    g1->nedge_owned = g->nedge_owned;
    g1->nface_owned = g->nface_owned;
    g1->nelem_owned = g->nelem_owned;

    g1->verts = phgAlloc(g1->nvert * sizeof(*g1->verts));
    if (g->types_vert != NULL)
	g1->types_vert = phgAlloc(g1->nvert * sizeof(*g1->types_vert));
#if USE_MPI
    if (g->L2Gmap_vert != NULL)
	g1->L2Gmap_vert = phgAlloc(g1->nvert * sizeof(*g1->L2Gmap_vert));
    if (g->owner_index_vert != NULL) {
	g1->owner_index_vert = phgAlloc(g1->nvert *
						sizeof(*g1->owner_index_vert));
	g1->owner_rank_vert = phgAlloc(g1->nvert *
						sizeof(*g1->owner_rank_vert));
    }
#endif	/* USE_MPI */

    for (i = j = 0; i < g->nvert; i++) {
	if (g->types_vert != NULL && g->types_vert[i] == UNREFERENCED)
	    continue;
	memcpy(g1->verts + j, g->verts + i, sizeof(COORD));
#if USE_MPI
	if (g->types_vert != NULL)
	    g1->types_vert[j] = g->types_vert[i];
	if (g1->L2Gmap_vert != NULL)
	    g1->L2Gmap_vert[j] = g->L2Gmap_vert[i];
	if (g1->owner_index_vert != NULL) {
	    g1->owner_index_vert[j] = g->owner_index_vert[i];
	    ii = g1->owner_rank_vert[j] = g->owner_rank_vert[i];
	    if (ii != -1 && ii != g->rank) {
		assert(scnts[ii] < INT_MAX);
		scnts[ii]++;
	    }
	}
#endif	/* USE_MPI */
	j++;
    }

    if (!(g->flags & EDGE_FLAG))
	goto face;
    g1->types_edge = phgAlloc(g1->nedge * sizeof(*g1->types_edge));
#if USE_MPI
    if (g->L2Gmap_edge != NULL)
	g1->L2Gmap_edge = phgAlloc(g1->nedge * sizeof(*g1->L2Gmap_edge));
    if (g->owner_index_edge != NULL) {
	g1->owner_index_edge = phgAlloc(g1->nedge *
						sizeof(*g1->owner_index_edge));
	g1->owner_rank_edge = phgAlloc(g1->nedge *
						sizeof(*g1->owner_rank_edge));
    }
#endif	/* USE_MPI */
    for (i = j = 0; i < g->nedge; i++) {
	if (g->types_edge[i] == UNREFERENCED)
	    continue;
	g1->types_edge[j] = g->types_edge[i];
#if USE_MPI
	if (g1->L2Gmap_edge != NULL)
	    g1->L2Gmap_edge[j] = g->L2Gmap_edge[i];
	if (g1->owner_index_edge != NULL) {
	    g1->owner_index_edge[j] = g->owner_index_edge[i];
	    ii = g1->owner_rank_edge[j] = g->owner_rank_edge[i];
	    if (ii != -1 && ii != g->rank) {
		assert(scnts[ii] < INT_MAX);
		scnts[ii]++;
	    }
	}
#endif	/* USE_MPI */
	j++;
    }

face:
    if (!(g->flags & FACE_FLAG))
	goto elem;
    g1->types_face = phgAlloc(g1->nface * sizeof(*g1->types_face));
#if USE_MPI
    if (g->L2Gmap_face != NULL)
	g1->L2Gmap_face = phgAlloc(g1->nface * sizeof(*g1->L2Gmap_face));
    if (g->owner_index_face != NULL) {
	g1->owner_index_face = phgAlloc(g1->nface *
						sizeof(*g1->owner_index_face));
	g1->owner_rank_face = phgAlloc(g1->nface *
						sizeof(*g1->owner_rank_face));
    }
#endif	/* USE_MPI */
    for (i = j = 0; i < g->nface; i++) {
	if (g->types_face[i] == UNREFERENCED)
	    continue;
	g1->types_face[j] = g->types_face[i];
#if USE_MPI
	if (g1->L2Gmap_face != NULL)
	    g1->L2Gmap_face[j] = g->L2Gmap_face[i];
	if (g1->owner_index_face != NULL) {
	    g1->owner_index_face[j] = g->owner_index_face[i];
	    ii = g1->owner_rank_face[j] = g->owner_rank_face[i];
	    if (ii != -1 && ii != g->rank) {
		assert(scnts[ii] < INT_MAX);
		scnts[ii]++;
	    }
	}
#endif	/* USE_MPI */
	j++;
    }

elem:
    if (!(g->flags & ELEM_FLAG))
	goto cont;
    g1->nroot = g1->nleaf = g1->ntree = g1->nelem;
    g1->types_elem = phgAlloc(g1->nelem * sizeof(*g1->types_elem));
#if USE_MPI
    if (g->L2Gmap_elem != NULL)
	g1->L2Gmap_elem = phgAlloc(g1->nelem * sizeof(*g1->L2Gmap_elem));
    if (g->owner_index_elem != NULL) {
	g1->owner_index_elem = phgAlloc(g1->nelem *
						sizeof(*g1->owner_index_elem));
	g1->owner_rank_elem = phgAlloc(g1->nelem *
						sizeof(*g1->owner_rank_elem));
    }
#endif	/* USE_MPI */
    for (i = j = 0; i < g->nelem; i++) {
	if (g->types_elem[i] == UNREFERENCED)
	    continue;
	g1->types_elem[j] = g->types_elem[i];
#if USE_MPI
	if (g1->L2Gmap_elem != NULL)
	    g1->L2Gmap_elem[j] = g->L2Gmap_elem[i];
	if (g1->owner_index_elem != NULL) {
	    g1->owner_index_elem[j] = g->owner_index_elem[i];
	    ii = g1->owner_rank_elem[j] = g->owner_rank_elem[i];
	    if (ii != -1 && ii != g->rank) {
		assert(scnts[ii] < INT_MAX);
		scnts[ii]++;
	    }
	}
#endif	/* USE_MPI */
	j++;
    }

#if USE_MPI
    if (g->nprocs > 1) {	/* update owner_index_xxxx[] arrays */
	MPI_Datatype type;
	int slen, rlen;
	struct {
	    INT		index;
	    GTYPE	type;
	} *sbuf, *rbuf, *pbuf;

	MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, g->comm);

	/* compute sdsps[] and rdsps[] */
	slen = rlen = 0;
	for (ii = 0; ii < g->nprocs; ii++) {
	    sdsps[ii] = slen;
	    rdsps[ii] = rlen;
	    assert(slen <= INT_MAX - scnts[ii]);
	    slen += scnts[ii];
	    assert(rlen <= INT_MAX - rcnts[ii]);
	    rlen += rcnts[ii];
	}
	sbuf = malloc((slen + rlen) * sizeof(*sbuf));
	rbuf = sbuf + slen;

	/* prepare sbuf[] */

	if (g1->owner_index_vert != NULL) {
	    for (j = 0; j < g1->nvert; j++) {
        	if ((ii = g1->owner_rank_vert[j]) == -1)
		    continue;
        	if (ii == g->rank) {
		    assert(g1->owner_index_vert[j] >= 0 &&
			   g1->owner_index_vert[j] < g->nvert);
		    g1->owner_index_vert[j] = map_vert[g1->owner_index_vert[j]];
		    continue;
		}
		pbuf = sbuf + (sdsps[ii]++);
		pbuf->index = g1->owner_index_vert[j];
		pbuf->type = VERTEX;
	    }
	}

	if (g1->owner_index_edge != NULL) {
	    for (j = 0; j < g1->nedge; j++) {
        	if ((ii = g1->owner_rank_edge[j]) == -1)
		    continue;
        	if (ii == g->rank) {
		    assert(g1->owner_index_edge[j] >= 0 &&
			   g1->owner_index_edge[j] < g->nedge);
		    g1->owner_index_edge[j] = map_edge[g1->owner_index_edge[j]];
		    continue;
		}
		pbuf = sbuf + (sdsps[ii]++);
		pbuf->index = g1->owner_index_edge[j];
		pbuf->type = EDGE;
	    }
	}

	if (g1->owner_index_face != NULL) {
	    for (j = 0; j < g1->nface; j++) {
        	if ((ii = g1->owner_rank_face[j]) == -1)
		    continue;
        	if (ii == g->rank) {
		    assert(g1->owner_index_face[j] >= 0 &&
			   g1->owner_index_face[j] < g->nface);
		    g1->owner_index_face[j] = map_face[g1->owner_index_face[j]];
		    continue;
		}
		pbuf = sbuf + (sdsps[ii]++);
		pbuf->index = g1->owner_index_face[j];
		pbuf->type = FACE;
	    }
	}

	if (g1->owner_index_elem != NULL) {
	    for (j = 0; j < g1->nelem; j++) {
        	if ((ii = g1->owner_rank_elem[j]) == -1)
		    continue;
        	if (ii == g->rank) {
		    assert(g1->owner_index_elem[j] >= 0 &&
			   g1->owner_index_elem[j] < g->nelem);
		    g1->owner_index_elem[j] = map_elem[g1->owner_index_elem[j]];
		    continue;
		}
		pbuf = sbuf + (sdsps[ii]++);
		pbuf->index = g1->owner_index_elem[j];
		pbuf->type = VOLUME;
	    }
	}

	/* restore sdsps[] */
	for (ii = 0; ii < g->nprocs; ii++)
	    sdsps[ii] -= scnts[ii];

	/* exchange indices */
	MPI_Type_contiguous(sizeof(*sbuf), MPI_BYTE, &type);
	MPI_Type_commit(&type);
	MPI_Alltoallv(sbuf, scnts, sdsps, type,
		      rbuf, rcnts, rdsps, type, g->comm);
	for (i = 0; i < rlen; i++) {
	    switch (rbuf[i].type) {
		case VERTEX:	rbuf[i].index = map_vert[rbuf[i].index]; break;
		case EDGE:	rbuf[i].index = map_edge[rbuf[i].index]; break;
		case FACE:	rbuf[i].index = map_face[rbuf[i].index]; break;
		case VOLUME:	rbuf[i].index = map_elem[rbuf[i].index]; break;
	    }
	}
	/* FIXME: only need to send back rbuf[].index */
	MPI_Alltoallv(rbuf, rcnts, rdsps, type,
		      sbuf, scnts, sdsps, type, g->comm);
	MPI_Type_free(&type);

	/* update indices */

	if (g1->owner_index_vert != NULL) {
	    for (j = 0; j < g1->nvert; j++) {
        	if ((ii = g1->owner_rank_vert[j]) == g->rank || ii == -1)
		    continue;
		pbuf = sbuf + (sdsps[ii]++);
		g1->owner_index_vert[j] = pbuf->index;
	    }
	}

	if (g1->owner_index_edge != NULL) {
	    for (j = 0; j < g1->nedge; j++) {
        	if ((ii = g1->owner_rank_edge[j]) == g->rank || ii == -1)
		    continue;
		pbuf = sbuf + (sdsps[ii]++);
		g1->owner_index_edge[j] = pbuf->index;
	    }
	}

	if (g1->owner_index_face != NULL) {
	    for (j = 0; j < g1->nface; j++) {
        	if ((ii = g1->owner_rank_face[j]) == g->rank || ii == -1)
		    continue;
		pbuf = sbuf + (sdsps[ii]++);
		g1->owner_index_face[j] = pbuf->index;
	    }
	}

	if (g1->owner_index_elem != NULL) {
	    for (j = 0; j < g1->nelem; j++) {
        	if ((ii = g1->owner_rank_elem[j]) == g->rank || ii == -1)
		    continue;
		pbuf = sbuf + (sdsps[ii]++);
		g1->owner_index_elem[j] = pbuf->index;
	    }
	}

	phgFree(scnts);
	phgFree(sbuf);
    }
#endif	/* USE_MPI */

cont:
    assert((g->flags & ELEM_FLAG));
    g1->roots = phgAlloc(g1->nroot * sizeof(*g1->roots));
    ForAllElements(g, e) {
	e1 = g1->roots + (map_elem == NULL ? e->index : map_elem[e->index]);
	*e1 = *e;
	e1->parent = NULL;
	e1->generation = 0;
	/* update neighbours */
	for (ii = 0; ii < NFace; ii++) {
	    if ((e->bound_type[ii] & REMOTE)) {
		e1->neighbours[ii] = NULL;
		continue;
	    }
	    else if ((e0 = e->neighbours[ii]) == NULL) {
		continue;
	    }
	    e1->neighbours[ii] = g1->roots + (map_elem == NULL ?
					e0->index : map_elem[e0->index]);
	}
	if (map_vert != NULL) {
	    for (ii = 0; ii < NVert; ii++)
		e1->verts[ii] = map_vert[e->verts[ii]];
	}
	if (map_edge != NULL && (g->flags & EDGE_FLAG)) {
	    for (ii = 0; ii < NEdge; ii++)
		e1->edges[ii] = map_edge[e->edges[ii]];
	}
	if (map_face != NULL && (g->flags & FACE_FLAG)) {
	    for (ii = 0; ii < NFace; ii++)
		e1->faces[ii] = map_face[e->faces[ii]];
	}
	if (map_elem != NULL && (g->flags & ELEM_FLAG))
	    e1->index = map_elem[e->index];
    }

#if USE_MPI
    if (g1->nprocs > 1) {
#if 0
	phgUpdateRNeighbours_(g1);
#else	/* 0 */
	ELEMENT **sbuf, **rbuf;
	MPI_Datatype type;

	g1->neighbours.counts = phgAlloc(g1->nprocs *
						sizeof(*g1->neighbours.counts));
	g1->neighbours.displs = phgAlloc(g1->nprocs *
						sizeof(*g1->neighbours.displs));
	memcpy(g1->neighbours.counts, g->neighbours.counts,
			g1->nprocs * sizeof(*g1->neighbours.counts));
	memcpy(g1->neighbours.displs, g->neighbours.displs,
			g1->nprocs * sizeof(*g1->neighbours.displs));

	g1->neighbours.allocated = g1->neighbours.count = g->neighbours.count;
	if (g1->neighbours.count > 0) {
	    g1->neighbours.list = phgAlloc(g1->neighbours.count *
						sizeof(*g1->neighbours.list));
	    memcpy(g1->neighbours.list, g->neighbours.list,
			g1->neighbours.count * sizeof(*g1->neighbours.list));
	    sbuf = phgAlloc(2 * g1->neighbours.count * sizeof(*sbuf));
	    rbuf = sbuf + g1->neighbours.count;
	    for (i = 0; i < g1->neighbours.count; i++)
		sbuf[i] = g1->neighbours.list[i].remote;
	    MPI_Type_contiguous(sizeof(*sbuf), MPI_BYTE, &type);
	    MPI_Type_commit(&type);
	    MPI_Alltoallv(
		sbuf, g1->neighbours.counts, g1->neighbours.displs, type,
		rbuf, g1->neighbours.counts, g1->neighbours.displs, type,
		g1->comm);
	    for (i = 0; i < g1->neighbours.count; i++) {
		j = rbuf[i]->index;
		rbuf[i] = g1->roots + (map_elem == NULL ? j : map_elem[j]);
	    }
	    MPI_Alltoallv(
		rbuf, g1->neighbours.counts, g1->neighbours.displs, type,
		sbuf, g1->neighbours.counts, g1->neighbours.displs, type,
		g1->comm);
	    MPI_Type_free(&type);
	    for (i = 0; i < g1->neighbours.count; i++) {
		RNEIGHBOUR *rn = g1->neighbours.list + i;
		rn->remote = sbuf[i];
		j = rn->local->index;
		rn->local = g1->roots + (map_elem == NULL ? j : map_elem[j]);
		rn->local->neighbours[rn->vertex] = (void *)((size_t)i);
	    }
	    phgFree(sbuf);
	}
#endif	/* 0 */
    }
#endif  /* USE_MPI */

    g1->elems = phgAlloc(g1->nelem * sizeof(*g1->elems));
    for (i = 0; i < g1->nroot; i++)
	g1->elems[i] = g1->roots + i;

    if ((g1->flags & GEOM_FLAG))
	phgGeomInit(g1);

    phgFree(map_vert);
    phgFree(map_edge);
    phgFree(map_face);
    phgFree(map_elem);

    return g1;
}
