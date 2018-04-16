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

/* $Id: option.c,v 1.173 2016/07/24 01:27:52 zlb Exp $ */

#include "phg.h"
#include "phg/option.h"
#include "phg/io.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static OPTION *options = NULL;
static int *options_index = NULL;
static size_t options_count = 0, options_allocated = 0;

/* 'key' is initialized to point to the last dummy entry in the option list
 * by phgOptionsParseCmdline. It is used for bsearch */
static OPTION *key = NULL;

static char *help_category = NULL;	/* will be set by '-help' option */

/* Before phgInit:
 * 	preset options are stored in argc_new/argv_new[].
 * After phgInit:
 *	the original cmdline options are saved in 'argc_save/argv_save',
 *	the '-oem_options' are saved in 'argc_new/argv_new[]'. */
static int argc_save = 0, argc_new = 0, argc_new_allocated = 0;
static char **argv_save = NULL, **argv_new = NULL;

/* Note: we use the address of the dummy variable 'null' to denote a string
 * argument whose value is NULL */
#define NULL_STRING	((void *)&dummy_null)
static char dummy_null;

static void save_option(int k, size_t size, void *ptr);
static BOOLEAN set_option0(OPTION *o, void *value);
static BOOLEAN set_option(const char *op_name, void *value, int type,
			  const char *func);
static void get_option0(OPTION *o, void **pvar);

/*---------------------------------------------------------------------------*/

/* list of obsoleted option names */
typedef struct {
    const char	*old, *new;
    BOOLEAN	warned;
} OBSOLETE_OPTION;
OBSOLETE_OPTION obsolete_options[] = {
    {"default_solver",			"solver",		FALSE},
    {"default_dof_type",		"dof_type",		FALSE},
    {"default_refine_strategy",		"refine_strategy",	FALSE},
    {"default_coarsen_strategy",	"coarsen_strategy",	FALSE},
    {"default_strategy",		"strategy",		FALSE},
    {"eigen_rtol",			"eigen_tol",		FALSE}
};
static int obsolete_options_n =
		sizeof(obsolete_options) / sizeof(obsolete_options[0]);

static int
obsolete_option_comp(const void *p0, const void *p1)
{
    return strcmp(((const OBSOLETE_OPTION *)p0)->old,
		  ((const OBSOLETE_OPTION *)p1)->old);
}

static const char *
obsolete_option_check(const char *opt)
/* check for obsolete options */
{
    OBSOLETE_OPTION key, *o;

    key.old = opt;
    o = bsearch(&key, obsolete_options, obsolete_options_n, sizeof(*o),
		obsolete_option_comp);
    if (o == NULL)
	return opt;
    if (!o->warned) {
	o->warned = TRUE;
	phgPrintf(
	    "\n************************************************************\n"
	    "Warning: the option \"-%s\" is obsolete,\n"
	    "the new option name \"-%s\" should be used instead.\n"
	    "*************************************************************\n\n",
		  opt, o->new);
    }
    return o->new;
}

/*---------------------------------------------------------------------------*/

static void
parse_options(int *argc, char ***argv, int *size, const char *optstr)
/* parses the space separated words in 'optstr' and appends them to the list
 * (*argc)/(*argv)[] */
{
    char quote = '\0', c, *p, *q;
    const char *optstr0 = optstr, *r;
    int ac = *argc;

    p = phgAlloc(strlen(optstr) + 1);
    while (TRUE) {
	while (isspace(*(const BYTE *)optstr))
	    optstr++;
	if (*optstr == '#') {
	    /* skip to eol */
	    while (*(++optstr) != '\n' && *optstr != '\0');
	    if (*optstr == '\n')
		continue;
	}
	if (*optstr == '\0')
	    break;
	q = p;
	while (TRUE) {
	    if ((c = *optstr) == '\0' || (quote == '\0' && isspace((BYTE)c)))
		break;
	    if (c == quote) {
		quote = '\0';
		optstr++;
		continue;
	    }
	    if (c == '\\') {
		if ((c = *(++optstr)) == '\0')
		    break;
		*(q++) = c;
		optstr++;
		continue;
	    }
	    if (quote == '\0' && (c == '\'' || c == '"')) {
		r = optstr;
		while (--r >= optstr0 && *r == '\\');
		if (((int)(optstr - r - 1)) % 2 == 0) {
		    quote = c;
		    optstr++;
		    continue;
		}
	    }
	    *(q++) = c;
	    optstr++;
	}
	if (quote != '\0') {
	    if (c == '\0')
		phgError(1, "invalid string: %s\n", optstr0);
	    optstr++;
	}
	*q = '\0';
	if (ac >= *size - 1) {
	    *argv = phgRealloc_(*argv, (*size + 16) * sizeof(**argv),
					(*size) * sizeof(**argv));
	    *size += 16;
	}
	(*argv)[ac++] = strdup(p);
    }
    phgFree(p);

    if (ac >= *size) {
	*argv = phgRealloc_(*argv, ((*size) + 1) * sizeof(**argv),
				   (*size) * sizeof(**argv));
	++(*size);
    }
    (*argv)[ac] = NULL;
    if (*argc < ac)
	*argc = ac;
}

static void
option_free(OPTION *o)
{
    char **p;

    phgFree(o->name);
    phgFree(o->help);
    o->name = o->help = NULL;
    if (o->type == VT_STRING || o->type == VT_FILENAME) {
	phgFree(*(char **)o->var);
	*(char **)o->var = NULL;
    }
    if (o->keys != NULL) {
	for (p = o->keys; *p != NULL; p++)
	    phgFree(*p);
	phgFree(o->keys);
	o->keys = NULL;
    }
}

static const char *
default_mapping_handler(OPTION *o, const char *arg, int prefix)
{
    int i;
    DOF_TYPE **t = NULL;

    Unused(prefix);

    if (arg == NULL) {
	/* No need to repeat the list of DOF_TYPES of "-dof_type" */
    }
    else if (arg[0] != '\0') {	/* Note the convention with arg[0]=='\0' */
	for (i = 0; phgDofTypesList[i].types_array != NULL; i++) {
	    t = phgDofTypesList[i].types_array;
	    for (; *t == NULL; t++);	/* skip leading NULLs */
	    for (; *t != NULL; t++) {
		if (!strcmp((*t)->name, arg))
		    break;
	    }
	    if (*t != NULL)
		break;
	}
	if (t == NULL || *t == NULL)
	    return NULL;
	MAPPING_DEFAULT = *t;
    }

    return MAPPING_DEFAULT->name;
}

static const char *
default_dof_handler(OPTION *o, const char *arg, int prefix)
{
    char *s;
    int i, m;
    DOF_TYPE **t = NULL, **t0 = NULL;

    Unused(prefix);

    if (arg == NULL) {
	/* print help message */
	phgOptionsPrintfHelp(o, "  List of available DOF types:");
	s = phgAlloc(128);
	for (i = 0; phgDofTypesList[i].types_array != NULL; i++) {
	    t0 = phgDofTypesList[i].types_array;
	    for (; *t0 == NULL; t0++);		/* skip leading NULLs */
	    for (t = t0; *t != NULL; t++);
	    m = t - t0 - 2;
	    sprintf(s, "    %s \"%s\"", phgDofTypesList[i].desc,
					(*(t0++))->name);
	    if (*t0 != NULL) {
		sprintf(s + strlen(s), ", \"%s\"", (*(t0++))->name);
		if (m <= 3) {
		    for (; *t0 != NULL;)
			sprintf(s + strlen(s), ", \"%s\"", (*(t0++))->name);
		}
		else {
		    sprintf(s + strlen(s), " ... \"%s\"", (*(--t))->name);
		}
	    }
	    phgOptionsPrintfHelp(o, s);
	}
	phgFree(s);
    }
    else if (arg[0] != '\0') {	/* Note the convention with arg[0]=='\0' */
	for (i = 0; phgDofTypesList[i].types_array != NULL; i++) {
	    t = phgDofTypesList[i].types_array;
	    for (; *t == NULL; t++);	/* skip leading NULLs */
	    for (; *t != NULL; t++) {
		if (!strcmp((*t)->name, arg))
		    break;
	    }
	    if (*t != NULL)
		break;
	}
	if (t == NULL || *t == NULL)
	    return NULL;
	DOF_DEFAULT = *t;
    }

    return DOF_DEFAULT->name;
}

static const char *
preset_only_handler(OPTION *o, const char *arg, int prefix)
/* This option handler will raise an error if the user tries to set the option
 * after phgInit. */
{
    Unused(prefix);

    if (arg == NULL) {
	/* print help (nothing to do) */
    }
    else if (arg[0] != '\0') {
	if (phgInitialized == TRUE) {
	    phgInfo(-1, "should not set \"-%s\" after phgInit.\n", o->name);
	    return NULL;
	}
    }

    return o->keys == NULL ? NULL : o->keys[0];
}

static void
parse_environment(void)
/* parses options in the environment variable PHG_OPTIONS */
{
    static BOOLEAN initialized = FALSE;
    const char *p;
#if USE_MPI
    char *tmp = NULL;
    int n, mpi_initialized;
#endif	/* USE_MPI */


    if (initialized)
	return;
    initialized = TRUE;

    p = getenv("PHG_OPTIONS");
#if USE_MPI
    /* bcast the value in proc 0 to ensure consistency across processes. */
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
	phgError(1, "%s:%d: unexpected error (shouldn't happen).\n",
			__FILE__, __LINE__);
    }
    else {
	if (phgRank == 0)
	    n = (p == NULL ? 0 : strlen(p));
	MPI_Bcast(&n, 1, MPI_INT, 0, phgComm);
	if (n == 0) {
	    p = NULL;
	}
	else {
	    tmp = malloc(n + 1);
	    if (phgRank == 0)
		strcpy(tmp, p);
	    MPI_Bcast(tmp, n + 1, MPI_BYTE, 0, phgComm);
	    p = tmp;
	}
    }
#endif	/* USE_MPI */
    if (p != NULL) {
	while (isspace(*(const BYTE *)p))
	    p++;
	if (*p != '\0') {
	    phgPrintf("PHG_OPTIONS = '%s'\n", p);
	    parse_options(&argc_new, &argv_new, &argc_new_allocated, p);
	}
    }
#if USE_MPI
    phgFree(tmp);
#endif	/* USE_MPI */
}

static void
option_register(const char *name, const char *help, const char **keys,
		void *var, VTYPE type, BOOLEAN append)
{
    OPTION *o;
    static BOOLEAN initialized = FALSE;

    if (name != NULL && phgInitialized && phgRank == 0) {
	phgPrintf("WARNING: phgOptionsRegisterXX() must be called before "
		  "phgInit().\n");
	phgPrintf("Option \"-%s\" not registered.\n", name);
	return;
    }

    if (!initialized && type != VT_INIT) {
	initialized = TRUE;
        /* Register title for user options */
	option_register("\nUser options:", "\n", NULL, "user", VT_TITLE, FALSE);
    }
    else if (type == VT_INIT) {
	extern char *phgLogFilename;	/* defined in utils.c */
	extern char *phgOutFilename;	/* defined in utils.c */
	extern BOOLEAN _phg_pause;	/* pause before execution flag */
	/* register some global options */
	initialized = TRUE;
	phgOptionsRegisterTitle("\nGeneric options:", "\n", "generic");
	phgOptionsRegisterString("-help", "Print options help then exit",
			&help_category);
	phgOptionsRegisterNoArg("-pause", "Pause before execution",
			&_phg_pause);
	phgOptionsRegisterInt("-verbosity", "Verbosity level", &phgVerbosity);
	option_register("-options_file", "Read options from given file",
			NULL, NULL, VT_INCLUDE, FALSE);
	phgOptionsRegisterFilename("-log_file",
			"Log filename (may be \"stderr\" or \"stdout\")",
			&phgLogFilename);
	phgOptionsRegisterFilename("-output_file",
			"Output filename (may be \"stdout\" or \"stderr\")",
			&phgOutFilename);
	phgOptionsRegisterHandler("-oem_options",
			"Options passed to external packages",
			preset_only_handler, TRUE);
	phgOptionsRegisterHandler("-mapping_type",
		"Default Mapping type (see list of available DOF types below)",
			default_mapping_handler, FALSE);
	phgOptionsRegisterHandler("-dof_type", "Default DOF type",
			default_dof_handler, FALSE);
	phgMapCreate(NULL, NULL);	/* register MAP options */
	phgImport(NULL, NULL, TRUE);	/* register import options */
	phgTrapSignals_(TRUE);		/* register -fpetrap */
	phgPerfInit();			/* register performance options */
	phgMatDumpMATLAB(NULL, NULL, NULL);	/* register -matlab_digits */
#if USE_MPI
	phgInitMPI(NULL, NULL);			/* register MPI options */
#endif	/* USE_MPI */
	phgIOOpen(NULL, NULL);			/* Register I/O options */
	phgExportVTK(NULL, NULL, NULL, NULL);	/* register VTK options */
	return;
    }

    if (name == NULL)
	return;

    if (*name == '-' || *name == '+')
	name++;

    if (options_index != NULL) {
	/* invalidate options_index */
	phgFree(options_index);
	options_index = NULL;
    }

    if (options_count >= options_allocated) {
	options = phgRealloc_(options,
			(options_allocated + 8) * sizeof(*options),
			options_allocated * sizeof(*options));
	options_allocated += 8;
    }

    o = options + (options_count++);

    o->name = strdup(name);
    o->help = help == NULL ? NULL : strdup(help);
    o->keys = NULL;
    o->var = var;
    o->type = type;
    o->used = FALSE;
    o->append = append;

    if (type == VT_STRING || type == VT_FILENAME) {
	if (*(char **)o->var != NULL) {
	    /* duplicate the string (it is then safe to free it) */
	    *((char **)o->var) = strdup(*(char **)o->var);
	}
    }
    else if (type == VT_KEYWORD) {
	/* make a copy of the keywords list */
	const char **p;
	char **q;

	if (keys == NULL) {
	    phgPrintf("phgOptionsRegisterKeyword(): keys should not be NULL "
			"(option \"-%s\").\n", name);
	    phgPrintf("Option not registered.\n");
	    option_free(options + options_count);
	    options_count--;
	    return;
	}

#if 0
	if (keys[0] == NULL) {
	    /* no keys defined, delete this option */
	    option_free(options + options_count);
	    options_count--;
	    return;
	}
#endif

	for (p = keys; *p != NULL; p++);
	o->keys = phgAlloc((p - keys + 1) * sizeof(*keys));
	for (p = keys, q = o->keys; *p != NULL; p++, q++) {
	    if ((*p)[0] == '\0')
		phgPrintf("WARNING: empty string in the keywords list for the "
			  "option \"-%s\".\n", o->name);
	    *q = strdup(*p);
	}
	*q = NULL;
    }
}

/* Wrapper functions for enforcing prototype checking */

void
phgOptionsRegisterInit_(void)
{
    option_register(NULL, NULL, NULL, NULL, VT_INIT, FALSE);
}

void
phgOptionsRegisterNoArg(const char *name, const char *help, BOOLEAN *var)
{
    option_register(name, help, NULL, var, VT_NONE, FALSE);
}

void
phgOptionsRegisterInt(const char *name, const char *help, INT *var)
{
    option_register(name, help, NULL, var, VT_INT, FALSE);
}

void
phgOptionsRegisterFloat(const char *name, const char *help, FLOAT *var)
{
    option_register(name, help, NULL, var, VT_FLOAT, FALSE);
}

void
phgOptionsRegisterString(const char *name, const char *help, char **var)
{
    option_register(name, help, NULL, var, VT_STRING, FALSE);
}

void
phgOptionsRegisterFilename(const char *name, const char *help, char **var)
{
    option_register(name, help, NULL, var, VT_FILENAME, FALSE);
}
void
phgOptionsRegisterKeyword(const char *name, const char *help,
			  const char **keys, int *var)
{
    option_register(name, help, keys, var, VT_KEYWORD, FALSE);
}

void
phgOptionsRegisterTitle(const char *str, const char *help, const char *category)
{
    /* Note: category will be stored in '->var' */
    option_register(str, help, NULL, (void *)category, VT_TITLE, FALSE);
}

void
phgOptionsRegisterHandler(const char *name, const char *help,
                                OPTION_HANDLER func, BOOLEAN append)
{
    option_register(name, help, NULL, func, VT_HANDLER, append);
}

void
phgOptionsPreset(const char *str)
/* Note: the preset options are saved in the list argc_new/argv_new[], the
 * latter will be processes before the system 'argc/argv[] */
{
    if (phgInitialized && phgRank == 0) {
	phgError(1, "phgOptionsPreset must be called before phgInit!\n");
    }

    parse_options(&argc_new, &argv_new, &argc_new_allocated, str);
}

static int
option_comp(const void *i0, const void *i1)
{
    OPTION *o0 = options + *(int *)i0, *o1 = options + *(int *)i1;

    /* This makes all separators to appear at the end of the list */
    if (o0->type == VT_TITLE || o1->type == VT_TITLE) {
	if (o0->type != VT_TITLE)
	    return -1;
	if (o1->type != VT_TITLE)
	    return 1;
	return (int)(o0 - o1);
    }

    return strcmp(o0->name, o1->name);
}

static void
option_sort(void)
{
    int i;

    if (options_index != NULL)
	return;

    /* append a dummy entry at the end of the list as the key for bsearch */
    option_register("dummy", NULL, NULL, NULL, VT_NONE, FALSE);
    options_index = phgAlloc(options_count * sizeof(*options_index));
    for (i = 0; i < options_count; i++)
	options_index[i] = i;
    options_count--;
    qsort(options_index, options_count, sizeof(*options_index), option_comp);
    phgFree(options[options_count].name);
    options[options_count].name = NULL;
}

#define STACK_MAX	16
static void **stack[STACK_MAX];
static BYTE *stack_used[STACK_MAX];
static int pstack = 0;
#define EIGHT		(sizeof(BYTE))

void
phgOptionsPush(void)
{
    if (pstack >= STACK_MAX)
	phgError(1, "options stack overflow.");
    stack[pstack] = phgCalloc(options_count, sizeof(*stack));
    stack_used[pstack] = phgCalloc((options_count + EIGHT - 1) / EIGHT,
				sizeof(*stack_used));
    pstack++;
}

void
phgOptionsPop(void)
{
    BOOLEAN append;
    int i, ps;
    OPTION *o;

    if (pstack <= 0)
	phgError(1, "options stack underflow.");
    ps = pstack - 1;
    pstack = 0;		/* Neutralize save_option() */
    for (i = 0, o = options; i < options_count; i++, o++) {
	if (stack[ps][i] != NULL) {
	    append = o->append;	/* save o->append */
	    o->append = FALSE;
	    set_option0(o, stack[ps][i]);
	    o->append = append;	/* restore o->append */
	    if (stack[ps][i] != NULL_STRING)
		phgFree(stack[ps][i]);
	    o->used = ((stack_used[ps][i / EIGHT] & (1 << (i % EIGHT))) != 0);
	}
    }
    phgFree(stack[ps]);
    phgFree(stack_used[ps]);
    pstack = ps;
}

static void
save_option(int k, size_t size, void *ptr)
{
    int i;

    assert(k >= 0 && k < options_count);
    if (pstack <= 0)
	return;
    i = pstack - 1;
    if (stack[i][k] != NULL)
	return;
    if (size > 0) {
	stack[i][k] = phgAlloc(size);
	memcpy(stack[i][k], ptr, size);
    }
    else {
	stack[i][k] = NULL_STRING;
    }

    /* save the used flag */
    if (options[k].used == TRUE)
	stack_used[i][k / EIGHT] |= (1 << (k % EIGHT));
}

void
phgOptionsReset(void)
{
    OPTION *o;
    int i;

    if (pstack > 0) {
	phgWarning("phgOptionsReset: non empty options stack.\n");
	while (pstack > 0)
	    phgOptionsPop();
    }

    if (options != NULL) {
	for (o = options; o < options + options_count; o++)
	    option_free(o);
	phgFree(options);
	phgFree(options_index);
	options = NULL;
	options_index = NULL;
	options_count = options_allocated = 0;
    }

    if (argv_new != NULL) {
#if 0
	if (!(USE_PETSC && phgInitialized) && argc_new != 1) {
	    for (i = 1; i < argc_new; i++)
		if (argv_new[i][0] == '-' || argv_new[i][0] == '+')
		    break;
	    if (i < argc_new) {
		phgPrintf("WARNING: the following options were not used "
			  "(misspelled?):\n   ");
		for (; i < argc_new; i++)
		    phgPrintf(" %s", argv_new[i]);
		phgPrintf("\n");
	    }
	}
#endif

	for (i = 0; i < argc_new; i++)
	    phgFree(argv_new[i]);
        phgFree(argv_new);
        argv_new = NULL;

	for (i = 0; i < argc_save; i++)
	    phgFree(argv_save[i]);
        phgFree(argv_save);
        argv_save = NULL;

        argc_new = argc_new_allocated = argc_save = 0;
    }
}

void
phgOptionsPrintfHelp(OPTION *o, const char *help)
/* format and print the help text of option 'o' */
{
    static char indent[] = "     ";
    char buffer[4096];
    char c, *p;
    const char *p0;
    size_t len;

    if (help == NULL)
	return;

    /* copy help text, expanding tabs */
    len = 0;
    p = buffer;
    p0 = help;
    while (*p0 != '\0' && p - buffer < sizeof(buffer) - 1) {
	c = *(p0++);
	if (c == '\t') {
	    *(p++) = ' ';
	    len++;
	    while (len % 8 != 0 && p - buffer < sizeof(buffer) - 1) {
		*(p++) = ' ';
		len++;
	    }
	    continue;
	}
	*(p++) = c;
	len++;
	if (c == '\n')
	    len = 0;
    }
    if (o->type == VT_KEYWORD) {
	char **pp;
	assert(p + 3 < buffer + sizeof(buffer));	/* +1 for '>' */
	memcpy(p, " <\"", 3);
	p += 3;
	for (pp = o->keys; *pp != NULL; pp++) {
	    if (pp > o->keys) {
		assert(p + 5 < buffer + sizeof(buffer));
		memcpy(p, "\", \"", 4);
		p += 4;
	    }
	    len = strlen(*pp);
	    assert(p + len + 1 < buffer + sizeof(buffer));
	    memcpy(p, *pp, len);
	    p += len;
	}
	*(p++) = '"';
	*(p++) = '>';
    }
    else if (o->type == VT_NONE) {
	p0 = " (the opposite is \"+%s\")";
	len = strlen(p0) - 2 + strlen(o->name);;
	assert(p + len + 1 < buffer + sizeof(buffer));
	sprintf(p, p0, o->name);
	p += len;
    }
    *p = '\0';

#if 0
    phgPrintf("%s%s\n", indent, buffer);
    return;
#endif
    len = sizeof(indent) - 1;
    phgPrintf(indent);
    p0 = p = buffer;
    while (*p != '\0') {
	while (*p != '\0' && *p != '\n' && isspace(*(BYTE *)p))
	    p++;
	while (*p != '\0' && *p != '\n' && !isspace(*(BYTE *)p))
	    p++;
	c = *p;
	if (*p != '\0')
	    *p = '\0';
	if (p == p0)
	    break;
	if (p - p0 + len + 1 >= 78 && len >= sizeof(indent)) {
	    phgPrintf("\n%s", indent);
	    len = sizeof(indent) - 1;
	    while (isspace(*(const BYTE *)p0))
		p0++;
	}
	phgPrintf("%s", p0);
	if (c != '\0') {
	    phgPrintf("%c", c);
	    *(p++) = c;
	}
	len += p - p0;
	if (c == '\n') {
	    phgPrintf("%s", indent);
	    len = sizeof(indent) - 1;
	    while (isspace(*(BYTE *)p))
		p++;
	}
	p0 = p;
    }
    phgPrintf("%s\n", p0);
}

void
phgOptionsShowCmdline(void)
{
    int i;

    if (phgRank > 0 || argc_save <= 1)
	return;

    phgPrintf("Command-line:");
    for (i = 0; i < argc_save; i++)
	phgPrintf(" %s", argv_save[i]);
    phgPrintf("\n");
}

void
phgOptionsShowUsed(void)
/* prints all options which have been called by the user */
{
    BOOLEAN flag = FALSE;
    OPTION *o;
    char **pp;
    const char *cat = NULL, *old_cat = NULL;

    if (!phgInitialized)
	phgError(1, "%s must be called after phgInit!\n", __func__);

    for (o = options; o < options + options_count; o++) {
	if (o->type == VT_TITLE)
	    cat = o->var;
	if (o->used == FALSE || o->type == VT_INCLUDE || o->type == VT_TITLE)
	    continue;
	if (!flag) {
	    phgPrintf("*-------------------- "
			"Parameter(s) set through PHG options "
			"--------------------\n");
	    flag = TRUE;
	}
	if (cat != NULL && cat != old_cat)
	    phgPrintf("* Category '%s':\n", old_cat = cat);
	switch (o->type) {
	    case VT_INIT:
		phgError(1, "unexpected.\n");
		break;
	    case VT_NONE:
		phgPrintf("*   %s: %s\n", o->help == NULL ? o->name : o->help,
				*(BOOLEAN *)o->var == TRUE ? "True" : "False");
		break;
	    case VT_INT:
		phgPrintf("*   %s: %"dFMT"\n",
				o->help == NULL ? o->name : o->help,
				*(INT *)o->var);
		break;
	    case VT_FLOAT:
		phgPrintf("*   %s: %lg\n", o->help == NULL ? o->name : o->help,
				(double)*(FLOAT *)o->var);
		break;
	    case VT_STRING:
		phgPrintf("*   %s: %s\n", o->help == NULL ? o->name : o->help,
				*(char **)o->var == NULL ? 
					"<null>" : *(char **)o->var);
		break;
	    case VT_FILENAME:
		phgPrintf("*   %s: %s\n", o->help == NULL ? o->name : o->help,
				*(char **)o->var == NULL ?
					"<null>" : *(char **)o->var);
		break;
	    case VT_KEYWORD:
		phgPrintf("*   %s: %s\n", o->help == NULL ? o->name : o->help,
				*(int *)o->var < 0 ?
					"<null>" : (o->keys)[*(int *)o->var]);
		break;
	    case VT_HANDLER:
		phgPrintf("*   %s", o->help == NULL ? o->name : o->help);
		if ((pp = o->keys) != NULL && *pp != NULL) {
		    phgPrintf(": %s", *pp);
		}
#if 0
		/* not necessary */
		else if (o->var != NULL) {
		    const char *p = ((OPTION_HANDLER)(o->var))(o, "", 0);
		    if (p != NULL)
			phgPrintf(": %s", p);
		}
#endif
		phgPrintf("\n");
		break;
	    case VT_INCLUDE:
	    case VT_TITLE:
		break;
	}
    }

    if (flag)
	phgPrintf("*-----------------------------------------------------"
		  "-------------------------\n");
}

BOOLEAN
phgOptionsIfUsed(const char *op_name)
/* returns TRUE if the option "opstr" has been called by the user. */
{
    int *k;
    OPTION *o;

    if (!phgInitialized)
	phgError(1, "%s must be called after phgInit!\n", __func__);

    if (op_name[0] == '-' || op_name[1] == '+')
	op_name++;

    key->name = (void *)obsolete_option_check(op_name);
    k = bsearch(options_index + options_count, options_index,
		options_count, sizeof(*options_index), option_comp);
    key->name = NULL;	/* reset key->name */
    if (k == NULL)
	phgError(1, "%s: unknown option \"-%s\"!\n", __func__, op_name);
    o = options + (*k);

    return o->used;
}

static int
comp_string(const void *p1, const void *p2)
{
    return strcmp(*(char **)p1, *(char **)p2);
}

void
phgOptionsHelp(void)
{
    OPTION *o;
    char **pp;
    const char *p;
    BOOLEAN all_flag, flag, matched;
    char **list = NULL;
    int i, list_count = 0, list_allocated = 0;

    if (help_category == NULL)
	return;

#if USE_MPI
    if (phgRank > 0)
	goto end;
#endif	/* USE_MPI */

    flag = TRUE;
    matched = FALSE;
    all_flag = !strcmp(help_category, "all");
    for (o = options; o < options + options_count; o++) {
	if (!all_flag && o->type == VT_TITLE) {
	    flag = (o->var == NULL || !strcmp(help_category, o->var));
	    if (o->var != NULL &&
		(list_count == 0 || strcmp(list[list_count - 1], o->var))) {
		if (list_count >= list_allocated) {
		    list = phgRealloc_(list,
					(list_allocated + 128) * sizeof(*list),
					list_allocated * sizeof(*list));
		    list_allocated += 128;
		}
		list[list_count++] = o->var;
	    }
	}
	if (!flag)
	    continue;
	/* skip option without help text */
	if (o->help == NULL)
	    continue;
	matched = TRUE;
	switch (o->type) {
	    case VT_INIT:
		phgError(1, "unexpected.\n");
		break;
	    case VT_NONE:
		phgPrintf("  -%s (%s)", o->name,
				*(BOOLEAN *)o->var ? "True" : "False");
		break;
	    case VT_INT:
		phgPrintf("  -%s <integer> (%"dFMT")", o->name, *(INT *)o->var);
		break;
	    case VT_FLOAT:
		phgPrintf("  -%s <real> (%lg)", o->name,
					(double)*(FLOAT *)o->var);
		break;
	    case VT_STRING:
		if (*(char **)o->var == NULL)
		    phgPrintf("  -%s <string> (<null>)", o->name);
		else
		    phgPrintf("  -%s <string> (\"%s\")",
					o->name, *(char **)o->var);
		break;
	    case VT_FILENAME:
		if (*(char **)o->var == NULL)
		    phgPrintf("  -%s <filename> (<null>)", o->name);
		else
		    phgPrintf("  -%s <filename> (\"%s\")",
					o->name, *(char **)o->var);
		break;
	    case VT_KEYWORD:
		if (*(int *)o->var < 0)
		    phgPrintf("  -%s <keyword> (<null>)", o->name);
		else
		    phgPrintf("  -%s <keyword> (\"%s\")",
					o->name, (o->keys)[*(int *)o->var]);
		break;
	    case VT_HANDLER:
		phgPrintf("  -%s <string>", o->name);
		if ((pp = o->keys) != NULL && *pp != NULL)
		    p = *pp;
		else if (o->var != NULL)
		    p = ((OPTION_HANDLER)(o->var))(o, "", 0);
		else
		    p = NULL;
		p == NULL ? phgPrintf(" (none)") : phgPrintf(" (\"%s\")", p);
		break;
	    case VT_INCLUDE:
		phgPrintf("  -%s <filename>", o->name);
		break;
	    case VT_TITLE:
		phgPrintf("%s", o->name);
		if (o->var != NULL)
		    phgPrintf(" (category \"%s\")", (char *)o->var);
		break;
	}
	phgPrintf("\n");
	if (o->help != NULL)
	    phgOptionsPrintfHelp(o, o->help);
	if (o->type == VT_HANDLER && o->var != NULL)
	    ((OPTION_HANDLER)o->var)(o, NULL, 0);
    }

    if (matched) {
	phgPrintf("\n");
    }
    else {
	if (strcmp(help_category, "help"))
	    phgPrintf("Unknown help category '%s'.\n", help_category);
	qsort(list, list_count, sizeof(*list), comp_string);
	phgPrintf("Usage:\n    %s -help <category>\n"
		  "where <category> should be one of:\n", argv_save[0]);
	phgPrintf("    all");
	list_allocated = 7;
	for (i = 0; i < list_count; i++) {
	    if ((list_allocated += strlen(list[i]) + 2) > 78) {
		list_allocated = strlen(list[i]) + 4;
		phgPrintf(",\n    %s", list[i]);
	    }
	    else {
		phgPrintf(", %s", list[i]);
	    }
	}
	phgPrintf("\n");
    }

    phgFree(list);

#if USE_MPI
end:
#endif	/* USE_MPI */

    phgOptionsReset();

#if USE_MPI
    MPI_Initialized(&i);
    if (i)
	MPI_Finalize();
#endif	/* USE_MPI */

    exit(0);
}

static BOOLEAN options_file_quiet_flag = FALSE;

static void
process_options_file(const char *fn)
/* processes options from file 'fn' */
{
    FILE *f;
    int i, argc = 0, allocated = 0;
    char *p, **argv = NULL, buffer[4096];

    if ((f = fopen(fn, "r")) == NULL) {
	phgPrintf("Cannot open options file \"%s\".\n", fn);
	phgAbort(1);
    }

    if (options_file_quiet_flag == FALSE)
	phgPrintf("Processing options file \"%s\".\n", fn);

    while (TRUE) {
	if (fgets(buffer, sizeof(buffer), f) == NULL)
	    break;
	p = buffer;
	while (isspace(*(BYTE *)p))
	    p++;
	if (*p == '#' || *p == '\0')
	    continue;
	parse_options(&argc, &argv, &allocated, p);
    }
    fclose(f);
    if (allocated == 0)
	return;
    argv[argc] = NULL;
    phgOptionsParseCmdline(&argc, &argv);
    for (i = 0; i < argc; i++)
	phgFree(argv[i]);
    phgFree(argv);
}

void
phgOptionsParseCmdline(int *argc, char ***argv)
/* parses cmdline parameters, processes and removes known options from
   the argument list */
{
    static int depth = 0;		/* recursion depth */
    static BOOLEAN firstcall = TRUE;	/* 1st time calling this function */
    OPTION *o;
    char **pp;
    char *p, *arg;
    int i, j;
    int *k = NULL;		/* points to sorted indices of options */
    size_t size;

    if (depth >= 8)
	phgError(1, "too many nested options files.\n");

    if (depth == 0 && firstcall) {
	int argc_tmp = 0;
	char **argv_tmp = NULL;

	/* Parse the environment variable PHG_OPTIONS */
	parse_environment();

	/* sort list of obsolete options */
	qsort(obsolete_options, obsolete_options_n,
		sizeof(obsolete_options[0]), obsolete_option_comp);

	if (options == NULL) {
	    /* this will register the options '-help', '-include', etc. */
	    option_register(NULL, NULL, NULL, NULL, VT_NONE, FALSE);
	}
	option_sort();
	key = options + options_count;
	/* check and warn duplicate options */
	j = 0;
	for (i = 1; i < options_count; i++) {
	    o = options + options_index[i];
	    if (o->type == VT_TITLE)
		break;
	    if (option_comp(options_index + i, options_index + j) == 0 &&
		phgRank == 0)
		phgWarning("duplicate option \"-%s\".\n",o->name);
	    j = i;
	}

	argc_save = *argc;
	argv_save = phgAlloc((argc_save + 1) * sizeof(*argv_save));
	for (i = 0; i < argc_save; i++)
	    argv_save[i] = strdup((*argv)[i]);
	argv_save[i] = NULL;

	/* save preset options */
	if (argc_new != 0) {
	    argc_tmp = argc_new;
	    argv_tmp = argv_new;
	}

	argv_new = phgAlloc((argc_new_allocated = 16) * sizeof(*argv_new));
	argv_new[0] = strdup(argv_save[0]);
	argc_new = 1;

	/* process preset options */
	if (argc_tmp != 0) {
	    depth++;
	    phgOptionsParseCmdline(&argc_tmp, &argv_tmp);
	    depth--;
	    for (i = 0; i < argc_tmp; i++)
	        phgFree(argv_tmp[i]);
	    phgFree(argv_tmp);
	}

	/* process default options file "basename(progname).options"
	 * FIXME: do we need to add an option to disable this feature? */
	if (TRUE) {
	    FILE *f;
	    arg = strrchr((*argv)[0], '/');
	    arg = (arg == NULL ? (*argv)[0] : arg + 1);
	    i = (int)strlen(arg);
	    p = phgAlloc(i + 8 + 1);
	    memcpy(p, arg, i);
	    memcpy(p + i, ".options", 8 + 1);
	    if ((f = fopen(p, "r")) != NULL) {
		fclose(f);
		depth++;
		process_options_file(p);
		depth--;
	    }
	}
	phgFree(p);
    }

    depth++;

    for (i = ((depth > 1 || !firstcall) ? 0 : 1); i < *argc; i++) {
	char *q;
	static char *prefixes = "+-*/%|&^";
	char *prefix = prefixes + strlen(prefixes);
	o = NULL;
	arg = NULL;
	k = NULL;
	if ((p = (*argv)[i])[0] == '-' || p[0] == '+') {
	    q = strdup(p[0] == '-' && p[1] == '-' ? p + 2 : p + 1);
	    if ((arg = strchr(p + 1, '=')) != NULL) {
		/* process "-option[prefix]=value" */
		if ((prefix = strchr(prefixes, *(arg - 1))) != NULL) {
		    /* the option is in the form "-option {prefix}=" */
		    q[arg - p - 1 - 1] = '\0';
		}
		else {
		    /* the option is in the form: "-option =" */
		    prefix = prefixes + strlen(prefixes);
		    q[arg - p - 1] = '\0';
		}
		arg++;
	    }
	    key->name = (void *)obsolete_option_check(q);
	    k = bsearch(options_index + options_count, options_index,
			options_count, sizeof(*options_index), option_comp);
	    phgFree(q);
	    key->name = NULL;
	}
	else {
	    /* append non-option arguments to OEM options */
	    key->name = "oem_options";
	    arg = p;
	    k = bsearch(options_index + options_count, options_index,
			options_count, sizeof(*options_index), option_comp);
	    assert(k != NULL);
	    key->name = NULL;
	}
	if (k == NULL)
	    phgError(1, "unknown option \"%s\"!\n", p);
	o = options + (*k);
	/* process option */
	if (o->type != VT_NONE) {
	    static char lbrace = '{', rbrace = '}';
	    if (arg == NULL && (arg = (*argv)[++i]) == NULL) {
		if (!strcmp(o->name, help_category = strdup("help"))) {
		    phgPrintf("Missing argument for option \"%s\".\n", p);
		    phgOptionsHelp();
		}
		phgError(1, "missing argument for option \"%s\".\n", p);
	    }
	    while (isspace((int)(*arg)) && *arg != '\0')
		arg++;
	    /* check whether there're curly braces in the argument */
	    if (strchr(arg, lbrace) != NULL) {
		static char *buffer = NULL;
		static int buffer_size = 0;
		int len = 0, nest = 0;
		FreeAtExit(buffer);
		while (TRUE) {
		    j = strlen(arg);
		    if (len + j + 2 > buffer_size) {
			p = phgAlloc((buffer_size = len + j + 2));
			if (len > 0)
			    memcpy(p, buffer, len);
			phgFree(buffer);
			buffer = p;
		    }
		    /* process curly braces in the argument */
		    for (q = arg; *q != '\0'; q++) {
			if (*q == lbrace) {
			    if (nest == 0) {
				/* drop top level braces */
				if (q > arg) {
				    memcpy(buffer + len, arg, q - arg);
				    len += q - arg;
				}
				arg = q + 1;
			    }
			    nest++;
			}
			else if (*q == rbrace) {
			    nest--;
			    if (nest < 0)
				phgError(1, "badly nested braces: %s\n", arg);
			    if (nest == 0) {
				/* drop top level braces */
				if (q > arg) {
				    memcpy(buffer + len, arg, q - arg);
				    len += q - arg;
				}
				arg = q + 1;
			    }
			}
		    }
		    if (q > arg) {
			memcpy(buffer + len, arg, q - arg);
			len += q - arg;
		    }
		    if (nest == 0)
			break;
		    /* insert a space between arguments */
		    buffer[len++] = ' ';
		    if ((arg = (*argv)[++i]) == NULL) {
			phgError(1, "invalid argument in the option \"-%s\".\n",
					o->name);
		    }
		}	/* while */
		buffer[len] = '\0';
		arg = buffer;
	    }
	}
	/* FIXME: use the same code for the loop below and the loop in
	 * set_option() */
	switch (o->type) {
	    case VT_INIT:
		phgError(1, "unexpected.\n");
		break;
	    case VT_NONE:
		save_option(*k, sizeof(BOOLEAN), o->var);
		if (*prefix != '\0')
		    phgError(1, "Invalid option \"-%s%c=\".\n",
						o->name, *prefix);
		*(BOOLEAN *)o->var = (p[0] == '-' ? TRUE : FALSE);
		o->used = TRUE;
		break;
	    case VT_INT:
		save_option(*k, sizeof(INT), o->var);
		switch (*prefix) {
		    case '\0':	*(INT *)o->var  = atoi(arg); break;
		    case '+':	*(INT *)o->var += atoi(arg); break;
		    case '-':	*(INT *)o->var -= atoi(arg); break;
		    case '*':	*(INT *)o->var *= atoi(arg); break;
		    case '/':	*(INT *)o->var /= atoi(arg); break;
		    case '%':	*(INT *)o->var %= atoi(arg); break;
		    case '|':	*(INT *)o->var |= atoi(arg); break;
		    case '&':	*(INT *)o->var &= atoi(arg); break;
		    case '^':	*(INT *)o->var ^= atoi(arg); break;
		    default:	phgError(1, "Invalid option \"-%s%c=\".\n",
							o->name, *prefix, arg);
		}
		o->used = TRUE;
		break;
	    case VT_FLOAT:
		save_option(*k, sizeof(FLOAT), o->var);
		switch (*prefix) {
		    case '\0':	*(FLOAT *)o->var  = atof(arg); break;
		    case '+':	*(FLOAT *)o->var += atof(arg); break;
		    case '-':	*(FLOAT *)o->var -= atof(arg); break;
		    case '*':	*(FLOAT *)o->var *= atof(arg); break;
		    case '/':	*(FLOAT *)o->var /= atof(arg); break;
		    default:	phgError(1, "Invalid option \"-%s%c=\".\n",
							o->name, *prefix, arg);
		}
		o->used = TRUE;
		break;
	    case VT_STRING:
	    case VT_FILENAME:
		p = *(char **)o->var;
		size = (p == NULL ? 0 : strlen(p) + 1);
		save_option(*k, size, p);
		if (*prefix != '\0') {
		    size_t size1;
		    if (*prefix != '+' || o->type == VT_FILENAME) 
			phgError(1, "Invalid option \"-%s%c=%s\".\n",
							o->name, *prefix, arg);
		    /* append 'arg' to existing string */
		    size1 = strlen(arg) + 1;
		    if (size == 1)
			size = 0;	/* empty string */
		    *(char **)o->var = p =
				phgRealloc_(p, size + size1, size);
		    if (size > 0)
			p[size - 1] = ' ';
		    /* remove multiple spaces */
		    while (size > 1 && p[size - 2] == ' ')
			size--;
		    memcpy(p + size, arg, size1);
		}
		else {
		    phgFree(*(char **)o->var);
		    *(char **)o->var = strdup(arg);
		}
		o->used = TRUE;
		break;
	    case VT_KEYWORD:
		p = (*(int *)o->var < 0 ? NULL : o->keys[*(int *)o->var]);
		save_option(*k, p == NULL ? 0 : strlen(p) + 1, p);
		if (*prefix != '\0')
		    phgError(1, "Invalid option \"-%s%c=\".\n",
						o->name, *prefix, arg);
	 	if (arg == NULL) {
		    *(int *)o->var = -1;
		    o->used = TRUE;
		    break;
		}
		for (pp = o->keys; *pp != NULL; pp++)
		    if (!strcmp(*pp, arg))
			break;
		if (*pp == NULL) {
		    if (phgRank == 0) {
			phgPrintf("Invalid argument \"%s\" "
				"for the option \"-%s\".\n", arg, o->name);
			phgPrintf("Valid keywords are: ");
	    		for (pp = o->keys; *pp != NULL; pp++)
			    phgPrintf("%s\"%s\"", pp == o->keys ? "":", ", *pp);
			phgPrintf("\n");
		    }
		    phgError(1, "abort.\n");
		    break;
		}
		*(int *)o->var = pp - o->keys;
		o->used = TRUE;
		break;
	    case VT_HANDLER:
		/* Note: the argument is saved in o->keys[0] */
		get_option0(o, (void *)&p);
		save_option(*k, p == NULL ? 0 : strlen(p) + 1, p);
		if (o->keys == NULL) {
		    o->keys = phgCalloc(2, sizeof(*o->keys));
		    j = 0;
		}
		else if (o->append || *prefix == '+') {
		    j = strlen(o->keys[0]);
		    o->keys[0][j++] = ' ';
		}
		else {
		    j = 0;
		    phgFree(o->keys[0]);
		    o->keys[0] = NULL;
		}
		o->keys[0] = phgRealloc_(o->keys[0], j + strlen(arg) + 1, j);
		strcpy(o->keys[0] + j, arg);
		/* call user supplied option handler */
		if (o->var != NULL) {
		    if (((OPTION_HANDLER)o->var)(o, o->keys[0], *prefix)
				== NULL) {
			phgInfo(-1, "invalid argument for \"-%s\" option.\n",
					o->name);
			/* show the help message for this options */
	    		((OPTION_HANDLER)o->var)(o, NULL, 0);
			phgError(1, "abort.\n");
		    }
		}
		o->used = TRUE;
		break;
	    case VT_INCLUDE:
		process_options_file(arg);
		break;
	    case VT_TITLE:
		/* cannot happen */
		break;
	}
    }

    depth--;

    if (depth == 0 && firstcall) { /* Process OEM options */
	firstcall = FALSE;
	key->name = "oem_options";
	k = bsearch(options_index + options_count, options_index,
		    options_count, sizeof(*options_index), option_comp);
	key->name = NULL;
	assert(k != NULL);
	o = options + (*k);
	for (i = 0; o->keys != NULL && (p = o->keys[i]) != NULL; i++)
	    parse_options(&argc_new, &argv_new, &argc_new_allocated, p);
	argv_new[argc_new] = NULL;
	*argc = argc_new;
	*argv = argv_new;
    }

    return;
}

static BOOLEAN
set_option0(OPTION *o, void *value)
{
    int j, k;
    char *p, **pp;

    if (value == NULL)
	return TRUE;

    if (value == NULL_STRING)
	value = NULL;

    k = (int)(o - options);

    /* FIXME: use the same code for the loop below and the loop in
     * phgOptionsParseCmdline() */
    switch (o->type) {
	case VT_NONE:
	    save_option(k, sizeof(BOOLEAN), o->var);
	    *(BOOLEAN *)o->var = *(BOOLEAN *)value;
	    o->used = TRUE;
	    break;
	case VT_INT:
	    save_option(k, sizeof(INT), o->var);
	    *(INT *)o->var = *(INT *)value;
	    o->used = TRUE;
	    break;
	case VT_FLOAT:
	    save_option(k, sizeof(FLOAT), o->var);
	    *(FLOAT *)o->var = *(FLOAT *)value;
	    o->used = TRUE;
	    break;
	case VT_STRING:
	case VT_FILENAME:
	    p = *(char **)o->var;
	    save_option(k, p == NULL ? 0 : strlen(p) + 1, p);
	    phgFree(*(char **)o->var);
	    if (value == NULL)
		*(char **)o->var = NULL;
	    else
		*(char **)o->var = strdup((const char *)value);
	    o->used = TRUE;
	    break;
	case VT_KEYWORD:
	    p = (*(int *)o->var < 0 ? NULL : o->keys[*(int *)o->var]);
	    save_option(k, p == NULL ? 0 : strlen(p) + 1, p);
	    if (value == NULL) {
		*(int *)o->var = -1;
		o->used = TRUE;
		break;
	    }
	    for (pp = o->keys; *pp != NULL; pp++)
		if (!strcmp(*pp, value))
		    break;
	    if (*pp == NULL) {
		if (phgRank == 0) {
		    phgPrintf("%s:%d, invalid argument \"%s\" "
				"for the option \"-%s\".\n", __FILE__, __LINE__,
				value, o->name);
		    phgPrintf("Valid keywords are: ");
	    	    for (pp = o->keys; *pp != NULL; pp++)
			phgPrintf("%s\"%s\"", pp == o->keys ? "":", ", *pp);
		    phgPrintf("\n");
		}
		phgError(1, "abort.\n");
		break;
	    }
	    *(int *)o->var = pp - o->keys;
	    o->used = TRUE;
	    break;
	case VT_HANDLER:
	    /* Note: the argument is saved in o->keys[0] */
	    get_option0(o, (void *)&p);
	    save_option(k, p == NULL ? 0 : strlen(p) + 1, p);
	    if (o->keys == NULL) {
		o->keys = phgCalloc(2, sizeof(*o->keys));
		j = 0;
	    }
	    else if (o->append) {
		j = strlen(o->keys[0]);
		o->keys[0][j++] = ' ';
	    }
	    else {
		j = 0;
		phgFree(o->keys[0]);
		o->keys[0] = NULL;
	    }
	    o->keys[0] = phgRealloc_(o->keys[0], j + strlen(value) + 1, j);
	    strcpy(o->keys[0] + j, value);
	    /* call user supplied option handler */
	    if (o->var != NULL) {
		if (((OPTION_HANDLER)o->var)(o, o->keys[0], 0) == NULL) {
		    phgInfo(-1, "invalid argument for \"-%s\" option.\n",
				o->name);
		    /* show the help message for this options */
	    	    ((OPTION_HANDLER)o->var)(o, NULL, 0);
		    phgError(1, "abort.\n");
		}
	    }
	    o->used = TRUE;
	    break;
	default:
	    /* not allowed */
	    phgError(1, "%s:%d: unsupported or unimplemented option type.\n",
			__FILE__, __LINE__);
    }

    return TRUE;
}

static BOOLEAN
set_option(const char *op_name, void *value, int type, const char *func)
{
    int *k;
    OPTION *o;

    if (!phgInitialized)
	phgError(1, "%s must be called after phgInit!\n", func);

    if (value == NULL)
	return TRUE;

    if (op_name[0] == '-' || op_name[1] == '+')
	op_name++;

    key->name = (void *)op_name;
    k = bsearch(options_index + options_count, options_index,
		options_count, sizeof(*options_index), option_comp);
    key->name = NULL;	/* reset key->name */
    if (k == NULL)
	phgError(1, "%s: unknown option \"-%s\"!\n", func, op_name);
    o = options + (*k);
    if (type >= 0 && o->type != type) {
	phgInfo(-1, "%s: wrong function type for \"-%s\".", func, op_name);
	switch (o->type) {
	    case VT_NONE:
		phgError(1, "Please use phgOptionsSetNoArg instead.\n");
	    case VT_INT:
		phgError(1, "Please use phgOptionsSetInt instead.\n");
	    case VT_FLOAT:
		phgError(1, "Please use phgOptionsSetFloat instead.\n");
	    case VT_STRING:
		phgError(1, "Please use phgOptionsSetString instead.\n");
	    case VT_FILENAME:
		phgError(1, "Please use phgOptionsSetFilename instead.\n");
	    case VT_KEYWORD:
		phgError(1, "Please use phgOptionsSetKeyword instead.\n");
	    case VT_HANDLER:
		phgError(1, "Please use phgOptionsSetHandler instead.\n");
	    default:
		phgError(1, "No phgOptionsSetXXXX function for \"-%s\".\n",
			op_name);
	}
    }

    return set_option0(o, value);
}

void
phgOptionsSetOptions(const char *str)
{
    int i, argc = 0, argc_allocated = 0;
    char **argv = NULL;

    if (!phgInitialized)
	phgError(1, "%s must be called after phgInit!\n", __func__);

    if (str == NULL)
	return;

    options_file_quiet_flag = TRUE;
    parse_options(&argc, &argv, &argc_allocated, str);
    phgOptionsParseCmdline(&argc, &argv);

    for (i = 0; i < argc; i++)
	phgFree(argv[i]);
    phgFree(argv);
}

BOOLEAN
phgOptionsSetNoArg(const char *op_name, BOOLEAN value)
{
    return set_option(op_name, &value, VT_NONE, __func__);
}

BOOLEAN
phgOptionsSetInt(const char *op_name, INT value)
{
    return set_option(op_name, &value, VT_INT, __func__);
}

BOOLEAN
phgOptionsSetFloat(const char *op_name, FLOAT value)
{
    return set_option(op_name, &value, VT_FLOAT, __func__);
}

BOOLEAN
phgOptionsSetKeyword(const char *op_name, const char *value)
{
    return set_option(op_name, (void *)value, VT_KEYWORD, __func__);
}

BOOLEAN
phgOptionsSetString(const char *op_name, const char *value)
{
    return set_option(op_name, (void *)value, VT_STRING, __func__);
}

BOOLEAN
phgOptionsSetFilename(const char *op_name, const char *value)
{
    return set_option(op_name, (void *)value, VT_FILENAME, __func__);
}

BOOLEAN
phgOptionsSetHandler(const char *op_name, const char *value)
{
    return set_option(op_name, (void *)value, VT_HANDLER, __func__);
}

static void
get_option0(OPTION *o, void **pvar)
{
    *pvar = NULL;
    switch (o->type) {
	case VT_NONE:
	case VT_INT:
	case VT_FLOAT:
	    *pvar = o->var;
	    break;
	case VT_STRING:
	case VT_FILENAME:
	    *pvar = *(char **)o->var;
	    break;
	case VT_KEYWORD:
	    *pvar = (*(int *)o->var < 0 ? "none" : o->keys[*(int *)o->var]);
	    break;
	case VT_HANDLER:
	    if (o->var == NULL)
		*pvar = (o->keys == NULL || o->keys[0] == NULL ?
				"none" : o->keys[0]);
	    else
		*pvar = (char *)((OPTION_HANDLER)o->var)(o, "", 0);
	    break;
	default:
	    /* not allowed */
	    phgError(1, "%s:%d: unsupported or unimplemented option type.\n",
			__FILE__, __LINE__);
    }
}

static BOOLEAN
get_option(const char *op_name, void **pvar, int type, const char *func)
{
    int *k;
    OPTION *o;

    if (!phgInitialized)
	phgError(1, "%s must be called after phgInit!\n", func);

    if (op_name[0] == '-' || op_name[1] == '+')
	op_name++;

    key->name = (void *)op_name;
    k = bsearch(options_index + options_count, options_index,
		options_count, sizeof(*options_index), option_comp);
    key->name = NULL;	/* reset key->name */
    if (k == NULL)
	phgError(1, "%s: unknown option \"-%s\"!\n", func, op_name);
    o = options + (*k);
    if (type >= 0 && o->type != type) {
	phgInfo(-1, "%s: wrong function type for \"-%s\".", func, op_name);
	switch (o->type) {
	    case VT_NONE:
		phgError(1, "Please use phgOptionsGetNoArg instead.\n");
	    case VT_INT:
		phgError(1, "Please use phgOptionsGetInt instead.\n");
	    case VT_FLOAT:
		phgError(1, "Please use phgOptionsGetFloat instead.\n");
	    case VT_STRING:
		phgError(1, "Please use phgOptionsGetString instead.\n");
	    case VT_FILENAME:
		phgError(1, "Please use phgOptionsGetFilename instead.\n");
	    case VT_KEYWORD:
		phgError(1, "Please use phgOptionsGetKeyword instead.\n");
	    case VT_HANDLER:
		phgError(1, "Please use phgOptionsGetHandler instead.\n");
	    default:
		phgError(1, "No phgOptionsGetXXXX function for \"-%s\".\n",
			op_name);
	}
    }

    get_option0(o, pvar);

    return *pvar == NULL ? FALSE : TRUE;
}

BOOLEAN
phgOptionsGetNoArg(const char *op_name)
{
    void *value;

    get_option(op_name, &value, VT_NONE, __func__);

    return *(BOOLEAN *)value;
}

INT
phgOptionsGetInt(const char *op_name)
{
    void *value;

    get_option(op_name, &value, VT_INT, __func__);

    return *(INT *)value;
}

FLOAT
phgOptionsGetFloat(const char *op_name)
{
    void *value;

    get_option(op_name, &value, VT_FLOAT, __func__);

    return *(FLOAT *)value;
}

const char *
phgOptionsGetKeyword(const char *op_name)
{
    void *value;

    get_option(op_name, &value, VT_KEYWORD, __func__);

    return value;
}

const char *
phgOptionsGetString(const char *op_name)
{
    void *value;

    get_option(op_name, &value, VT_STRING, __func__);

    return value;
}

const char *
phgOptionsGetFilename(const char *op_name)
{
    void *value;

    get_option(op_name, &value, VT_FILENAME, __func__);

    return value;
}

const char *
phgOptionsGetHandler(const char *op_name)
{
    void *value;

    get_option(op_name, &value, VT_HANDLER, __func__);

    return value;
}

#if 0	/*===================== test code =====================*/

static FLOAT v1 = 0;

static BOOLEAN v2 = FALSE;

static int v3 = 0;
static const char *keys[] = {"P1", "P2", "P3", "P4", NULL};

static char *v4 = "string argument";

int main(int argc, char *argv[])
{
    int i;

    phgOptionsRegisterFloat("option1", "option with a real argument", &v1);
    phgOptionsRegisterNoArg("option2", "option without argument", &v2);
    phgOptionsRegisterKeyword("option3", "option with a keyword argument",
		keys, &v3);
    phgOptionsRegisterString("option4", "option with a string argument", &v4);

#if 0
    phgOptionsParseCmdline(&argc, &argv);
    phgOptionsHelp();
#else
    phgInit(&argc, &argv);
#endif

    phgVerbosity = 3;

    phgInfo(-1, "Options remained:\n");
    for (i = 0; i < argc; i++)
	phgInfo(-1, "argv[%d]=%s\n", i, argv[i]);

    phgOptionsPush();

    phgOptionsSetFloat("option1", (FLOAT)25.);
    phgOptionsSetNoArg("option2", TRUE);
    phgOptionsSetKeyword("option3", "P4");
    phgOptionsSetString("option4", "xyz");

    phgInfo(-1, "option1: %le\n", (double)phgOptionsGetFloat("option1"));
    phgInfo(-1, "option2: %s\n", phgOptionsGetNoArg("option2") ? "TRUE":"FALSE");
    phgInfo(-1, "option3: %s\n", phgOptionsGetKeyword("option3"));
    phgInfo(-1, "option4: %s\n", phgOptionsGetString("option4"));

    phgOptionsPop();

    phgInfo(-1, "option1: %le\n", (double)phgOptionsGetFloat("option1"));
    phgInfo(-1, "option2: %s\n", phgOptionsGetNoArg("option2") ? "TRUE":"FALSE");
    phgInfo(-1, "option3: %s\n", phgOptionsGetKeyword("option3"));
    phgInfo(-1, "option4: %s\n", phgOptionsGetString("option4"));

    phgOptionsReset();
    phgFinalize();

    return 0;
}
#endif	/*===================== test code =====================*/
