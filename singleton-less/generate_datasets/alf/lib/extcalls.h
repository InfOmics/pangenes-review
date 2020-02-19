#ifndef EXTCALLSH
#define EXTCALLSH
#define EXTPARAM_LIST 0
#define EXTPARAM_CHAR 1
#define EXTPARAM_INT 2
#define EXTPARAM_DBL 3
#define EXTID_INVAL (-2)

typedef struct extref {
  int type;
  int len;
  union { char c[1]; int n[1]; double x[1]; 
	  struct extref *p[1]; } data;
} extref;
  
typedef struct extblock {
  int inid;
  int size;
  int end;
  int nrparam;  /* -1: p[0]->data.c contains error message,
                   -2: p[1]->data.c contains additional message */
  extref *p[1];
} extblock;

#ifdef _DEBUG
/* Macros to access parameters with range checks 
   if compiled with -DDARWIN_DEBUG */
int extindexerr (char *file, int line, char *name, int i, int len);
#define EXTC(r,i) ((r)->data.c[(i)>=0 && (i)<((r)->len)? (i) :\
		extindexerr(__FILE__,__LINE__, "r",i,(r)->len)])
#define EXTN(r,i) ((r)->data.n[(i)>=0 && (i)<((r)->len)? (i) :\
		extindexerr(__FILE__,__LINE__, "r",i,(r)->len)])
#define EXTX(r,i) ((r)->data.x[(i)>=0 && (i)<((r)->len)? (i) :\
		extindexerr(__FILE__,__LINE__, "r",i,(r)->len)])
#define EXTP(r,i) ((r)->data.p[(i)>=0 && (i)<((r)->len)? (i) :\
		extindexerr(__FILE__,__LINE__, "r",i,(r)->len)])
#else
/* Macros to access parameters if compiled without -DDARWIN_DEBUG */
#define EXTC(r,i) ((r)->data.c[(i)])
#define EXTN(r,i) ((r)->data.n[(i)])
#define EXTX(r,i) ((r)->data.x[(i)])
#define EXTP(r,i) ((r)->data.p[(i)])
#endif

int
  extenter (char *arg1, int maxoutsize, int maxoutparams);
  /* To be called at the beginning of an externally called program.
     Returns the number of parameters passed from Darwin.
     arg1: first argument passed to the main program (argv[1])
     maxoutsize: maximum size of all output parameters (in bytes)
     maxoutparams: maximum number of return parameters. */

void
  extexit (),
  /* To be called at the end of an externally called program.
     [Normal termination]. */
  extuserror (char *msg),
  /* To be called when a user error occurs (results in error call 
     in Darwin) [Error termination]. */
  extterminate (char *msg);
  /* To be called when a fatal error occurs (results in msg being
     written to standard error and Darwin not reading any data from
     shared memory). [Fatal error termination]. */

  /* The following functions provide access to the parameters
     passed from Darwin.

     For a Darwin call like: CallExternal( C_program, arg1, arg2, .. )
     arg1, arg2, .. can be accessed by the appropriate extget function.

     The arguments can be integers, numbers, strings or lists of
     any of the previous, including lists or lists.

     In general
     	handle1 = extget...( 1, ... )
     will create a handle that accesses the 1st argument.  This
     handle provides access to:
     	handle1->len	- the length of the array/list
			  the length of a string
			  1 in the case of an integer or double
	handle1->type	- the type
				EXTPARAM_LIST	list
				EXTPARAM_CHAR	string
				EXTPARAM_INT	{integer,array(integer)}
				EXTPARAM_DBL	{numeric,array(numeric)}
	direct usage		macro		definition
	handle1->data.c[i]  EXTC(handle1,i)	ith character
	handle1->data.n[i]  EXTN(handle1,i)	ith integer
	handle1->data.x[i]  EXTX(handle1,i)	ith double
	handle1->data.p[i]  EXTP(handle1,i)	ith extref (recursive structure)
	(they all go from i=0 to handle1->len-1, i.e. C-style)
									*/

extref
  *extget (int nr),
  /* Returns reference to parameter nr (first: nr = 1). Calls
     extuserror if no such parameter exists. */
  *extgetlist (int nr, int dim),
  /* Returns reference to parameter nr and verifies it being a
     dim-dimensional list.  Scalars are 0-dimensional, lists are
     1-dimensional, matrices are 2-dimensional, etc.  Calls extuserror
     if the parameter does not exist or has a different type. */
  *extgettext (int nr, int dim),
  /* Returns reference to parameter nr and verifies it being a
     dim-dimensional array of characters. A 1-dimensional array is
     a string in Darwin. Calls extuserror if the parameter does not
     exist or has a different type. */
  *extgetint (int nr, int dim),
  /* Returns reference to parameter nr and verifies it being a
     dim-dimensional array of integers. Calls extuserror if the
     parameter does not exist or has a different type. */
  *extgetdbl (int nr, int dim);
  /* Returns reference to parameter nr and verifies it being a
     dim-dimensional array of doubles. Calls extuserror if the
     parameter does not exist or has a different type. */

void
  extaddretval (extref *p);
  /* Appends the data at p to the list of return values. Normally
     called as
     extaddretval (extalloc...)   or  extaddretval (extput...). */

extref
  *extalloclist (int len, extref ***dest),
  /* Allocates a list of length len, returns the reference to it,
     and assigns the pointer to the list (i.e. to e->data.p) 
     to *dest. */
  *extalloctext (int len, char **dest),
  /* Allocates a text of length len, returns the reference to it,
     and assigns the pointer to the text (i.e. to e->data.c) 
     to *dest. */
  *extallocint (int len, int **dest),
  /* Allocates an array of integers of length len, returns the
     reference to it, and assigns the pointer to the array (i.e. 
     to e->data.n) to *dest. */
  *extallocdbl (int len, double **dest),
  /* Allocates an array of doubles of length len, returns the
     reference to it, and assigns the pointer to the array (i.e. to
     e->data.x) to *dest. */
  *extputtext (int len, char *c),
  /* Allocates a text of length len, returns the reference to it,
     and fills it with the text at c. */
  *extputint (int len, int *n),
  /* Allocates an array of integers of length len, returns the
     reference to it, and fills it with the integers at n. */
  *extputdbl (int len, double *x);
  /* Allocates an array of doubles of length len, returns the
     reference to it, and fills it with the doubles at x. */
#endif
