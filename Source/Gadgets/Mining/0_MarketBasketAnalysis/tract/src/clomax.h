/*----------------------------------------------------------------------
  File    : clomax.h
  Contents: prefix tree management for closed and maximal item sets
  Author  : Christian Borgelt
  History : 2009.10.08 file created as pfxtree.c
            2009.10.26 function cmt_dir() added (get item order dir.)
            2010.02.10 function cmt_proj() added (project a c/m tree)
            2010.03.12 function cmt_xproj() added (project a c/m tree)
            2010.06.21 generalized to support type RSUPP (int/double)
            2010.07.09 function cmt_check() added (with given support)
            2010.07.22 closed/maximal item set filter functions added
            2011.05.10 bug for RSUPP=double fixed (CMTREE.max)
----------------------------------------------------------------------*/
/* This version uses a top-down structure for the repository trees  */
/* and their processing. A frequent pattern tree structure was also */
/* tried, but turned out to be slower while needing more memory.    */

#ifndef __CLOMAX__
#define __CLOMAX__
#include "memsys.h"
#include "tract.h"

#ifdef _MSC_VER
#define INFINITY    (DBL_MAX+DBL_MAX)
#endif                          /* MSC still does not support C99 */

/*--------------------------------------------------------------------*/

#ifndef RSUPP
#define RSUPP       SUPP        /* support type for repository */
#define RSUPP_MAX   SUPP_MAX    /* maximum support value */
#define RSUPP_FMT   SUPP_FMT    /* printf format code for SUPP_T */

#else
#define int         1           /* to check definition of RSUPP_T */
#define long        2           /* for certain types */
#define double      3

#if   RSUPP==double
#ifndef RSUPP_MAX
#define RSUPP_MAX   INFINITY    /* maximum support value */
#endif
#ifndef RSUPP_FMT
#define RSUPP_FMT   "g"         /* printf format code for double */
#endif

#elif RSUPP==int
#ifndef RSUPP_MAX
#define RSUPP_MAX   INT_MAX     /* maximum support value */
#endif
#ifndef RSUPP_FMT
#define RSUPP_FMT   "d"         /* printf format code for int */
#endif

#elif RSUPP==long
#ifndef RSUPP_MAX
#define RSUPP_MAX   LONG_MAX    /* maximum support value */
#endif
#ifndef RSUPP_FMT
#define RSUPP_FMT   "ld"        /* printf format code for long */
#endif

#else                           /* assuming ptrdiff_t */
#ifndef RSUPP_MAX
#define RSUPP_MAX   PTRDIFF_MAX /* maximum support value */
#endif
#ifndef RSUPP_FMT
  #ifdef _MSC_VER
  #define RSUPP_FMT "Id"        /* printf format code for ptrdiff_t */
  #else
  #define RSUPP_FMT "zd"        /* printf format code for ptrdiff_t */
  #endif                        /* MSC still does not support C99 */
#endif
#endif

#undef int                      /* remove preprocessor definitions */
#undef long                     /* needed for the type checking */
#undef double
#endif

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct cmnode {         /* --- c/m prefix tree node --- */
  RSUPP         supp;           /* support of represented item set */
  ITEM          item;           /* associated item (last item in set) */
  struct cmnode *sibling;       /* successor node in sibling list */
  struct cmnode *children;      /* list of child nodes */
} CMNODE;                       /* (c/m prefix tree node) */

typedef struct {                /* --- c/m prefix tree --- */
  MEMSYS  *mem;                 /* memory management system */
  int     dir;                  /* direction of item order */
  ITEM    size;                 /* (maximum) number of items */
  ITEM    item;                 /* associated prefix item */
  RSUPP   max;                  /* maximum support for prefix */
  CMNODE  root;                 /* root node of the tree */
  int     keep[1];              /* flags for cmt_xproj() calls */
} CMTREE;                       /* (c/m prefix tree) */

typedef struct {                /* --- closed/maximal filter --- */
  int     dir;                  /* direction of item order */
  ITEM    size;                 /* maximum number of prefix trees */
  ITEM    cnt;                  /* current number of prefix trees */
  CMTREE  *trees[1];            /* conditional prefix trees */
} CLOMAX;                       /* (closed/maximal filter) */

/*----------------------------------------------------------------------
  Prefix Tree Functions
----------------------------------------------------------------------*/
extern CMTREE* cmt_create (MEMSYS *mem, int dir, ITEM size);
extern void    cmt_clear  (CMTREE *cmt);
extern void    cmt_delete (CMTREE *cmt, int delms);
extern MEMSYS* cmt_memsys (CMTREE *cmt);
extern ITEM    cmt_cnt    (CMTREE *cmt);
extern int     cmt_dir    (CMTREE *cmt);
extern RSUPP   cmt_supp   (CMTREE *cmt);
extern RSUPP   cmt_max    (CMTREE *cmt);
extern int     cmt_valid  (CMTREE *cmt);

extern int     cmt_add    (CMTREE *cmt, const ITEM *items, ITEM n,
                           RSUPP supp);
extern RSUPP   cmt_get    (CMTREE *cmt, const ITEM *items, ITEM n);
extern RSUPP   cmt_check  (CMTREE *cmt, const ITEM *items, ITEM n,
                           RSUPP supp);
extern void    cmt_prune  (CMTREE *cmt, ITEM item);
extern CMTREE* cmt_proj   (CMTREE *dst, CMTREE *src, ITEM item);
extern CMTREE* cmt_xproj  (CMTREE *dst, CMTREE *src, ITEM item,
                           const ITEM *keep, ITEM n);

#ifndef NDEBUG
extern void    cmt_show   (CMTREE *cmt, ITEMBASE *base, int ind);
#endif

/*----------------------------------------------------------------------
  Closed/Maximal Filter Functions
----------------------------------------------------------------------*/
extern CLOMAX* cm_create  (int dir, ITEM size);
extern void    cm_delete  (CLOMAX *cm);
extern ITEM    cm_cnt     (CLOMAX *cm);
extern int     cm_dir     (CLOMAX *cm);
extern RSUPP   cm_supp    (CLOMAX *cm);
extern CMTREE* cm_tree    (CLOMAX *cm, ITEM i);

extern int     cm_add     (CLOMAX *cm, ITEM item, RSUPP supp);
extern int     cm_addnc   (CLOMAX *cm, ITEM item, RSUPP supp);
extern void    cm_remove  (CLOMAX *cm, ITEM n);
extern RSUPP   cm_tail    (CLOMAX *cm, const ITEM *items, ITEM n);
extern int     cm_update  (CLOMAX *cm, const ITEM *items, ITEM n,
                           RSUPP supp);
#ifndef NDEBUG
extern void    cm_show    (CLOMAX *cm, ITEMBASE *base, int ind);
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define cmt_memsys(t)    ((t)->mem)
#define cmt_cnt(t)       ((t)->cnt)
#define cmt_dir(t)       ((t)->dir)
#define cmt_supp(t)      ((t)->root.supp)
#define cmt_max(t)       ((t)->max)
#define cmt_valid(t)     ((t)->max >= -1)

#define cm_cnt(f)        ((f)->cnt)
#define cm_dir(f)        ((f)->dir)
#define cm_tree(f,i)     ((f)->trees[i])

#endif
