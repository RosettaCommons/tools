/* Output from p2c, the Pascal-to-C translator */
/* From input file "maxd.pas" */


#include <p2c/p2c.h>


#define MAXD_G
#include "maxd.h"
/* p2c: maxd.pas, line 143: 
 * Warning: Expected IMPLEMENT, found 'FORTRAN' [227] */
/* p2c: maxd.pas, line 143: Warning: Expected END, found 'FORTRAN' [227] */


void _maxd_init()
{
  /* ----------------------------------------------------------------- */
  /* ----------------------------------------------------------------- */

  /* 9-Jan-89 build-range modification: chainstart */
  long curfragment, longest, curjump;

  static int _was_initialized = 0;
  if (_was_initialized++)
    return;
  if (start > 0) {
/* p2c: maxd.pas, line 231: Warning: Symbol 'START' is not defined [221] */
    for (i = start; i <= goal; i++) {
/* p2c: maxd.pas, line 232: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 232: Warning: Symbol 'START' is not defined [221] */
/* p2c: maxd.pas, line 232: Warning: Symbol 'GOAL' is not defined [221] */
      if (fragmentlist[i] != NULL) {
/* p2c: maxd.pas, line 232:
 * Warning: Symbol 'FRAGMENTLIST' is not defined [221] */
/* p2c: maxd.pas, line 232: Warning: Symbol 'I' is not defined [221] */
	curfragment = fragmentlist[i];
/* p2c: maxd.pas, line 234:
 * Warning: Symbol 'FRAGMENTLIST' is not defined [221] */
/* p2c: maxd.pas, line 234: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 234:
 * Warning: Symbol 'CURFRAGMENT' is not defined [221] */
	if (i == chainstart) {
/* p2c: maxd.pas, line 235: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 235:
 * Warning: Symbol 'CHAINSTART' is not defined [221] */
	  curfragment.longeur = chainstart;
/* p2c: maxd.pas, line 235: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 235:
 * Warning: No field called LONGEUR in that record [288] */
/* p2c: maxd.pas, line 235:
 * Warning: Symbol 'CHAINSTART' is not defined [221] */
	}
	if (i > chainstart) {
/* p2c: maxd.pas, line 236: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 236:
 * Warning: Symbol 'CHAINSTART' is not defined [221] */
	  while (curfragment != NULL) {
	    longest = 0;
/* p2c: maxd.pas, line 239:
 * Warning: Symbol 'LONGEST' is not defined [221] */
	    curjump = curfragment.firstjump;
/* p2c: maxd.pas, line 240: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 240:
 * Warning: No field called FIRSTJUMP in that record [288] */
/* p2c: maxd.pas, line 240:
 * Warning: Symbol 'CURJUMP' is not defined [221] */
	    while (curjump != NULL) {
	      if (curjump.jumpfromfragment.longeur > longest) {
/* p2c: maxd.pas, line 243: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 243:
 * Warning: No field called JUMPFROMFRAGMENT in that record [288] */
/* p2c: maxd.pas, line 243: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 243:
 * Warning: No field called LONGEUR in that record [288] */
		longest = curjump.jumpfromfragment.longeur;
/* p2c: maxd.pas, line 244: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 244:
 * Warning: No field called JUMPFROMFRAGMENT in that record [288] */
/* p2c: maxd.pas, line 244: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 244:
 * Warning: No field called LONGEUR in that record [288] */
	      }
	      curjump = curjump.nextjump;
/* p2c: maxd.pas, line 245: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 245:
 * Warning: No field called NEXTJUMP in that record [288] */
	    }
	    curfragment.longeur = longest + 1;
/* p2c: maxd.pas, line 247: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 247:
 * Warning: No field called LONGEUR in that record [288] */
	    curfragment = curfragment.nextfragment;
/* p2c: maxd.pas, line 248: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 248:
 * Warning: No field called NEXTFRAGMENT in that record [288] */
	  }
	}
      }
    }
  }
  if (start <= 0)
    return;
/* p2c: maxd.pas, line 251: Warning: Symbol 'START' is not defined [221] */
  for (i = start; i <= goal; i++) {
/* p2c: maxd.pas, line 251: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 251: Warning: Symbol 'START' is not defined [221] */
/* p2c: maxd.pas, line 251: Warning: Symbol 'GOAL' is not defined [221] */
    curfragment = fragmentlist[i];
/* p2c: maxd.pas, line 253:
 * Warning: Symbol 'FRAGMENTLIST' is not defined [221] */
/* p2c: maxd.pas, line 253: Warning: Symbol 'I' is not defined [221] */
    if (curfragment != NULL) {
      gap[i] = (curfragment.longeur != i);
/* p2c: maxd.pas, line 254: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 254: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 254:
 * Warning: No field called LONGEUR in that record [288] */
/* p2c: maxd.pas, line 254: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 254: Warning: Symbol 'GAP' is not defined [221] */
    } else {
      gap[i] = true;
/* p2c: maxd.pas, line 254: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 254: Warning: Symbol 'GAP' is not defined [221] */
    }
    while (curfragment != NULL && gap[i]) {
/* p2c: maxd.pas, line 255: Warning: Symbol 'GAP' is not defined [221] */
/* p2c: maxd.pas, line 255: Warning: Symbol 'I' is not defined [221] */
      if (curfragment.longeur != i) {
/* p2c: maxd.pas, line 256: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 256:
 * Warning: No field called LONGEUR in that record [288] */
/* p2c: maxd.pas, line 256: Warning: Symbol 'I' is not defined [221] */
	curfragment = curfragment.nextfragment;
/* p2c: maxd.pas, line 256: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 256:
 * Warning: No field called NEXTFRAGMENT in that record [288] */
	continue;
      }
      gap[i] = false;
/* p2c: maxd.pas, line 258: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 258: Warning: Symbol 'GAP' is not defined [221] */
      if (curfragment.prev != NULL) {
/* p2c: maxd.pas, line 259: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 259:
 * Warning: No field called PREV in that record [288] */
	curfragment.prev.nextfragment = curfragment.nextfragment;
/* p2c: maxd.pas, line 259: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 259:
 * Warning: No field called PREV in that record [288] */
/* p2c: maxd.pas, line 259: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 259:
 * Warning: No field called NEXTFRAGMENT in that record [288] */
/* p2c: maxd.pas, line 259: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 259:
 * Warning: No field called NEXTFRAGMENT in that record [288] */
      }
      if (curfragment.nextfragment != NULL) {
/* p2c: maxd.pas, line 260: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 260:
 * Warning: No field called NEXTFRAGMENT in that record [288] */
	curfragment.nextfragment.prev = curfragment.prev;
/* p2c: maxd.pas, line 260: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 260:
 * Warning: No field called NEXTFRAGMENT in that record [288] */
/* p2c: maxd.pas, line 260: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 260:
 * Warning: No field called PREV in that record [288] */
/* p2c: maxd.pas, line 260: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 260:
 * Warning: No field called PREV in that record [288] */
      }
      fragmentlist[i]->prev = curfragment;
/* p2c: maxd.pas, line 261: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 261:
 * Warning: Symbol 'FRAGMENTLIST' is not defined [221] */
      curfragment.nextfragment = fragmentlist[i];
/* p2c: maxd.pas, line 262: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 262:
 * Warning: No field called NEXTFRAGMENT in that record [288] */
/* p2c: maxd.pas, line 262:
 * Warning: Symbol 'FRAGMENTLIST' is not defined [221] */
/* p2c: maxd.pas, line 262: Warning: Symbol 'I' is not defined [221] */
      fragmentlist[i] = curfragment;
/* p2c: maxd.pas, line 263: Warning: Symbol 'I' is not defined [221] */
/* p2c: maxd.pas, line 263:
 * Warning: Symbol 'FRAGMENTLIST' is not defined [221] */
      curfragment.prev = NULL;
/* p2c: maxd.pas, line 264: Warning: bad pointer dereference [165] */
/* p2c: maxd.pas, line 264:
 * Warning: No field called PREV in that record [288] */
    }
  }
}
/* p2c: maxd.pas, line 270: 
 * Warning: Expected a '.', found PROCEDURE [227] */
/* p2c: maxd.pas, line 270: 
 * Warning: Junk at end of input file ignored [277] */
/* p2c: Note: Remember to call _maxd_init() in main program [215] */





/* End. */
