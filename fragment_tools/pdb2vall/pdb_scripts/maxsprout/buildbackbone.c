/* Output from p2c, the Pascal-to-C translator */
/* From input file "buildbackbone.pas" */


#include "p2c.h"
#include "fbind.h"

/* #define DEBUG */ /* Remove the comments if you want (some!) debug info */
/*
__________________________________________________________________________

This set of programs constructs full atomic coordinates of a protein
from a given C(alpha) trace and optimizes side chain geometry.

Copyright by Liisa Holm and Chris Sander, 1989-1991.

No redistribution, no program changes, no commercial use.

For details, see J.Mol.Biol. 218, 183-194 (1991).

___________________________________________________________________________
*/

/* array dimensioning constants */

#define maxlen          30
#define maxres          10000
#define scaler          100
#define maxprot         100

/* Command line argumane stuff. Added by Chris Dodge - EBI - Oct 96 */
#define MAXARGLEN       81

char DBName[MAXARGLEN];    /* Protein list database file name */
char FileName[MAXARGLEN];  /* Input filename */
int  Default;              /* Use default parameters */

#define DEF_NOTDEF  -1     /* Use default not defined on command line */
#define DEF_NO       0     /* Don't use the defualt values */
#define DEF_YES      1     /* Use the default values */

char DefaultDir[] = "/data/research/pdb";  /* Defualt input dir */
char DefaultOut[] = "./";                  /* Defualt output dir */


/* string type declarations */

typedef Char string[81];

typedef Char a1;

typedef Char a2[2];
typedef Char a3[3];
typedef Char a4[4];
typedef Char a5[5];
typedef Char a6[6];
typedef Char a10[10];

/* openfile types */

typedef enum {
  infile, outfile
} inout;

/* CA distances */

typedef long dist_array[maxlen + 1];
typedef dist_array distancematrix[maxres];

typedef struct protlist_record {
  string protname;
  long start, length;
} protlist_record;

typedef protlist_record protlist[maxprot];


typedef struct score_record {
  long begresidue, matchresidue, matchlength, matchprotlistindex;
  string matchprotname;
  double score, carms;
  struct score_record *prevscore, *nextscore;
} score_record;

typedef enum {
  fast, slow
} slidemode;
typedef boolean gaparray[maxres];
typedef long redundancyarray[maxres];

/* protein list type declarations */

typedef double xyz[3];

typedef struct main_record {
  xyz n, ca, c, o;
} main_record;


typedef struct atom_record {
  struct atom_record *nextatom, *prevatom;
  xyz coord;
  a4 atomname;
  long icode;   /* free to use */
} atom_record;

typedef struct backbone_record {
  a4 resname;
  a5 resno;
  main_record main;
  a10 pdbresid;
} backbone_record;

typedef struct residue_record {
  struct residue_record *nextresidue, *prevresidue;
  backbone_record backbone;
  atom_record *firstside;
} residue_record;

typedef struct protein_record {
  struct protein_record *prevprotein, *nextprotein;
  string protname;
  long length, extra;   /* user's fancy */
  residue_record *firstresidue;
} protein_record;

typedef struct dist_record {
  dist_array d;
  residue_record *respo;
} dist_record;

typedef dist_record distm[maxres];

/* -------------------------------------------------------------------------- */


/*u3b*/
typedef double t1[maxres * 4];
typedef xyz t2[maxres * 4];
typedef double rotationmatrix[3][3];
typedef double translationvector[3];

/* -------------------------------------------------------------------------- */


typedef backbone_record fragmentcoord[maxlen + 1];


typedef struct fragment_record {
  double rmssum, maxdelta;
  long sourceprotno, sourceresno, longeur, njumpto;
  struct fragment_record *prevfragment, *prev, *nextfragment;
  backbone_record main0, main1, main2;
  struct jump_record *firstjump;
} fragment_record;

typedef struct jump_record {
  struct jump_record *nextjump, *prev;
  fragment_record *jumpfromfragment;
  double biggestdelta, rmsdelta;
} jump_record;

typedef fragment_record *fragmentlistarray[maxres];

/* lookup table */


#define maxl            8
#define maxdist         4000   /* distance unit is 0.01 A */


typedef struct hit_record {
  struct hit_record *nexthit;
  long protno, resno;
} hit_record;

typedef hit_record lookuptable[maxl][maxdist + 1];

/* ----------------------------------------------------------------- */
/* external procedures */
/*procedure u3b(w:t1;x,y:t2;n,mode: integer;var rms:real;
                var u:rotationmatrix;var t:translationvector;var ier:integer); FORTRAN;*/


extern Char *str_ PP((Char *Result, long i));

extern Void init_fragmentlist PP((fragment_record **fragmentlist, long len));

extern Void make_new_protein PP((protein_record **firstprotein));

extern Void make_new_residue PP((protein_record *protein,
				 residue_record **firstresidue,
				 residue_record **lastresidue));

extern Void backbone_to_list PP((protein_record *lastprotein,
				 residue_record **firstresidue,
				 Char *infilename, long *length));

extern Void fill_protein_record PP((protein_record *protein,
				    residue_record *first, Char *proteinname,
				    long l));

extern Void generate_cb PP((double *n, double *ca, double *c0, double *cb));

extern Void add_sidechains PP((protein_record *protein,
			       protein_record **firstprotein));

extern Void makeside PP((residue_record *residue, Char *atnm, double *cd));

extern Void mainrecord_to_brk PP((main_record main, backbone_record backbone,
				  FILE **tapeout));

extern Void oldmainrecord_to_brk PP((main_record main, Char *resname,
				     Char *resno, FILE **tapeout));

extern Void write_brk PP((Char *atomname, double *x, Char *pdbname,
			  Char *title, FILE **brkout));

extern Void oldwrite_brk PP((Char *atomname, double *x, Char *resname,
			     int chainid, Char *resno, Char *title,
			     FILE **brkout));

extern Void matchca PP((residue_record *residue1, residue_record *residue2,
			long length, double (*u)[3], double *t, double *rms));

extern Void transrot PP((double *r, double *rnew, double (*u)[3], double *t));

extern Void update_jump_fromoldtonew PP((fragment_record **fragmentlist,
    long k, long begresidue, long testlen, fragment_record **newfragment,
    double tolerance));

extern Void update_jump_fromnewtoold PP((fragment_record **fragmentlist,
    long k, long begresidue, long testlen, fragment_record **newfragment,
    double tolerance));

extern Void extend_to_right PP((long j, long begresidue, long testlen,
				long length, long cutoff1, long comparestart,
				long templatestart, dist_record *dist,
				long *nright, double *s));

extern Void extend_to_left PP((long j, long begresidue, long testlen,
			       long length, long cutoff1, long comparestart,
			       long templatestart, dist_record *dist,
			       long *nleft, double *s, long nright));

extern Void update_pathsums PP((long start, long goal,
				fragment_record **fragmentlist));

extern double distance PP((double *a, double *b));

extern Void load_distances PP((long *nprot, Char *master,
			       protlist_record *proteins, dist_record *dist,
			       protein_record **firstprotein,
			       protein_record **prot1));

extern Void fill_lookup PP((hit_record (*lookup)[maxdist + 1],
			    dist_record *dist, long nprot,
			    protlist_record *proteins));

extern Void openfile PP((FILE **filvar, Char *filename, inout mode));


#define maxmaxscores    50
/* Up'd this from 1000 --> 6000. CJD - 17.10.96 */
#define maxprotlength   6000


typedef score_record candarray[maxmaxscores];
typedef candarray scorearray[maxprotlength];
typedef long ilist[maxres];


#define superscale      1000


#define superscaler     1000


#define datadir         "csdata:"


/* Local variables for build: */
struct LOC_build {
  long l1, c1, l2, c2, maxscores, minlooplen;
  distm dist;
  protein_record *prot1;
  residue_record *firstresidue;
  gaparray gap, filled;
  long nscore[maxres], fragindex[maxres], redundancy[maxres];
  scorearray candidates;
  long chainstart, chainend, leftzone, rightzone;
  double tolerance;
  string master, statfilename, coorfilename, logfilename;
  FILE *statfile, *coorfile, *logfile;
  fragmentlistarray fragmentlist;
  protlist proteins;
  long testlen;
  lookuptable lookup;
  long looplen, cutoff1;
  boolean cb_only;
} ;

Local Void diff(x, y, z, LINK)
double *x, *y, *z;
struct LOC_build *LINK;
{
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  z[2] = x[2] - y[2];
}

Local Void cross(x, y, z, LINK)
double *x, *y, *z;
struct LOC_build *LINK;
{
  z[0] = x[1] * y[2] - y[1] * x[2];
  z[1] = x[2] * y[0] - y[2] * x[0];
  z[2] = x[0] * y[1] - y[0] * x[1];
}

Local Void generate_cb_(n, ca, c0, cb, LINK)
double *n, *ca, *c0, *cb;
struct LOC_build *LINK;
{
  long i;
  xyz a, b, c, d;
  double lc, ld, TEMP;

  /* a = ca to n; b = ca to c; c = right-vector; d = up-vector; cb = ca + 1.32 * right-vector + 0.765 * up-vector */
  diff(n, ca, a, LINK);
  diff(c0, ca, b, LINK);
  cross(a, b, c, LINK);   /* c = a x b */
  lc = 0.0;
  ld = 0.0;
  for (i = 0; i <= 2; i++) {
    TEMP = c[i];
    lc += TEMP * TEMP;
    d[i] = a[i] + b[i];
    TEMP = d[i];
    ld += TEMP * TEMP;
  }
  for (i = 0; i <= 2; i++) {
    c[i] /= sqrt(lc);
    d[i] = -(d[i] / sqrt(ld));
  }
  for (i = 0; i <= 2; i++)
    cb[i] = ca[i] + 1.32 * c[i] + 0.765 * d[i];
}


Local Void fill_protein_record_(protein, first, proteinname, l, LINK)
protein_record *protein;
residue_record *first;
Char *proteinname;
long l;
struct LOC_build *LINK;
{
  strcpy(protein->protname, proteinname);
  protein->length = l;
  protein->firstresidue = first;
}


Local Void make_new_residue_(protein, firstresidue, lastresidue, LINK)
protein_record *protein;
residue_record **firstresidue, **lastresidue;
struct LOC_build *LINK;
{
  residue_record *newresidue;

  newresidue = (residue_record *)Malloc(sizeof(residue_record));
  if (*firstresidue == NULL)
    *firstresidue = newresidue;
  if (*lastresidue != NULL)
    (*lastresidue)->nextresidue = newresidue;
  newresidue->prevresidue = *lastresidue;
  newresidue->nextresidue = NULL;
  newresidue->firstside = NULL;
  *lastresidue = newresidue;
  if (protein->firstresidue == NULL)
    protein->firstresidue = newresidue;
}

Local Void makeside_(residue, atnm, cd, LINK)
residue_record *residue;
Char *atnm;
double *cd;
struct LOC_build *LINK;
{
  atom_record *lastatom, *newatom;

  newatom = (atom_record *)Malloc(sizeof(atom_record));
  lastatom = residue->firstside;
  /* go to end of sidechain */
  if (residue->firstside == NULL)
    residue->firstside = newatom;
  else {
    while (lastatom->nextatom != NULL)
      lastatom = lastatom->nextatom;
  }
  if (lastatom != NULL)
    lastatom->nextatom = newatom;
  newatom->nextatom = NULL;
  memcpy(newatom->coord, cd, sizeof(xyz));
  memcpy(newatom->atomname, atnm, sizeof(a4));
}

Local Void backbone_to_list_(lastprotein, firstresidue, infilename, length,
			     LINK)
protein_record *lastprotein;
residue_record **firstresidue;
Char *infilename;
long *length;
struct LOC_build *LINK;
{
  FILE *tapein;
  residue_record *lastresidue;

  tapein = fopen(infilename, "rb");
  if (tapein != NULL)
    rewind(tapein);
  else {
    fprintf(stderr, "ERROR: Unable to open file: %s\n", infilename);
    fprintf(stderr, "ERROR: in function: backbone_to_list_()\n");
    fprintf(stderr, "Bye for now folks........\n\n");
    exit(1);
  }
  *firstresidue = NULL;
  lastresidue = NULL;
  *length = 0;
  while (!P_eof(tapein)) {
    make_new_residue_(lastprotein, firstresidue, &lastresidue, LINK);
    (*length)++;
    fread(&lastresidue->backbone, sizeof(backbone_record), 1, tapein);
  }
  if (tapein != NULL)
    fclose(tapein);
}

Local Void init_fragmentlist_(fragmentlist, len, LINK)
fragment_record **fragmentlist;
long len;
struct LOC_build *LINK;
{
  long i;

  for (i = 0; i < len; i++)
    fragmentlist[i] = NULL;
}

Local Void update_pathsums_(start, goal, fragmentlist, LINK)
long start, goal;
fragment_record **fragmentlist;
struct LOC_build *LINK;
{
  long i;
  double smallest, t;
  fragment_record *curfragment;
  jump_record *curjump, *remjump;

  for (i = start + 1; i <= goal; i++) {
    curfragment = fragmentlist[i - 1];
    while (curfragment != NULL) {
      smallest = 999999.0;
      curjump = curfragment->firstjump;
      while (curjump != NULL) {
	if (i == start + 1)
	  t = curjump->rmsdelta;
	else
	  t = curjump->rmsdelta + curjump->jumpfromfragment->rmssum;
	if (t <= smallest) {
	  smallest = t;
	  remjump = curjump;
	}
	curjump = curjump->nextjump;
      }
      curfragment->prevfragment = remjump->jumpfromfragment;
      curfragment->rmssum = smallest;
      curfragment->maxdelta = remjump->biggestdelta;
      curfragment = curfragment->nextfragment;
    }
  }
}

Local Void make_new_protein_(firstprotein, LINK)
protein_record **firstprotein;
struct LOC_build *LINK;
{
  protein_record *newprotein;

  newprotein = (protein_record *)Malloc(sizeof(protein_record));
  /* insert to 1st */
  newprotein->nextprotein = *firstprotein;
  if (*firstprotein != NULL)
    (*firstprotein)->prevprotein = newprotein;
  newprotein->prevprotein = NULL;
  *firstprotein = newprotein;
  newprotein->firstresidue = NULL;
}

Local Void update_gap(chainstart, start, goal, LINK)
long chainstart, start, goal;
struct LOC_build *LINK;
{
  /* 9-Jan-89 build-range modification: chainstart */
  fragment_record *curfragment;
  long i, longest;
  jump_record *curjump;

  if (start > 0) {
    for (i = start; i <= goal; i++) {
      if (LINK->fragmentlist[i - 1] != NULL) {
	curfragment = LINK->fragmentlist[i - 1];
	if (i == chainstart)
	  curfragment->longeur = chainstart;
	if (i > chainstart) {
	  while (curfragment != NULL) {
	    longest = 0;
	    curjump = curfragment->firstjump;
	    while (curjump != NULL) {
	      if (curjump->jumpfromfragment->longeur > longest)
		longest = curjump->jumpfromfragment->longeur;
	      curjump = curjump->nextjump;
	    }
	    curfragment->longeur = longest + 1;
	    curfragment = curfragment->nextfragment;
	  }
	}
      }
    }
  }
  if (start <= 0)
    return;
  for (i = start; i <= goal; i++) {
    curfragment = LINK->fragmentlist[i - 1];
    if (curfragment != NULL)
      LINK->gap[i - 1] = (curfragment->longeur != i);
    else
      LINK->gap[i - 1] = true;
    while (curfragment != NULL && LINK->gap[i - 1]) {
      if (curfragment->longeur != i) {
	curfragment = curfragment->nextfragment;
	continue;
      }
      LINK->gap[i - 1] = false;
      if (curfragment->prev != NULL)
	curfragment->prev->nextfragment = curfragment->nextfragment;
      if (curfragment->nextfragment != NULL)
	curfragment->nextfragment->prev = curfragment->prev;
      LINK->fragmentlist[i - 1]->prev = curfragment;
      curfragment->nextfragment = LINK->fragmentlist[i - 1];
      LINK->fragmentlist[i - 1] = curfragment;
      curfragment->prev = NULL;
    }
  }
}

Local Void transrot_(r, rnew, u, t, LINK)
double *r, *rnew;
double (*u)[3];
double *t;
struct LOC_build *LINK;
{
  long i, j;

  /* rnew = u*r+t */
  for (i = 0; i <= 2; i++) {
    rnew[i] = t[i];
    for (j = 0; j <= 2; j++)
      rnew[i] += u[j][i] * r[j];
  }
}

Local protein_record *findproteinpointer(master, firstprotein, LINK)
Char *master;
protein_record *firstprotein;
struct LOC_build *LINK;
{
  protein_record *prot;
  boolean found;

  prot = firstprotein;
  found = false;
  while (!found && prot != NULL) {
    if (strcmp(prot->protname, master))
      prot = prot->nextprotein;
    else
      found = true;
  }
  return prot;
}

Local residue_record *findresnoresname(protein, resno, resname, LINK)
protein_record *protein;
Char *resno;
Char *resname;
struct LOC_build *LINK;
{
  residue_record *curresidue;
  boolean found;

  curresidue = protein->firstresidue;
  if (curresidue == NULL)
    printf(" firstresidue = nil !!! \n");
  found = false;
  while (curresidue != NULL && !found) {
    if (strncmp(curresidue->backbone.resno, resno, sizeof(a5)) ||
	strncmp(curresidue->backbone.resname, resname, sizeof(a4)))
      curresidue = curresidue->nextresidue;
    else
      found = true;
  }
  return curresidue;
}

Local Char *str__(Result, v, LINK)
Char *Result;
long v;
struct LOC_build *LINK;
{
  long i;
  Char aa[11];
  a4 s;

  if (v > 9999) {
    printf("too large: %12ld\n", v);
    _Escape(0);
  }
  memcpy(aa, "0123456789 ", 11L);
  i = (long)(v / 1000.0);
  s[0] = aa[i];
  v -= i * 1000;
  i = (long)(v / 100.0);
  s[1] = aa[i];
  v -= i * 100;
  i = (long)(v / 10.0);
  s[2] = aa[i];
  v -= i * 10;
  s[3] = aa[v];
  i = 1;
  while (s[i - 1] == '0' && i < 4) {
    s[i - 1] = ' ';
    i++;
  }
  return memcpy(Result, s, sizeof(a4));
}

Local Char *codetoresname(Result, code, LINK)
Char *Result;
long code;
struct LOC_build *LINK;
{
  Char aa[150];
  Char s[5][30];
  long i, l, k;
  a4 X;

  memcpy(s[0], "ALAARGASNASPCYSGLUGLNGLYHISILE", 30L);
  memcpy(s[1], "LEULYSMETPHEPROSERTHRTRPTYRVAL", 30L);
  memcpy(s[2], "ASXGLXACDALBALIABUAROBASBETHSE", 30L);
  memcpy(s[3], "HYPHYLORNPCASARTAUTHYUNKACEFOR", 30L);
  memcpy(s[4], "CYHCSHCSSCYXILUPRZPROCPRTRYHOH", 30L);
  l = 0;
  for (k = 0; k <= 4; k++) {
    for (i = 0; i <= 29; i++) {
      l++;
      aa[l - 1] = s[k][i];
    }
  }
  for (i = 1; i <= 3; i++)
    X[i - 1] = aa[code + i - 2];
  X[3] = ' ';
  return memcpy(Result, X, sizeof(a4));
}


Local Char *codetoatomname(Result, code, LINK)
Char *Result;
long code;
struct LOC_build *LINK;
{
  Char aa[210];
  Char s[7][30];
  long i, l, k;
  a4 X;

  memcpy(s[0], "N  CA C  O  CB CG SG CG1CG2OG ", 30L);
  memcpy(s[1], "OG1OG2CD CD1CD2ND ND1ND2OD1OD2", 30L);
  memcpy(s[2], "SD CE NE NE1NE2OE1OE2CE1CE2CE3", 30L);
  memcpy(s[3], "CZ NZ CZ2CZ3CZ1OH OXTNH1NH2CH2", 30L);
  memcpy(s[4], "AE1AE2CH1DH D  HA HB HG HD DE ", 30L);
  memcpy(s[5], "HD1HD2H  DG1DG2HG1HG2DH1DH2DD1", 30L);
  memcpy(s[6], "DD2HE1HE2DE1DE2HE HZ DZ    ???", 30L);
  l = 0;
  for (k = 0; k <= 6; k++) {
    for (i = 0; i <= 29; i++) {
      l++;
      aa[l - 1] = s[k][i];
    }
  }
  for (i = 1; i <= 3; i++)
    X[i - 1] = aa[code + i - 2];
  X[3] = ' ';
  return memcpy(Result, X, sizeof(a4));
}


Local Void side_to_list(protein, infilename, LINK)
protein_record *protein;
Char *infilename;
struct LOC_build *LINK;
{
  FILE *tapein;
  residue_record *curresidue;
  long i, j, n, code, a, b, c, length;
  xyz cd;
  a4 atnm, resname;
  a5 resno;

  tapein = fopen(infilename, "rb");
  if (tapein != NULL)
    rewind(tapein);
  else {
    fprintf(stderr, "ERROR: Unable to open file: %s\n", infilename);
    fprintf(stderr, "ERROR: in function: side_to_list_()\n");
    fprintf(stderr, "Bye for now folks........\n\n");
    exit(1);
  }

  /* assume .MAIN has been read & protein-list initialized */
  fread(&length, sizeof(long), 1, tapein);
  for (n = 1; n <= length; n++) {   /* residue-counter */
    fread(&code, sizeof(long), 1, tapein);
    codetoresname(resname, code, LINK);
    fread(&code, sizeof(long), 1, tapein);
    str__(resno, code, LINK);
    curresidue = findresnoresname(protein, resno, resname, LINK);
    if (curresidue == NULL)
      _Escape(0);
    fread(&i, sizeof(long), 1, tapein);
    for (j = 1; j <= i; j++) {   /* atom-counter */
      fread(&code, sizeof(long), 1, tapein);
      fread(&a, sizeof(long), 1, tapein);
      fread(&b, sizeof(long), 1, tapein);
      fread(&c, sizeof(long), 1, tapein);
      codetoatomname(atnm, code, LINK);
      cd[0] = (double)a / superscaler;
      cd[1] = (double)b / superscaler;
      cd[2] = (double)c / superscaler;
      makeside_(curresidue, atnm, cd, LINK);
    }
  }
  if (tapein != NULL)
    fclose(tapein);
}


Local Void getprotein(master, sidechain, firstprotein, protein, LINK)
Char *master;
boolean sidechain;
protein_record **firstprotein, **protein;
struct LOC_build *LINK;
{
  /* return pointer to protein master - read coordinates if it is not in list */
  residue_record *firstresidue;
  long l;
  Char STR1[86];

  *protein = findproteinpointer(master, *firstprotein, LINK);
  if (*protein != NULL)
    return;
  make_new_protein_(firstprotein, LINK);
  sprintf(STR1, "%s.main", master);
  backbone_to_list_(*firstprotein, &firstresidue, STR1, &l, LINK);
  if (sidechain) {
    sprintf(STR1, "%s.side", master);
    side_to_list(*firstprotein, STR1, LINK);
  }
  fill_protein_record_(*firstprotein, firstresidue, master, l, LINK);
  *protein = *firstprotein;
}


Local Void generate_sidechains(ca_, cb_, c_, n_, residue, firstprotein, LINK)
double *ca_, *cb_, *c_, *n_;
residue_record *residue;
protein_record **firstprotein;
struct LOC_build *LINK;
{

  /* add sidechain atoms to residuepointer^.firstside^... */
  xyz ca, cb, c, n;
  protein_record *prot1;
  t1 w;
  t2 x, y;
  int mode, ilen, ier=0;
  double rms=0;
  rotationmatrix u;
  translationvector t;
  atom_record *curatom;
  xyz cbeta, rnew;
  Char STR1[12];
  main_record *WITH;

  memcpy(ca, ca_, sizeof(xyz));
  memcpy(cb, cb_, sizeof(xyz));
  memcpy(c, c_, sizeof(xyz));
  memcpy(n, n_, sizeof(xyz));
  if (!strncmp(residue->backbone.resname, "GLY ", sizeof(a4)))
    goto _Lquit;
  sprintf(STR1, "%s%.4s", datadir, residue->backbone.resname);
  getprotein(STR1, true, firstprotein, &prot1, LINK);
  WITH = &prot1->firstresidue->backbone.main;
  /* get transrot matrices for n-ca-cb-c */
  generate_cb_(WITH->n, WITH->ca, WITH->c, cbeta, LINK);
  memcpy(x[0], WITH->n, sizeof(xyz));
  memcpy(x[1], WITH->ca, sizeof(xyz));
  memcpy(x[2], cbeta, sizeof(xyz));
  memcpy(x[3], WITH->c, sizeof(xyz));
  memcpy(y[0], n, sizeof(xyz));
  memcpy(y[1], ca, sizeof(xyz));
  memcpy(y[2], cb, sizeof(xyz));
  memcpy(y[3], c, sizeof(xyz));
  w[0] = 1.0;
  w[1] = 10.0;
  w[2] = 10.0;
  w[3] = 1.0;

  mode = 1; ilen = 4; 
#ifdef DEBUG
  t[0] = 12.09876;
  printf("Calling u3b from generate_sidechains().\n");
  printf("C: %f %f %d %d %f %f %d\n", w[0],w[ilen-1],ilen,mode,rms,t[0],ier);
  printf("C: x=%f %f %f y=%f %f %f\n", x[0][0], x[0][1], x[0][2], y[0][0], y[0][1], y[0][2]);
#endif
  u3b(w, x, y, &ilen, &mode, &rms, u, t, &ier);
  /* copy transrotated sidechain atoms to residue^.firstside... */
  curatom = prot1->firstresidue->firstside;
  while (curatom != NULL) {
    transrot_(curatom->coord, rnew, u, t, LINK);
    makeside_(residue, curatom->atomname, rnew, LINK);
    curatom = curatom->nextatom;
  }
_Lquit: ;
}

Local Void add_sidechains_(protein, firstprotein, LINK)
protein_record *protein, **firstprotein;
struct LOC_build *LINK;
{
  residue_record *curresidue;
  xyz cb;
  backbone_record *WITH;

  curresidue = protein->firstresidue;
  while (curresidue != NULL) {
    WITH = &curresidue->backbone;
    if (curresidue->firstside == NULL)
      generate_cb_(WITH->main.n, WITH->main.ca, WITH->main.c, cb, LINK);
    else
      memcpy(cb, curresidue->firstside->coord, sizeof(xyz));
    generate_sidechains(WITH->main.ca, cb, WITH->main.c, WITH->main.n,
			curresidue, firstprotein, LINK);
    curresidue = curresidue->nextresidue;
  }
}

Local Void write_brk_(atomname, x, pdbname, segid, brkout, LINK)
Char *atomname;
double *x;
Char *pdbname;
Char *segid;
FILE **brkout;
struct LOC_build *LINK;
{
  fprintf(*brkout,
	  "ATOM         %.4s%.10s   %8.3f%8.3f%8.3f                  %.4s\n",
	  atomname, pdbname, x[0], x[1], x[2], segid);
}


Local Void mainrecord_to_brk_(main, backbone, tapeout, LINK)
main_record main;
backbone_record backbone;
FILE **tapeout;
struct LOC_build *LINK;
{
  write_brk_("N   ", main.n, backbone.pdbresid, "    ", tapeout, LINK);
  write_brk_("CA  ", main.ca, backbone.pdbresid, "    ", tapeout, LINK);
  write_brk_("C   ", main.c, backbone.pdbresid, "    ", tapeout, LINK);
  write_brk_("O   ", main.o, backbone.pdbresid, "    ", tapeout, LINK);
}



Local Void path(start, goal, LINK)
long start, goal;
struct LOC_build *LINK;
{
  long i;
  double smallest;
  fragment_record *goalfragment, *curfragment;
  fragmentlistarray trace;
  protein_record *prot1, *firstprotein;
  residue_record *curresidue;
  atom_record *curratom;
  xyz cb;
  Char STR1[86];
  backbone_record *WITH;

  firstprotein = NULL;
  init_fragmentlist_(trace, goal + 1, LINK);
  update_pathsums_(start, goal, LINK->fragmentlist, LINK);
  /* take smallest pathsum */
  goalfragment = LINK->fragmentlist[goal - 1];
  smallest = goalfragment->rmssum;
  while (goalfragment != NULL) {
    if (goalfragment->rmssum <= smallest) {
      curfragment = goalfragment;
      smallest = goalfragment->rmssum;
    }
    goalfragment = goalfragment->nextfragment;
  }
  /* copy temporarily to trace */
  trace[goal - 1] = curfragment;
  for (i = goal - 1; i >= start; i--)
    trace[i - 1] = trace[i]->prevfragment;
  /* make a copy of the C(alpha) trace */
  prot1 = NULL;
  make_new_protein_(&prot1, LINK);
  sprintf(STR1, "%s.main", LINK->master);
  backbone_to_list_(prot1, &LINK->firstresidue, STR1, &i, LINK);
  fill_protein_record_(prot1, LINK->firstresidue, LINK->master, i, LINK);
  curresidue = prot1->firstresidue;
  for (i = 1; i < start; i++)
    curresidue = curresidue->nextresidue;
  for (i = start - 1; i < goal; i++) {
    curresidue->backbone.main = trace[i]->main1.main;
    curresidue = curresidue->nextresidue;
  }
  if (goal < prot1->length)
    curresidue->backbone.main = trace[goal - 1]->main2.main;
  /* wind protein pointer's first residue to start of fragment */
  curresidue = prot1->firstresidue;
  for (i = 2; i <= start; i++)
    curresidue = curresidue->nextresidue;
  prot1->firstresidue = curresidue;
  /* add c(beta) */
  curresidue = prot1->firstresidue;
  for (i = start; i <= goal; i++) {
    WITH = &curresidue->backbone;
    if (strncmp(WITH->resname, "GLY ", sizeof(a4))) {
      if (LINK->cb_only) {   /* add cbetas or whole sidechain */
	generate_cb_(WITH->main.n, WITH->main.ca, WITH->main.c, cb, LINK);
	makeside_(curresidue, "CB  ", cb, LINK);
      } else
	add_sidechains_(prot1, &firstprotein, LINK);
    }
    curresidue = curresidue->nextresidue;
  }
  /* write only coordinates to first gap */
  curresidue = prot1->firstresidue;
      /* this points now to start of fragment ! */
  for (i = start + 1; i <= goal; i++)
    curresidue = curresidue->nextresidue;
  curresidue->nextresidue = NULL;
  curresidue = prot1->firstresidue;
  while (curresidue != NULL) {
    mainrecord_to_brk_(curresidue->backbone.main, curresidue->backbone,
		      &LINK->coorfile, LINK);
    curratom = curresidue->firstside;
    while (curratom != NULL) {
      WITH = &curresidue->backbone;
      write_brk_(curratom->atomname, curratom->coord,
		curresidue->backbone.pdbresid, "    ", &LINK->coorfile);
      curratom = curratom->nextatom;
    }
    curresidue = curresidue->nextresidue;
  }
  /* write fragment source to logfile */
  for (i = start - 1; i < goal; i++) {
    fprintf(LINK->logfile, "%5ld%5ld %11s%5ld%5ld%5ld%5.1f",
	    i + 1, trace[i]->sourceresno,
	    LINK->proteins[trace[i]->sourceprotno - 1].protname,
	    LINK->nscore[i], LINK->fragindex[i], LINK->redundancy[i],
	    trace[i]->rmssum);
    if (LINK->fragindex[i] > 0 && LINK->fragindex[i] <= LINK->maxscores)
      fprintf(LINK->logfile, "%8.1f%8.1f\n",
	      LINK->candidates[i][LINK->fragindex[i] - 1].score,
	      LINK->candidates[i][LINK->fragindex[i] - 1].carms);
    else
      putc('\n', LINK->logfile);
  }
}  /* path */

/* Local variables for score_to_fragmentlist: */
struct LOC_score_to_fragmentlist {
  struct LOC_build *LINK;
} ;

Local Void matchca_(residue1, residue2, length, u, t, rms, LINK)
residue_record *residue1, *residue2;
long length;
double (*u)[3];
double *t;
double *rms;
struct LOC_build *LINK;
{
  t2 x, y;   /* coordinates */
  long j;
  residue_record *res1, *res2;
  t1 w;   /* atom pair weight factors */
  int ilen, mode, ier=0;

  /* 1 moves, 2 is static */
  /* copy CA coordinates to onedimensional array */
  res1 = residue1;
  res2 = residue2;
  for (j = 0; j < length; j++) {
    memcpy(x[j], res1->backbone.main.ca, sizeof(xyz));
    memcpy(y[j], res2->backbone.main.ca, sizeof(xyz));
    res1 = res1->nextresidue;
    res2 = res2->nextresidue;
    w[j] = 1.0;
  }

  mode = 1; ilen = length;
#ifdef DEBUG
  ier = 2334, t[0] = 988.2;
  printf("Calling u3b from matchca_().\n");
  printf("C: %f %f %d %d %f %f %d\n", w[0],w[ilen-1],ilen,mode,*rms,t[0],ier);
  printf("C: x=%f %f %f y=%f %f %f\n", x[0][0], x[0][1], x[0][2], y[0][0], y[0][1], y[0][2]);
#endif
  u3b(w, x, y, &ilen, &mode, rms, u, t, &ier);
#ifdef DEBUG
  printf("Results after u3b.\n");
  printf("C: %f %f %d %d %f %f %d\n", w[0],w[ilen-1],ilen,mode,*rms,t[0],ier);
  printf("C: x=%f %f %f y=%f %f %f\n", x[0][0], x[0][1], x[0][2], y[0][0], y[0][1], y[0][2]);
#endif
}


Local Void get_transrot_fragment(fragment, beg1, beg2, length, p, LINK)
backbone_record *fragment;
long beg1, beg2, length, p;
struct LOC_score_to_fragmentlist *LINK;
{
  residue_record *res1, *res2, *curresidue;
  long i;
  rotationmatrix u;
  translationvector t;
  double rms;
  backbone_record *WITH;

  /* 1 moves, 2 is static */
  res1 = LINK->LINK->dist[LINK->LINK->proteins[p - 1].start + beg1 - 1].respo;
  curresidue = res1;   /* save for later use */
  res2 = LINK->LINK->dist[beg2 - 1].respo;
  /* match c(alphas) */
  matchca_(res1, res2, length, u, t, &rms);
  /* now transrotate, store result in fragment[] */
  for (i = 0; i < length; i++) {
    WITH = &curresidue->backbone;
    transrot_(WITH->main.n, fragment[i].main.n, u, t);
    transrot_(WITH->main.ca, fragment[i].main.ca, u, t);
    transrot_(WITH->main.c, fragment[i].main.c, u, t);
    transrot_(WITH->main.o, fragment[i].main.o, u, t);
    curresidue = curresidue->nextresidue;
    /* write all fragment coordinates to a file */
    /*with fragment[i].main do
                            begin
                                    oldwrite_brk('CA  ',ca,'UNK ',' ',str(i),
                                            substr(score.matchprotname,4,4),fragfile);
                                    oldwrite_brk('N   ',n,'UNK ',' ',str(i),
                                            substr(score.matchprotname,4,4),fragfile);
                                    oldwrite_brk('O   ',o,'UNK ',' ',str(i),
                                            substr(score.matchprotname,4,4),fragfile);
                                    oldwrite_brk('C   ',c,'UNK ',' ',str(i),
                                            substr(score.matchprotname,4,4),fragfile);
                            end;
    */
  }
}  /* get_transrot_fragment */

Local double distance_(a_, b_, LINK)
double *a_, *b_;
struct LOC_build *LINK;
{
  double Result;
  xyz a, b, c;
  double d;
  long i;

  memcpy(a, a_, sizeof(xyz));
  memcpy(b, b_, sizeof(xyz));
  diff(a, b, c, LINK);
  for (i = 0; i <= 2; i++) {
    if (c[i] > 1e10)
      goto _Lwarn;
  }
  for (i = 0; i <= 2; i++) {
    if (c[i] < -1e10)
      goto _Lwarn;
  }
  d = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
  Result = sqrt(d);
  goto _Lok;
_Lwarn:
  /*writeln('WARNING: distance out of range ',c[1], c[2], c[3]);
*/
  Result = 99999.999;
_Lok:
  return Result;
}


Local Void jump(jumpfrom, jumpto, rms, suurin, LINK)
fragment_record *jumpfrom, *jumpto;
double *rms, *suurin;
struct LOC_build *LINK;
{
  double delta[5];
  double sum;
  long i;
  double TEMP;

  delta[0] = distance_(jumpfrom->main2.main.ca, jumpto->main1.main.ca, LINK);
  delta[1] = distance_(jumpfrom->main2.main.n, jumpto->main1.main.n, LINK);
  delta[2] = distance_(jumpfrom->main1.main.c, jumpto->main0.main.c, LINK);
  delta[3] = distance_(jumpfrom->main1.main.ca, jumpto->main0.main.ca, LINK);
  delta[4] = distance_(jumpfrom->main1.main.o, jumpto->main0.main.o, LINK);
  /*delta[5] := 0.5 * distance(jumpfrom^.main1.main.o, jumpto^.main0.main.o);
*/
  *suurin = delta[0];
  for (i = 1; i <= 4; i++) {
    if (delta[i] > *suurin)
      *suurin = delta[i];
  }
  sum = 0.0;
  for (i = 0; i <= 4; i++) {
    TEMP = delta[i];
    sum += TEMP * TEMP;
  }
  *rms = sqrt(sum / 4);
}

Local Void update_jump_fromnewtoold_(fragmentlist, k, begresidue, testlen,
				     newfragment, tolerance, LINK)
fragment_record **fragmentlist;
long k, begresidue, testlen;
fragment_record **newfragment;
double tolerance;
struct LOC_build *LINK;
{
  fragment_record *jumpto;
  jump_record *newjump, *lastjump;
  double rms, suurin;

  if (k + begresidue > testlen)
    return;
  jumpto = fragmentlist[k + begresidue - 1];
  while (jumpto != NULL) {
    jump(*newfragment, jumpto, &rms, &suurin, LINK);
    if (suurin <= tolerance) {
      newjump = (jump_record *)Malloc(sizeof(jump_record));
      newjump->nextjump = NULL;
      newjump->jumpfromfragment = *newfragment;
      /* newjump^.*delta */
      newjump->rmsdelta = rms;
      newjump->biggestdelta = suurin;
      (*newfragment)->njumpto++;
      lastjump = jumpto->firstjump;
      newjump->nextjump = lastjump;
      if (lastjump != NULL)
	lastjump->prev = newjump;
      newjump->prev = NULL;
      jumpto->firstjump = newjump;
    }
    jumpto = jumpto->nextfragment;
  }
}

Local Void update_jump_fromoldtonew_(fragmentlist, k, begresidue, testlen,
				     newfragment, tolerance, LINK)
fragment_record **fragmentlist;
long k, begresidue, testlen;
fragment_record **newfragment;
double tolerance;
struct LOC_build *LINK;
{
  fragment_record *jumpfrom;
  jump_record *newjump, *lastjump;
  double rms, suurin;

  if (k + begresidue <= 2)
    return;
  jumpfrom = fragmentlist[k + begresidue - 3];
  while (jumpfrom != NULL) {
    jump(jumpfrom, *newfragment, &rms, &suurin, LINK);
    if (suurin <= tolerance) {
      newjump = (jump_record *)Malloc(sizeof(jump_record));
      newjump->nextjump = NULL;
      newjump->jumpfromfragment = jumpfrom;
      newjump->rmsdelta = rms;
      newjump->biggestdelta = suurin;
      lastjump = (*newfragment)->firstjump;
      newjump->nextjump = lastjump;
      if (lastjump != NULL)
	lastjump->prev = newjump;
      newjump->prev = NULL;
      jumpfrom->njumpto++;
      (*newfragment)->firstjump = newjump;
    }
    jumpfrom = jumpfrom->nextfragment;
  }
}



Local Void score_to_fragmentlist(score, LINK)
score_record *score;
struct LOC_build *LINK;
{
  struct LOC_score_to_fragmentlist V;
  long k;
  fragment_record *lastfragment, *newfragment;
  fragmentcoord tempfragment;
  boolean unique[maxlen + 1];   /* max matchlength on maxlen+1 ! */
  long FORLIM;

  V.LINK = LINK;
  get_transrot_fragment(tempfragment, score->matchresidue, score->begresidue,
			score->matchlength, score->matchprotlistindex, &V);
  FORLIM = score->matchlength;
  /* checkaa uniikit */
  for (k = 1; k <= FORLIM; k++) {
    unique[k - 1] = true;
    lastfragment = LINK->fragmentlist[k + score->begresidue - 2];
    while (lastfragment != NULL && unique[k - 1]) {
      if (distance_(lastfragment->main1.main.c, tempfragment[k - 1].main.c) <=
	  0.2) {
	if (distance_(lastfragment->main2.main.c, tempfragment[k].main.c) <= 0.2) {
	  if (distance_(lastfragment->main1.main.n, tempfragment[k - 1].main.n) <=
	      0.2) {
	    if (distance_(lastfragment->main2.main.n, tempfragment[k].main.n) <=
		0.2) {
	      if (distance_(lastfragment->main1.main.ca,
			   tempfragment[k - 1].main.ca) <= 0.2) {
		if (distance_(lastfragment->main2.main.ca,
			     tempfragment[k].main.ca) <= 0.2)
		  unique[k - 1] = false;
	      }
	    }
	  }
	}
      }
      lastfragment = lastfragment->nextfragment;
    }
  }
  FORLIM = score->matchlength;
  for (k = 1; k <= FORLIM; k++) {
    if (unique[k - 1])
    {  /* make new fragment, update links to other fragments */
      newfragment = (fragment_record *)Malloc(sizeof(fragment_record));
      newfragment->nextfragment = NULL;
      newfragment->prevfragment = NULL;
      newfragment->sourceprotno = score->matchprotlistindex;
      newfragment->sourceresno = score->matchresidue + k - 1;
      newfragment->rmssum = 0.0;
      newfragment->njumpto = 0;
      newfragment->maxdelta = 0.0;
      newfragment->main1 = tempfragment[k - 1];
      if (k < score->matchlength)
	newfragment->main2 = tempfragment[k];
      if (k > 1)
	newfragment->main0 = tempfragment[k - 2];
      newfragment->firstjump = NULL;
      newfragment->prev = NULL;
      update_jump_fromoldtonew_(LINK->fragmentlist, k, score->begresidue,
			       LINK->prot1->length, &newfragment,
			       LINK->tolerance);
      update_jump_fromnewtoold_(LINK->fragmentlist, k, score->begresidue,
			       LINK->prot1->length, &newfragment,
			       LINK->tolerance);
      newfragment->nextfragment = LINK->fragmentlist[k + score->begresidue - 2];
      if (LINK->fragmentlist[k + score->begresidue - 2] != NULL)
	LINK->fragmentlist[k + score->begresidue - 2]->prev = newfragment;
      newfragment->prev = NULL;
      LINK->fragmentlist[k + score->begresidue - 2] = newfragment;
      LINK->redundancy[k + score->begresidue - 2]++;
    }
  }
}  /* score_to_fragmentlist */


Local long firstgap(chainstart, testlen, LINK)
long chainstart, testlen;
struct LOC_build *LINK;
{
  long goal;

  goal = chainstart - 1;
  do {
    goal++;
  } while (!(LINK->gap[goal - 1] || goal > testlen));
  return goal;
}

Local Void change_parameters(LINK)
struct LOC_build *LINK;
{
  Char c, *TEMP;

  /* set range of residues to build */
  printf("build range Y/[N] ?  whole protein %12ld%12ld\n",
	 LINK->chainstart, LINK->chainend);
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c == 'y' || c == 'Y') {
    printf("enter new range \n");
    scanf("%ld%ld%*[^\n]", &LINK->chainstart, &LINK->chainend);
    getchar();
  }
  /* set width of gap repair zone */
  printf("change repair zone Y/[N] ?  currently %12ld%12ld\n",
	 LINK->leftzone, LINK->rightzone);
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c == 'y' || c == 'Y') {
    printf("enter new leftzone, rightzone\n");
    scanf("%ld%ld%*[^\n]", &LINK->leftzone, &LINK->rightzone);
    getchar();
  }
  /* set max no of fragments per residue to fetch */
  printf("change maxscores Y/[N] ?  currently %12ld\n", LINK->maxscores);
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c == 'y' || c == 'Y') {
    printf("enter new maxscores\n");
    scanf("%ld%*[^\n]", &LINK->maxscores);
    getchar();
    if (LINK->maxscores > maxmaxscores) {
      printf("out of range: setting to %12ld\n", (long)maxmaxscores);
      LINK->maxscores = maxmaxscores;
    }
  }
  /* change fragment search parameters */
  printf("change distance deviation cutoffs Y/[N] ? \n");
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c == 'y' || c == 'Y') {
    printf("QuickSearch looplen currently %12ld enter new ? \n", LINK->l1);
    scanf("%ld%*[^\n]", &LINK->l1);
    getchar();
    printf("QuickSearch cutoff1 currently %12ld enter new ? \n", LINK->c1);
    scanf("%ld%*[^\n]", &LINK->c1);
    getchar();
    printf("FullSearch looplen currently %12ld enter new ? \n", LINK->l2);
    scanf("%ld%*[^\n]", &LINK->l2);
    getchar();
    printf("FullSearch cutoff1 currently %12ld enter new ? \n", LINK->c2);
    scanf("%ld%*[^\n]", &LINK->c2);
    getchar();
  }
  /* relax looplen if fewer than maxscores fragments ? */
  printf("change minlooplen Y/[N] ?  currently %12ld\n", LINK->minlooplen);
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c == 'y' || c == 'Y') {
    printf("enter new minlooplen\n");
    scanf("%ld%*[^\n]", &LINK->minlooplen);
    getchar();
    if (LINK->maxscores < 3) {
      printf("out of range: setting to %12d\n", 3);
      LINK->minlooplen = 3;
    }
  }
  /* set tolerance for joining */
  printf("change tolerance Y/[N] ?  currently % .5E\n", LINK->tolerance);
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c == 'y' || c == 'Y') {
    printf("enter new tolerance\n");
    scanf("%lg%*[^\n]", &LINK->tolerance);
    getchar();
  }
  /* output file names */
  printf("coordinate output to %s OK ? [Y]/N \n", LINK->coorfilename);
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c == 'n' || c == 'N') {
    printf("enter new filename \n");
    fgets(LINK->coorfilename, 81, stdin);
    TEMP = strchr(LINK->coorfilename, '\n');
    if (TEMP != NULL)
      *TEMP = 0;
  }
  printf("log to %s OK ? [Y]/N \n", LINK->logfilename);
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c == 'n' || c == 'N') {
    printf("enter new filename \n");
    fgets(LINK->logfilename, 81, stdin);
    TEMP = strchr(LINK->logfilename, '\n');
    if (TEMP != NULL)
      *TEMP = 0;
  }
  printf("fragment statistics to %s OK ? [Y]/N \n", LINK->statfilename);
  scanf("%c%*[^\n]", &c);
  getchar();
  if (c == '\n')
    c = ' ';
  if (c != 'n' && c != 'N')
    return;
  printf("enter new filename \n");
  fgets(LINK->statfilename, 81, stdin);
  TEMP = strchr(LINK->statfilename, '\n');
  if (TEMP != NULL)
    *TEMP = 0;
}

/* Local variables for fetch: */
struct LOC_fetch {
  struct LOC_build *LINK;
  long begres, looplen, cutoff1;
  boolean quick;
  long *nscore;
  boolean enough;
  score_record candide[maxres];
  ilist di;
} ;

/* Local variables for scanfragments: */
struct LOC_scanfragments {
  struct LOC_fetch *LINK;
  long i, j;
} ;

Local Void extend_to_right_(j, begresidue, testlen, length, cutoff1,
			    comparestart, templatestart, dist, nright, s,
			    LINK)
long j, begresidue, testlen, length, cutoff1, comparestart, templatestart;
dist_record *dist;
long *nright;
double *s;
struct LOC_build *LINK;
{
  boolean diverge;
  long delta, k, l;
  double sc;

  diverge = false;
  k = 0;
  do {   /* eteen */
    k++;
    sc = 0.0;
    for (l = 0; l < k; l++) {
      delta = dist[j + l - 1].d[k - l] - dist[begresidue + l - 1].d[k - l];
      if (labs(delta) > cutoff1)
	diverge = true;
      sc += delta * delta;
    }
    if (!diverge)   /* tilaa nleftille! */
      *s += sc;
  } while (!(diverge || k == testlen - 1 || k == maxlen - 2 ||
	     j - comparestart + k == length ||
	     begresidue - templatestart + k == testlen));
  if (diverge)   /* begresidue lasketaan mukaan */
    *nright = k;
  else
    *nright = k + 1;
}

Local Void extend_to_left_(j, begresidue, testlen, length, cutoff1,
			   comparestart, templatestart, dist, nleft, s,
			   nright, LINK)
long j, begresidue, testlen, length, cutoff1, comparestart, templatestart;
dist_record *dist;
long *nleft;
double *s;
long nright;
struct LOC_build *LINK;
{
  boolean diverge;
  long delta, k, l, m;
  double sc;

  /* taakse */
  diverge = false;
  k = 0;
  l = nright;
  do {
    k++;
    sc = 0.0;
    for (m = 1; m < l + k; m++) {
      if (m <= maxlen && j - comparestart - k > 0 &&
	  begresidue - templatestart - k > 0) {
	delta = dist[j - k - 1].d[m] - dist[begresidue - k - 1].d[m];
	if (labs(delta) > cutoff1)
	  diverge = true;
	sc += delta * delta;
      } else
	diverge = true;
    }
    if (!diverge)
      *s += sc;
  } while (!(j - comparestart - k == 1 ||
	     begresidue - templatestart - k == 1 || l + k - 1 == testlen ||
	     diverge || l + k == maxlen));
  if (diverge)
    *nleft = k - 1;
  else
    *nleft = k;
}


Local Void test(LINK)
struct LOC_scanfragments *LINK;
{   /* test */
  long dchir, nright, nleft, ndist;
  double dssq;
  protlist_record *WITH;
  score_record *WITH1;

  WITH = &LINK->LINK->LINK->proteins[LINK->i - 1];
  dssq = 0.0;
  extend_to_right_(WITH->start + LINK->j,
    LINK->LINK->LINK->proteins[0].start + LINK->LINK->begres,
    LINK->LINK->LINK->testlen, WITH->length, LINK->LINK->cutoff1, WITH->start,
    LINK->LINK->LINK->proteins[0].start, LINK->LINK->LINK->dist, &nright,
    &dssq);
  if (nright < LINK->LINK->looplen)
    return;
  extend_to_left_(WITH->start + LINK->j,
		 LINK->LINK->LINK->proteins[0].start + LINK->LINK->begres,
		 LINK->LINK->LINK->testlen, WITH->length, LINK->LINK->cutoff1,
		 WITH->start, LINK->LINK->LINK->proteins[0].start,
		 LINK->LINK->LINK->dist, &nleft, &dssq, nright);
  dchir = labs(LINK->LINK->LINK->dist[LINK->LINK->begres - nleft - 1].d[0] -
	       LINK->LINK->LINK->dist[WITH->start + LINK->j - nleft - 1].d[0]);
  if (dchir >= 90 && dchir <= 270)  /* save score */
    return;
  (*LINK->LINK->nscore)++;   /* new score */
  WITH1 = &LINK->LINK->candide[*LINK->LINK->nscore - 1];
  strcpy(WITH1->matchprotname,
	 LINK->LINK->LINK->proteins[LINK->i - 1].protname);
  WITH1->matchprotlistindex = LINK->i;
  WITH1->matchlength = nright + nleft;
  WITH1->matchresidue = LINK->j - nleft;
  WITH1->begresidue = LINK->LINK->begres - nleft;
  ndist = WITH1->matchlength * (WITH1->matchlength + 1) / 2;
  WITH1->score = dssq / ndist;   /* square of d-rms */
  if (LINK->LINK->quick)   /* pick first only in quickfill */
    LINK->LINK->enough = true;
  /* writeln(begresidue:4,matchresidue:4,matchlength:4,score:8:1,matchprotname);
*/
}  /* test */

Local Void scanfragments(e, row, LINK)
long e, row;
struct LOC_fetch *LINK;
{
  struct LOC_scanfragments V;
  hit_record *curhit;

  V.LINK = LINK;
  if (e > maxdist)   /* array dimension check: skip chain break */
    goto _Lquit;
  curhit = LINK->LINK->lookup[row - 1][e].nexthit;
  while (curhit != NULL) {
    V.j = curhit->resno;
    V.i = curhit->protno;
    if (V.j > 1)
      test(&V);
    if (LINK->enough)
      goto _Lquit;
    curhit = curhit->nexthit;
  }
_Lquit: ;
}  /* scanfragments */

Local Void drmssort(i, j, d, LINK)
long i, j;
long *d;
struct LOC_fetch *LINK;
{
  long m, left, right, e;
  double x;

  left = i;
  right = j;
  m = (left + right) / 2;
  x = LINK->candide[d[m - 1] - 1].score;
  do {
    while (LINK->candide[d[i - 1] - 1].score < x)
      i++;
    while (LINK->candide[d[j - 1] - 1].score > x)
      j--;
    if (i <= j) {
      e = d[j - 1];
      d[j - 1] = d[i - 1];
      d[i - 1] = e;
      i++;
      j--;
    }
  } while (i <= j);
  if (left < j)
    drmssort(left, j, d, LINK);
  if (i < right)
    drmssort(i, right, d, LINK);
}  /* drmssort */

Local Void carmssort(i, j, d, LINK)
long i, j;
long *d;
struct LOC_fetch *LINK;
{
  long m, left, right, e;
  double x;

  left = i;
  right = j;
  m = (left + right) / 2;
  x = LINK->candide[d[m - 1] - 1].carms;
  do {
    while (LINK->candide[d[i - 1] - 1].carms < x)
      i++;
    while (LINK->candide[d[j - 1] - 1].carms > x)
      j--;
    if (i <= j) {
      e = d[j - 1];
      d[j - 1] = d[i - 1];
      d[i - 1] = e;
      i++;
      j--;
    }
  } while (i <= j);
  if (left < j)
    carmssort(left, j, d, LINK);
  if (i < right)
    carmssort(i, right, d, LINK);
}  /* carmssort */

Local Void calculate_carmses(LINK)
struct LOC_fetch *LINK;
{
  int i, j, ier, mode, ilen;
  residue_record *res1, *res2;
  t2 x, y;
  t1 w;
  rotationmatrix u;
  translationvector t;
  double cassq;
  long FORLIM;
  score_record *WITH;
  long FORLIM1;

  FORLIM = P_min(*LINK->nscore, LINK->LINK->maxscores);
  for (i = 0; i < FORLIM; i++) {
    WITH = &LINK->candide[LINK->di[i] - 1];
    res1 = LINK->LINK->dist[LINK->LINK->proteins[WITH->matchprotlistindex - 1].
			    start + WITH->matchresidue - 1].respo;
    res2 = LINK->LINK->dist[WITH->begresidue - 1].respo;
    FORLIM1 = WITH->matchlength;
    for (j = 0; j < FORLIM1; j++) {
      memcpy(x[j], res1->backbone.main.ca, sizeof(xyz));
      memcpy(y[j], res2->backbone.main.ca, sizeof(xyz));
      w[j] = 1.0;
      res1 = res1->nextresidue;
      res2 = res2->nextresidue;
    }

    mode = 0; ilen = WITH->matchlength;
#ifdef DEBUG
    printf("Calling u3b from calculate_carmses().\n");
    printf("C: %f %f %d %d %f %f %d\n", w[0],w[ilen-1],ilen,mode,cassq,t[0],ier);
    printf("C: x=%f %f %f y=%f %f %f\n", x[0][0], x[0][1], x[0][2], y[0][0], y[0][1], y[0][2]);
#endif
    u3b(w, x, y, &ilen, &mode, &cassq, u, t, &ier);
#ifdef DEBUG
    printf("Results after u3b.\n");
    printf("C: %f %f %d %d %f %f %d\n", w[0],w[ilen-1],ilen,mode,cassq,t[0],ier);
    printf("C: x=%f %f %f y=%f %f %f\n", x[0][0], x[0][1], x[0][2], y[0][0], y[0][1], y[0][2]);
#endif

	/* 0: get sum of squares only */
    WITH->carms = cassq / WITH->matchlength;   /* C(alpha) rms squared */
  }
}  /* calculate_carmses */

/* ------------------------------------------------------------------ */

Local Void fetch(begres_, looplen_, cutoff1_, quick_, dosort, nscore_, LINK)
long begres_, looplen_, cutoff1_;
boolean quick_, dosort;
long *nscore_;
struct LOC_build *LINK;
{
  struct LOC_fetch V;
  long i, d, e, FORLIM;

  V.LINK = LINK;
  V.begres = begres_;
  V.looplen = looplen_;
  V.cutoff1 = cutoff1_;
  V.quick = quick_;
  V.nscore = nscore_;
  *V.nscore = 0;
  V.enough = false;
  while (!V.enough && V.looplen >= LINK->minlooplen &&
	 *V.nscore < LINK->maxscores) {
    d = LINK->dist[V.begres - 1].d[V.looplen - 1];
    scanfragments(d, V.looplen - 1, &V);
    FORLIM = V.cutoff1;
    for (e = 1; e <= FORLIM; e++) {
      if (!V.enough) {
	scanfragments(d + e, V.looplen - 1, &V);
	if (!V.enough) {
	  if (d - e > 0)
	    scanfragments(d - e, V.looplen - 1, &V);
	}
      }
    }
    V.looplen--;
  }
  FORLIM = *V.nscore;
  for (i = 1; i <= FORLIM; i++)
    V.di[i - 1] = i;
  if (*V.nscore > 1) {
    if (dosort) {
      drmssort(1L, *V.nscore, V.di, &V);
      calculate_carmses(&V);
      carmssort(1L, (long)P_min(*V.nscore, LINK->maxscores), V.di, &V);
    }
  }
  FORLIM = P_min(*V.nscore, LINK->maxscores);
  for (i = 0; i < FORLIM; i++) {
    if (V.begres==maxprotlength) {
      fprintf(stderr, "\n\nERROR!! Is it the year 2046?\n");
      fprintf(stderr, "Oh no... this shouldn't happen!: It looks like the input\n");
      fprintf(stderr, "is too big. Don't dispair, it can probably be fixed quite\n" );
      fprintf(stderr, "easily by increasing the size of maxprotlength!\n");
      fprintf(stderr, "\nNow, where the heck is the source code?????\n");
      fprintf(stderr, "If Chris Dodge is still around, ask him.\n\n");
      exit(1);
    }
    LINK->candidates[V.begres - 1][i] = V.candide[V.di[i] - 1];
  }
  *V.nscore = P_min(*V.nscore, LINK->maxscores);
}  /* fetch */

/* ------------------------------------------------------------------ */

Local Void quickfill(LINK)
struct LOC_build *LINK;
{
  long i, FORLIM;

  printf("quickfill %5ld%5ld\n", LINK->chainstart, LINK->chainend);
  FORLIM = LINK->chainend - LINK->looplen;
  for (i = LINK->chainstart - 1; i <= FORLIM; i++) {
    fetch(i + 1, LINK->looplen, LINK->cutoff1, true, false, &LINK->nscore[i],
	  LINK);
	/* pick first fragment */
    if (LINK->nscore[i] > 0)
      score_to_fragmentlist(LINK->candidates[i], LINK);
  }
}

Local Void seal(chainstart, chainend, goal, LINK)
long chainstart, chainend, *goal;
struct LOC_build *LINK;
{
  long i, neof, l, l0, left, right;
  boolean allfetched;

  *goal = firstgap(chainstart, chainend, LINK);
  printf("seal %12ld%12ld%12ld\n", *goal, chainstart, chainend);
  while (*goal <= chainend) {
    allfetched = false;
    left = P_max(chainstart, *goal - LINK->leftzone);
    right = P_min(*goal + LINK->rightzone, chainend);
    do {   /* add fragments until gap moves forward or out of fragments */
      if (!allfetched) {
	for (i = left; i <= right; i++) {
	  if (!LINK->filled[i - 1])
	  {  /* fetch fragments for residues in gap zone */
	    fetch(i, LINK->looplen, LINK->cutoff1, false, true,
		  &LINK->nscore[i - 1], LINK);
	    LINK->filled[i - 1] = true;
	  }
	}
      }
      allfetched = true;
      /* add more fragments to fragmentlist */
      l0 = chainend;   /* lowest residue covered by new fragments */
      neof = 0;
      for (i = left - 1; i < right; i++) {
	LINK->fragindex[i]++;
	if (LINK->fragindex[i] <= LINK->nscore[i]) {
	  score_to_fragmentlist(&LINK->candidates[i]
				[LINK->fragindex[i] - 1], LINK);
	  l = LINK->candidates[i][LINK->fragindex[i] - 1].begresidue;
	  if (l < l0)
	    l0 = l;
	} else
	  neof++;
      }
      update_gap(chainstart, l0, *goal, LINK);
      if (LINK->gap[*goal - 1]) {
	if (neof > right - left)
	  goto _Lgiveup;
      }
      /* all fragments in gap zone have been used */
    } while (LINK->gap[*goal - 1]);   /* gap has been removed */
    /* how far did the gap move ? */
    printf("gap moved from %12ld", *goal);
    update_gap(chainstart, *goal - 1, chainend, LINK);
    *goal = firstgap(chainstart, chainend, LINK);
    printf(" to %12ld\n", *goal);
    /*if goal>8 then freespace(cs, goal-8, fragmentlist);
*/
  }
_Lgiveup: ;
}  /* seal */

Local Void write_headers(LINK)
struct LOC_build *LINK;
{
  fprintf(LINK->coorfile,
    "REMARK  Backbone construction from CA trace by the program MAX-SPROUT, \n");
  fprintf(LINK->coorfile, "REMARK  ver. 1 Oct 1990.  \n");
  fprintf(LINK->coorfile,
    "REMARK  Reference L. Holm and C. Sander, J.Mol.Biol. (1991) 218, 183-194. \n");
  fprintf(LINK->coorfile, "REMARK  %s residues %5ld to %5ld\n",
	  LINK->master, LINK->chainstart, LINK->chainend);
  fprintf(LINK->statfile,
    "*  Backbone construction from CA trace by the program MAX-SPROUT, \n");
  fprintf(LINK->statfile, "*  ver. 1 Oct 1990.  \n");
  fprintf(LINK->statfile,
    "*  Reference L. Holm and C. Sander, J.Mol.Biol. (1991) 218, 183-194. \n");
  fprintf(LINK->statfile, "*  %s residues %5ld to %5ld\n",
	  LINK->master, LINK->chainstart, LINK->chainend);
  fprintf(LINK->logfile,
    "*  Backbone construction from CA trace by the program MAX-SPROUT, \n");
  fprintf(LINK->logfile, "*  ver. 1 Oct 1990.  \n");
  fprintf(LINK->logfile,
    "*  Reference L. Holm and C. Sander, J.Mol.Biol. (1991) 218, 183-194. \n");
  fprintf(LINK->logfile, "*   %s residues %5ld to %5ld\n",
	  LINK->master, LINK->chainstart, LINK->chainend);
  fprintf(LINK->logfile, "*\n");
  fprintf(LINK->logfile, "*   residue number in reconstructed protein\n");
  fprintf(LINK->logfile,
	  "*   !    residue number in fragment origin protein\n");
  fprintf(LINK->logfile, "*   !    !    fragment origin protein\n");
  fprintf(LINK->logfile,
    "*   !    !    !          number of fragments fetched from the database\n");
  fprintf(LINK->logfile,
	  "*   !    !    !          !    number of test cycles\n");
  fprintf(LINK->logfile,
    "*   !    !    !          !    !    number of tripeptides to choose from\n");
  fprintf(LINK->logfile,
    "*   !    !    !          !    !    !   SSQ(distance deviations) col 5 fragment\n");
  fprintf(LINK->logfile,
    "*   !    !    !          !    !    !       !     C(alpha) rms of fragment\n");
  fprintf(LINK->logfile,
    "*   !    !    !          !    !    !       !       !  indicated in 5th column\n");
}

Local Void logstatistics(LINK)
struct LOC_build *LINK;
{
  /* append fragment statistics to logfile */
  long i, ngap, l, maxl_, resno, oldresno;
  Char oldprotname[11], protname[11];
  long hist[maxres + 1];
  Char c;
  boolean start;

  LINK->logfile = fopen(LINK->logfilename, "r");
  if (LINK->logfile != NULL)
    rewind(LINK->logfile);
  else {
    fprintf(stderr, "ERROR: Unable to open file: %s\n", LINK->logfilename);
    fprintf(stderr, "ERROR: in function: logstatistics()\n");
    fprintf(stderr, "Bye for now folks........\n\n");
    exit(1);
  }
  /* read three first lines with REMARK ! */
  for (i = 1; i <= 3; i++) {
    fscanf(LINK->logfile, "%*[^\n]");
    getc(LINK->logfile);
  }
  start = true;
  l = 0;
  maxl_ = 0;
  for (i = 1; i <= maxres; i++)
    hist[i] = 0;
  ngap = 0;
  oldresno = 0;
  memcpy(oldprotname, "!\"#$%&/()= ", 11L);
  while (!P_eof(LINK->logfile)) {
    c = getc(LINK->logfile);
    if (c == '\n')
      c = ' ';
    if (c == '*') {   /* skip header */
      fscanf(LINK->logfile, "%*[^\n]");
      getc(LINK->logfile);
      continue;
    }
    if (c == 'R') {
      ngap++;
      start = true;
      if (l > maxl_)
	maxl_ = l;
      if (l > 0)
	hist[l]++;
      l = 0;
      fscanf(LINK->logfile, "%*[^\n]");
      getc(LINK->logfile);
      memcpy(oldprotname, "!@#$%^&*()=", 11L);
      continue;
    }
    fscanf(LINK->logfile, "%ld%ld", &i, &resno);
    P_readlnpaoc(LINK->logfile, protname, 11);
    /* writeln(statfile, resno:5,oldresno:5,protname:20, oldprotname:20, l:5);
*/
    if (start || !strncmp(protname, oldprotname, 11) && resno == oldresno + 1) {
      l++;
      memcpy(oldprotname, protname, 11L);
      oldresno = resno;
      start = false;
      continue;
    }
    if (l > maxl_)
      maxl_ = l;
    hist[l]++;
    l = 1;
    oldresno = resno;
    memcpy(oldprotname, protname, 11L);
    start = false;
  }
  if (l > maxl_)
    maxl_ = l;
  hist[l]++;
  /* write histogram */
  fprintf(LINK->statfile,
	  "fragment statistics: longest: %5ld number of gaps: %5ld\n",
	  maxl_ + 1, ngap);
  fprintf(LINK->statfile, "length     count\n");
  for (i = 1; i <= maxl_; i++)
    fprintf(LINK->statfile, "%5ld%10ld\n", i + 1, hist[i]);

  /* read fragment data */
}  /*logstatistics */

Local Void x(LINK)
struct LOC_build *LINK;
{
}

Local Void init_residuelist(firstresidue, lastresidue, LINK)
residue_record **firstresidue, **lastresidue;
struct LOC_build *LINK;
{
  *firstresidue = NULL;
  *lastresidue = NULL;
}

Local long val(s, LINK)
Char *s;
struct LOC_build *LINK;
{
  long i, v;
  Char aa[11];

  memcpy(aa, "0123456789 ", 11L);
  v = 0;
  i = 0;
  do {
    i++;
  } while (aa[i - 1] != s[0] && i != 11);
  if (i < 11)
    v = (i - 1) * 1000;
  i = 0;
  do {
    i++;
  } while (aa[i - 1] != s[1] && i != 11);
  if (i < 11)
    v += (i - 1) * 100;
  i = 0;
  do {
    i++;
  } while (aa[i - 1] != s[2] && i != 11);
  if (i < 11)
    v += (i - 1) * 10;
  i = 0;
  do {
    i++;
  } while (aa[i - 1] != s[3] && i != 11);
  if (i < 11)
    v += i - 1;
  return v;
}

Local long resnametocode(resnm, LINK)
Char *resnm;
struct LOC_build *LINK;
{
  Char aa[150];
  Char s[5][30];
  long i, l, k;

  memcpy(s[0], "ALAARGASNASPCYSGLUGLNGLYHISILE", 30L);
  memcpy(s[1], "LEULYSMETPHEPROSERTHRTRPTYRVAL", 30L);
  memcpy(s[2], "ASXGLXACDALBALIABUAROBASBETHSE", 30L);
  memcpy(s[3], "HYPHYLORNPCASARTAUTHYUNKACEFOR", 30L);
  memcpy(s[4], "CYHCSHCSSCYXILUPRZPROCPRTRYHOH", 30L);
  l = 0;
  for (k = 0; k <= 4; k++) {
    for (i = 0; i <= 29; i++) {
      l++;
      aa[l - 1] = s[k][i];
    }
  }
  l = 0;
  k = 1;
  i = 1;
  while (k < 51 && l == 0) {
    if (aa[i - 1] == resnm[0]) {
      if (aa[i] == resnm[1]) {
	if (aa[i + 1] == resnm[2])
	  l = i;
      }
    }
    i += 3;
    k++;
  }
  return l;
}

Local long atomnametocode(atnm, LINK)
Char *atnm;
struct LOC_build *LINK;
{
  Char aa[210];
  Char s[7][30];
  long i, l, k;

  memcpy(s[0], "N  CA C  O  CB CG SG CG1CG2OG ", 30L);
  memcpy(s[1], "OG1OG2CD CD1CD2ND ND1ND2OD1OD2", 30L);
  memcpy(s[2], "SD CE NE NE1NE2OE1OE2CE1CE2CE3", 30L);
  memcpy(s[3], "CZ NZ CZ2CZ3CZ1OH OXTNH1NH2CH2", 30L);
  memcpy(s[4], "AE1AE2CH1DH D  HA HB HG HD DE ", 30L);
  memcpy(s[5], "HD1HD2H  DG1DG2HG1HG2DH1DH2DD1", 30L);
  memcpy(s[6], "DD2HE1HE2DE1DE2HE HZ DZ    ???", 30L);
  l = 0;
  for (k = 0; k <= 6; k++) {
    for (i = 0; i <= 29; i++) {
      l++;
      aa[l - 1] = s[k][i];
    }
  }
  l = 0;
  k = 1;
  i = 1;
  while (k < 71 && l == 0) {
    if (aa[i - 1] == atnm[0]) {
      if (aa[i] == atnm[1]) {
	if (aa[i + 1] == atnm[2])
	  l = i;
      }
    }
    i += 3;
    k++;
  }
  if (k == 71) {
    printf("unknown atomname *%.4s*\n", atnm);
    l = 68;
  }
  return l;
}


Local Void save_backbone(outfilename, protein, LINK)
Char *outfilename;
protein_record *protein;
struct LOC_build *LINK;
{
  FILE *tapeout;
  residue_record *curresidue;

  tapeout = fopen(outfilename, "w");
  if (tapeout != NULL)
    rewind(tapeout);
  else {
    fprintf(stderr, "ERROR: Unable to open file: %s\n", outfilename);
    fprintf(stderr, "ERROR: in function: save_backbone()\n");
    fprintf(stderr, "Bye for now folks........\n\n");
    exit(1);
  }
  curresidue = protein->firstresidue;
  while (curresidue != NULL) {
    fwrite(&curresidue->backbone, sizeof(backbone_record), 1, tapeout);
    curresidue = curresidue->nextresidue;
  }
  if (tapeout != NULL)
    fclose(tapeout);
  tapeout = NULL;
}



Local Void oldwrite_brk_(atomname, x, resname, chainid, resno, segid, brkout,
			 LINK)
Char *atomname;
double *x;
Char *resname;
Char chainid;
Char *resno;
Char *segid;
FILE **brkout;
struct LOC_build *LINK;
{
  fprintf(*brkout,
    "ATOM         %.4s%.4s%c%.5s   %8.3f%8.3f%8.3f                  %.4s\n",
    atomname, resname, chainid, resno, x[0], x[1], x[2], segid);
}

Local Void oldmainrecord_to_brk_(main, resname, resno, tapeout, LINK)
main_record main;
Char *resname;
Char *resno;
FILE **tapeout;
struct LOC_build *LINK;
{
  oldwrite_brk_("N   ", main.n, resname, ' ', resno, "    ", tapeout, LINK);
  oldwrite_brk_("CA  ", main.ca, resname, ' ', resno, "    ", tapeout, LINK);
  oldwrite_brk_("C   ", main.c, resname, ' ', resno, "    ", tapeout, LINK);
  oldwrite_brk_("O   ", main.o, resname, ' ', resno, "    ", tapeout, LINK);
}


Local double atan2_xx(y, x, LINK)
double y, x;
struct LOC_build *LINK;
{
  double z;

  if (x != 0.0)
    z = atan(y / x);
  else if (y > 0.0)
    z = 1.570796;
  else if (y < 0.0)
    z = -1.570796;
  else
    z = 6.283185;
  if (x >= 0.0)
    return z;
  if (y > 0.0)
    z += 3.141593;
  else
    z -= 3.141593;
  return z;
}

Local double dot(x, y, LINK)
double *x, *y;
struct LOC_build *LINK;
{
  return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}


Local double dihedralangle(v1, v2, v3, v4, LINK)
double *v1, *v2, *v3, *v4;
struct LOC_build *LINK;
{
  double Result, u, v;
  xyz v12, v43, x, y, z, p;

  diff(v1, v2, v12, LINK);
  diff(v4, v3, v43, LINK);
  diff(v2, v3, z, LINK);
  cross(z, v12, p, LINK);
  cross(z, v43, x, LINK);
  cross(z, x, y, LINK);
  u = dot(x, x, LINK);
  v = dot(y, y, LINK);
  Result = 360.0;
  if (u <= 0.0 || v <= 0.0)
    return Result;
  u = dot(p, x, LINK) / sqrt(u);
  v = dot(p, y, LINK) / sqrt(v);
  if (u != 0.0 || v != 0.0)
    return (atan2_xx(v, u, LINK) * 57.29578);
  return Result;
}

int openfile_(filvar, filename, mode, LINK)
FILE **filvar;
Char *filename;
inout mode;
struct LOC_build *LINK;
{
  long x;

  switch (mode) {

  case infile:
    /* If open is unsuccessful, then ignore & continue */
    *filvar = fopen(filename, "rb");
    if (*filvar != NULL) {
      rewind(*filvar);
      printf("reading file %s\n", filename);
    } else {
      printf("file %s not found. \n", filename);
      return -1;
    }
    break;

  case outfile:
    *filvar = fopen(filename, "wb");
    if (*filvar != NULL) {
      rewind(*filvar);
      printf("writing to %s\n", filename);
    } else {
      fprintf(stderr, "ERROR: Unable to open output file: %s\n", filename);
      fprintf(stderr, "ERROR: in function: openfile_()\n");
      fprintf(stderr, "Bye for now folks........\n\n");
      exit(1);
    }

    break;
  }
  return 0;
}


/* make list - procedures */


Local Void brk_to_list(lastprotein, firstresidue, infilename, headerfile,
		       length, LINK)
protein_record *lastprotein;
residue_record **firstresidue;
Char *infilename, *headerfile;
long *length;
struct LOC_build *LINK;
{
  /* NEW modification: header file, but missing from this - discovered 11-sep-89 */
  a6 header, skip6;
  a4 atnm, resnm;
  a5 resno, oldresno;
  xyz cd;
  a1 d, d1;
  FILE *tapein;
  residue_record *lastresidue;
  a10 pdbname;
  main_record *WITH;

  /* Brookhaven format: a6, i5, 1x, a4, a1, a3, 1x, a1, i4, a1, 3x, 3f8.3
   header=ATOM  ,atom serial number, atom name, alternate location indicator for
   atoms, residue name, chain identifier, residue sequence number, code for insertions
   of residue, x, y, z */

  /* update proteinlist */

  tapein = NULL;
  openfile_(&tapein, infilename, infile, LINK);
  init_residuelist(firstresidue, &lastresidue, LINK);
  *length = 0;

  /* read brk file */
  memcpy(oldresno, "*****", sizeof(a5));
  do {
    P_readpaoc(tapein, header, (int)(sizeof(a6)));
    if (strncmp(header, "ATOM  ", sizeof(a6))) {
      fscanf(tapein, "%*[^\n]");
      getc(tapein);
    }
  } while (strncmp(header, "ATOM  ", sizeof(a6)));
  while (!P_eof(tapein) && !strncmp(header, "ATOM  ", sizeof(a6))) {
    P_readpaoc(tapein, skip6, (int)(sizeof(a6)));
    d = getc(tapein);
    if (d == '\n')
      d = ' ';
    P_readpaoc(tapein, atnm, (int)(sizeof(a4)));
    P_readpaoc(tapein, resnm, (int)(sizeof(a4)));
    d1 = getc(tapein);
    if (d1 == '\n')
      d1 = ' ';
    P_readpaoc(tapein, resno, (int)(sizeof(a5)));
    fscanf(tapein, "%lg%lg%lg%*[^\n]", cd, &cd[1], &cd[2]);
    getc(tapein);
    if (strncmp(oldresno, resno, sizeof(a5))) {
      memcpy(oldresno, resno, sizeof(a5));
      /*
      if (findresnoresname(lastprotein,resno,resnm) = nil) then
      begin */
      /*if lastresidue<>nil then writeln('new residue: ',atnm,resnm,resno,d,'****',
lastresidue^.backbone.resname,lastresidue^.backbone.resno);
*/
      make_new_residue_(lastprotein, firstresidue, &lastresidue, LINK);
      memcpy(lastresidue->backbone.resno, resno, sizeof(a5));
      memcpy(lastresidue->backbone.resname, resnm, sizeof(a4));
      pdbname[0] = resnm[0];
      pdbname[1] = resnm[1];
      pdbname[2] = resnm[2];
      pdbname[3] = resnm[3];
      pdbname[4] = d1;
      pdbname[5] = resno[0];
      pdbname[6] = resno[1];
      pdbname[7] = resno[2];
      pdbname[8] = resno[3];
      pdbname[9] = resno[4];
      memcpy(lastresidue->backbone.pdbresid, pdbname, sizeof(a10));
      (*length)++;
    }

    WITH = &lastresidue->backbone.main;
    if (!strncmp(atnm, "CA  ", sizeof(a4)))
      memcpy(WITH->ca, cd, sizeof(xyz));
    else {
      if (!strncmp(atnm, "C   ", sizeof(a4)))
	memcpy(WITH->c, cd, sizeof(xyz));
      else {
	if (!strncmp(atnm, "N   ", sizeof(a4)))
	  memcpy(WITH->n, cd, sizeof(xyz));
	else {
	  if (!strncmp(atnm, "O   ", sizeof(a4)))
	    memcpy(WITH->o, cd, sizeof(xyz));
	  else
	    makeside_(lastresidue, atnm, cd, LINK);
	}
      }
    }
    if (!P_eof(tapein))
      P_readpaoc(tapein, header, (int)(sizeof(a6)));
  }
  if (tapein != NULL)
    fclose(tapein);
  tapein = NULL;
  if (tapein != NULL)
    fclose(tapein);
}

Local Void list_to_side(protein, outfilename, LINK)
protein_record *protein;
Char *outfilename;
struct LOC_build *LINK;
{
  FILE *tapeout;
  residue_record *curresidue;
  atom_record *currside;
  long atoms[100][4];
  long i, j, n;
  a4 resn;
  backbone_record *WITH;
  long TEMP;

  tapeout = fopen(outfilename, "rb");
  if (tapeout != NULL)
    rewind(tapeout);
  else {
    fprintf(stderr, "ERROR: Unable to open file: %s\n", outfilename);
    fprintf(stderr, "ERROR: in function: list_to_side()\n");
    fprintf(stderr, "Bye for now folks........\n\n");
    exit(1);
  }

  fwrite(&protein->length, sizeof(long), 1, tapeout);
  curresidue = protein->firstresidue;
  while (curresidue != NULL) {
    WITH = &curresidue->backbone;
    resn[0] = WITH->resno[0];
    resn[1] = WITH->resno[1];
    resn[2] = WITH->resno[2];
    resn[3] = WITH->resno[3];
    TEMP = resnametocode(WITH->resname, LINK);
    fwrite(&TEMP, sizeof(long), 1, tapeout);
    TEMP = val(resn, LINK);
    fwrite(&TEMP, sizeof(long), 1, tapeout);
    /* copy coordinates to atoms */
    n = 0;   /* atomcounter */
    currside = curresidue->firstside;
    while (currside != NULL) {
      n++;
      if (n > 100) {
	printf("n>100\n");
	_Escape(0);
      }
      for (i = 1; i <= 3; i++)
	atoms[n - 1]
	  [i] = (long)floor(superscale * currside->coord[i - 1] + 0.5);
      atoms[n - 1][0] = atomnametocode(currside->atomname, LINK);
      if (atoms[n - 1][0] == 68)   /* unknown atomname */
	n--;
      currside = currside->nextatom;
    }
    fwrite(&n, sizeof(long), 1, tapeout);
    for (i = 0; i < n; i++) {
      for (j = 0; j <= 3; j++)
	fwrite(&atoms[i][j], sizeof(long), 1, tapeout);
    }
    curresidue = curresidue->nextresidue;
  }
  if (tapeout != NULL)
    fclose(tapeout);
}

#undef superscale


#undef superscaler


Local Void fill_lookup_(lookup, dist, nprot, proteins, LINK)
hit_record (*lookup)[maxdist + 1];
dist_record *dist;
long nprot;
protlist_record *proteins;
struct LOC_build *LINK;
{
  /* tabulate i+2, i+3, .., i+maxl distances */
  long i, j, d, l, q, s;
  hit_record *newhit, *WITH;
  protlist_record *WITH1;
  long FORLIM2;

  /* init lookup table */
  for (d = 0; d <= maxdist; d++) {
    for (l = 1; l < maxl; l++) {
      WITH = &lookup[l][d];
      WITH->nexthit = NULL;
    }
  }
  /* loop through dist[] */
  for (i = 2; i <= nprot; i++) {   /* i = protein */
    WITH1 = &proteins[i - 1];
    for (l = 2; l <= maxl; l++) {
      FORLIM2 = WITH1->length - l;
      for (j = 1; j <= FORLIM2; j++) {
	s = WITH1->start + j;
	q = dist[s - 1].d[l];
	if (q < maxdist) {
	  WITH = &lookup[l - 1][q];
	  newhit = (hit_record *)Malloc(sizeof(hit_record));
	  newhit->nexthit = WITH->nexthit;
	  newhit->protno = i;
	  newhit->resno = j;
	  WITH->nexthit = newhit;
	}
      }
    }
  }
}

Local Void maintodist(master, LINK)
Char *master;
struct LOC_build *LINK;
{
  FILE *tapeout;
  protein_record *firstprotein, *prot1;
  residue_record *curresidue;
  long i, j, l;
  dist_array TEMPDIST;
  double x[maxres], y[maxres], z[maxres];
  xyz r[4];
  long k;
  main_record *WITH;
  double TEMP, TEMP1, TEMP2;
  char outfile[FILENAME_MAX];

  tapeout = NULL;
  firstprotein = NULL;
  if (*master != '\0') {
    getprotein(master, false, &firstprotein, &prot1, LINK);
    l = prot1->length;
    /* copy ca coordinates to array */
    curresidue = prot1->firstresidue;
    for (i = 0; i < l; i++) {
      WITH = &curresidue->backbone.main;
      x[i] = WITH->ca[0];
      y[i] = WITH->ca[1];
      z[i] = WITH->ca[2];
      curresidue = curresidue->nextresidue;
    }
    /* write protlength, calculated distances to tapeout */
    /* dihedralangle to tempdist[0] */
    sprintf(outfile, "%s.dist", master);
    tapeout = fopen(outfile, "w");
    if (tapeout != NULL)
      rewind(tapeout);
    else {
      fprintf(stderr, "ERROR: Unable to open output file: %s\n", outfile);
      fprintf(stderr, "ERROR: in function: maintodist()\n");
      fprintf(stderr, "Bye for now folks........\n\n");
      exit(1);
    }

    TEMPDIST[0] = 0;
    TEMPDIST[1] = l;
    for (i = 2; i <= maxlen; i++)
      TEMPDIST[i] = 0;
    fwrite(TEMPDIST, sizeof(dist_array), 1, tapeout);
    /* CALCULATE DISTANCES */
    for (i = 1; i < l; i++) {
      for (j = 1; j <= maxlen; j++) {
	if (i + j <= l) {
	  TEMP = x[i - 1] - x[i + j - 1];
	  TEMP1 = y[i - 1] - y[i + j - 1];
	  TEMP2 = z[i - 1] - z[i + j - 1];
	  TEMPDIST[j] = (long)floor(
		scaler * sqrt(TEMP * TEMP + TEMP1 * TEMP1 + TEMP2 * TEMP2) + 0.5);
	} else
	  TEMPDIST[j] = 0;
      }
      if (TEMPDIST[1] > 500 || TEMPDIST[1] < 200)
	printf("WARNING: weird distance %5ld%10ld%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f\n",
	       i, TEMPDIST[1], x[i - 1], y[i - 1], z[i - 1], x[i], y[i],
	       z[i]);
      if (i <= l - 3) {
	for (k = 1; k <= 4; k++) {
	  r[k - 1][0] = x[i + k - 2];
	  r[k - 1][1] = y[i + k - 2];
	  r[k - 1][2] = z[i + k - 2];
	}
	TEMPDIST[0] = (long)floor(dihedralangle(r[0], r[1], r[2], r[3], LINK) + 0.5);
      } else
	TEMPDIST[0] = 0;
      fwrite(TEMPDIST, sizeof(dist_array), 1, tapeout);
    }
    if (tapeout != NULL)
      fclose(tapeout);
    tapeout = NULL;
  }
  if (tapeout != NULL)
    fclose(tapeout);
}

/* Local variables for load_distances_: */
struct LOC_load_distances_ {
  struct LOC_build *LINK;
  long *nprot;
  protlist_record *proteins;
  dist_record *dist;
  protein_record **firstprotein;
  long where;
} ;

Local Void readdists(master_, LINK)
Char *master_;
struct LOC_load_distances_ *LINK;
{
  string master;
  protein_record *prot;
  residue_record *curresidue;
  FILE *tapein2;
  FILE *tapein1;
  long i;
  dist_array temp;
  protlist_record *WITH;
  long FORLIM;
  main_record *WITH1;
  char infile[FILENAME_MAX];

#ifdef DEBUG
  printf("In readdists()\n");
#endif

  strcpy(master, master_);
  sprintf(infile, "%s.main", master_);
#ifdef DEBUG
  printf("About to open file: %s\n", infile);
#endif
  tapein1 = fopen(infile, "rb");
  if (tapein1 == NULL) {
    printf("file %s not found.  Skipping protein.\n", infile);
    goto _Lskip;
  }
  rewind(tapein1);
  getprotein(master, false, LINK->firstprotein, &prot, LINK->LINK);
  (*LINK->nprot)++;
_Ll1:

  sprintf(infile, "%s.dist", master_);
  tapein2 = fopen(infile, "rb");
  if (tapein2 == NULL) {   /* create missing .dist file */
    printf("creating .dist file for %s\n", master_);
    maintodist(master, LINK->LINK);
    goto _Ll1;
  }
  rewind(tapein2);
  curresidue = prot->firstresidue;
  fread(temp, sizeof(dist_array), 1, tapein2);
  WITH = &LINK->proteins[*LINK->nprot - 1];
  WITH->length = temp[1];
  if (WITH->length != prot->length)
    printf(
      "WARNING: %sdistance and coordinate files have different no of residues !!!\n",
      master);
  WITH->start = LINK->where;
  strcpy(WITH->protname, master);
  FORLIM = WITH->length;
  for (i = 1; i < FORLIM; i++) {
    LINK->where++;
    fread(LINK->dist[LINK->where - 1].d, sizeof(dist_array), 1, tapein2);
    LINK->dist[LINK->where - 1].respo = curresidue;
    curresidue = curresidue->nextresidue;
  }
  LINK->where++;
  LINK->dist[LINK->where - 1].respo = curresidue;
  if (tapein2 != NULL)
    fclose(tapein2);
  tapein2 = NULL;
  if (tapein1 != NULL)
    fclose(tapein1);
  tapein1 = NULL;

  /* check chain continuity */
  curresidue = prot->firstresidue;
  i = LINK->proteins[*LINK->nprot - 1].start;
  while (curresidue != NULL) {
    i++;
    WITH1 = &curresidue->backbone.main;
    if (WITH1->ca[0] == 0) {
      if (WITH1->ca[1] == 0) {
	if (WITH1->ca[2] == 0)
	  printf("WARNING: missing CA coordinates for residue %12ld%.5s%.4s\n",
		 i, curresidue->backbone.resno, curresidue->backbone.resname);
      }
    }
    if (i < prot->length) {
      if (LINK->dist[i - 1].d[1] < 200 || LINK->dist[i - 1].d[1] > 480)
	printf(
	  "LOOK OUT !  C(alpha)-C(alpha-next) distance out of range for residue %12ld%.5s%.4s%12ld\n",
	  i, curresidue->backbone.resno, curresidue->backbone.resname,
	  LINK->dist[i - 1].d[1]);
    }
    curresidue = curresidue->nextresidue;
  }
_Lskip:
  if (tapein2 != NULL)
    fclose(tapein2);
  if (tapein1 != NULL)
    fclose(tapein1);

#ifdef DEBUG
  printf("In readdists()\n");
#endif

}
/* Process command line args. This part has been added so that the 
   program can be called as part of an automatic update procesdure
   on unix. The old VMS pascal version used a command file to pass
   the arguments.

   Chris Dodge - EBI - Oct 96

   */
void Print_Usage()
{
  printf("\nUSAGE: buildbackbone <args>\n\n");
  printf("              -h, -help   : Print this message.\n");
  printf("              -pdb <file> : PDB file name.\n");
  printf("              -pl  <file> : Protein list database file name.\n");
  printf("              -d   <Y/N>) : Do/don't use default parameters.\n");
  printf("\n");
  exit(1);
}

void Sort_Args(int argc, char *argv[])
{
  int i;

  DBName[0]   = '\0';
  FileName[0] = '\0';
  Default     = DEF_NOTDEF;  /* Default values not defined */

  for (i=1; i<argc; i++) {
    
    /* HELP */
    if ((!strcmp(argv[i], "-h"))||(!strcmp(argv[i],"-help"))
	||(!strcmp(argv[i], "-H"))||(!strcmp(argv[i],"-HELP"))) {
      Print_Usage();
    } 

    /* FILE NAME */
    /* OK, how this works is that if the '-pdb' is given, it looks to
       see if there is a following argument. If not, set the string to
       the default value. If there is a follwing arg, look to see if it
       is another flag, ie. starting with '-', in which case assume that
       we wanted the default value.
       
       If, on the other hand the value doesn't start with '-', then use 
       this value.

       Is that now clear?
       */
    else if ((!strcmp(argv[i], "-pdb"))||(!strcmp(argv[i], "-PDB"))) {
      if ((i+1)<argc) {
	if (strlen(argv[i+1])>=MAXARGLEN) {
	  fprintf(stderr, "ERROR: Sorry old chap. File name too long.\n");
	  exit(1);
	}
	if (argv[i+1][0]!='-') {
	  strcpy(FileName, argv[i+1]);
	  i++;
	} 
      }
    }

    /* Database of file names */
    else if ((!strcmp(argv[i], "-pl"))||(!strcmp(argv[i], "-PL"))) {
      if ((i+1)<argc) {
	if (strlen(argv[i+1])>=MAXARGLEN) {
	  fprintf(stderr, "ERROR: Sorry old chap. File name too long.\n");
	  exit(1);
	}
	if (argv[i+1][0]!='-') {
	  strcpy(DBName, argv[i+1]);
	  i++;
	} 
      }
    }

    /* Default values? */
    else if ((!strcmp(argv[i], "-d"))||(!strcmp(argv[i], "-D"))) {
      if ((i+1)<argc) {
	if ((argv[i+1][0]=='y')||(argv[i+1][0]=='Y')) {
	  Default = DEF_YES;
	} else if ((argv[i+1][0]=='n')||(argv[i+1][0]=='N')) {
	  Default = DEF_NO; 
	}
	i++;
	/* Otherwise Default remains not defined */
      }
    }

    /* .... ummm, pass. */
    else {
      printf("\n\n**** SORRY, DON'T UNDERSTAND %s\n", argv[i]);
      if (!strcmp(argv[i], "-hoo-ha"))
	printf("**** Maybe you meant: \"-ha-hoo\"?\n");
      else
	printf("**** Maybe you meant: \"-hoo-ha\"?\n");
      Print_Usage();
    }
  }
}

Local Void load_distances_(nprot_, master, proteins_, dist_, firstprotein_,
			   prot1, LINK)
long *nprot_;
Char *master;
protlist_record *proteins_;
dist_record *dist_;
protein_record **firstprotein_, **prot1;
struct LOC_build *LINK;
{
  struct LOC_load_distances_ V;
  string prot;
  FILE *tapein;
  Char *TEMP;

  V.LINK = LINK;
  V.nprot = nprot_;
  V.proteins = proteins_;
  V.dist = dist_;
  V.firstprotein = firstprotein_;
  tapein = NULL;
  *V.nprot = 0;
  V.where = 0;
  /* read test-protein to 'slot 1' */
  if (strlen(FileName)==0) {
    printf("enter name of testprotein (e.g. 1sn3) \n");
    fgets(master, 81, stdin);
    TEMP = strchr(master, '\n');
    if (TEMP != NULL)
      *TEMP = 0;
#ifdef DEBUG
    printf("Input protein name is: %s\n", master);
#endif
  } else {
    strcpy(master, FileName);
    printf("Input protein name is: %s\n", master);
  }
  readdists(master, &V);
  printf("%12ld residues in test-protein\n", V.where);
  *V.firstprotein = NULL;
  getprotein(master, false, V.firstprotein, prot1, LINK);
  /* read database distances */
  strcpy(prot, master);
  if (strlen(DBName)==0) {
    printf("enter name of database protein names file (e.g. dglp.list) \n");
    fgets(prot, 81, stdin);
    TEMP = strchr(prot, '\n');
    if (TEMP != NULL)
      *TEMP = 0;
  } else {
    strcpy(prot, DBName);
    printf("Input protein list file name is: %s\n", prot);
  }
  if (openfile_(&tapein, prot, infile, LINK)==-1) {
    fprintf(stderr, "Can't open protein list file %s.\n", prot);
    fprintf(stderr, "Exiting......\n\n");
    exit(1);
  }
  while (fgets(prot, 81, tapein) != NULL) {
    TEMP = strchr(prot, '\n');
    if (TEMP != NULL)
      *TEMP = 0;
    printf("loading %s to %12ld%12ld\n", prot, V.where, *V.nprot + 1);
    readdists(prot, &V);
  }
  if (tapein != NULL)
    fclose(tapein);
  tapein = NULL;
  printf("%12ld proteins read from distance library\n", *V.nprot - 1);
  printf("%12ld residues in distance library\n",
	 V.where - V.proteins[0].length);
  if (tapein != NULL)
    fclose(tapein);
}


#undef datadir

/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */

Static Void build()
{
  /* struct LOC_build V; */

  /* Allocate V on heap instead of on stack.  Changed to avoid crashes
     on Unix machines with restrictions on stack size (LOC_build is rather large).
     
     Wouter Boomsma, 2006
  */
  struct LOC_build *V = (struct LOC_build *)malloc(sizeof(struct LOC_build));
  long i, cs;
  protein_record *firstprotein;
  Char c;
  long goal, nprot, FORLIM;


  V->logfile = NULL;
  V->coorfile = NULL;
  V->statfile = NULL;
  /* init all variables */
  firstprotein = NULL;
  load_distances_(&nprot, V->master, V->proteins, V->dist, &firstprotein,
		  &V->prot1, V);
  fill_lookup_(V->lookup, V->dist, nprot, V->proteins, V);
  V->chainstart = 1;
  V->chainend = V->prot1->length;
  FORLIM = V->chainend;
  for (i = V->chainstart - 1; i < FORLIM; i++) {
    V->filled[i] = false;
    V->gap[i] = true;
    V->fragindex[i] = 0;
    V->redundancy[i] = 0;
    V->fragmentlist[i] = NULL;
    V->nscore[i] = 0;
  }
  V->cb_only = true;   /* option to build cbetas or whole sidechain */
  V->leftzone = 6;
  V->rightzone = 2;
  V->minlooplen = 4;
      /* looplen can be relaxed to fetch at least maxscores fragments */
  V->maxscores = maxmaxscores;
  V->tolerance = 0.4;
  V->l1 = 5;
  V->c1 = 20;
  V->l2 = 4;
  V->c2 = 80;
  sprintf(V->coorfilename, "%s.brk_mod", V->master);
  sprintf(V->logfilename, "%s.log", V->master);
  sprintf(V->statfilename, "%s.stat", V->master);
  if (Default == DEF_NOTDEF) {
    printf("use default parameters [Y]/N ?\n");
    scanf("%c%*[^\n]", &c);
    getchar();
  }
  if (c == '\n')
    c = ' ';
  if (c == 'n' || c == 'N' || Default == DEF_NO)
    change_parameters(V);
  openfile_(&V->coorfile, V->coorfilename, outfile, V);
  openfile_(&V->logfile, V->logfilename, outfile, V);
  openfile_(&V->statfile, V->statfilename, outfile, V);
  /*openfile(fragfile, 'frag.pdb', outfile);*/
  write_headers(V);

  V->looplen = V->l1;
  V->cutoff1 = V->c2;
  quickfill(V);
  V->looplen = V->l2;
  V->cutoff1 = V->c2;
  update_gap(V->chainstart, V->chainstart, V->chainend, V);
  cs = V->chainstart;
  do {
    seal(cs, V->chainend, &goal, V);
    if (goal <= V->chainend) {
      printf("Cannot seal gap %12ld%12ld\n", goal - 1, goal);
      if (goal - 1 > cs)
	path(cs, goal - 1, V);
      fprintf(V->logfile, "REMARK      Cannot seal gap %12ld%12ld\n",
	      goal - 1, goal);
      if (goal == cs)
	goal++;
      cs = goal;
      FORLIM = V->chainend;
      for (i = V->chainstart - 1; i < FORLIM; i++) {
	V->filled[i] = false;
	V->gap[i] = true;
	V->fragindex[i] = 0;
	V->redundancy[i] = 0;
	V->fragmentlist[i] = NULL;
	V->nscore[i] = 0;
      }
      update_gap(cs, cs, V->chainend, V);
      goal = firstgap(cs, V->chainend, V);
      printf("Starting new chain from %12ld\n", cs);
    }
  } while (goal <= V->chainend && V->chainend - cs >= V->looplen - 1);
      /* won't find fragments for very short chains */
  if (goal - 1 > cs)
    path(cs, goal - 1, V);

  if (V->coorfile != NULL)
    fclose(V->coorfile);
  V->coorfile = NULL;
  if (V->logfile != NULL)
    fclose(V->logfile);
  V->logfile = NULL;
  logstatistics(V);
  if (V->statfile != NULL)
    fclose(V->statfile);
  V->statfile = NULL;
  if (V->statfile != NULL)
    fclose(V->statfile);
  if (V->coorfile != NULL)
    fclose(V->coorfile);
  if (V->logfile != NULL)
    fclose(V->logfile);
  free(V);  
}  /*build */

#undef maxmaxscores
#undef maxprotlength


main(argc, argv)
int argc;
Char *argv[];
{
  printf("__________________________________________________________________________ \n");
  printf("                                                                           \n");
  printf("This set of programs constructs full atomic coordinates of a protein       \n");
  printf("from a given C(alpha) trace and optimizes side chain geometry.             \n");
  printf("                                                                           \n");
  printf("Copyright by Liisa Holm and Chris Sander, 1989-1991.                       \n");
  printf("                                                                           \n");
  printf("No redistribution, no program changes, no commercial use.                  \n");
  printf("                                                                           \n");
  printf("For details, see J.Mol.Biol. 218, 183-194 (1991).                          \n");
  printf("                                                                           \n");
  printf("___________________________________________________________________________\n");

  PASCAL_MAIN(argc, argv);
  Sort_Args(argc, argv);
  build();
  exit(EXIT_SUCCESS);
}
/* p2c: buildbackbone.pas, line 2036: 
 * Warning: Junk at end of input file ignored [277] */
/* module */




/* End. */
