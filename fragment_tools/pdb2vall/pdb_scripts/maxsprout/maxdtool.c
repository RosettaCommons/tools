/* Output from p2c, the Pascal-to-C translator */
/* From input file "maxdtool.pas" */


#include <p2c/p2c.h>
/* p2c: maxdtool.pas, line 1: 
 * Warning: Expected a colon, found a ')' [227] */


#define maxlen          30
#define maxres          10000
#define scaler          100
#define maxprot         100


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


typedef struct scorelist_record *scorepointersarray[maxres];

typedef struct scorelist_record {
  long begresidue, matchresidue, matchlength, matchprotlistindex;
  string matchprotname;
  double score;
  struct scorelist_record *prevscore, *nextscore;
} scorelist_record;

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


#define maxscores       1000


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

/*procedure u3b(w: t1; x, y: t2; n, mode: integer; var rms: real;
 var u:rotationmatrix; var t: translationvector; var ier: integer);
 FORTRAN; */
/* p2c: maxdtool.pas, line 159: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 166: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 181: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 202: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 228: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 246: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 259: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 286: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 294: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 305: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 313: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 324: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 342: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 362: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 383: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 392: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 395: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 398: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 405: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 423: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 449: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 464: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 480: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 496: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 512: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 583: 
 * Warning: Type attribute GLOBAL ignored [128] */


#define superscale      1000
/* p2c: maxdtool.pas, line 624: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 643: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 663: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 672: 
 * Warning: Type attribute GLOBAL ignored [128] */


#define superscaler     1000
/* p2c: maxdtool.pas, line 711: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 723: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 742: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 768: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 827: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 929: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 956: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 989: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 1023: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 1041: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 1076: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 1110: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 1115: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 1137: 
 * Warning: Type attribute GLOBAL ignored [128] */


#define datadir         "csdata:"
/* p2c: maxdtool.pas, line 1177: 
 * Warning: Type attribute GLOBAL ignored [128] */
/* p2c: maxdtool.pas, line 1195: 
 * Warning: Type attribute GLOBAL ignored [128] */


Local Void x()
{
}

Local Void init_residuelist(firstresidue, lastresidue)
residue_record **firstresidue, **lastresidue;
{
  *firstresidue = NULL;
  *lastresidue = NULL;
}

Local long val(s)
Char *s;
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

Local long resnametocode(resnm)
Char *resnm;
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

Local long atomnametocode(atnm)
Char *atnm;
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


Local Void save_backbone(outfilename, protein)
Char *outfilename;
protein_record *protein;
{
  FILE *tapeout;
  residue_record *curresidue;

/* p2c: maxdtool.pas, line 235:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 235: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 235:
 * Note: Ignoring block size in binary file [182] */
  tapeout = NULL;
  if (tapeout != NULL)
    rewind(tapeout);
  else
    tapeout = tmpfile();
  if (tapeout == NULL)
    _EscIO(FileNotFound);
  if (tapeout != NULL) {
/* p2c: maxdtool.pas, line 235:
 * Warning: Expected a ')', found a ':=' [227] */
/* p2c: maxdtool.pas, line 236:
 * Note: REWRITE does not specify a name [181] */
    rewind(tapeout);
  } else
    tapeout = tmpfile();
  if (tapeout == NULL)
    _EscIO(FileNotFound);
  curresidue = protein->firstresidue;
  while (curresidue != NULL) {
    fwrite(&curresidue->backbone, sizeof(backbone_record), 1, tapeout);
    curresidue = curresidue->nextresidue;
  }
  if (tapeout != NULL)
    fclose(tapeout);
  tapeout = NULL;
  if (tapeout != NULL)
    fclose(tapeout);
}

Local Void transrot(r, rnew, u, t)
double *r, *rnew;
double (*u)[3];
double *t;
{
  long i, j;

  /* rnew = u*r+t */
  for (i = 0; i <= 2; i++) {
    rnew[i] = t[i];
    for (j = 0; j <= 2; j++)
      rnew[i] += u[j][i] * r[j];
  }
}

Local Void matchca(residue1, residue2, length, u, t, rms)
residue_record *residue1, *residue2;
long length;
double (*u)[3];
double *t;
double *rms;
{
  t2 x, y;   /* coordinates */
  long j, ier;
  residue_record *res1, *res2;
  t1 w;   /* atom pair weight factors */

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

  u3b(w, x, y, length, 1, *rms, u, t, ier);
/* p2c: maxdtool.pas, line 281:
 * Warning: Symbol 'U3B' is not defined [221] */

}


Local Void write_brk(atomname, x, pdbname, segid, brkout)
Char *atomname;
double *x;
Char *pdbname;
Char *segid;
FILE **brkout;
{
  fprintf(*brkout,
	  "ATOM         %.4s%.10s   %8.3f%8.3f%8.3f                  %.4s\n",
	  atomname, pdbname, x[0], x[1], x[2], segid);
}

Local Void mainrecord_to_brk(main, backbone, tapeout)
main_record main;
backbone_record backbone;
FILE **tapeout;
{
  write_brk("N   ", main.n, backbone.pdbresid, "    ", tapeout);
  write_brk("CA  ", main.ca, backbone.pdbresid, "    ", tapeout);
  write_brk("C   ", main.c, backbone.pdbresid, "    ", tapeout);
  write_brk("O   ", main.o, backbone.pdbresid, "    ", tapeout);
}


Local Void oldwrite_brk(atomname, x, resname, chainid, resno, segid, brkout)
Char *atomname;
double *x;
Char *resname;
Char chainid;
Char *resno;
Char *segid;
FILE **brkout;
{
  fprintf(*brkout,
    "ATOM         %.4s%.4s%c%.5s   %8.3f%8.3f%8.3f                  %.4s\n",
    atomname, resname, chainid, resno, x[0], x[1], x[2], segid);
}

Local Void oldmainrecord_to_brk(main, resname, resno, tapeout)
main_record main;
Char *resname;
Char *resno;
FILE **tapeout;
{
  oldwrite_brk("N   ", main.n, resname, ' ', resno, "    ", tapeout);
  oldwrite_brk("CA  ", main.ca, resname, ' ', resno, "    ", tapeout);
  oldwrite_brk("C   ", main.c, resname, ' ', resno, "    ", tapeout);
  oldwrite_brk("O   ", main.o, resname, ' ', resno, "    ", tapeout);
}


Local Char *str_(Result, v)
Char *Result;
long v;
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


Local Char *codetoresname(Result, code)
Char *Result;
long code;
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


Local Char *codetoatomname(Result, code)
Char *Result;
long code;
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

Local double atan2(y, x)
double y, x;
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

Local Void diff(x, y, z)
double *x, *y, *z;
{
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  z[2] = x[2] - y[2];
}

Local double dot(x, y)
double *x, *y;
{
  return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

Local Void cross(x, y, z)
double *x, *y, *z;
{
  z[0] = x[1] * y[2] - y[1] * x[2];
  z[1] = x[2] * y[0] - y[2] * x[0];
  z[2] = x[0] * y[1] - y[0] * x[1];
}

Local double dihedralangle(v1, v2, v3, v4)
double *v1, *v2, *v3, *v4;
{
  double Result, u, v;
  xyz v12, v43, x, y, z, p;

  diff(v1, v2, v12);
  diff(v4, v3, v43);
  diff(v2, v3, z);
  cross(z, v12, p);
  cross(z, v43, x);
  cross(z, x, y);
  u = dot(x, x);
  v = dot(y, y);
  Result = 360.0;
  if (u <= 0.0 || v <= 0.0)
    return Result;
  u = dot(p, x) / sqrt(u);
  v = dot(p, y) / sqrt(v);
  if (u != 0.0 || v != 0.0)
    return (atan2(v, u) * 57.29578);
  return Result;
}

Local Void openfile(filvar, filename, mode)
FILE **filvar;
Char *filename;
inout mode;
{
  long x;

  switch (mode) {

  case infile:
    if (*filvar != NULL)
      rewind(*filvar);
    else
      *filvar = tmpfile();
    if (*filvar == NULL)
      _EscIO(FileNotFound);
    rewind(*filvar);   /**, error:= continue); **/
    x = status(*filvar);
/* p2c: maxdtool.pas, line 432:
 * Warning: Symbol 'STATUS' is not defined [221] */
    printf("status opening %s = %12ld\n", filename, x);
    if (status(*filvar) != 0)
      printf("file %s not found. \n", filename);
    else
      printf("reading file %s\n", filename);
/* p2c: maxdtool.pas, line 434:
 * Warning: Symbol 'STATUS' is not defined [221] */
    break;
/* p2c: maxdtool.pas, line 430:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 430: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 430:
 * Note: Ignoring block size in binary file [182] */

  case outfile:
    if (*filvar != NULL)
      rewind(*filvar);
    else
      *filvar = tmpfile();
    if (*filvar == NULL)
      _EscIO(FileNotFound);
    if (*filvar != NULL)
      rewind(*filvar);
    else
      *filvar = tmpfile();
    if (*filvar == NULL)
      _EscIO(FileNotFound);
    printf("writing to %s\n", filename);
    break;
/* p2c: maxdtool.pas, line 439:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 439: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 439:
 * Note: Ignoring block size in binary file [182] */
  }

/* p2c: maxdtool.pas, line 430:
 * Warning: Expected a ')', found a ':=' [227] */
/* p2c: maxdtool.pas, line 439:
 * Warning: Expected a ')', found a ':=' [227] */
/* p2c: maxdtool.pas, line 440:
 * Note: REWRITE does not specify a name [181] */
}


/* make list - procedures */

Local Void make_new_protein(firstprotein)
protein_record **firstprotein;
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

Local Void make_new_residue(protein, firstresidue, lastresidue)
protein_record *protein;
residue_record **firstresidue, **lastresidue;
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

Local residue_record *findresnoresname(protein, resno, resname)
protein_record *protein;
Char *resno;
Char *resname;
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

Local Void makeside(residue, atnm, cd)
residue_record *residue;
Char *atnm;
double *cd;
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

Local Void brk_to_list(lastprotein, firstresidue, infilename, headerfile,
		       length)
protein_record *lastprotein;
residue_record **firstresidue;
Char *infilename, *headerfile;
long *length;
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
  openfile(&tapein, infilename, infile);
  init_residuelist(firstresidue, &lastresidue);
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
      make_new_residue(lastprotein, firstresidue, &lastresidue);
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
	    makeside(lastresidue, atnm, cd);
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

Local Void list_to_side(protein, outfilename)
protein_record *protein;
Char *outfilename;
{
  FILE *tapeout;
  residue_record *curresidue;
  atom_record *currside;
  long atoms[100][4];
  long i, j, n;
  a4 resn;
  backbone_record *WITH;
  long TEMP;

/* p2c: maxdtool.pas, line 594:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 594: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 594:
 * Note: Ignoring block size in binary file [182] */
  tapeout = NULL;
  if (tapeout != NULL)
    rewind(tapeout);
  else
    tapeout = tmpfile();
  if (tapeout == NULL)
    _EscIO(FileNotFound);
  if (tapeout != NULL) {
/* p2c: maxdtool.pas, line 594:
 * Warning: Expected a ')', found a ':=' [227] */
/* p2c: maxdtool.pas, line 595:
 * Note: REWRITE does not specify a name [181] */
    rewind(tapeout);
  } else
    tapeout = tmpfile();
  if (tapeout == NULL)
    _EscIO(FileNotFound);
  fwrite(&protein->length, sizeof(long), 1, tapeout);
  curresidue = protein->firstresidue;
  while (curresidue != NULL) {
    WITH = &curresidue->backbone;
    resn[0] = WITH->resno[0];
    resn[1] = WITH->resno[1];
    resn[2] = WITH->resno[2];
    resn[3] = WITH->resno[3];
    TEMP = resnametocode(WITH->resname);
    fwrite(&TEMP, sizeof(long), 1, tapeout);
    TEMP = val(resn);
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
      atoms[n - 1][0] = atomnametocode(currside->atomname);
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

Local double distance(a_, b_)
double *a_, *b_;
{
  double Result;
  xyz a, b, c;
  double d;
  long i;

  memcpy(a, a_, sizeof(xyz));
  memcpy(b, b_, sizeof(xyz));
  diff(a, b, c);
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

Local Void backbone_to_list(lastprotein, firstresidue, infilename, length)
protein_record *lastprotein;
residue_record **firstresidue;
Char *infilename;
long *length;
{
  FILE *tapein;
  residue_record *lastresidue;

/* p2c: maxdtool.pas, line 650:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 650: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 650:
 * Note: Ignoring block size in binary file [182] */
  tapein = NULL;
  if (tapein != NULL)
    rewind(tapein);
  else
    tapein = tmpfile();
  if (tapein == NULL)
    _EscIO(FileNotFound);
/* p2c: maxdtool.pas, line 650:
 * Warning: Expected a ')', found a ':=' [227] */
  rewind(tapein);
  *firstresidue = NULL;
  lastresidue = NULL;
  *length = 0;
  while (!P_eof(tapein)) {
    make_new_residue(lastprotein, firstresidue, &lastresidue);
    (*length)++;
    fread(&lastresidue->backbone, sizeof(backbone_record), 1, tapein);
  }
  if (tapein != NULL)
    fclose(tapein);
}

Local Void fill_protein_record(protein, first, proteinname, l)
protein_record *protein;
residue_record *first;
Char *proteinname;
long l;
{
  strcpy(protein->protname, proteinname);
  protein->length = l;
  protein->firstresidue = first;
}


Local Void side_to_list(protein, infilename)
protein_record *protein;
Char *infilename;
{
  FILE *tapein;
  residue_record *curresidue;
  long i, j, n, code, a, b, c, length;
  xyz cd;
  a4 atnm, resname;
  a5 resno;

/* p2c: maxdtool.pas, line 685:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 685: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 685:
 * Note: Ignoring block size in binary file [182] */
  tapein = NULL;
  if (tapein != NULL)
    rewind(tapein);
  else
    tapein = tmpfile();
  if (tapein == NULL)
    _EscIO(FileNotFound);
/* p2c: maxdtool.pas, line 685:
 * Warning: Expected a ')', found a ':=' [227] */
  rewind(tapein);
  /* assume .MAIN has been read & protein-list initialized */
  fread(&length, sizeof(long), 1, tapein);
  for (n = 1; n <= length; n++) {   /* residue-counter */
    fread(&code, sizeof(long), 1, tapein);
    codetoresname(resname, code);
    fread(&code, sizeof(long), 1, tapein);
    str_(resno, code);
    curresidue = findresnoresname(protein, resno, resname);
    if (curresidue == NULL)
      _Escape(0);
    fread(&i, sizeof(long), 1, tapein);
    for (j = 1; j <= i; j++) {   /* atom-counter */
      fread(&code, sizeof(long), 1, tapein);
      fread(&a, sizeof(long), 1, tapein);
      fread(&b, sizeof(long), 1, tapein);
      fread(&c, sizeof(long), 1, tapein);
      codetoatomname(atnm, code);
      cd[0] = (double)a / superscaler;
      cd[1] = (double)b / superscaler;
      cd[2] = (double)c / superscaler;
      makeside(curresidue, atnm, cd);
    }
  }
  if (tapein != NULL)
    fclose(tapein);
}

#undef superscaler


Local protein_record *findproteinpointer(master, firstprotein)
Char *master;
protein_record *firstprotein;
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

Local Void getprotein(master, sidechain, firstprotein, protein)
Char *master;
boolean sidechain;
protein_record **firstprotein, **protein;
{
  /* return pointer to protein master - read coordinates if it is not in list */
  residue_record *firstresidue;
  long l;
  Char STR1[86];

  *protein = findproteinpointer(master, *firstprotein);
  if (*protein != NULL)
    return;
  make_new_protein(firstprotein);
  sprintf(STR1, "%s.main", master);
  backbone_to_list(*firstprotein, &firstresidue, STR1, &l);
  if (sidechain) {
    sprintf(STR1, "%s.side", master);
    side_to_list(*firstprotein, STR1);
  }
  fill_protein_record(*firstprotein, firstresidue, master, l);
  *protein = *firstprotein;
}

Local Void fill_lookup(lookup, dist, nprot, proteins)
hit_record (*lookup)[maxdist + 1];
dist_record *dist;
long nprot;
protlist_record *proteins;
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

Local Void maintodist(master)
Char *master;
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

  tapeout = NULL;
  firstprotein = NULL;
  if (*master != '\0') {
    getprotein(master, false, &firstprotein, &prot1);
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
/* p2c: maxdtool.pas, line 795:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 795: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 795:
 * Note: Ignoring block size in binary file [182] */
    if (tapeout != NULL)
      rewind(tapeout);
    else
      tapeout = tmpfile();
    if (tapeout == NULL)
      _EscIO(FileNotFound);
    if (tapeout != NULL) {
/* p2c: maxdtool.pas, line 795:
 * Warning: Expected a ')', found a ':=' [227] */
/* p2c: maxdtool.pas, line 796:
 * Note: REWRITE does not specify a name [181] */
      rewind(tapeout);
    } else
      tapeout = tmpfile();
    if (tapeout == NULL)
      _EscIO(FileNotFound);
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
	TEMPDIST[0] = (long)floor(dihedralangle(r[0], r[1], r[2], r[3]) + 0.5);
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

/* Local variables for load_distances: */
struct LOC_load_distances {
  long *nprot;
  protlist_record *proteins;
  dist_record *dist;
  protein_record **firstprotein;
  long where;
} ;

Local Void readdists(master_, LINK)
Char *master_;
struct LOC_load_distances *LINK;
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

  /*open(tapein1, file_name := master+'.main', history := readonly,
                  error := continue);
  */
/* p2c: maxdtool.pas, line 850:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 850: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 850:
 * Note: Ignoring block size in binary file [182] */
  strcpy(master, master_);
  tapein1 = NULL;
  tapein2 = NULL;
  if (tapein1 != NULL)
    rewind(tapein1);
  else
    tapein1 = tmpfile();
  if (tapein1 == NULL)
    _EscIO(FileNotFound);
/* p2c: maxdtool.pas, line 850:
 * Warning: Expected a ')', found a ':=' [227] */
  if (status(tapein1) != 0) {
    printf("%s.main not found.  Skipping protein. Status= %12ld\n",
	   master, status(tapein1));
    goto _Lskip;
  }
/* p2c: maxdtool.pas, line 851:
 * Warning: Symbol 'STATUS' is not defined [221] */
/* p2c: maxdtool.pas, line 852:
 * Warning: Symbol 'STATUS' is not defined [221] */
/* p2c: maxdtool.pas, line 854:
 * Warning: Symbol 'ERROR' is not defined [221] */
/* p2c: maxdtool.pas, line 854:
 * Note: Can't interpret name argument in RESET [180] */
/* p2c: maxdtool.pas, line 854:
 * Note: Ignoring block size in binary file [182] */
  rewind(tapein1);
/* p2c: maxdtool.pas, line 854:
 * Warning: Expected a ')', found a ':=' [227] */
  getprotein(master, false, LINK->firstprotein, &prot);
  (*LINK->nprot)++;
_Ll1:
/* p2c: maxdtool.pas, line 857:
 * Warning: Symbol 'FILE_NAME' is not defined [221] */
/* p2c: maxdtool.pas, line 857: Note: OPEN does not specify a name [181] */
/* p2c: maxdtool.pas, line 857:
 * Note: Ignoring block size in binary file [182] */
  if (tapein2 != NULL)
    rewind(tapein2);
  else
    tapein2 = tmpfile();
  if (tapein2 == NULL)
    _EscIO(FileNotFound);
/* p2c: maxdtool.pas, line 857:
 * Warning: Expected a ')', found a ':=' [227] */
  if (status(tapein2) != 0) {   /* create missing .dist file */
    printf("creating .dist file for %sstatus=%12ld\n",
	   master, status(tapein2));
    maintodist(master);
    goto _Ll1;
  }
/* p2c: maxdtool.pas, line 859:
 * Warning: Symbol 'STATUS' is not defined [221] */
/* p2c: maxdtool.pas, line 861:
 * Warning: Symbol 'STATUS' is not defined [221] */
/* p2c: maxdtool.pas, line 865:
 * Warning: Symbol 'ERROR' is not defined [221] */
/* p2c: maxdtool.pas, line 865:
 * Note: Can't interpret name argument in RESET [180] */
/* p2c: maxdtool.pas, line 865:
 * Note: Ignoring block size in binary file [182] */
  rewind(tapein2);
/* p2c: maxdtool.pas, line 865:
 * Warning: Expected a ')', found a ':=' [227] */
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
}

Local Void load_distances(nprot_, master, proteins_, dist_, firstprotein_,
			  prot1)
long *nprot_;
Char *master;
protlist_record *proteins_;
dist_record *dist_;
protein_record **firstprotein_, **prot1;
{
  struct LOC_load_distances V;
  string prot;
  FILE *tapein;
  Char *TEMP;

  V.nprot = nprot_;
  V.proteins = proteins_;
  V.dist = dist_;
  V.firstprotein = firstprotein_;
  tapein = NULL;
  *V.nprot = 0;
  V.where = 0;
  /* read test-protein to 'slot 1' */
  printf("enter name of testprotein (e.g. 1sn3) \n");
  fgets(master, 81, stdin);
  TEMP = strchr(master, '\n');
  if (TEMP != NULL)
    *TEMP = 0;
  readdists(master, &V);
  printf("%12ld residues in test-protein\n", V.where);
  *V.firstprotein = NULL;
  getprotein(master, false, V.firstprotein, prot1);
  /* read database distances */
  strcpy(prot, master);
  printf("enter name of database protein names file (e.g. dglp.list) \n");
  fgets(prot, 81, stdin);
  TEMP = strchr(prot, '\n');
  if (TEMP != NULL)
    *TEMP = 0;
  openfile(&tapein, prot, infile);
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

Local Void extend_to_rigth(j, begresidue, testlen, length, cutoff1,
			   comparestart, templatestart, dist, nright, s)
long j, begresidue, testlen, length, cutoff1, comparestart, templatestart;
dist_record *dist;
long *nright;
double *s;
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

Local Void extend_to_left(j, begresidue, testlen, length, cutoff1,
			  comparestart, templatestart, dist, nleft, s, nright)
long j, begresidue, testlen, length, cutoff1, comparestart, templatestart;
dist_record *dist;
long *nleft;
double *s;
long nright;
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
	     begresidue - templatestart - k == 1 ||
	     l + k - 1 == testlen || diverge || l + k == maxlen));
  if (diverge)
    *nleft = k - 1;
  else
    *nleft = k;
}


Local Void update_pathsums(start, goal, fragmentlist)
long start, goal;
fragment_record **fragmentlist;
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

Local Void jump(jumpfrom, jumpto, rms, suurin)
fragment_record *jumpfrom, *jumpto;
double *rms, *suurin;
{
  double delta[5];
  double sum;
  long i;
  double TEMP;

  delta[0] = distance(jumpfrom->main2.main.ca, jumpto->main1.main.ca);
  delta[1] = distance(jumpfrom->main2.main.n, jumpto->main1.main.n);
  delta[2] = distance(jumpfrom->main1.main.c, jumpto->main0.main.c);
  delta[3] = distance(jumpfrom->main1.main.ca, jumpto->main0.main.ca);
  delta[4] = distance(jumpfrom->main1.main.o, jumpto->main0.main.o);
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

Local Void update_jump_fromnewtoold(fragmentlist, k, begresidue, testlen,
				    newfragment, tolerance)
fragment_record **fragmentlist;
long k, begresidue, testlen;
fragment_record **newfragment;
double tolerance;
{
  fragment_record *jumpto;
  jump_record *newjump, *lastjump;
  double rms, suurin;

  if (k + begresidue > testlen)
    return;
  jumpto = fragmentlist[k + begresidue - 1];
  while (jumpto != NULL) {
    jump(*newfragment, jumpto, &rms, &suurin);
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

Local Void update_jump_fromoldtonew(fragmentlist, k, begresidue, testlen,
				    newfragment, tolerance)
fragment_record **fragmentlist;
long k, begresidue, testlen;
fragment_record **newfragment;
double tolerance;
{
  fragment_record *jumpfrom;
  jump_record *newjump, *lastjump;
  double rms, suurin;

  if (k + begresidue <= 2)
    return;
  jumpfrom = fragmentlist[k + begresidue - 3];
  while (jumpfrom != NULL) {
    jump(jumpfrom, *newfragment, &rms, &suurin);
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

Local Void init_fragmentlist(fragmentlist, len)
fragment_record **fragmentlist;
long len;
{
  long i;

  for (i = 0; i < len; i++)
    fragmentlist[i] = NULL;
}

Local Void generate_cb(n, ca, c0, cb)
double *n, *ca, *c0, *cb;
{
  long i;
  xyz a, b, c, d;
  double lc, ld, TEMP;

  /* a = ca to n; b = ca to c; c = right-vector; d = up-vector; cb = ca + 1.32 * right-vector + 0.765 * up-vector */
  diff(n, ca, a);
  diff(c0, ca, b);
  cross(a, b, c);   /* c = a x b */
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


Local Void generate_sidechains(ca_, cb_, c_, n_, residue, firstprotein)
double *ca_, *cb_, *c_, *n_;
residue_record *residue;
protein_record **firstprotein;
{

  /* add sidechain atoms to residuepointer^.firstside^... */
  xyz ca, cb, c, n;
  protein_record *prot1;
  t1 w;
  t2 x, y;
  long ier;
  double rms;
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
  getprotein(STR1, true, firstprotein, &prot1);
  WITH = &prot1->firstresidue->backbone.main;
  /* get transrot matrices for n-ca-cb-c */
  generate_cb(WITH->n, WITH->ca, WITH->c, cbeta);
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
  u3b(w, x, y, 4, 1, rms, u, t, ier);
/* p2c: maxdtool.pas, line 1165:
 * Warning: Symbol 'U3B' is not defined [221] */
  /* copy transrotated sidechain atoms to residue^.firstside... */
  curatom = prot1->firstresidue->firstside;
  while (curatom != NULL) {
    transrot(curatom->coord, rnew, u, t);
    makeside(residue, curatom->atomname, rnew);
    curatom = curatom->nextatom;
  }
_Lquit: ;
}

#undef datadir

Local Void add_sidechains(protein, firstprotein)
protein_record *protein, **firstprotein;
{
  residue_record *curresidue;
  xyz cb;
  backbone_record *WITH;

  curresidue = protein->firstresidue;
  while (curresidue != NULL) {
    WITH = &curresidue->backbone;
    if (curresidue->firstside == NULL)
      generate_cb(WITH->main.n, WITH->main.ca, WITH->main.c, cb);
    else
      memcpy(cb, curresidue->firstside->coord, sizeof(xyz));
    generate_sidechains(WITH->main.ca, cb, WITH->main.c, WITH->main.n,
			curresidue, firstprotein);
    curresidue = curresidue->nextresidue;
  }
}

Local Void read_brk()
{
  /* read backbone file to list, do list-to-brk */
  protein_record *firstprotein;
  residue_record *firstresidue;
  string readdir, writedir, master, ext;
  long protlength, i, j;
  Char STR1[256];
  Char *TEMP;
  long FORLIM;
  Char STR3[242];
  Char STR4[170];
  Char STR5[168];

  printf("enter readdirectory (e.g. cs:) \n");
  fgets(readdir, 81, stdin);
  TEMP = strchr(readdir, '\n');
  if (TEMP != NULL)
    *TEMP = 0;
  printf("enter writedirectory (e.g. work$disk:) \n");
  fgets(writedir, 81, stdin);
  TEMP = strchr(writedir, '\n');
  if (TEMP != NULL)
    *TEMP = 0;
  firstprotein = NULL;
  do {
    printf("enter filename with extension [.brk] (e.g. 1sn3.brk_mod) ");
    fgets(master, 81, stdin);
    TEMP = strchr(master, '\n');
    if (TEMP != NULL)
      *TEMP = 0;
    if (*master != '\0') {
      /* parse filename > extension */
      j = strlen(master) + 1;
      FORLIM = strlen(master);
      for (i = 1; i <= FORLIM; i++) {
	if (master[i - 1] == '.')
	  j = i;
      }
      if (j <= strlen(master))
	sprintf(ext, "%.*s", (int)(strlen(master) - j + 1), master + j - 1);
      else
	strcpy(ext, ".brk");
      sprintf(master, "%.*s", (int)(j - 1), strcpy(STR1, master));
      make_new_protein(&firstprotein);
      sprintf(STR3, "%s%s%s", readdir, master, ext);
      sprintf(STR4, "%s%s.brk_hdr", writedir, master);
      brk_to_list(firstprotein, &firstresidue, STR3, STR4, &protlength);
      fill_protein_record(firstprotein, firstresidue, master, protlength);
      sprintf(STR5, "%s%s.main", writedir, master);
      save_backbone(STR5, firstprotein);
      sprintf(STR5, "%s%s.side", writedir, master);
      list_to_side(firstprotein, STR5);
    }
  } while (*master != '\0');
}


Static Void maxd_toolbox(input, output)
Anyptr input, output;
{
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
}

#undef maxlen
#undef maxres
#undef scaler
#undef maxprot
#undef maxscores
#undef maxl
#undef maxdist
/* p2c: maxdtool.pas, line 1231: 
 * Warning: Expected a semicolon, found a '.' [227] */


/* module */





/* End. */
