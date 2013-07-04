program readbrk(input, output); 

(*
__________________________________________________________________________

This set of programs constructs full atomic coordinates of a protein
from a given C(alpha) trace and optimizes side chain geometry.

Copyright by Liisa Holm and Chris Sander, 1989-1991.

No redistribution, no program changes, no commercial use.

For details, see J.Mol.Biol. 218, 183-194 (1991).

___________________________________________________________________________
*)
(* array dimensioning constants *)

const
   maxlen  = 30;
   maxres  = 10000;
   scaler  = 100;
   maxprot = 100;

(* string type declarations *)

type
   string = varying[80] of char;
   a1	  = char;
   a2	  = packed array[1..2] of char;
   a3	  = packed array[1..3] of char;
   a4	  = packed array[1..4] of char;
   a5	  = packed array[1..5] of char;
   a6	  = packed array[1..6] of char;
   a10	  = packed array[1..10] of char;

(* openfile types *)

inout = (infile, outfile);

(* CA distances *)

dist_array = array[0..maxlen] of integer;
   distancematrix = array[1..maxres] of dist_array;
   protlist_record = record protname  : string;
			start, length : integer;
		     end;	      
   protlist = array[1..maxprot] of protlist_record;

type 
   scorepointer	      = ^scorelist_record;
   scorepointersarray = array[1..maxres] of scorepointer;
   scorelist_record   = record
			   begresidue, matchresidue, matchlength,
			   matchprotlistindex	: integer;
			   matchprotname	: string;
			   score		: real;
			   prevscore, nextscore	: scorepointer;
			end;			
   slidemode = (fast, slow);
   gaparray = array[1..maxres] of boolean;

redundancyarray = array[1..maxres] of integer;

(* protein list type declarations *)

   xyz = array[1..3] of real;
   main_record = record n,ca,c,o: xyz; end;

   proteinpointer = ^protein_record;
   residuepointer = ^residue_record;
   atompointer = ^atom_record;
   
atom_record =
   record
      nextatom, prevatom : atompointer;
      coord		 : xyz;
      atomname		 : a4;
      icode		 : integer; (* free to use *)
   end;			 

backbone_record =
   record
      resname  : a4;
      resno    : a5;
      main     : main_record;
      pdbresid : a10;
   end;	       

residue_record =
   record
      nextresidue, prevresidue : residuepointer;
      backbone		       : backbone_record;
      firstside		       : atompointer;
   end;			       

protein_record = 
   record
      prevprotein, nextprotein : proteinpointer;
      protname		       : string;
      length		       : integer;
      extra		       : integer; (* user's fancy *)
      firstresidue	       : residuepointer;
   end;			       

dist_record = record d: dist_array; respo: residuepointer; end;
   distm = array[1..maxres] of dist_record;
 

const maxscores = 1000;

(* -------------------------------------------------------------------------- *)
type                                            (*u3b*)
   t1		     = array[1..4*maxres] of real;
   t2		     = array[1..4*MAXRES] of xyz;
   rotationmatrix    = array[1..3,1..3] of real;
   translationvector = array[1..3] of real;
     
(* -------------------------------------------------------------------------- *)


type 
   fragmentcoord = array[1..maxlen+1] of backbone_record;
type
   fragmentpointer   = ^fragment_record;
   jumppointer	     = ^jump_record;
   fragment_record   = record
			  rmssum, maxdelta			      : real;			   
			  sourceprotno, sourceresno, longeur, njumpto : integer;
			  prevfragment, prev, nextfragment	      : fragmentpointer;
			  main0,main1, main2			      : backbone_record;
			  firstjump				      : jumppointer;
		       end;					      
   jump_record	     = record
			  nextjump, prev	 : jumppointer;
			  jumpfromfragment	 : fragmentpointer;
			  biggestdelta, rmsdelta : real;
		       end;			 
   fragmentlistarray = array[1..maxres] of fragmentpointer;        

(* lookup table *)
const
   maxl        = 8; maxdist = 4000; (* distance unit is 0.01 A *)
type				 
   hitpointer  = ^hit_record;
		    hit_record  = record nexthit																		 : hitpointer; protno, resno: integer; end;
   lookuptable = array[1..maxl, 0..maxdist] of hit_record;


(* ----------------------------------------------------------------- *)

(*procedure u3b(w: t1; x, y: t2; n, mode: integer; var rms: real;
 var u:rotationmatrix; var t: translationvector; var ier: integer); 
 FORTRAN; *)

procedure x();
begin
end;

 procedure init_residuelist(var firstresidue, lastresidue: residuepointer);

begin
   firstresidue := nil;
   lastresidue := nil;
end;

 function val(s: a4): integer;
var
	i,v: integer;
	aa: packed array[1..11] of char;

begin
	aa:='0123456789 ';
V := 0;
	i:=0; repeat i:=i+1 until (aa[i]=s[1]) OR (I=11); IF I<11 THEN  v:=(I-1)*1000;
	i:=0; repeat i:=i+1 until (aa[i]=s[2]) OR (I=11); IF I<11 THEN  v:=v+(I-1)*100;
	i:=0; repeat i:=i+1 until (aa[i]=s[3]) OR (I=11); IF I<11 THEN  v:=v+(I-1)*10;
	i:=0; repeat i:=i+1 until (aa[i]=s[4]) OR (I=11); IF I<11 THEN  v:=v+(I-1);
	val:=v;
	end;

 function resnametocode(resnm: a4): integer;
var
	aa: packed array[1..150] of char;
	s: array[1..5] of packed array[1..30] of char;
	i,l,k: integer;

begin
	S[1] := 'ALAARGASNASPCYSGLUGLNGLYHISILE';
	S[2] := 'LEULYSMETPHEPROSERTHRTRPTYRVAL';
	S[3] := 'ASXGLXACDALBALIABUAROBASBETHSE';
	S[4] := 'HYPHYLORNPCASARTAUTHYUNKACEFOR';
	S[5] := 'CYHCSHCSSCYXILUPRZPROCPRTRYHOH';
	L := 0; FOR K := 1 TO 5 DO FOR I := 1 TO 30 DO
	BEGIN L:=L+1;AA[L]:=s[K,I];END;
	L:=0; K:=1; I:=1;
	WHILE (K<51) AND (L=0) DO
	BEGIN IF AA[I]=RESNM[1] THEN IF AA[I+1]=RESNM[2] 
		THEN IF AA[I+2]=RESNM[3] THEN L:=I; I:=I+3;K:=K+1;END;
	RESNAMETOCODE:=L;
	END;
	
 function atomnametocode(atnm: a4): integer;
var
	aa: packed arRay[1..210] of char;
	s: array[1..7] of packed array[1..30] of char;
	i,l,k: integer;
	X: A4;

begin
	S[1] := 'N  CA C  O  CB CG SG CG1CG2OG ';
	S[2] := 'OG1OG2CD CD1CD2ND ND1ND2OD1OD2';
	S[3] := 'SD CE NE NE1NE2OE1OE2CE1CE2CE3';
	S[4] := 'CZ NZ CZ2CZ3CZ1OH OXTNH1NH2CH2';
	s[5] := 'AE1AE2CH1DH D  HA HB HG HD DE ';
	S[6] := 'HD1HD2H  DG1DG2HG1HG2DH1DH2DD1';
	S[7] := 'DD2HE1HE2DE1DE2HE HZ DZ    ???';
	L := 0; FOR K := 1 TO 7 DO FOR I := 1 TO 30 DO
	BEGIN L:=L+1;AA[L]:=s[K,I];END;
	L:=0; K:=1; I:=1;
	WHILE (K<71) AND (L=0) DO
	BEGIN IF AA[I]=atnm[1] THEN IF AA[I+1]=atnm[2] 
		THEN IF AA[I+2]=atnm[3] THEN L:=I; I:=I+3;K:=K+1;END;
	if k=71 then begin writeln('unknown atomname *',atnm,'*'); l:=68; end;
	ATOMNAMETOCODE:=L;
	END;


 procedure save_backbone(outfilename: string; protein: proteinpointer);

var
	tapeout: file of backbone_record;
	curresidue: residuepointer;

begin
	open(tapeout, file_name := outfilename);
	rewrite(tapeout);
	curresidue := protein^.firstresidue;
	while curresidue <> nil do
	begin
		write(tapeout, curresidue^.backbone);
		curresidue := curresidue^.nextresidue;
		end;
	close(tapeout);
	end;

 procedure transrot(r: xyz; var rnew: xyz; u: rotationmatrix; t: translationvector);
var    
        i,j: integer;
     
begin
        (* rnew = u*r+t *)
        for i:=1 to 3 do
        begin
                rnew[i]:=t[i];
                for j:=1 to 3 do rnew[i]:=rnew[i]+u[j,i]*r[j];
                end;
        end;
     
 procedure matchca(residue1, residue2: residuepointer; length: integer;
		var u: rotationmatrix; var t: translationvector; var rms: real);
     
var     x,y: t2;                (* coordinates *)
var     i,j,ier: integer;
     	res1,res2: residuepointer;
var     w: t1;                  (* atom pair weight factors *)

begin
	(* 1 moves, 2 is static *)
        (* copy CA coordinates to onedimensional array *)
	res1 := residue1;
	res2 := residue2;
	for j := 1 to length do
	begin
		x[j] := res1^.backbone.main.ca;
		y[j] := res2^.backbone.main.ca;
		res1 := res1^.nextresidue;
		res2 := res2^.nextresidue;
		w[j] := 1;
                end;
     
        u3b(w,x,y,length,1,rms,u,t,ier);
     
        end;
     

 procedure write_brk(atomname: a4; x: xyz; pdbname: a10; 
	segid: a4; var brkout: text);

begin
	writeln(brkout,'ATOM         ',atomname,pdbname,
	'   ',x[1]:8:3,x[2]:8:3,x[3]:8:3,'                  ',segid);
end;

 procedure mainrecord_to_brk(main: main_record; 
	backbone: backbone_record; var tapeout: text);
 
begin
	write_brk('N   ',main.n,backbone.pdbresid,'    ',tapeout);
	write_brk('CA  ',main.ca,backbone.pdbresid,'    ',tapeout);
	write_brk('C   ',main.c,backbone.pdbresid,'    ',tapeout);
	write_brk('O   ',main.o,backbone.pdbresid,'    ',tapeout);
end;


 procedure oldwrite_brk(atomname: a4; x: xyz; resname: a4; chainid: a1;
			resno: a5; segid: a4; var brkout: text);

begin
	writeln(brkout,'ATOM         ',atomname,resname,chainid,
		resno,'   ',x[1]:8:3,x[2]:8:3,x[3]:8:3,'                  ',segid);
end;

 procedure oldmainrecord_to_brk(main: main_record; resname: a4;
	resno: a5; var tapeout: text);
 
begin
	oldwrite_brk('N   ',main.n,resname,' ',resno,'    ',tapeout);
	oldwrite_brk('CA  ',main.ca,resname,' ',resno,'    ',tapeout);
	oldwrite_brk('C   ',main.c,resname,' ',resno,'    ',tapeout);
	oldwrite_brk('O   ',main.o,resname,' ',resno,'    ',tapeout);
end;


 function str(v: integer): a4;
var
        i: integer;
        aa: packed array[1..11] of char;
        s: a4;

begin
        if v > 9999 then begin writeln('too large: ',v); HALT; end;
        aa := '0123456789 ';
        i := trunc(v/1000); s[1] := aa[i+1];
        v := v-i*1000; i := trunc(v/100); s[2] := aa[i+1];
        v := v-i*100; i:= trunc(v/10); s[3] := aa[i+1];
        v := v-i*10; s[4] := aa[v+1];
        i:=1; while (s[i]='0') and (i<4) do begin s[i] := ' ';i := i+1; end;
        str := s;
        end;
        
        
 function codetoresname(code: integer): a4;
var
        aa: packed array[1..150] of char;
        s: array[1..5] of packed array[1..30] of char;
        i,l,k: integer;
        X: A4;

begin
        S[1] := 'ALAARGASNASPCYSGLUGLNGLYHISILE';
        S[2] := 'LEULYSMETPHEPROSERTHRTRPTYRVAL';
        S[3] := 'ASXGLXACDALBALIABUAROBASBETHSE';
        S[4] := 'HYPHYLORNPCASARTAUTHYUNKACEFOR';
        S[5] := 'CYHCSHCSSCYXILUPRZPROCPRTRYHOH';
        L := 0; FOR K := 1 TO 5 DO FOR I := 1 TO 30 DO
        BEGIN L:=L+1;AA[L]:=s[K,I];END;
        FOR I:=1 TO 3 DO X[I]:=AA[CODE+I-1]; X[4]:=' ';
        CODETORESNAME:=X;
        END;
 

 function codetoatomname(code: integer): a4;
var
        aa: packed arRay[1..210] of char;
        s: array[1..7] of packed array[1..30] of char;
        i,l,k: integer;
        X: A4;

begin
        S[1] := 'N  CA C  O  CB CG SG CG1CG2OG ';
        S[2] := 'OG1OG2CD CD1CD2ND ND1ND2OD1OD2';
        S[3] := 'SD CE NE NE1NE2OE1OE2CE1CE2CE3';
        S[4] := 'CZ NZ CZ2CZ3CZ1OH OXTNH1NH2CH2';
        s[5] := 'AE1AE2CH1DH D  HA HB HG HD DE ';
        S[6] := 'HD1HD2H  DG1DG2HG1HG2DH1DH2DD1';
        S[7] := 'DD2HE1HE2DE1DE2HE HZ DZ    ???';
        L := 0; FOR K := 1 TO 7 DO FOR I := 1 TO 30 DO
        BEGIN L:=L+1;AA[L]:=s[K,I];END;
        FOR I:=1 TO 3 DO X[I]:=AA[CODE+I-1]; X[4]:=' ';
        codetoatomname:=X;
        END;

 function atan2(y,x: real): real;
var     z: real;
begin
        if x<>0.0 then z := arctan(y/x) else if y>0.0 then z := 1.570796
        else if y<0.0 then z := -1.570796 else z := 6.283185;
        if x<0.0 then if y>0.0 then z := z+3.141593 else z := z-3.141593;
        atan2 := z;
        end;

 procedure diff(x,y: xyz; var z: xyz);
begin z[1] := x[1]-y[1]; z[2] := x[2]-y[2]; z[3] := x[3]-y[3]; end;

 function dot(x,y: xyz): real;
begin dot := x[1]*y[1]+x[2]*y[2]+x[3]*y[3]; end;

 procedure cross(x,y: xyz; var z: xyz);
begin
        z[1] := x[2]*y[3]-y[2]*x[3];
        z[2] := x[3]*y[1]-y[3]*x[1];
        z[3] := x[1]*y[2]-y[1]*x[2];
        end;

 function dihedralangle(v1, v2, v3, v4: xyz): real;
var
        i: integer;
        u, v: real;
        v12, v43, x, y, z, p: xyz;

begin
        diff(v1, v2, v12); diff(v4, v3, v43); diff(v2,v3,z);
        cross(z, v12, p); cross(z, v43, x); cross(z, x, y);
        u := dot(x,x); v := dot(y, y);
        dihedralangle := 360.0;
        if (u>0.0) and (v>0.0) then
        begin
                u := dot(p,x)/sqrt(u); v := dot(p,y)/sqrt(v);
                if (u<>0.0) or (v<>0.0) then dihedralangle := atan2(v,u)*57.29578;
                end;
        end;

 procedure openfile(var filvar: text; filename: string; mode: inout);
var
	x: integer;     

begin
        case mode of
                infile: begin
                                open(filvar,file_name:=filename, history:= readonly, error:= continue);
                                reset(filvar); (**, error:= continue); **)
				x := status(filvar);
                		writeln('status opening ',filename,' = ',x);
		                if status(filvar)<>0 then
                                writeln('file ',filename,' not found. ')
                                else writeln('reading file ',filename);
                                end;
                outfile:begin
                                open(filvar,file_name:=filename);
                                rewrite(filvar);
                                writeln('writing to ',filename);
                                end;
                end;
        end;
     

(* make list - procedures *)

 procedure make_new_protein(var firstprotein: proteinpointer);

var
	newprotein: proteinpointer;

begin
	new(newprotein);
	(* insert to 1st *)
	newprotein^.nextprotein := firstprotein;
	if firstprotein <> nil then firstprotein^.prevprotein := newprotein;
	newprotein^.prevprotein := nil;
	firstprotein := newprotein;
	newprotein^.firstresidue := nil;
	end;

 procedure make_new_residue(protein: proteinpointer; var firstresidue, lastresidue: residuepointer); 

var
	newresidue: residuepointer;

begin
	new(newresidue);
	if firstresidue = nil then firstresidue := newresidue;
	if lastresidue <> nil then lastresidue^.nextresidue := newresidue;
	newresidue^.prevresidue := lastresidue;
	newresidue^.nextresidue := nil;
	newresidue^.firstside := nil;
	lastresidue := newresidue;
	if protein^.firstresidue = nil then protein^.firstresidue := newresidue;
	end;

 function findresnoresname(protein: proteinpointer; resno: a5; resname: a4): residuepointer;

var
	curresidue: residuepointer;
	found: boolean;

begin
	curresidue:= protein^.firstresidue;
if curresidue = nil then writeln(' firstresidue = nil !!! ');
	found := false;
	while (curresidue <> nil) and not found do 
	if (curresidue^.backbone.resno <> resno) or 
	   (curresidue^.backbone.resname <> resname) then curresidue := curresidue^.nextresidue else found := true;
	findresnoresname := curresidue;
	end;
		
 procedure makeside(residue: residuepointer; atnm: a4; cd: xyz); 
var
	lastatom, newatom: atompointer;

begin 
	new(newatom);
	lastatom := residue^.firstside;
	(* go to end of sidechain *)
	if residue^.firstside = nil then residue^.firstside := newatom
	else while lastatom^.nextatom <> nil do lastatom := lastatom^.nextatom;	
	if lastatom <> nil then lastatom^.nextatom := newatom;
	newatom^.nextatom := nil;
	newatom^.coord := cd;
	newatom^.atomname := atnm;
	end;

 procedure brk_to_list(lastprotein: proteinpointer; 
var firstresidue: residuepointer; infilename, headerfile: string; var length: integer);
(* NEW modification: header file, but missing from this - discovered 11-sep-89 *) 
var  
	header,skip6: a6;   
	atnm, resnm, top: a4;
	resno: a5;
	oldresno: a5;
     cd: xyz;
     d,d1,d2: a1;
     i: integer;
	tapein: text;
 	lastresidue: residuepointer;
	pdbname: a10;

(* Brookhaven format: a6, i5, 1x, a4, a1, a3, 1x, a1, i4, a1, 3x, 3f8.3
 header=ATOM  ,atom serial number, atom name, alternate location indicator for
 atoms, residue name, chain identifier, residue sequence number, code for insertions
 of residue, x, y, z *)

begin
	(* update proteinlist *)
	
	openfile(tapein, infilename, infile);
	init_residuelist(firstresidue, lastresidue);
	length := 0;

	(* read brk file *)
	oldresno :='*****';
	repeat 
		read(tapein, header);
		if header<>'ATOM  ' then readln(tapein);
		until header='ATOM  ';
	while (not eof(tapein)) and (header='ATOM  ') do
	begin 
		readln(tapein, skip6, d, atnm, resnm, d1, resno, cd[1], cd[2], cd[3]);
		if oldresno <> resno then begin oldresno:=resno;  {
		if (findresnoresname(lastprotein,resno,resnm) = nil) then
		begin }
(*if lastresidue<>nil then writeln('new residue: ',atnm,resnm,resno,d,'****',
lastresidue^.backbone.resname,lastresidue^.backbone.resno);
*)			make_new_residue(lastprotein,firstresidue, lastresidue);
			lastresidue^.backbone.resno := resno;
			lastresidue^.backbone.resname := resnm;
			pdbname[1] := resnm[1];
			pdbname[2] := resnm[2];
			pdbname[3] := resnm[3];
			pdbname[4] := resnm[4];
			pdbname[5] := d1;
			pdbname[6] := resno[1];
			pdbname[7] := resno[2];
			pdbname[8] := resno[3];
			pdbname[9] := resno[4];
			pdbname[10] := resno[5];
			lastresidue^.backbone.pdbresid := pdbname;
			length := length + 1;
			end;
		
		with lastresidue^.backbone.main do
		begin
			if atnm = 'CA  ' then ca := cd else
			if atnm = 'C   ' then c  := cd else
			if atnm = 'N   ' then n  := cd else
			if atnm = 'O   ' then o  := cd else
			makeside(lastresidue, atnm, cd);
			end;
		if not eof(tapein) then read(tapein, header);
		end;
	close(tapein);
	end;

 procedure list_to_side(protein: proteinpointer; outfilename: string);
var
	tapeout: file of integer;
	curresidue: residuepointer;
	currside: atompointer;
	atoms: array[1..100,0..3] of integer;
	i,j,n: integer;
	resn: a4;
const	superscale = 1000;

begin
	open(tapeout, file_name := outfilename);
	rewrite(tapeout);
	write(tapeout, protein^.length);
	curresidue := protein^.firstresidue;
	while curresidue <> nil do 
	begin
		with curresidue^.backbone do 
		begin
			resn[1]:=resno[1];resn[2]:=resno[2];
			resn[3]:=resno[3];resn[4]:=resno[4];
			write(tapeout, resnametocode(resname), val(resn));
		end;
		(* copy coordinates to atoms *)
		n := 0; (* atomcounter *)
		currside:=curresidue^.firstside;
		while currside<>nil do
		begin
			n:=n+1;
if n>100 then begin writeln('n>100'); HALT; end;
			for i:=1 to 3 do atoms[n,i]:=round(superscale*currside^.coord[i]);
			atoms[n,0]:=atomnametocode(currside^.atomname);
			if atoms[n,0] = 68 then n := n-1; (* unknown atomname *)
			currside := currside^.nextatom;
			end;
		write(tapeout, n);
		for i:=1 to n do for j:=0 to 3 do write(tapeout, atoms[i,j]);
		curresidue:=curresidue^.nextresidue;
		end;
	end;

 function distance(a,b: xyz): real;
var
	c: xyz;
	d: real;
	i: integer;
label warn, ok;

begin
	diff(a,b,c);
	for i := 1 to 3 do if c[i] > 1e10 then goto warn;
	for i := 1 to 3 do if c[i] < -1e10 then goto warn;
	d := (c[1]*c[1]+c[2]*c[2]+c[3]*c[3]);
	distance := sqrt(d);
	goto ok;
warn:
(*	writeln('WARNING: distance out of range ',c[1], c[2], c[3]);
*)	distance := 99999.999;
ok:	end;

 procedure backbone_to_list(lastprotein: proteinpointer; var firstresidue: residuepointer; infilename: string; var length: integer);

var
	tapein: file of backbone_record;
	lastresidue: residuepointer;

begin
	open(tapein, file_name := infilename, history := readonly);
	reset(tapein);
	firstresidue := nil;
	lastresidue := nil;
	length := 0;
	while not eof(tapein) do
	begin
		make_new_residue(lastprotein, firstresidue, lastresidue);
		length := length + 1;
		read(tapein, lastresidue^.backbone);
		end;
	end;

 procedure fill_protein_record(protein: proteinpointer; first: residuepointer; proteinname: string; l: integer);

begin
	protein^.protname := proteinname;
	protein^.length := l;
	protein^.firstresidue := first;
	end;


 procedure side_to_list(protein: proteinpointer; infilename: string);

var
	tapein: file of integer;
	curresidue: residuepointer;
	i, j, n, code, a, b, c, length: integer;
	cd: xyz;
	atnm, resname: a4;
	resno: a5;

const	superscaler = 1000;

begin
	open(tapein, file_name := infilename, history := readonly);
	reset(tapein);
	(* assume .MAIN has been read & protein-list initialized *)
	read(tapein, length);
	for n := 1 to length do (* residue-counter *)
	begin
		read(tapein, code);
		resname := codetoresname(code);
		read(tapein, code);
		resno := str(code);
		curresidue := findresnoresname(protein, resno, resname);
		if curresidue = nil then HALT;
		read(tapein, i); 
		for j := 1 to i do (* atom-counter *)
		begin
			read(tapein, code, a, b, c);
			atnm := codetoatomname(code);
			cd[1] := a/superscaler; 
			cd[2] := b/superscaler;
			cd[3] := c/superscaler;
			makeside(curresidue, atnm, cd);
			end;
		end;
	end;


 function findproteinpointer(master: string; firstprotein: proteinpointer): proteinpointer;
var	prot: proteinpointer;
	found: boolean;

begin
	prot := firstprotein;
	found := false;
	while not found and (prot<>nil) do if prot^.protname<>master then
		prot := prot^.nextprotein else found := true;	
	findproteinpointer := prot;
	end;

 procedure getprotein(master: string; sidechain: boolean; 
			var firstprotein, protein: proteinpointer);
	(* return pointer to protein master - read coordinates if it is not in list *)
var
	firstresidue: residuepointer;
	l: integer;

begin
	protein := findproteinpointer(master, firstprotein);
	if protein = nil then
     	begin
		make_new_protein(firstprotein);
		backbone_to_list(firstprotein, firstresidue, master+'.main',l);
		if sidechain then side_to_list(firstprotein, master+'.side');
		fill_protein_record(firstprotein, firstresidue, master, l);
		protein := firstprotein;
		end;
	end;

 procedure fill_lookup(var lookup: lookuptable; var dist: distm; nprot: integer;
	var proteins: protlist);
(* tabulate i+2, i+3, .., i+maxl distances *)
var
	i,j,d,l,q,s: integer;
	newhit: hitpointer;
begin
	(* init lookup table *)
	for d := 0 to maxdist do for l := 2 to maxl do with lookup[l,d] do nexthit := nil; 
	(* loop through dist[] *)
	for i := 2 to nprot do with proteins[i] do (* i = protein *)
	for l := 2 to maxl do for j := 1 to length-l do
	begin
		s := start+j;
		q := dist[s].d[l];
		if q<maxdist then with lookup[l,q] do
		begin
			new(newhit);
			newhit^.nexthit := nexthit; 
			newhit^.protno := i;
			newhit^.resno := j;
			nexthit := newhit;
			end;
		end;
	end;

 procedure maintodist(master: string);

var
	tapeout: file of DIST_ARRAY;
	firstprotein, prot1: proteinpointer;
	curresidue,firstresidue: residuepointer;
	i,j,l: integer;
	TEMPDIST: DIST_ARRAY;
        x,y,z: array[1..maxres] of real;
	r: array[1..4] of xyz;
	k: integer;

begin
	firstprotein := nil;
        if master<>'' then
        begin
		getprotein(master, false, firstprotein, prot1);
		l := prot1^.length;
		(* copy ca coordinates to array *)
		curresidue := prot1^.firstresidue;
		for i:=1 to l do with curresidue^.backbone.main do
		begin 
			x[i] := ca[1]; y[i] := ca[2]; z[i] := ca[3]; 
			curresidue := curresidue^.nextresidue;
			end;
		(* write protlength, calculated distances to tapeout *)
		(* dihedralangle to tempdist[0] *)
	        open(tapeout, file_name:=master+'.dist');
		rewrite(tapeout);
		tempdist[0] := 0;
		TEMPDIST[1]:=L; 
		FOR I := 2 TO MAXLEN DO TEMPDIST[I] := 0;
		write(tapeout, TEMPDIST);
		(* CALCULATE DISTANCES *)
		for i := 1 to l-1 do 
		BEGIN	for j := 1 to maxlen do
			if i+j <= l then TEMPDIST[J] := 
			round(scaler*(sqrt(sqr(x[i]-x[i+j])+sqr(y[i]-y[i+j])+sqr(z[i]-z[i+j])))) 
			ELSE TEMPDIST[J]:=0;
			if (tempdist[1] > 500) OR (tempdist[1] < 200)
			then writeln('WARNING: weird distance ',i:5,
					tempdist[1]:10, x[i]:6:1,y[i]:6:1,z[i]:6:1,x[i+1]:6:1,y[i+1]:6:1,z[i+1]:6:1);
			if i <= l-3 then 
			begin
				for k := 1 to 4 do
				begin
					r[k,1] := x[i+k-1];
					r[k,2] := y[i+k-1];
					r[k,3] := z[i+k-1];
					end;
				tempdist[0] := round(dihedralangle(r[1], r[2], r[3], r[4]));
				end
			else tempdist[0] := 0;
			WRITE(TAPEOUT, TEMPDIST);
			END;
		close(tapeout);
		end;
	end;
          
 procedure load_distances(var nprot: integer; var master: string;
          var proteins: protlist; var dist: distm;
          var firstprotein, prot1: proteinpointer);

var
   where,i,j: integer;
   prot: string;
   tapein: text;

	procedure readdists(master: string);

	var
		prot: proteinpointer;
		curresidue: residuepointer;
		tapein2: file of dist_array;
		tapein1: file of backbone_record;
		i: integer;
		temp: dist_array;
	label skip, l1;
	
	begin
(*		open(tapein1, file_name := master+'.main', history := readonly,
		error := continue); 
*)		open(tapein1, file_name := master+'.main', history := readonly); 
		if status(tapein1)<>0 then
		begin writeln(master,'.main not found.  Skipping protein. Status= ',status(tapein1));
			goto skip; end;
		reset(tapein1, error := continue);
		getprotein(master, false, firstprotein, prot);
		nprot := nprot+1;
l1:		open(tapein2, file_name := master+'.dist', history := readonly,
			error := continue);
		if status(tapein2)<>0 then (* create missing .dist file *)
		begin
			writeln('creating .dist file for ',master,'status=',status(tapein2));
			maintodist(master);
			goto l1;
			end;
		reset(tapein2, error := continue);
		curresidue := prot^.firstresidue;
		read(tapein2, temp);
		with proteins[nprot] do
		begin
			length := temp[1];
			if length <> prot^.length then
writeln('WARNING: ',master,'distance and coordinate files have different no of residues !!!');
			start := where; protname := master;
			for i := 1 to length-1 do
			begin
				where := where+1;
				read(tapein2, dist[where].d);
				dist[where].respo := curresidue;
				curresidue := curresidue^.nextresidue;
				end;
			where := where+1; dist[where].respo := curresidue;
			end;
		close(tapein2);
		close(tapein1);

        (* check chain continuity *)
        curresidue := prot^.firstresidue; 
	i := proteins[nprot].start;
        while curresidue <> nil do
        begin
             i := i+1;
             with curresidue^.backbone.main do if ca[1] = 0 then if ca[2] = 0 then if ca[3] = 0 then
             begin writeln('WARNING: missing CA coordinates for residue ',i,curresidue^.backbone.resno, 
curresidue^.backbone.resname); end;
		if i < prot^.length then if (dist[i].d[1] < 200) or (dist[i].d[1] > 480) then
             begin writeln('LOOK OUT !  C(alpha)-C(alpha-next) distance out of range for residue ',
i,curresidue^.backbone.resno, curresidue^.backbone.resname,dist[i].d[1]); end;
             curresidue := curresidue^.nextresidue;
             end;
skip:
		end;
		
begin
        nprot := 0;
        where := 0;
        (* read test-protein to 'slot 1' *)
        writeln('enter name of testprotein (e.g. 1sn3) ');
        readln(master);
        readdists(master);
        writeln(where,' residues in test-protein');
        firstprotein := nil;
        getprotein(master, false, firstprotein, prot1);
        (* read database distances *)
        prot := master;
                writeln('enter name of database protein names file (e.g. dglp.list) ');
                readln(prot);
                openfile(tapein, prot, infile);
                while not eof(tapein) do
                begin
                        readln(tapein, prot);
                  	writeln('loading ',prot,' to ',where,nprot+1);
		        readdists(prot);
                        end;
                close(tapein);
        writeln(nprot-1,' proteins read from distance library');
        writeln(where-proteins[1].length,' residues in distance library');
     end;

 procedure extend_to_rigth(j, begresidue, testlen, length, cutoff1, comparestart, templatestart: integer;
			  var dist: distm;
			  var nright: integer;
			  var s: real);
var
	diverge: boolean;
	delta, k, l: integer;
	sc: real;

begin
	diverge := false;
	k := 0;
	repeat (* eteen *)
		k := k+1;
		sc := 0;		
		for l := 0 to k-1 do 
		begin
			delta := dist[j+l].d[k-l] - dist[begresidue+l].d[k-l];
			if ABS(DELTA) > cutoff1 then diverge := true;
			sc := sc+sqr(delta);
			end;
		if not diverge then s := s+sc;
		until diverge or(k=testlen-1) or (k=maxlen-2) (* tilaa nleftille! *)
		or (j-comparestart+k = length) or (begresidue-templatestart+k = testlen);
	if diverge then nright := k else nright := k+1;  (* begresidue lasketaan mukaan *)
	end;

 procedure extend_to_left(j, begresidue, testlen, length, cutoff1, comparestart, templatestart: integer;
			  var dist: distm;
			  var nleft: integer;
			  var s: real;
			  nright: integer);

var
	diverge: boolean;
	delta, k, l, m: integer;
	sc: real;

begin
	(* taakse *)
	diverge := false;
	k := 0;
	l := nright;
	repeat
		k := k+1;
		sc := 0;
		for m := 1 to l+k-1 do if  (m<=maxlen) and (j-comparestart-k>0)
 and (begresidue-templatestart-k>0) then
		begin
			delta := dist[j-k].d[m] - dist[begresidue-k].d[m];
			if ABS(delta)>cutoff1 then diverge := true;
			sc := sc+sqr(delta);
			end else diverge := true;
		if not diverge then s := s+sc;
		until (j-comparestart-k=1) or (begresidue-templatestart-k=1) or (l+k-1=testlen) or
			 diverge or (l+k=maxlen);
	if diverge then nleft := k-1 else nleft := k;
	end;


 procedure update_pathsums(start, goal: integer; var fragmentlist: fragmentlistarray);
var
	i: integer;
	smallest, t: real;
	curfragment: fragmentpointer;
	curjump, remjump: jumppointer;

begin
	for i := start+1 to goal do
	begin
		curfragment := fragmentlist[i];
		while curfragment <> nil do
		begin
			smallest := 999999;
			curjump := curfragment^.firstjump;
			while curjump <> nil do
			begin
				if i=start+1 then t := curjump^.rmsdelta else
				t := curjump^.rmsdelta+curjump^.jumpfromfragment^.rmssum;
				if t <= smallest then
				begin
					smallest := t;
					remjump := curjump;
					end;
				curjump := curjump^.nextjump;
				end;
			curfragment^.prevfragment := remjump^.jumpfromfragment;
			curfragment^.rmssum := smallest;
			curfragment^.maxdelta := remjump^.biggestdelta;
			curfragment := curfragment^.nextfragment;
			end;
		end;
	end;

 procedure jump(jumpfrom, jumpto: fragmentpointer; var rms, suurin: real);
var
	delta: array[1..5] of real;
	sum: real;
	i: integer;

begin
	delta[1] := distance(jumpfrom^.main2.main.ca, jumpto^.main1.main.ca);
	delta[2] := distance(jumpfrom^.main2.main.n, jumpto^.main1.main.n);
	delta[3] := distance(jumpfrom^.main1.main.c, jumpto^.main0.main.c);
	delta[4] := distance(jumpfrom^.main1.main.ca, jumpto^.main0.main.ca);
	delta[5] := distance(jumpfrom^.main1.main.o, jumpto^.main0.main.o);
(*	delta[5] := 0.5 * distance(jumpfrom^.main1.main.o, jumpto^.main0.main.o);
*)	suurin := delta[1]; for i := 2 to 5 do if delta[i]>suurin then suurin := delta[i];
	sum := 0; for i := 1 to 5 do sum := sum+sqr(delta[i]);
	rms := sqrt(sum/4);
	end;

 procedure update_jump_fromnewtoold(var fragmentlist: fragmentlistarray; 
				k, begresidue, testlen: integer;
				var newfragment: fragmentpointer; tolerance: real);
var
	jumpto: fragmentpointer;
	newjump, lastjump: jumppointer;
	rms, suurin: real;

begin
	if k+begresidue <= testlen then
	begin
		jumpto := fragmentlist[k+begresidue];
		while jumpto <> nil do
		begin
			jump(newfragment, jumpto, rms, suurin);
if suurin <= tolerance then
		begin
			new(newjump);
			newjump^.nextjump := nil;
			newjump^.jumpfromfragment := newfragment;
			(* newjump^.*delta *)
			newjump^.rmsdelta := rms;
			newjump^.biggestdelta := suurin;
			newfragment^.njumpto := newfragment^.njumpto+1;
			lastjump := jumpto^.firstjump;
			newjump^.nextjump := lastjump;
			if lastjump <> nil then lastjump^.prev := newjump;
			newjump^.prev := nil;
			jumpto^.firstjump := newjump;
			end;
			jumpto := jumpto^.nextfragment;
			end;
		end;
	end;

 procedure update_jump_fromoldtonew(var fragmentlist: fragmentlistarray; 
				k, begresidue, testlen: integer;
				var newfragment: fragmentpointer; tolerance: real);
var
	jumpfrom: fragmentpointer;
	newjump, lastjump: jumppointer;
	rms, suurin: real;

begin
	if k+begresidue > 2 then
	begin
		jumpfrom := fragmentlist[k+begresidue-2];
		while jumpfrom <> nil do
		begin
			jump(jumpfrom, newfragment, rms, suurin);
if suurin <= tolerance then
		begin
			new(newjump);
			newjump^.nextjump := nil;
			newjump^.jumpfromfragment := jumpfrom;
			newjump^.rmsdelta := rms;
			newjump^.biggestdelta := suurin;
			lastjump := newfragment^.firstjump;
			newjump^.nextjump := lastjump;
			if lastjump <> nil then lastjump^.prev := newjump;
			newjump^.prev := nil;
			jumpfrom^.njumpto := jumpfrom^.njumpto+1;
			newfragment^.firstjump := newjump;
			end;
			jumpfrom := jumpfrom^.nextfragment;
			end;
		end;
	end;

 procedure init_fragmentlist(var fragmentlist: fragmentlistarray; len: integer);
var 	i: integer;
begin	for i := 1 to len do fragmentlist[i] := nil; 
	end;

 procedure generate_cb(n, ca, c0: xyz; var cb: xyz);

var	i: integer;
	a,b,c,d: xyz;
	lc,ld: real;

begin
	(* a = ca to n; b = ca to c; c = right-vector; d = up-vector; cb = ca + 1.32 * right-vector + 0.765 * up-vector *)
	diff(n,ca,a); diff(c0,ca,b);
	cross(a,b,c); (* c = a x b *)
	lc := 0; ld := 0;
	for i := 1 to 3 do
	begin
		lc := lc+sqr(c[i]);
		d[i] := a[i]+b[i];
		ld := ld+sqr(d[i]);
		end;
	for i := 1 to 3 do begin c[i] := c[i]/sqrt(lc); d[i] := -d[i]/sqrt(ld);end;
	for i := 1 to 3 do cb[i] := ca[i]+1.32*c[i]+0.765*d[i];
	end;


 procedure generate_sidechains(ca,cb,c,n: xyz; residue: residuepointer;
	 var firstprotein: proteinpointer);

(* add sidechain atoms to residuepointer^.firstside^... *)
var
	prot1: proteinpointer;
	w: t1;
	x,y: t2;
	ier,i: integer;
	rms: real;
	u: rotationmatrix;
	t: translationvector;
	curatom: atompointer;
	cbeta, rnew: xyz;
const   datadir = 'csdata:';
label	quit;

begin
	if (residue^.backbone.resname='GLY ') then goto quit;
	getprotein(datadir+residue^.backbone.resname, true, firstprotein, prot1);
	(* get transrot matrices for n-ca-cb-c *)
	with prot1^.firstresidue^.backbone.main do
	begin
		generate_cb(n, ca, c, cbeta);
		x[1] := n; x[2] := ca; x[3] := cbeta; x[4] := c;
		end;
	y[1] := n; y[2] := ca; y[3] := cb; y[4] := c;
	w[1] := 1; w[2] := 10; w[3] := 10; w[4] := 1;
	u3b(w,x,y,4,1,rms,u,t,ier);
	(* copy transrotated sidechain atoms to residue^.firstside... *)
	curatom := prot1^.firstresidue^.firstside;
	while curatom <> nil do
	begin
		transrot(curatom^.coord, rnew, u, t);
		makeside(residue, curatom^.atomname, rnew);
		curatom := curatom^.nextatom;
		end;
quit:
	end;

 procedure add_sidechains(protein: proteinpointer; var firstprotein: proteinpointer);
var
	curresidue: residuepointer;
	cb: xyz;

begin
	curresidue := protein^.firstresidue;
	while curresidue <> nil do
	begin
		with curresidue^.backbone do 
		begin
			if curresidue^.firstside=nil then generate_cb(main.n,main.ca,main.c,cb) else cb := curresidue^.firstside^.coord;
			generate_sidechains(main.ca,cb,main.c,main.n,curresidue, firstprotein);
			end;
		curresidue := curresidue^.nextresidue;
		end;
	end;

 procedure read_brk;
(* read backbone file to list, do list-to-brk *)

var
	firstprotein: proteinpointer;
	firstresidue: residuepointer;
	readdir, writedir, master,ext: string;
	protlength,i,j: integer;

begin
	writeln('enter readdirectory (e.g. cs:) '); readln(readdir);
	writeln('enter writedirectory (e.g. work$disk:) '); readln(writedir);
	firstprotein := nil;
	repeat
		write('enter filename with extension [.brk] (e.g. 1sn3.brk_mod) ');
		readln(master);
		if master <> '' then
		begin
			(* parse filename > extension *)
			j :=length(master)+1;
			for i :=1 to length(master) do if master[i]='.' then j :=i;
			if j<=length(master) then ext := substr(master, j, length(master)-j+1) 
				else ext:='.brk'; 
			master := substr(master, 1, j-1);
			make_new_protein(firstprotein);
			brk_to_list(firstprotein, firstresidue, 
				readdir+master+ext,writedir+master+'.brk_hdr',
				protlength);
			fill_protein_record(firstprotein,firstresidue,
				master, protlength);
			save_backbone(writedir+master+'.main', firstprotein);
			list_to_side(firstprotein,writedir+master+'.side');
			end;
		until master='';
	end;

end. (* module *)




begin 
       (* print MaxSprout header *)
writeln('__________________________________________________________________________ ');
writeln('                                                                           ');
writeln('This set of programs constructs full atomic coordinates of a protein       ');
writeln('from a given C(alpha) trace and optimizes side chain geometry.             ');
writeln('                                                                           ');
writeln('Copyright by Liisa Holm and Chris Sander, 1989-1991.                       ');
writeln('                                                                           ');
writeln('No redistribution, no program changes, no commercial use.                  ');
writeln('                                                                           ');
writeln('For details, see J.Mol.Biol. 218, 183-194 (1991).                          ');
writeln('                                                                           ');
writeln('___________________________________________________________________________');
	read_brk;
end.

