module maxd(input, output);
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
       maxlen = 30;
       maxres = 10000;
       scaler = 100;
       maxprot = 100;

(* string type declarations *)

type
        string = varying[80] of char;
        a1 = char;
        a2 = packed array[1..2] of char;
        a3 = packed array[1..3] of char;
        a4 = packed array[1..4] of char;
        a5 = packed array[1..5] of char;
        a6 = packed array[1..6] of char;
        a10= packed array[1..10] of char;

(* openfile types *)

inout = (infile, outfile);

(* CA distances *)

dist_array = array[0..maxlen] of integer;
distancematrix = array[1..maxres] of dist_array;
protlist_record = record protname: string;
                        start, length: integer;
                        end;
protlist = array[1..maxprot] of protlist_record;

type 
        scorepointer = ^score_record;
        score_record = record
                                begresidue, matchresidue, matchlength,
                                matchprotlistindex: integer;
                                matchprotname: string;
                                score, carms: real;
                                prevscore, nextscore: scorepointer;
                                end;
gaparray = array[1..maxres] of boolean;

(* protein list type declarations *)

xyz = array[1..3] of real;
main_record = record n,ca,c,o: xyz; end;

proteinpointer = ^protein_record;
residuepointer = ^residue_record;
atompointer = ^atom_record;

atom_record =
        record
                nextatom, prevatom: atompointer;
                coord: xyz;
                atomname: a4;
                icode: integer; (* free to use *)
                end;

backbone_record =
        record
                resname: a4;
		resno: a5;
                main: main_record;
		pdbresid: a10;
                end;

residue_record =
        record
                nextresidue, prevresidue: residuepointer;
                backbone: backbone_record;
                firstside: atompointer;
                end;

protein_record = 
        record
                prevprotein, nextprotein: proteinpointer;
                protname: string;
                length: integer;
                extra: integer; (* user's fancy *)
                firstresidue: residuepointer;
                end;

dist_record = record d: dist_array; respo: residuepointer; end;
distm = array[1..maxres] of dist_record;
 
(* -------------------------------------------------------------------------- *)
type                                            (*u3b*)
        t1 = array[1..4*maxres] of real;
        t2 = array[1..4*MAXRES] of xyz;
        rotationmatrix = array[1..3,1..3] of real;
        translationvector = array[1..3] of real;
     
(* -------------------------------------------------------------------------- *)


type 
        fragmentcoord = array[1..maxlen+1] of backbone_record;
type
        fragmentpointer = ^fragment_record;
        jumppointer = ^jump_record;
        fragment_record = record
                                rmssum, maxdelta: real;
                                sourceprotno, sourceresno, longeur, njumpto: integer;
                                prevfragment, prev, nextfragment: fragmentpointer;
                                main0,main1, main2: backbone_record;
                                firstjump: jumppointer;
                                end;
        jump_record = record
                                nextjump, prev: jumppointer;
                                jumpfromfragment: fragmentpointer;
                                biggestdelta, rmsdelta: real;
                                end;
        fragmentlistarray = array[1..maxres] of fragmentpointer;        

(* lookup table *)
const
	maxl = 8; maxdist = 4000; (* distance unit is 0.01 A *)
type
	hitpointer = ^hit_record;
	hit_record = record nexthit: hitpointer; protno, resno: integer; end;
	lookuptable = array[1..maxl, 0..maxdist] of hit_record;

(* ----------------------------------------------------------------- *)
(* external procedures *)
procedure u3b(w:t1;x,y:t2;n,mode: integer;var rms:real;
                var u:rotationmatrix;var t:translationvector;var ier:integer); FORTRAN;
(* ----------------------------------------------------------------- *)
function str(i: integer): a4; external;
procedure init_fragmentlist(var fragmentlist: fragmentlistarray; len: integer); external;
procedure make_new_protein(var firstprotein: proteinpointer); external;
procedure make_new_residue(protein: proteinpointer; 
          var firstresidue, lastresidue: residuepointer);  external;
procedure backbone_to_list(lastprotein: proteinpointer; 
          var firstresidue: residuepointer; infilename: string; 
          var length: integer); external;
procedure fill_protein_record(protein: proteinpointer; 
          first: residuepointer; proteinname: string; l: integer); external;
procedure generate_cb(n, ca, c0: xyz; var cb: xyz); external;
procedure add_sidechains(protein: proteinpointer; 
          var firstprotein: proteinpointer); external;
procedure makeside(residue: residuepointer; atnm: a4; cd: xyz); external;
procedure mainrecord_to_brk(main: main_record; backbone: backbone_record;
        var tapeout: text); external;
procedure oldmainrecord_to_brk(main: main_record; resname: a4; resno: a5;
        var tapeout: text); external;
procedure write_brk(atomname: a4; x: xyz; pdbname: a10; title: a4; 
	var brkout: text); external;
procedure oldwrite_brk(atomname: a4; x: xyz; resname: a4; chainid: a1;
                        resno: a5; title: a4; var brkout: text); external;
procedure matchca(residue1, residue2: residuepointer; length: integer;
                var u: rotationmatrix; var t: translationvector; var rms: real); external;
procedure transrot(r: xyz; var rnew: xyz; u: rotationmatrix; t: translationvector); external;
procedure update_jump_fromoldtonew(var fragmentlist: fragmentlistarray;
          k, begresidue, testlen: integer; 
          var newfragment: fragmentpointer; tolerance: real); external;
procedure update_jump_fromnewtoold(var fragmentlist: fragmentlistarray;
          k, begresidue, testlen: integer; 
          var newfragment: fragmentpointer; tolerance: real); external;
procedure extend_to_rigth(j, begresidue, testlen, length, cutoff1, comparestart, templatestart: integer;
                          var dist: distm;
                          var nright: integer;
                          var s: real); external;
procedure extend_to_left(j, begresidue, testlen, length, cutoff1, comparestart, templatestart: integer;
                          var dist: distm;
                          var nleft: integer;
                          var s: real;
                          nright: integer); external;
procedure update_pathsums(start, goal: integer; var fragmentlist: fragmentlistarray); external;
function distance(a,b: xyz): real; external;
procedure load_distances(var nprot: integer; var master: string; 
          var proteins: protlist; var dist: distm; 
          var firstprotein, prot1: proteinpointer); external;
procedure fill_lookup(var lookup: lookuptable; var dist: distm; 
      nprot: integer; var proteins: protlist); external;
procedure openfile(var filvar: text; filename: string; mode: inout); external;
(* ----------------------------------------------------------------- *)
[global] procedure build;
const
	maxmaxscores = 50;
	maxprotlength = 1000;
type
	candarray = array[1..maxmaxscores] of score_record;
	scorearray = array[1..maxprotlength] of candarray;
	ilist = array[1..maxres] of integer;
var
	l1,c1,l2,c2,i, cs, maxscores, minlooplen: integer;
	dist: distm;
	firstprotein, prot1: proteinpointer;
	firstresidue: residuepointer;
	gap, filled: gaparray;
	nscore, fragindex, redundancy: array[1..maxres] of integer;
	candidates: scorearray;
	chainstart, chainend, leftzone, rightzone: integer;
	tolerance: real;
	c: char;
	master, statfilename, coorfilename, logfilename: string;
	statfile, coorfile, logfile, fragfile: text; 
	fragmentlist: fragmentlistarray;
	proteins: protlist;
	testlen, goal, nprot, startti: integer;
	lookup: lookuptable;
	looplen, cutoff1: integer;
	cb_only: boolean;

procedure update_gap(chainstart, start, goal: integer);
(* 9-Jan-89 build-range modification: chainstart *)
var
        curfragment: fragmentpointer;
        i: integer;
	longest: integer;
	curjump: jumppointer;
	
begin
	if (start>0) then
	for i := start to goal do if fragmentlist[i] <> nil then
	begin
		curfragment := fragmentlist[i];
		if i=chainstart then curfragment^.longeur := chainstart;
		if i>chainstart then
		while curfragment <> nil do
		begin
			longest := 0;
			curjump := curfragment^.firstjump;
			while curjump <> nil do
			begin
				if curjump^.jumpfromfragment^.longeur > longest
				then longest := curjump^.jumpfromfragment^.longeur;
				curjump := curjump^.nextjump;
			end;
			curfragment^.longeur := longest+1;
			curfragment := curfragment^.nextfragment;
		end;
	end;
	if start>0 then for i := start to goal do
	begin
		curfragment := fragmentlist[i];
		if curfragment <> nil then gap[i] := (curfragment^.longeur <> i) else gap[i] := true;
		while (curfragment<>nil) AND gap[i] do
		if curfragment^.longeur <> i then curfragment := curfragment^.nextfragment else
		begin
			gap[i] := false;
			if curfragment^.prev <> nil then curfragment^.prev^.nextfragment := curfragment^.nextfragment;
			if curfragment^.nextfragment <> nil then curfragment^.nextfragment^.prev := curfragment^.prev;
			fragmentlist[i]^.prev := curfragment;
			curfragment^.nextfragment := fragmentlist[i];
			fragmentlist[i] := curfragment;
			curfragment^.prev := nil;
		end;
	end;
end;


procedure path(start, goal: integer);

var
	i: integer;
	smallest, t: real;
	goalfragment, curfragment: fragmentpointer;
	curjump, remjump: jumppointer;
	trace: fragmentlistarray;
	prot1, firstprotein: proteinpointer;
	curresidue: residuepointer;
	curratom: atompointer;
	cb: xyz;

begin (* path *)
	firstprotein := nil;
	init_fragmentlist(trace, goal+1);
	update_pathsums(start, goal, fragmentlist);
	(* take smallest pathsum *)
	goalfragment := fragmentlist[goal];
	smallest := goalfragment^.rmssum;
	while goalfragment <> nil do
	begin
		if goalfragment^.rmssum<=smallest then
		begin
			curfragment := goalfragment;
			smallest := goalfragment^.rmssum;
		end;
		goalfragment := goalfragment^.nextfragment;
	end;
	(* copy temporarily to trace *)
	trace[goal] := curfragment;
	for i := goal-1 downto start do trace[i] := trace[i+1]^.prevfragment;
	(* make a copy of the C(alpha) trace *)
	prot1 := nil;
	make_new_protein(prot1);
	backbone_to_list(prot1, firstresidue, master+'.main', i);
	fill_protein_record(prot1, firstresidue, master, i);
	curresidue := prot1^.firstresidue;
	for i := 1 to start-1 do curresidue := curresidue^.nextresidue;
	for i := start to goal do
	begin
		curresidue^.backbone.main := trace[i]^.main1.main;
		curresidue := curresidue^.nextresidue;
	end;
	if goal<prot1^.length then curresidue^.backbone.main := trace[goal]^.main2.main;
	(* wind protein pointer's first residue to start of fragment *)
	curresidue := prot1^.firstresidue;
	for i := 2 to start do curresidue := curresidue^.nextresidue;
	prot1^.firstresidue := curresidue;
	(* add c(beta) *)
	curresidue := prot1^.firstresidue;
	for i := start to goal do
	begin
		with curresidue^.backbone do if resname <> 'GLY ' then
		begin
			if cb_only then (* add cbetas or whole sidechain *)
			begin 
				generate_cb(main.n, main.ca, main.c, cb);
				makeside(curresidue, 'CB  ', cb);
			end
			else add_sidechains(prot1, firstprotein);
		end;
		curresidue := curresidue^.nextresidue;
	end;
	(* write only coordinates to first gap *)
	curresidue := prot1^.firstresidue; (* this points now to start of fragment ! *)
	for i := start+1 to goal do curresidue := curresidue^.nextresidue;
	curresidue^.nextresidue := nil;
	curresidue := prot1^.firstresidue;
	while curresidue <> nil do
	begin
		mainrecord_to_brk(curresidue^.backbone.main, 
			curresidue^.backbone, coorfile);
		curratom := curresidue^.firstside;
		while curratom <> nil do with curresidue^.backbone do
		begin
			write_brk(curratom^.atomname,curratom^.coord,
				curresidue^.backbone.pdbresid,
				'    ',coorfile);
			curratom := curratom^.nextatom;
		end;
		curresidue := curresidue^.nextresidue;
	end;
	(* write fragment source to logfile *)
	for i := start to goal do 
	begin
		write(logfile, i:5, trace[i]^.sourceresno:5,' ',
		proteins[trace[i]^.sourceprotno].protname:11, nscore[i]:5,
		fragindex[i]:5, redundancy[i]:5, trace[i]^.rmssum:5:1);
		if (fragindex[i]>0) and (fragindex[i]<=maxscores) 
		then writeln(logfile, 
		candidates[i, fragindex[i]].score:8:1,
		candidates[i, fragindex[i]].carms:8:1) else writeln(logfile);
	end;
end; (* path *)

procedure score_to_fragmentlist(var score: score_record);
var
	k: integer;
	lastfragment, newfragment: fragmentpointer;
	tempfragment: fragmentcoord;
	prot2: proteinpointer;
	unique: array[1..maxlen+1] of boolean; (* max matchlength on maxlen+1 ! *)

	procedure get_transrot_fragment(var fragment: fragmentcoord;
		beg1, beg2, length,p: integer);
	var
		res1, res2, curresidue: residuepointer;
		i: integer;
		u: rotationmatrix; t: translationvector; rms: real;
		x, y: t2; j, ier: integer; w: t1;
	
	begin (* get_transrot_fragment *)
		(* 1 moves, 2 is static *)
		res1 := dist[proteins[p].start+beg1].respo;
		curresidue := res1; (* save for later use *)
		res2 := dist[beg2].respo;
		(* match c(alphas) *)
		matchca(res1, res2, length, u, t, rms);
		(* now transrotate, store result in fragment[] *)
		for i := 1 to length do with curresidue^.backbone do
		begin
			transrot(main.n, fragment[i].main.n, u, t);
			transrot(main.ca, fragment[i].main.ca, u, t);
			transrot(main.c, fragment[i].main.c, u, t);
			transrot(main.o, fragment[i].main.o, u, t);
			curresidue := curresidue^.nextresidue;
			(* write all fragment coordinates to a file *)
(*			with fragment[i].main do
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
*)		end;
	end; (* get_transrot_fragment *)

begin (* score_to_fragmentlist *)
	with score do
	begin
		get_transrot_fragment(tempfragment, matchresidue, begresidue, 
			matchlength, matchprotlistindex);
		(* 	checkaa uniikit *)
		for k := 1 to matchlength do 
		begin
			unique[k] := true;
			lastfragment := fragmentlist[k+begresidue-1];
			while (lastfragment <> nil) AND unique[k] do
			begin
				if (distance(lastfragment^.main1.main.c, 
					tempfragment[k].main.c) <= 0.2)
				then if (distance(lastfragment^.main2.main.c,
					tempfragment[k+1].main.c) <= 0.2)
				then if (distance(lastfragment^.main1.main.n,
					tempfragment[k].main.n) <= 0.2)
				then if (distance(lastfragment^.main2.main.n,
					tempfragment[k+1].main.n) <= 0.2)
				then if (distance(lastfragment^.main1.main.ca,
					tempfragment[k].main.ca) <= 0.2)
				then if (distance(lastfragment^.main2.main.ca,
					tempfragment[k+1].main.ca) <= 0.2)
				then unique[k] := false;
				lastfragment := lastfragment^.nextfragment;
			end;
		end;
		for k := 1 to matchlength do if unique[k] then
		begin	(* make new fragment, update links to other fragments *)
			new(newfragment);
			newfragment^.nextfragment := nil;
			newfragment^.prevfragment := nil;
			newfragment^.sourceprotno := matchprotlistindex;
			newfragment^.sourceresno := matchresidue+k-1;
			newfragment^.rmssum := 0;
			newfragment^.njumpto := 0;
			newfragment^.maxdelta := 0; 
			newfragment^.main1 := tempfragment[k];
			if k<matchlength then newfragment^.main2 := tempfragment[k+1];
			if k>1 then newfragment^.main0 := tempfragment[k-1];
			newfragment^.firstjump := nil;
			newfragment^.prev := nil;
			update_jump_fromoldtonew(fragmentlist, k, begresidue, prot1^.length, newfragment,tolerance);
			update_jump_fromnewtoold(fragmentlist, k, begresidue, prot1^.length, newfragment,tolerance);
			newfragment^.nextfragment := fragmentlist[k+begresidue-1];
			if fragmentlist[k+begresidue-1]<>nil then fragmentlist[k+begresidue-1]^.prev := newfragment;
			newfragment^.prev := nil;
			fragmentlist[k+begresidue-1] := newfragment;
			redundancy[k+begresidue-1] := redundancy[k+begresidue-1]+1;
		end;
	end;
end; (* score_to_fragmentlist *)


function firstgap(chainstart, testlen: integer): integer;
var      goal: integer;
begin
      goal := chainstart-1;
      repeat goal := goal+1 until gap[goal] or (goal>testlen);
      firstgap := goal;
end;

procedure change_parameters;
var
	c: char;
begin
	(* set range of residues to build *)
	writeln('build range Y/[N] ?  whole protein ',chainstart, chainend);
	readln(c);
	if c in ['Y','y'] then 
	begin
		writeln('enter new range ');
		readln(chainstart, chainend);
	end;
	(* set width of gap repair zone *)
	writeln('change repair zone Y/[N] ?  currently ',leftzone, rightzone);
	readln(c);
	if c in ['Y','y'] then 
	begin
		writeln('enter new leftzone, rightzone');
		readln(leftzone, rightzone);
	end;
	(* set max no of fragments per residue to fetch *)
	writeln('change maxscores Y/[N] ?  currently ',maxscores);
	readln(c);
	if c in ['Y','y'] then 
	begin
		writeln('enter new maxscores');
		readln(maxscores);
		if maxscores > maxmaxscores then
		begin
			writeln('out of range: setting to ',maxmaxscores);
			maxscores := maxmaxscores;
		end;
	end;
	(* change fragment search parameters *)
	writeln('change distance deviation cutoffs Y/[N] ? '); 
	readln(c);
	if c in ['Y', 'y'] then
	begin
		writeln('QuickSearch looplen currently ',l1,' enter new ? ');
		readln(l1);
		writeln('QuickSearch cutoff1 currently ',c1,' enter new ? ');
		readln(c1);
		writeln('FullSearch looplen currently ',l2,' enter new ? ');
		readln(l2);
		writeln('FullSearch cutoff1 currently ',c2,' enter new ? ');
		readln(c2);
	end;
	(* relax looplen if fewer than maxscores fragments ? *)
	writeln('change minlooplen Y/[N] ?  currently ',minlooplen);
	readln(c);
	if c in ['Y','y'] then 
	begin
		writeln('enter new minlooplen');
		readln(minlooplen);
		if maxscores < 3 then
		begin
			writeln('out of range: setting to ', 3);
			minlooplen := 3;
		end;
	end;
	(* set tolerance for joining *)
	writeln('change tolerance Y/[N] ?  currently ',tolerance);
	readln(c);
	if c in ['Y','y'] then 
	begin
		writeln('enter new tolerance');
		readln(tolerance);
	end;
	(* output file names *)
	writeln('coordinate output to ',coorfilename,' OK ? [Y]/N ');
	readln(c);
	if c in ['N','n'] then 
	begin
		writeln('enter new filename ');
		readln(coorfilename);
	end;
	writeln('log to ', logfilename,' OK ? [Y]/N ');
	readln(c);
	if c in ['N','n'] then 
	begin
		writeln('enter new filename ');
		readln(logfilename);
	end;
	writeln('fragment statistics to ',statfilename,' OK ? [Y]/N ');
	readln(c);
	if c in ['N','n'] then
	begin
		writeln('enter new filename ');
		readln(statfilename);
	end;
end;

(* ------------------------------------------------------------------ *)

procedure fetch(begres, looplen, cutoff1: integer;
                quick, dosort: boolean;
            var nscore: integer);

var
	enough: boolean;
	i, d, e: integer;
	candide: array[1..maxres] of score_record;
	di: ilist;

	procedure scanfragments(e, row: integer);
	var
		curhit: hitpointer;
		i,j: integer;
	label quit;
	
		procedure test;
		var
			dchir, nright, nleft, ndist: integer;
			dssq: real;

		begin with proteins[i] do
		begin (* test *)
			dssq := 0;
			extend_to_rigth(start+j, proteins[1].start+begres, 
testlen, length, cutoff1, start, proteins[1].start, dist, nright, dssq);
			if nright >= looplen then
			begin
				extend_to_left(start+j, proteins[1].start+begres, 
testlen, length, cutoff1, start, proteins[1].start, dist, nleft, dssq, nright);
				dchir := ABS(dist[begres-nleft].d[0] - dist[start+j-nleft].d[0]);
				if ((dchir < 90) or (dchir>270)) then
				begin (* save score *)
					nscore := nscore+1;
					with candide[nscore] do (* new score *)
					begin
						matchprotname := proteins[i].protname;
						matchprotlistindex := i;
						matchlength := nright+nleft;
						matchresidue := j-nleft;
						begresidue := begres-nleft;
						ndist := matchlength*(matchlength+1) div 2;
						score := dssq/ndist; (* square of d-rms *)
						if quick then enough := true; (* pick first only in quickfill *)
(* writeln(begresidue:4,matchresidue:4,matchlength:4,score:8:1,matchprotname);
*)					end;
				end;
			end;
		end;
		end; (* test *)

	begin (* scanfragments *)
	  if (e>maxdist) then goto quit; (* array dimension check: skip chain break *)
		curhit := lookup[row,e].nexthit;
		while curhit <> nil do
		begin
			j := curhit^.resno;
			i := curhit^.protno;
			if j>1 then test;
			if enough then goto quit;
			curhit := curhit^.nexthit;
		end;
	quit:
	end; (* scanfragments *)

	procedure drmssort(i,j: integer; var d: ilist);
	var
		m, left, right, e: integer;
		x: real;

	begin (* drmssort *)
		left := i;
		right := j;
		m := (left+right) div 2;
		x := candide[d[m]].score;
		repeat
			while candide[d[i]].score<x do i := i+1;
			while candide[d[j]].score>x do j := j-1;
			if i<=j then
			begin
				e := d[j];
				d[j] := d[i];
				d[i] := e;
				i := i+1;
				j := j-1;
			end;
		until i>j;
		if left<j then drmssort(left, j, d);
		if i<right then drmssort(i, right, d);
	end; (* drmssort *)

	procedure carmssort(i,j: integer; var d: ilist);
	var
		m, left, right, e: integer;
		x: real;

	begin (* carmssort *)
		left := i;
		right := j;
		m := (left+right) div 2;
		x := candide[d[m]].carms;
		repeat
			while candide[d[i]].carms<x do i := i+1;
			while candide[d[j]].carms>x do j := j-1;
			if i<=j then
			begin
				e := d[j];
				d[j] := d[i];
				d[i] := e;
				i := i+1;
				j := j-1;
			end;
		until i>j;
		if left<j then carmssort(left, j, d);
		if i<right then carmssort(i, right, d);
	end; (* carmssort *)

	procedure calculate_carmses;
	var
		i, j, ier: integer;
		res1, res2: residuepointer;
		x,y: t2; w: t1; u: rotationmatrix; t: translationvector;
		cassq: real;

	begin (* calculate_carmses *)
		for i := 1 to min(nscore, maxscores) do with candide[di[i]] do
		begin
			res1 := dist[proteins[matchprotlistindex].start + matchresidue].respo;
			res2 := dist[begresidue].respo;
			for j := 1 to matchlength do
			begin
				x[j] := res1^.backbone.main.ca;
				y[j] := res2^.backbone.main.ca;
				w[j] := 1;
				res1 := res1^.nextresidue;
				res2 := res2^.nextresidue;
			end;
			u3b(w, x, y, matchlength, 0, cassq, u, t, ier); (* 0: get sum of squares only *)
			carms := cassq/matchlength; (* C(alpha) rms squared *)
		end;
	end; (* calculate_carmses *)

begin (* fetch *)
	nscore := 0;
	enough := false;
	while (not enough) and (looplen>=minlooplen) and (nscore<maxscores) do
	begin
		d := dist[begres].d[looplen-1];
		scanfragments(d, looplen-1);
		for e := 1 to cutoff1 do if not enough then
		begin
			scanfragments(d+e, looplen-1);
			if not enough then if d-e>0 then
			scanfragments(d-e, looplen-1);
		end;
		looplen := looplen-1;
	end;
	for i := 1 to nscore do di[i] := i;
	if nscore>1 then if dosort then 
	begin
		drmssort(1, nscore, di);
		calculate_carmses;
		carmssort(1, min(nscore, maxscores), di);
	end;
	for i := 1 to min(nscore, maxscores) do 
		candidates[begres, i] := candide[di[i]];
	nscore := min(nscore, maxscores);
end;  (* fetch *)

(* ------------------------------------------------------------------ *)

procedure quickfill;
var
	i: integer;

begin
	writeln('quickfill ',chainstart:5, chainend:5);
	for i := chainstart to chainend-looplen+1 do
	begin
		fetch(i, looplen, cutoff1, true, false, nscore[i]); (* pick first fragment *)
		if nscore[i]>0 then score_to_fragmentlist(candidates[i,1]);
	end;
end;

procedure seal(chainstart, chainend: integer; var goal: integer);
label giveup;
var
	i, neof, l, l0, left, right: integer;
	allfetched: boolean;

begin (* seal *)
	goal := firstgap(chainstart, chainend);
	writeln('seal ',goal, chainstart, chainend);
	while goal<=chainend do
	begin
		allfetched := false;
		left := max(chainstart, goal-leftzone);
		right := min(goal+rightzone, chainend);
		repeat (* add fragments until gap moves forward or out of fragments *)
			if not allfetched then 
			for i := left to right do if not filled[i] then 
			begin (* fetch fragments for residues in gap zone *)
				fetch(i, looplen, cutoff1, false, true, nscore[i]);
				filled[i] := true;
			end;
			allfetched := true;
			(* add more fragments to fragmentlist *)
			l0 := chainend; (* lowest residue covered by new fragments *)
			neof := 0;
			for i := left to right do
			begin
				fragindex[i] := fragindex[i]+1;
				if fragindex[i]<=nscore[i] then 
				begin
					score_to_fragmentlist(candidates[i, fragindex[i]]);
					l := candidates[i, fragindex[i]].begresidue;
					if l<l0 then l0 := l;
				end
				else neof := neof+1;
			end;
			update_gap(chainstart, l0, goal); 
			if gap[goal] then 
			if neof>( right - left) then goto giveup; 
				(* all fragments in gap zone have been used *)
		until (not gap[goal]); (* gap has been removed *)
		(* how far did the gap move ? *)
		write('gap moved from ',goal);
		update_gap(chainstart, goal-1, chainend); 
		goal := firstgap(chainstart, chainend);
		writeln(' to ',goal);
(*		if goal>8 then freespace(cs, goal-8, fragmentlist);
*)	end;
giveup:
end; (* seal *)

	procedure write_headers;
	begin
		writeln(coorfile,'REMARK  Backbone construction from CA trace by the program MAX-SPROUT, ');
		writeln(coorfile,'REMARK  ver. 1 Oct 1990.  ');
		writeln(coorfile,'REMARK  Reference L. Holm and C. Sander, J.Mol.Biol. (1991) 218, 183-194. ');
		writeln(coorfile,'REMARK  ',master,' residues ',chainstart:5,' to ',chainend:5);
		writeln(statfile,'*  Backbone construction from CA trace by the program MAX-SPROUT, ');
		writeln(statfile,'*  ver. 1 Oct 1990.  ');
		writeln(statfile,'*  Reference L. Holm and C. Sander, J.Mol.Biol. (1991) 218, 183-194. ');
		writeln(statfile,'*  ',master,' residues ',chainstart:5,' to ',chainend:5);
		writeln(logfile,'*  Backbone construction from CA trace by the program MAX-SPROUT, ');
		writeln(logfile,'*  ver. 1 Oct 1990.  ');
		writeln(logfile,'*  Reference L. Holm and C. Sander, J.Mol.Biol. (1991) 218, 183-194. ');
		writeln(logfile,
'*   ',master,' residues ',chainstart:5,' to ',chainend:5);
		writeln(logfile,'*');
		writeln(logfile,
'*   residue number in reconstructed protein');
		writeln(logfile,
'*   !    residue number in fragment origin protein');
		writeln(logfile,
'*   !    !    fragment origin protein');
		writeln(logfile,
'*   !    !    !          number of fragments fetched from the database');
		writeln(logfile,
'*   !    !    !          !    number of test cycles');
		writeln(logfile,
'*   !    !    !          !    !    number of tripeptides to choose from');
		writeln(logfile,
'*   !    !    !          !    !    !   SSQ(distance deviations) col 5 fragment');
		writeln(logfile,
'*   !    !    !          !    !    !       !     C(alpha) rms of fragment');
		writeln(logfile,
'*   !    !    !          !    !    !       !       !  indicated in 5th column');
	end;

procedure logstatistics; (* append fragment statistics to logfile *)
var
	i, ngap, l, maxl, resno, oldresno: integer;
	oldprotname, protname: packed array[1..11] of char;
	hist: array[0..maxres] of integer;
	c: char;
	start: boolean;

begin (* logstatistics *)
	open(logfile, file_name := logfilename, history := readonly);
	reset(logfile);
	(* read three first lines with REMARK ! *)
	for i := 1 to 3 do readln(logfile);
	start := true;
	l := 0;
	maxl := 0;
	for i := 1 to maxres do hist[i] := 0;
	ngap := 0;
	oldresno := 0;
	oldprotname := '!"#$%&/()=';
	while not eof(logfile) do
	begin
		read(logfile, c);
		if c='*' (* skip header *) then readln(logfile) else 
		if c='R' then 
		begin ngap := ngap+1; start := true;  
			if l > maxl then maxl := l;
			if l>0 then hist[l] := hist[l]+1; l := 0;
			readln(logfile); oldprotname := '!@#$%^&*()=' end else
		begin (* read fragment data *)
			readln(logfile, i, resno, protname);
(* writeln(statfile, resno:5,oldresno:5,protname:20, oldprotname:20, l:5);
*)			if (start OR 
			((protname = oldprotname) and (resno = oldresno+1))) 
			then begin l := l+1; oldprotname := protname; 
					oldresno := resno; start := false; end 
			else
			begin
				if l > maxl then maxl := l;
				hist[l] := hist[l]+1;
				l := 1;
				oldresno := resno;
				oldprotname  := protname;
				start := false;
			end;
		end;
	end;
	if l > maxl then maxl := l;
	hist[l] := hist[l]+1;
	(* write histogram *)
	writeln(statfile, 'fragment statistics: longest: ',maxl+1:5,
			 ' number of gaps: ',ngap:5);
	writeln(statfile, 'length     count');
	for i := 1 to maxl do writeln(statfile, i+1:5, hist[i]:10);
end; (*logstatistics *)

begin (* build *)
	(* init all variables *)
	firstprotein := nil;
	load_distances(nprot, master, proteins, dist, firstprotein, prot1);
	fill_lookup(lookup, dist, nprot, proteins);
	chainstart := 1;
	chainend := prot1^.length;
	for i := chainstart to chainend do
	begin
		filled[i] := false;
		gap[i] := true;
		fragindex[i] := 0;
		redundancy[i] := 0;
		fragmentlist[i] := nil;
		nscore[i] := 0;
	end;
	cb_only := true; (* option to build cbetas or whole sidechain *)
	leftzone := 6;
	rightzone := 2;
	minlooplen := 4; (* looplen can be relaxed to fetch at least maxscores fragments *)
	maxscores := maxmaxscores;
	tolerance := 0.4;
	l1 := 5; c1 := 20; l2 := 4; c2 := 80;
	coorfilename := master+'.brk_mod';
	logfilename := master+'.log';
	statfilename := master+'.stat';
	writeln('use default parameters [Y]/N ?'); readln(c);
	if c in ['N','n'] then change_parameters;
	openfile(coorfile, coorfilename, outfile);
	openfile(logfile, logfilename, outfile);
	openfile(statfile, statfilename, outfile);
(*	openfile(fragfile, 'frag.pdb', outfile);
*)	write_headers;

	looplen := l1; cutoff1 := c2;
	quickfill;
	looplen := l2; cutoff1 := c2;
	update_gap(chainstart, chainstart, chainend);
	cs := chainstart;
	repeat
		seal(cs, chainend, goal);
		if goal<=chainend then 
		begin
			writeln('Cannot seal gap ', goal-1, goal);
			if goal-1>cs then path(cs, goal-1);
			writeln(logfile, 'REMARK      Cannot seal gap ',goal-1, goal);
			if goal=cs then goal := goal+1;
			cs := goal;
			for i := chainstart to chainend do
			begin
				filled[i] := false;
				gap[i] := true;
				fragindex[i] := 0;
				redundancy[i] := 0;
				fragmentlist[i] := nil;
				nscore[i] := 0;
			end;
			update_gap(cs, cs, chainend);
			goal := firstgap(cs, chainend);
			writeln('Starting new chain from ', cs);
		end;
	until (goal>chainend) or (chainend-cs<looplen-1); (* won't find fragments for very short chains *)
	if goal-1>cs then path(cs, goal-1);

	close(coorfile);
	close(logfile);
	logstatistics;
	close(statfile);
(*	close(fragfile);
*)end; (* build *)

end. (* module *)



