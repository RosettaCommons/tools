c___________________________________________________________________________
c
c	This set of programs constructs full atomic coordinates of a protein
c	from a given C(alpha) trace and optimizes side chain geometry.
c
c 	Copyright by Liisa Holm and Chris Sander, 1989-1991.
c
c	No redistribution, no program changes, no commercial use.
c
c	For details, see J.Mol.Biol. 218, 183-194 (1991).
c
c___________________________________________________________________________
c	
c     **********************************************************
c     **                                                      **
c     **  torso.for  --  main program			      **
c     **                                                      **
c     **********************************************************
c
c	VAX/VMS version: 	ostype = 'VAX '
c	Unix versions:		ostype = 'DECs','Iris','Sun4'
c
c     **********************************************************
c
c	read a Brookhaven file
c	setup a Monte Carlo optimization
c	do the optimization
c	logical solvate controls use of solvation term (this version is DRY)
c
c	input channels:	torso.lib -- for071
c			mc_param.dat -- for081
c	input/output	*.tmp -- for080
c
	program torso

	implicit none
	include 'sizes.for'
	
c	general:
	integer inp, iout, ierr, irotlib, ireslib, idglib, itmpfile,
     $		ientropy
	character*80 infile, outfile, parse_outfile
	character*20 brkdfltfile, dsspdfltfile
	character*4 ostype
c	getcoor:
	integer natom, nres, atmres(maxatm)
	real xyz(3, maxatm), bvalue(maxatm)
c	getdssp
	integer acc(maxres)
	character seq(maxres), struc(maxres)
c	getstatic:
	logical staticside(maxres), staticmain(maxres)
	integer staticatm(maxatm), nstaticatm, firstatm(maxres)
	character*3 resnam(maxres), atmnam(maxatm)
c	getlibrary:	
	integer nrotlib
	real rotkhi1(maxrot), rotkhi2(maxrot)
	character*3 rotresnam(maxrot), rotatmnam(maxrot, 14)
	real xyzrot(3, maxrot, 14), wala(3,4)
	logical found(maxrot, 14)
c	getrotamers
	integer resrot(maxatm), nresrot, rottores(maxres), rottolib(maxres)
	integer atmtype(maxatm), firstrotatm(maxres)
c	fillgrid
	integer grid(-maxgrid:maxgrid,-maxgrid:maxgrid,-maxgrid:maxgrid)
	integer overlay(maxatm), nocc, occ(3,maxatm)
c	search:
	integer search(3,1000), nsearch
	integer searchbox(3,2000), nbox
c	energy
	real mcmene(8,8,0:3600)	
	real estatic(maxres)
	logical activeatm(maxatm), activerot(maxres), jactive(maxpair),
     $		rejected(maxres)
	integer ipair(maxpair), jpair(maxpair), 
     $		ni(maxres)
	integer npair, j_index(maxpair), i_first(maxres)
	real epair(maxpair), eres(maxres), etot
c	more pointers
	integer ires_to_irot(maxres), nirot(maxres), rotamer(maxres)
c	Monte Carlo
	real bestemin
c	solvation:
c	occstatic -- iatm's occupancy due to all static atoms
c	occiatm, occjrot, occs have the same index
c	occiatm -- iatm
c	occjrot -- iatm's contacting rotamer
c	occs -- iatm's occupancy due to jrot
c	niocc(iatm)= number of stored occupancies for iatm
c	firstocc(iatm)= index to occ_index
c	
	real occstatic(maxatm), occs(maxnoccs)
	real volume(1000), minocc(1000), maxocc(1000),gsol(1000), gtrf(1000)
	integer noccs, niocc(maxatm), njrot(maxres),firstjrot(maxres), i4,j4
	integer solvatomtype(maxatm), occiatm(maxnoccs), occjrot(maxnoccs)
	integer occ_index(maxnoccs),firstocc(maxatm),occjrot_index(maxnoccs)
c	locals:
	character c
	integer i,j,n,ires, k, irot, jrot
	real property(maxatm), ecutoff
	logical solvate
	integer pdbresno(maxres)
c	khi deviations
	real khi1(maxres), khi2(maxres)
c	testing:
	real testxyz(3,maxatm)
	character*3 testatmnam(maxatm),aa
	integer testatmres(maxatm)
c----------------------------------------------------------------------
c	VAX/Unix specifics
c
c	VAX/VMS version: 	ostype = 'VAX '
c	Unix versions:		ostype = 'DECs','Iris','Sun4'
c
c	Sun4 does not use tmp file !
c
c
c
c	>>>> set ostype according to your system <<<<
c
	ostype = 'Linux'
c
c
c
	if(ostype.eq.'VAX ') then
		brkdfltfile = '.brk'
		dsspdfltfile = '.dssp'
	else
		brkdfltfile = ' '
		dsspdfltfile = ' '
	end if
c----------------------------------------------------------------------
c	print MaxSprout header
c
	write(*,*)
     $  '________________________________________________________'
	write(*,*)
     $ ' '
	write(*,*)
     $ 'This set of programs constructs full atomic coordinates'
	write(*,*)
     $ 'of a protein from a given C(alpha) trace and optimizes '
	write(*,*) 'side chain geometry.'
	write(*,*) ' '
	write(*,*) 'Copyright by Liisa Holm and Chris Sander, 1989-1991.'
	write(*,*) ' '
	write(*,*)
     $ 'No redistribution, no program changes, no commercial use.'
	write(*,*) ' '
	write(*,*) 'For details, see J.Mol.Biol. 218, 183-194 (1991) '
	write(*,*) 'and Proteins 14, 213-223 (1992).'
	write(*,*) ' '
	write(*,*)
     $ '________________________________________________________'

c----------------------------------------------------------------------
c	define input, output units
c
	inp = 5
	iout=6
	irotlib = 71
	ireslib = 72
	idglib = 21
	itmpfile = 80
	ientropy = 11
c----------------------------------------------------------------------
c	input Brookhaven file
c
	write(*,*) ' give input brk file name (default=',brkdfltfile,') '
	read(*,13) infile
	outfile=parse_outfile(infile, ostype)
	write(*,*) outfile,' is the output file'
13	format(a80)
c----------------------------------------------------------------------
c	start directly from tmp file ?
c
c	write(*,*) ' start directly from tmp file Y/N ? '
c	read(*,12) c
12	format(a1)
c	if((c.eq.'Y').or.(c.eq.'y')) goto 890
c	if((c.eq.'E').or.(c.eq.'e')) goto 910
c----------------------------------------------------------------------
c	use solvation ? YES
c
	solvate=.true.
c	write(*,*) ' wet (Y) or dry (N - default) ? '
c	read(*,12) c
c	if((c.eq.'Y').or.(c.eq.'y')) solvate=.true.
c----------------------------------------------------------------------
c	input DSSP file
c
	call getcoor(10, iout, infile, xyz,natom,atmnam,
     $			atmres, nres,resnam, .true., ierr,bvalue, 
     $			brkdfltfile,pdbresno)
c	write(*,*) ' give input DSSP file name (default=',dsspdfltfile,') '
c	read(*,13) infile
c	call getdssp(infile, seq, struc, acc, ca, i, j, dsspdfltfile)
c        if(j.eq.-1) write(*, 500) infile
c        if(j.gt.0) write(*, 510) infile, j
500     format(' !!! DSSP file not found: ', a9)
510     format(' !!! error reading DSSP file: ',a9,' after line ',i5)
c	convert acc to relative acc in %
c	IF(J.EQ.0) THEN
c	  do i=1,nres
c		acc(i)=1000*acc(i)/reference(seq(i))
c	  end do
c	END IF
c----------------------------------------------------------------------
c	get static coordinates = all main chain + CB atoms by atmnam
c	
	call checkCbeta(natom, nres, atmnam, resnam, atmres, xyz)
	call getkhi(xyz, atmnam, resnam, atmres, natom, nres, khi1, khi2)
c
c	read fixed sidechains 
c
c	default: no fixation
	do i=1,nres
		staticside(i)=.false.
		staticmain(i)=.true.
	end do	
c	write(*,*) ' do you want to fix sidechains ? (Y/N) '
c	read(*,80) c
	c='n'
80	format(a1)
	if((c.eq.'Y').or.(c.eq.'y')) then
	  write(*,*) ' 0 - free all, fix by resno; 1 - fix all, free by resno '
	  read(*,*) i
	  if(i.eq.0) then
301		write(*,*) ' enter resno to fix (0 to quit) '
		read(*,*) i
		if(i.gt.0) then
			staticside(i)=.true.
			write(*,*) i, resnam(i), ' fixed'
			goto 301
		end if	
	  else
		do i=1,nres
			staticside(i)=.true.
		end do
310		write(*,*) ' enter resno to free (0 to quit) '
		read(*,*) i
		if(i.gt.0) then
			staticside(i)=.false.
			write(*,*) i, resnam(i), ' free'
			goto 310
		end if
	  end if
	end if
c----------------------------------------------------------------------
c	write(*,*) ' do you want to mutate the sequence ? (Y/N) '
c	read(*,80) c
	c='n'
	if((c.eq.'Y').or.(c.eq.'y')) then
410		write(*,*) 
     $		' enter resno to mutate, new residue (i4,1x,a3) (0 to quit) '
		read(*,570) i, aa
		if(i.gt.0) then
			write(*,*) i, resnam(i), ' mutated to ', aa
			staticside(i)=.false.
			resnam(i)=aa
			goto 410
		end if
	end if	
570	format(i4,1x,a3)
c----------------------------------------------------------------------
c
	call getstatic(natom, atmnam, atmres, staticside, staticmain, 
     $			nstaticatm, staticatm)
	write (iout, 120) nstaticatm
120	format (/' number of static atoms :', i10)
c
c	compress x,y,z,natom to include only static atoms
c
	call compress_static(xyz,natom, nstaticatm, staticatm, atmnam,
     $		atmres, firstatm)
c----------------------------------------------------------------------
c	get rotamer library
c	superpose all non-glycine rotamers on alanine
c	
	call getlibrary(iout,irotlib,nrotlib,rotkhi1, rotkhi2, rotatmnam, 
     $		rotresnam, xyzrot, found, ierr, brkdfltfile)
c        write (*, *) "getlibrary finished"
	call standardizelibrary(nrotlib, rotresnam, xyzrot, found, wala)
c       write (*, *) "standardizelibrary finished"
c	wala is now the common backbone of all rotamers
c----------------------------------------------------------------------
c	load rotamers for varying residues (staticside(ires)=.false.)
c	default rotamers are only native sequence
c	add coordinates to xyz
c	CBETAs were already built in static, not added here
c
c	index to 'rotamers'
	nresrot=0
	call getrotamers(firstatm, atmnam, natom, xyz, xyzrot, found,
     $		rotatmnam, rotresnam, nresrot, resrot, staticside, nrotlib,
     $		nres, resnam, atmres, rottores, rottolib, wala, firstrotatm)
c        write (*, *) "getrotamers finished"
	write(*,*) ' nresrot ', nresrot
	if (nresrot.gt.maxres)
     $ 	stop ' Too many rotamers ! - increase maxres '
c----------------------------------------------------------------------
c	get atom types
c	define energy function
c	create prototype grid-box shell
c
	call atomtype(iout, natom, atmres, atmnam, resnam, atmtype, nres)
c        write (*, *) "atomtype finished"
	call init_mcmenergy(mcmene)
c        write (*, *) "init_mcmenergy finished"
	call init_search(search, nsearch)
c        write (*, *) "init_search finished"
c	call init_solvation(ireslib,idglib,volume,minocc,maxocc,gtrf,gsol)
c	
c	get solvation library atomtypes 
c
c	do i=1,natom
c		solvatomtype(i)=atomid(atmnam(i), resnam(atmres(i)))
c		write(*,*) i, solvatomtype(i), atmnam(i), resnam(atmres(i))
c	end do
c----------------------------------------------------------------------
c	"local" rotamer-static energy:
c	ires-1,ires,ires+1 static atoms ... irot varying atoms
c
	call geteself(estatic, nres, nresrot, rottores, firstrotatm,
     $		resrot, firstatm, atmres, atmnam, resnam, xyz, mcmene,atmtype)
c        write (*, *) "geteself finished"
c
c	dirty trick: set PRO eself to zero
c
	do i=1,nresrot
		if(resnam(rottores(i)).eq.'PRO') estatic(i)=0.0
	end do
c
c	"long-range" rotamer-static energy:
c	1..ires-2,ires+2..nres static atoms ... irot varying atoms
c
c	put all atoms to grid 
c
	do i=1,natom
		activeatm(i)=.true.
		occstatic(i)=0.0
	end do
	nocc=0
	call fillgrid(iout, grid, overlay, xyz, natom, nocc, occ,
     $		activeatm)
c        write (*, *) "fillgrid finished"
c
	call getestaticsolv(estatic, nstaticatm, search, nsearch, xyz,
     $		grid, overlay, resrot, atmres, mcmene, atmnam, resnam,
     $		atmtype, occstatic, solvatomtype, volume, solvate)
c        write (*, *) "getestaticsolv finished"
c----------------------------------------------------------------------
c	use estatic to reject clashing rotamers
c	reject all rotamers whose estatic > max(emin+10.0, 10.0)
c	by flagging activeatm as .false.
c
	call reject(estatic, activeatm, activerot, nres, nresrot,
     $		nstaticatm, natom, rottores, resrot, rejected)
c        write (*, *) "reject finished"
c
c	count rejected rotamers
c
	j=0
	do i=1,nresrot
		if(.not.activerot(i)) j=j+1
	end do
	write(*,*) j, ' rotamers of ',nresrot,' rejected by estatic '
c
c	make static atoms active
c
	do i=1,nstaticatm
		activeatm(i)=.true.
	end do
c
c	clear grid; put in only active atoms
c
	call fillgrid(iout, grid, overlay, xyz, natom, nocc, occ,
     $		activeatm)
c        write (*, *) "fillgrid finished again"
c----------------------------------------------------------------------
c	rotamer-rotamer energies
c	active varying ... active varying
c
	npair=0
	noccs=0
	ecutoff=0.1
	do i=1,nres
		call getsearchboxsolv(i, searchbox, natom, nstaticatm,
     $			firstrotatm, atmres, xyz, grid, overlay,
     $			nsearch, search, nbox, estatic, resrot, 
     $			mcmene, atmnam, resnam, atmtype, nresrot, rottolib,
     $			rottores, activeatm, activerot,
     $			npair, ipair, jpair, epair, ecutoff,
     $			occiatm, occjrot, occs, solvatomtype, volume, noccs,
     $			solvate)
	end do
c        write (*, *) "getsearchboxsolv finished"
	write(iout, *) npair, ' rotamer-rotamer interactions saved'
	write(iout, *) noccs,  ' partial atomic occupancies saved '
	j4=0
	do i4=1,npair
		if(epair(i4).gt.0.0) j4=j4+1
        	write(*,*) ' epair=', epair(i4), ' i4=', i4
	end do
	write(*,*) j4, ' positives'
c----------------------------------------------------------------------
c	construct pointers for i-j-epair arrays
c
	call getpairpointers(npair, ipair, jpair, i_first, ni, 
     $		j_index, nresrot, ires_to_irot, nirot, rottores)
c        write (*, *) "getpairpointers finished"
c
c	construct pointers for occ* array
c
	do i=1,natom
		niocc(i)=0
	end do
	call getoccpointers(noccs, niocc, njrot, firstjrot, occ_index,
     $		firstocc, occjrot_index, occiatm, occjrot)
c        write (*, *) "getoccpointers finished"
c----------------------------------------------------------------------
c	what is the max # of pairs per irot ?
c
	n=0
	do i=1,nresrot
		if(ni(i).gt.n) then
			n=ni(i)
			j=i
		end if
	end do
	write(*,*) ' max(ni) is ',n, j
c----------------------------------------------------------------------
	if(noccs.gt.maxnoccs) write(*,*) ' WARNING: noccs exceeds maxnoccs ',
     $		noccs, maxnoccs
c----------------------------------------------------------------------
c	write tables for Monte Carlo to a temp file
c
c
c	Sun4 skips temp file
c
c	if(ostype.eq.'Linux') goto 900
c
	write(itmpfile) npair, nres, nresrot, noccs, nstaticatm, natom
	write(itmpfile) (j_index(i4),i4=1,npair)
	write(itmpfile) (rejected(i), i=1,nresrot)
	write(itmpfile) (i_first(i), ni(i), rottores(i), i=1,nresrot)
	write(itmpfile) (ipair(i4), jpair(i4), epair(i4), i4=1,npair)
	write(itmpfile) (nirot(i), ires_to_irot(i), i=1,nres)
	write(itmpfile) (estatic(i), i=1,nresrot)
	write(itmpfile) (firstocc(i), occstatic(i), niocc(i), i=1,natom)
	write(itmpfile) (occjrot(i4), occs(i4), occ_index(i4), i4=1,noccs)
	write(itmpfile) (maxocc(i), minocc(i), gsol(i), gtrf(i), i=1,1000)
	write(itmpfile) (solvatomtype(i), atmres(i), resrot(i), i=1,natom)
	write(itmpfile) (firstrotatm(i), i=1,nres)
	write(itmpfile) (njrot(i), firstjrot(i), i=1,nresrot)
	write(itmpfile) (occiatm(i4), occjrot_index(i4), i4=1,noccs)
	write(itmpfile) (xyz(1,i), xyz(2,i), xyz(3,i), atmnam(i), i=1,natom)
	write(itmpfile) (resnam(i), i=1,nres)
	write(itmpfile) (seq(i), struc(i), acc(i), i=1,nres)
	write(itmpfile) (khi1(i), khi2(i), i=1,nres)
	write(itmpfile) nrotlib, (rotkhi1(i), rotkhi2(i), i=1,nrotlib)
	write(itmpfile) (rottolib(i), i=1, nresrot)
	write(itmpfile) solvate
	close(itmpfile)
	goto 900
c----------------------------------------------------------------------
c	read tables for Monte Carlo from a temp file
c
890	read(itmpfile) npair, nres, nresrot, noccs, nstaticatm, natom
	read(itmpfile) (j_index(i4),i4=1,npair)
	read(itmpfile) (rejected(i), i=1,nresrot)
	read(itmpfile) (i_first(i), ni(i), rottores(i), i=1,nresrot)
	read(itmpfile) (ipair(i4), jpair(i4), epair(i4), i4=1,npair)
	read(itmpfile) (nirot(i), ires_to_irot(i), i=1,nres)
	read(itmpfile) (estatic(i), i=1,nresrot)
	read(itmpfile) (firstocc(i), occstatic(i), niocc(i), i=1,natom)
	read(itmpfile) (occjrot(i4), occs(i4), occ_index(i4), i4=1,noccs)
	read(itmpfile) (maxocc(i), minocc(i), gsol(i), gtrf(i), i=1,1000)
	read(itmpfile) (solvatomtype(i), atmres(i), resrot(i), i=1,natom)
	read(itmpfile) (firstrotatm(i), i=1,nres)
	read(itmpfile) (njrot(i), firstjrot(i), i=1,nresrot)
	read(itmpfile) (occiatm(i4), occjrot_index(i4), i4=1,noccs)
	read(itmpfile) (xyz(1,i), xyz(2,i), xyz(3,i), atmnam(i), i=1,natom)
	read(itmpfile) (resnam(i), i=1,nres)
	read(itmpfile) (seq(i), struc(i), acc(i), i=1,nres)
	read(itmpfile) (khi1(i), khi2(i), i=1,nres)
	read(itmpfile) nrotlib, (rotkhi1(i), rotkhi2(i), i=1,nrotlib)
	read(itmpfile) (rottolib(i), i=1, nresrot)
	read(itmpfile) solvate
	close(itmpfile)
c----------------------------------------------------------------------
c	how many YWF in core ?
c
	do i=1, nres
	  n=0
	  if(acc(i).le.200) then
	    if((seq(i).eq.'Y').or.(seq(i).eq.'W').or.(seq(i).eq.'F')) then
		do j=1,nresrot
		  if((rottores(j).eq.i).and.(.not.rejected(j))) n=n+1
		end do
		write(*,*) seq(i), i, n, acc(i)
	    end if
	  end if
	end do
c----------------------------------------------------------------------
c	Monte Carlo 
c
900	if(noccs.eq.0) noccs=1
c	 write (*,*) 'npair='        , npair
c	 write (*,*) 'j_index='      , j_index
c	 write (*,*) 'nres='         , nres
c	 write (*,*) 'nresrot='      , nresrot
c	 write (*,*) 'rotamer='      , rotamer
c	 write (*,*) 'rejected='     , rejected
c	 write (*,*) 'i_first='      , i_first
c	 write (*,*) 'ni='           , ni
c	 write (*,*) 'rottores='     , rottores
c	 write (*,*) 'ipair='        , ipair
c	 write (*,*) 'jpair='        , jpair
c	 write (*,*) 'nirot='        , nirot
c	 write (*,*) 'bestemin='     , bestemin
c	 write (*,*) 'eres='         , eres
c	 write (*,*) 'etot='         , etot
c	 write (*,*) 'epair='        , epair
c	 write (*,*) 'estatic='      , estatic
c	 write (*,*) 'jactive='      , jactive
c	 write (*,*) 'ires_to_irot'  , ires_to_irot
c	 write (*,*) 'activeatm='    , activeatm
c	 write (*,*) 'activerot='    , activerot
c	 write (*,*) 'nstaticatm='   , nstaticatm
c	 write (*,*) 'natom='        , natom
c	 write (*,*) 'firstocc='     , firstocc
c	 write (*,*) 'occstatic='    , occstatic
c	 write (*,*) 'niocc='        , niocc
c	 write (*,*) 'occjrot='      , occjrot
c	 write (*,*) 'occs='         , occs
c	 write (*,*) 'occ_index='    , occ_index
c	 write (*,*) 'maxocc='       , maxocc
c	 write (*,*) 'minocc='       , minocc
c	 write (*,*) 'solvatomtype=' , solvatomtype
c	 write (*,*) 'gsol='         , gsol
c	 write (*,*) 'gtrf='         , gtrf
c	 write (*,*) 'atmres='       , atmres
c	 write (*,*) 'resrot='       , resrot
c	 write (*,*) 'firstrotatm='  , firstrotatm
c	 write (*,*) 'occiatm='      , occiatm
c	 write (*,*) 'occjrot_index=', occjrot_index
c	 write (*,*) 'njrot='        , njrot
c	 write (*,*) 'firstjrot='    , firstjrot
c	 write (*,*) 'noccs='        , noccs
c	 write (*,*) 'solvate='      , solvate
c	 write (*,*) 'khi1='         , khi1
c	 write (*,*) 'khi2='         , khi2
c	 write (*,*) 'rotkhi1='      , rotkhi1
c	 write (*,*) 'rotkhi2='      , rotkhi2
c	 write (*,*) 'rottolib='     , rottolib
	call protomcmsolv(npair, j_index, nres, nresrot, rotamer, rejected,
     $		i_first, ni, rottores, ipair, jpair, nirot, bestemin,
     $		eres, etot, epair, estatic, jactive, ires_to_irot,
     $		activeatm, activerot, nstaticatm, natom, firstocc,
     $		occstatic, niocc, occjrot, occs, occ_index, maxocc, minocc,
     $		solvatomtype, gsol, gtrf, atmres, resrot, firstrotatm,
     $		occiatm, occjrot_index, njrot, firstjrot, noccs, solvate,
     $		khi1, khi2, rotkhi1, rotkhi2, rottolib)
c	write(*,*) 'protomcmsolv finished'        
c----------------------------------------------------------------------
c	write out entropies to ientropy
c
c910	call entropy(84, s, freq, seq, acc)
c----------------------------------------------------------------------
c	energy check 
c
	write(*,*) ' bestemin ', bestemin
	do i4=1,npair
		jactive(i4)=.false.
	end do
	do i=1,nres
		jrot=rotamer(i)
		if(jrot.gt.0) then
			do i4=i_first(jrot),i_first(jrot)+ni(jrot)-1
				jactive(j_index(i4))=.true.
			end do
		end if
	end do
	etot=0
	do ires=1,nres
		irot=rotamer(ires)
		eres(ires)=estatic(irot)
		do i=0,ni(irot)-1
			i4=i_first(irot)+i
			if(jactive(i4)) eres(ires)=eres(ires)+epair(i4)
		end do
		write(*,*) ires, eres(ires), rotamer(ires)
		etot=etot+eres(ires)
	end do
	write(*,*) ' etot ', etot
c----------------------------------------------------------------------
c	write out optimized coordinates
c
	j=0
	do i=1,nstaticatm
		j=j+1
		do k=1,3
			testxyz(k,j)=xyz(k,j)
		end do
		testatmnam(j)=atmnam(i)
		testatmres(j)=atmres(i)
		property(j)=eres(atmres(i))
	end do
	do i=nstaticatm+1, natom
	  if(resrot(i).eq.rotamer(atmres(i))) then
		j=j+1
		do k=1,3
			testxyz(k,j)=xyz(k,i)
		end do
		testatmnam(j)=atmnam(i)
c restore original numbering
		testatmres(j)=pdbresno(atmres(i))
		property(j)=eres(atmres(i))
	  end if
	end do
	open(17, file=outfile)
	call putcoor(17, testxyz, j, testatmnam, testatmres, resnam, property)
	close(17)
	write(*,*) j, ' atoms in optimized structure ', outfile
c----------------------------------------------------------------------
c	
c
      end

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c----------------------------------------------------------------------

