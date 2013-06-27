c__________________________________________________________________________
c
c	This set of programs constructs full atomic coordinates of a protein
c	from a given C(alpha) trace and optimizes side chain geometry.
c
c	Copyright by Liisa Holm and Chris Sander, 1989-1991.
c
c	No redistribution, no program changes, no commercial use.
c
c	For details, see J.Mol.Biol. 218, 183-194 (1991).
c
c___________________________________________________________________________
c     **********************************************************
c	torsotool.for
c     **********************************************************
c
	subroutine checkCbeta(natom, nres, atmnam, resnam, atmres, xyz)
c
c	build missing Cbeta coordinates
c
	implicit none
	include 'sizes.for'
	integer natom, nres, atmres(maxatm)
	real xyz(3, maxatm)
	character*3 atmnam(maxatm), resnam(nres)
c
	integer ires, iatm, i, j
	real cbeta(3), ca(3,maxres), c(3,maxres), n(3,maxres), 
     $		ca0(3), c0(3), n0(3)
	logical cbfound(maxres), cafound(maxres), cfound(maxres), 
     $		nfound(maxres), ofound(maxres)
	character*3 atnm
c	sort:
	integer v(maxatm), tempatmres(maxatm), code
	integer n4, d(maxatm)
	character*3 tempatmnam(maxatm)
	real temp(3, maxatm)
c
c	put backbone coordinates in table, flag logicals *found
c
	do ires=1,nres
		cafound(ires)=.false.
		cfound(ires)=.false.
		nfound(ires)=.false.
		ofound(ires)=.false.
		cbfound(ires)=.false.
	end do
	do iatm=1, natom
		ires=atmres(iatm)
		atnm=atmnam(iatm)
		if(atnm.eq.'CA ') then
			cafound(ires)=.true.
			do i=1,3
				ca(i,ires)=xyz(i, iatm)
			end do
		else if(atnm.eq.'C  ') then
			cfound(ires)=.true.
			do i=1,3
				c(i,ires)=xyz(i, iatm)
			end do
		else if(atnm.eq.'N  ') then
			nfound(ires)=.true.
			do i=1,3
				n(i,ires)=xyz(i, iatm)
			end do
		else if(atnm.eq.'O  ') then
			ofound(ires)=.true.
		else if(atnm.eq.'CB ') then
			cbfound(ires)=.true.
		end if
	end do
	do ires=1, nres
		if(.not.cafound(ires)) write(*, 150) ires, 'CA '
		if(.not.cfound(ires)) write(*, 150) ires, 'C  '
		if(.not.ofound(ires)) write(*, 150) ires, 'O  '
		if(.not.nfound(ires)) write(*, 150) ires, 'N  '
	  	if(resnam(ires).eq.'GLY') goto 99
		if(.not.cbfound(ires)) then
			write(*, 150) ires, 'CB '
			write(*,*) 'generating CB coordinates '
c	  		generate CBETA from n-ca-c
			do i=1,3
				ca0(i)=ca(i,ires)
				c0(i)=c(i,ires)
				n0(i)=n(i,ires)
			end do
			call generate_cb(n0, ca0, c0, cbeta)
			natom=natom+1
			do i=1,3
				xyz(i, natom)=cbeta(i)
				atmnam(natom)='CB '
				atmres(natom)=ires
			end do
		end if
99	end do
150	format(' Warning: incomplete backbone ',i5,2x,a3,' missing')
c
c	sort
c
	do i=1,natom
		d(i)=i
		v(i)=atmres(i)*14+code(atmnam(i))
		do j=1,3
			temp(j,i)=xyz(j,i)
		end do
		tempatmres(i)=atmres(i)
		tempatmnam(i)=atmnam(i)
	end do
	n4=natom
	call j_index_Qsort(1, n4, n4, d, v)
	do i=1, natom
		do j=1,3
			xyz(j,i)=temp(j, d(i))
		end do
		atmres(i)=tempatmres(d(i))
		atmnam(i)=tempatmnam(d(i))
	end do

	return
	end
c
c     **********************************************************
c
	subroutine getkhi(xyz, atmnam, resnam, atmres, natom, nres, 
     $		khi1, khi2)
c
c	find khi1 and khi2 angles for all residues in xyz
c
	implicit none
	include 'sizes.for'
	integer natom, nres, atmres(natom)
	real xyz(3, natom), khi1(nres), khi2(nres)
	character*3 atmnam(natom), resnam(nres)
c
	integer i, j, k
	real table(3, 0:7, maxres), v1(3), v2(3), v3(3), v4(3)
	real dihedralangle

	do i=1,natom
		j=0
		if(atmnam(i).eq.'C  ') j=1
		if(atmnam(i).eq.'O  ') j=2
		if(atmnam(i).eq.'N  ') j=3
		if(atmnam(i).eq.'CA ') j=4
		if(atmnam(i).eq.'CB ') j=5
		if((atmnam(i).eq.'CG ').or.(atmnam(i).eq.'SG ').or.
     $			   (atmnam(i).eq.'OG ').or.
     $			   (atmnam(i).eq.'CG1').or.(atmnam(i).eq.'OG1')) j=6
		if((atmnam(i).eq.'CD ').or.(atmnam(i).eq.'OD1').or.
     $			   (atmnam(i).eq.'ND1').or.(atmnam(i).eq.'CD1').or.
     $			   (atmnam(i).eq.'SD ')) j=7
c
c		  khi1=n-ca-cb-... khi2=ca-cb-cg-... khi3=cb-cg-cd-...
c		  table ix: n-3 ca-4 cb-5 'cg'-6 'cd'-7 '...'- up to 14
		do k=1, 3
			table(k, j, atmres(i))=xyz(k,i)
		end do
	end do	
	do i=1, nres
		if((resnam(i).eq.'GLY').or.(resnam(i).eq.'ALA')) then
			khi1(i)=0.0
			khi2(i)=0.0
		else if((resnam(i).eq.'VAL').or.(resnam(i).eq.'SER').or.
     $			(resnam(i).eq.'THR').or.(resnam(i).eq.'CYS')) then
			do j=1,3
				v1(j)=table(j,3,i)
				v2(j)=table(j,4,i)
				v3(j)=table(j,5,i)
				v4(j)=table(j,6,i)
			end do
			khi1(i)=dihedralangle(v1, v2, v3, v4)
			khi2(i)=0.0
		else
			do j=1,3
				v1(j)=table(j,3,i)
				v2(j)=table(j,4,i)
				v3(j)=table(j,5,i)
				v4(j)=table(j,6,i)
			end do
			khi1(i)=dihedralangle(v1,v2,v3,v4)
			do j=1,3
				v1(j)=table(j,7,i)
			end do
			khi2(i)=dihedralangle(v2,v3,v4,v1)
		end if
	end do

	return
	end
c
c     **********************************************************
c
	subroutine entropy(irotfile, s, freq, seq, acc)
c
c	input:  irotfile - unit with nrun * [energy, nres, rotamer()]
c			   (output by MC)
c		seq, acc
c	output: s	 - entropy per residue, SUM freq(irot)*ln(freq(irot))
c		freq	 - frequency of rotamers
c
	implicit none
	include 'sizes.for'
	integer irotfile
	integer acc(maxres)
	real s(maxres), freq(maxres)
	character seq(maxres)
c
	integer nrun, nres, i, rotamer(maxres), rottores(maxres), irot,
     $		nresrot
	real emin

c	read file
	rewind(irotfile)
	nrun=0
	nresrot=0
10	read(irotfile, end=19) emin, nres, (rotamer(i), i=1,nres)
	nrun=nrun+1
	do i=1, nres
		irot=rotamer(i)
		if(irot.gt.0) then
			freq(irot)=freq(irot)+1
			rottores(irot)=i
			if(irot.gt.nresrot) nresrot=irot
		end if
	end do
	goto 10
19	close(irotfile)
	write(*,*) nrun,' runs read from ',irotfile, nresrot, nres
	do i=1,nres
		s(i)=0
	end do
	if(nrun.eq.0) nrun=1
	do irot=1,nresrot
		freq(irot)=freq(irot)/nrun
c		write(*,*) irot, freq(irot)
		if(freq(irot).gt.0) s(rottores(irot))=s(rottores(irot))-
     $			freq(irot)*log(freq(irot))
	end do
        do i=1, nres
                write(11, 191) i, char(9),
     $                  seq(i), char(9), acc(i), char(9), s(i)
        end do
191     format(i8,3a1,i8,a1,f10.3)

	return
	end
c
c     **********************************************************
c
	subroutine getdssp(filnam, seq, struc, acc, ca, nres, ier, dfltfile)
c
c	read dssp file to get seq(), struc(), acc(), ca()
c
	implicit none
	include 'sizes.for'
	character*80 filnam
	character*20 dfltfile
	character seq(maxres), struc(maxres), c
	real ca(3,maxres)
	integer nres, acc(maxres), ier
c
	integer ires, i

c	empty string -> file not found
	if(len(filnam).eq.0) goto 120
c	until line with # 
	ier=0
	i=0
	open(20, file=filnam, status='old',
     $		err=120)
	goto 101
c	file not found handler
120		ier=-1
		goto 2999
c	end file not found handler
101	read(20, 1501) c
	if(c.ne.'#') goto 101
1501	format(2x,a1)
c	read residue, sec-struc data
201	read(20, 1601, err=1700, end=301) ires, seq(ires), struc(ires), 
     $		acc(ires), ca(1, ires), ca(2, ires), ca(3, ires)
1601	format(i5, 8x, a1, 2x, a1, 18x, i3, 69x, 3f7.1)
	i=ires
	goto 300
c	error handle reading DSSP file format
1700		ier=i
		goto 301
c	end error handle reading DSSP file format
300	goto 201
301	close(20)

2999	return
	end
c
c     **********************************************************
c
	function parse_outfile(infile, ostype)
	implicit none
	character*80 infile, parse_outfile
	character*4 ostype
c
	integer i,j,k,l

c	parse infile to get outfile: j is start of code, k is end of code
c
c	ostype.eq.'VAX '
c	VAX syntax {dir:}{[.subdir{.subdir}]}code{.brk}
c                      j   k       k       j      k  l
c
c	ostype.ne.'VAX ' 
c	Unix syntax {.,/dir/subdir/}code{.suffixes}
c                                 j      k       l
c
	if (ostype.eq.'VAX ') then
		do i=1,80
			parse_outfile(i:i)=' '
		end do
		j=1
		i=1
		do while(infile(i:i).eq.' ')
			i=i+1
		end do
		l=i-1
			do i=1,len(infile)
			if(infile(i:i).eq.':') j=i+1
			if(infile(i:i).eq.']') j=i+1
			if(infile(i:i).eq.'.') k=i-1
			if(infile(i:i).ne.' ') l=l+1
		end do
		if(k.eq.0) k=l
		l=0
		do i=j,k
			l=l+1
			parse_outfile(l:l)=infile(i:i)
		end do
		parse_outfile(l+1:l+7)='_sc.brk'
	else
		do i=1,80
			parse_outfile(i:i)=' '
		end do
		j=1
		i=1
		do while(infile(i:i).eq.' ')
			i=i+1
		end do
		l=i-1
		k=0
		do i=1,len(infile)
			if(infile(i:i).eq.'/') j=i+1
			if(infile(i:i).eq.'.') k=i-1
			if(infile(i:i).ne.' ') l=l+1
		end do
		if(k.eq.0) k=l
		l=0
		do i=j,k
			l=l+1
			parse_outfile(l:l)=infile(i:i)
		end do
		parse_outfile(l+1:l+7)='_sc.brk'
	end if

	return
	end
c----------------------------------------------------------------------
c
c
	subroutine init_solvation(ireslib, idglib, volume,
     $		minocc, maxocc, gtrf, gsol)
	implicit none
	integer ireslib, idglib, ierr
	character c
	integer i
	real volume(1000),minocc(1000),maxocc(1000),gsol(1000),gtrf(1000),cc

c	read solvation parameters from reslib
c	the index is atomid(atmnam, resnam)
c
	ierr=-1
	open(ireslib, status='old', err=45)
	ierr=0
43	read(ireslib,132, err=45) c
	if(c.ne.'#') goto 43
44	read(ireslib,131,err=45) i,volume(i),minocc(i),maxocc(i),cc,cc
c	write(*,*) i, volume(i),minocc(i),maxocc(i),gsol(i),gtrf(i)
	goto 44
131	format(i5,10x,5f8.3)
132 	format(a1)
45	close(ireslib)
	if(ierr.eq.-1) write(*,*) ' unable to open occupancy parameters file '	
c
c	read deltaG from dglib
c
	ierr=-1
	open(idglib, err=153, status='old')
	ierr=0
53	read(idglib,152,end=153, err=153) i, gtrf(i), gsol(i)
152	format(i5,2f6.2)
c	gsol(i)=gsol(i)*10
	goto 53
153	close(idglib)
	if(ierr.eq.-1) write(*,*) ' unable to open solvation parameters file '	
	
	return
	end
c----------------------------------------------------------------------
c
c
	subroutine getoccpointers(noccs, niocc, njrot, firstjrot, occ_index,
     $		firstocc, occjrot_index, occiatm, occjrot)
	implicit none
	include 'sizes.for'
	integer noccs, niocc(maxatm), njrot(maxres), firstjrot(maxres),
     $		occ_index(maxnoccs), firstocc(maxatm), occjrot_index(maxnoccs)
	integer occiatm(maxnoccs), occjrot(maxnoccs)
c
	integer noccs4, i4, one
	integer i

	do i4=1,noccs
		occ_index(i4)=i4
		occjrot_index(i4)=i4
	end do
	noccs4=noccs
	i4=1
	if(noccs4.gt.0) then
 	  call j_index_Qsort(i4,noccs4,noccs4, occ_index, occiatm)
	  call j_index_Qsort(i4,noccs4,noccs4, occjrot_index, occjrot)
	end if
c	niocc(iatm)= number of stored occupancies for iatm
c	firstocc(iatm)= index to occ_index
c	occjrot_index is sorted by occjrot
c	njrot(irot)=number of partial-occupancy-entries for irot
c	firstjrot(irot)= index to occjrot_index
	one=1
	do i4=1,noccs
		i=occiatm(occ_index(i4))
		niocc(i)=niocc(i)+one
		if(firstocc(i).eq.0) firstocc(i)=i4
		i=occjrot(occjrot_index(i4))
		njrot(i)=njrot(i)+one
		if(firstjrot(i).eq.0) firstjrot(i)=i4
	end do

	return
	end
c	
c     **********************************************************
c
	subroutine getsearchboxsolv(ires, searchbox, natom, nstaticatm, 
     $		firstrotatm, atmres, xyz, grid, overlay, nsearch, search,
     $		nbox, e_x, resrot, mcmene, atmnam, resnam, atmtype, nresrot,
     $		rottolib, rottores, activeatm, activerot, 
     $		npair, ipair, jpair, epair, ecutoff, occiatm, occjrot, occs,
     $		solvatomtype, volume, noccs, solvate)
c	
c	npair=counter for rotamer-rotamer interactions
c	ipair=irot [resrot#]
c	jpair=jrot [resrot#]
c	epair=energy(irot,jrot) -- abs(e)>ecutoff
c
c	solvation:
c	solvate = logical flag for including solvation term
c	occiatm = iatm
c	occjrot = contacting rotamer jrot
c	occs    = occupancy of iatm due to jrot
c	noccs    = number of saved partial occupancies
c	solvatomtype, volume from torso.reslib
c	
	implicit none
	include 'sizes.for'
	integer noccs, npair
	integer ires, searchbox(3,2000), natom, nstaticatm, nresrot,
     $		firstrotatm(maxres), atmres(maxatm), overlay(maxatm),
     $		grid(-maxgrid:maxgrid,-maxgrid:maxgrid,-maxgrid:maxgrid),
     $		nsearch, search(3,1000), nbox, resrot(maxatm),atmtype(maxatm),
     $		rottolib(maxres), rottores(maxres), ipair(maxpair),
     $		jpair(maxpair), solvatomtype(1000), 
     $		occiatm(maxnoccs), occjrot(maxnoccs)
	real xyz(3,maxatm), e_x(maxres), mcmene(8,8,0:3600), efun,
     $		epair(maxpair), ecutoff, volume(1000), occs(maxnoccs)
	character*3 atmnam(natom), resnam(maxres)
	logical activeatm(natom), activerot(maxres), solvate
c
	integer iatmlist(maxatm), i, n, gx, gy, gz, ncell,
     $		j, iatm, icellist(3,1000), nx, jatm, jres, nrot,
     $		haamu(-maxgrid:maxgrid,-maxgrid:maxgrid,-maxgrid:maxgrid),
     $		irot, jrot, irotlist(400), tsiltori(maxres)
	integer fung
	real enaybor(400, maxres), e
	logical hasnaybors(maxrot)

c
c	skip GLY, ALA
c
	if(resnam(ires).eq.'ALA') goto 99
	if(resnam(ires).eq.'GLY') goto 99
c	make list of resrot# --> list index [irotlist] and reverse [tsiltori]
c	can handle max 400 rotamers at one residue -- warn if more !
	nrot=0
	do i=1,nresrot
		if((activerot(i)).and.(rottores(i).eq.ires)) then
			nrot=nrot+1
			irotlist(nrot)=i
			tsiltori(i)=nrot
		end if
	end do
c	write(*,*) nrot, ' active rotamers for ', ires
c	write(*,*) (irotlist(i),i=1,nrot)
	if(nrot.gt.400) 
     $	stop ' Too many rotamers per site -- modify getsearchbox'
c	make list of varying iatms in ires
	i=firstrotatm(ires)
	n=0
	nx=0
	do while (atmres(i).eq.ires)
	  if(activeatm(i)) then
		n=n+1
		iatmlist(n)=i
	  else
		nx=nx+1
	  end if
	  i=i+1
	end do
c	write(*,*) ires, n, ' varying atoms', nx, ' inactives'
c	make list of occupied cells to icellist() and flag haamu()
c	put self-occupied cells to box-list
	ncell=0
	nbox=0
	do i=1,n
		iatm=iatmlist(i)
		gx=fung(xyz(1,iatm))
		gy=fung(xyz(2,iatm))
		gz=fung(xyz(3,iatm))
		if(haamu(gx,gy,gz).ne.ires) then
			ncell=ncell+1
			icellist(1,ncell)=gx
			icellist(2,ncell)=gy
			icellist(3,ncell)=gz
			haamu(gx,gy,gz)=ires
			jatm=grid(gx,gy,gz)
c			add to list -- including only ires cells !!!
			nbox=nbox+1
			searchbox(1,nbox)=gx
			searchbox(2,nbox)=gy
			searchbox(3,nbox)=gz
		end if
	end do
c	write(*,*) ncell, ' occupied cells'
c	create list of search cells to searchbox() with absolute grid 
c		coordinates
	nx=0
	do i=1,ncell
c	  new cells
	  do j=1,nsearch
		gx=icellist(1,i)+search(1,j)
		gy=icellist(2,i)+search(2,j)
		gz=icellist(3,i)+search(3,j)
c		reject if empty
		if(grid(gx,gy,gz).eq.0) goto 20
c		check if new cell is already in list
		if(haamu(gx,gy,gz).eq.ires) goto 20
		haamu(gx,gy,gz)=ires
c		add to list -- cell is new and occupied
		nbox=nbox+1
		searchbox(1,nbox)=gx
		searchbox(2,nbox)=gy
		searchbox(3,nbox)=gz
20	  end do
	end do
c	write(*,*) ires, nbox, ' boxes in searchbox'
c	calculate estatic:
c	all boxes in searchbox
c		all jatms !!! only varying left in grid !!! -- NOT ANY MORE !!!
c			all iatms in iatmlist
c
c	zero naybor counter 
	do i=1,nrot
		hasnaybors(i)=.false.
		do j=1,nresrot
			enaybor(i,j)=0.0
		end do
	end do
	do i=1,nbox
		jatm=grid(searchbox(1,i), searchbox(2,i), searchbox(3,i))
110		if(jatm.eq.0) goto 120
		jres=atmres(jatm)
		jrot=resrot(jatm)
c		solvation block
c		if(solvate) then
c		sum up contacts of jatm in occ(irot)
c
c		do irot=1,nrot
c			occ(irot)=0.0
c		end do
c		do j=1,n
c			iatm=iatmlist(j)
c			irot=tsiltori(resrot(iatm))
c			occ(irot)=occ(irot)+volume(solvatomtype(iatm))*
c    $			envelope(distance(xyz(1,iatm),xyz(2,iatm),xyz(3,iatm),
c     $				          xyz(1,jatm),xyz(2,jatm),xyz(3,jatm)))
c		end do
c		save contacts in occiatm, occjrot, occs
c		do irot=1,nrot
c		  do not save contacts irot1-irot2 which belong to the same ires
c		  if(occ(irot).gt.0.0) then
c		    	if((jres.eq.ires).and.(jatm.gt.nstaticatm).and.
c    $			   (irotlist(irot).ne.jrot)) goto 111
c			noccs=noccs+1
c			if(mod(noccs,100).eq.0) write(*,*) ' noccs',noccs
c			occiatm(noccs)=jatm
c			occjrot(noccs)=irotlist(irot)
c			occs(noccs)=occ(irot)
c111		  end if
c		end do
c		end if
c		end of solvation block
c		no ires in vdw calculation !
c		and exclude static atoms, too !
		if((jatm.gt.nstaticatm).and.(jres.ne.ires)) then
			jrot=resrot(jatm)
     			do j=1,n
			    	iatm=iatmlist(j)
				irot=tsiltori(resrot(iatm))
				hasnaybors(irot)=.true.
	e=efun(iatm,atmtype(iatm),jatm,atmtype(jatm),xyz,mcmene)
	enaybor(irot, jrot)=enaybor(irot, jrot)+e
		 	end do
		end if
		jatm=overlay(jatm)
		goto 110
120	end do

c	compress enaybor to erotarota(irot, [naybor-no]), nayblist
	do irot=1,nrot
	  if(hasnaybors(irot)) then
	    do j=1,nresrot
c	      if(enaybor(irot,j).gt.10.0) then
c		write(*,*) ires, resnam(ires), irotlist(irot), j, rottores(j),
c     $			enaybor(irot,j)
c	      end if
c	reject interactions weaker than ecutoff
	      if(abs(enaybor(irot,j)).gt.ecutoff) then
		npair=npair+1
		if(npair.gt.maxpair) write(*,*) npair, maxpair
		ipair(npair)=irotlist(irot)
		jpair(npair)=j
		epair(npair)=enaybor(irot,j)
	      end if
	    end do
c	    write(*,*) irot, rottores(irot), npair
	  end if
	end do
c	check for too many pairs
	if(npair.gt.maxpair)
     $	stop ' Too many rotamer-rotamer pairs -- increase MAXPAIR'

99	return
	end
c
c     **********************************************************
c
	subroutine getestaticsolv(estatic, nstaticatm, search, nsearch,
     $		 xyz, grid, overlay, resrot, atmres, mcmene, atmnam, resnam,
     $		 atmtype, occstatic, solvatomtype, volume, solvate)
c
c	get interaction energy between varying rotamers and static atoms
c	interactions between ires and ires-1..ires+1 have been calculated
c	in eself
c	now include only guaranteedly nonbonded interations between
c	ires and 1..ires-2, ires+2..nres
c	(use in next step to reject clashing rotamers !)
c
c	occstatic = occupancy for all jatms due to static atoms
c	solvatomtype is used in torso.reslib (different from mcmene!)
c	volume is corresponding volume
c
c	logical solvate controls inclusion of solvation term
c
	implicit none
	include 'sizes.for'
	integer nstaticatm, nsearch, search(3,nsearch), overlay(maxatm),
     $		grid(-maxgrid:maxgrid,-maxgrid:maxgrid,-maxgrid:maxgrid),
     $		resrot(maxatm), atmres(maxatm), atmtype(maxatm), 
     $		solvatomtype(maxatm)
	real estatic(maxres), xyz(3, maxatm), mcmene(8,8,0:3600),
     $		occstatic(maxatm), volume(10000)
	character*3 atmnam(maxatm), resnam(maxres)
	logical solvate
c
	integer iatm, jatm, icell, gx,gy,gz, ires, jres
	integer fung
	real efun

c	loop over all static atoms; save energy with irot>0 in estatic(irot)
	do iatm=1,nstaticatm
	  ires=atmres(iatm)
c	  all search grid cells
	  gx=fung(xyz(1,iatm))
	  gy=fung(xyz(2,iatm))
	  gz=fung(xyz(3,iatm))
	  do icell=1,nsearch
	    jatm=grid(gx+search(1,icell), gy+search(2,icell), 
     $		   gz+search(3,icell))
10	    if(jatm.ne.0) then
c	      if(solvate) then
c	        get contact strength 
c	        occstatic(jatm)=occstatic(jatm)+volume(solvatomtype(iatm))*
c     $		  envelope(distance(xyz(1,iatm),xyz(2,iatm),xyz(3,iatm),
c     $			          xyz(1,jatm),xyz(2,jatm),xyz(3,jatm)))
c	      end if
c	      do pairwise energy with varying nonbonded atoms in cell+overlay
	      if(jatm.gt.nstaticatm) then
		      	jres=atmres(jatm)
			if((jres.lt.ires-1).or.(jres.gt.ires+1)) then
c	write(*,*) ires, iatm, atmnam(iatm), resnam(ires),
c     $		   jres, jatm, atmnam(jatm), resnam(jres), resrot(jatm)
     			  estatic(resrot(jatm))= estatic(resrot(jatm))
     $			  +efun(iatm, atmtype(iatm), jatm, atmtype(jatm),
     $				xyz, mcmene)
			end if
	      end if
	      jatm=overlay(jatm)
	      goto 10
	    end if	
	  end do
	end do

	return
	end
c	
c     **********************************************************
c
c	
c     **********************************************************
c
	subroutine j_index_Qsort(i,j,n,d,v)
c	
c	non-recursive Quicksort
c	sorts integer array v which has n elements from smallest to largest
c	d contains indices to v
c
	implicit none
	integer i,j,n,d(n)
	integer v(n)
c
	integer p,top,bottom,stack(100), nstack, partition

	stack(1)=j
	stack(2)=i
	nstack=2
	do while(nstack.ne.0)
		top=stack(nstack)
		bottom=stack(nstack-1)
		nstack=nstack-2
		do while(top.lt.bottom)
			p=partition(top, bottom, n, d, v)
			if((p-top).gt.(bottom-p)) then
				stack(nstack+1)=p-1
				stack(nstack+2)=top
				top=p+1
				nstack=nstack+2
				if(nstack.gt.100) 
     $	stop ' Stack overflow in Qsort '
			else
				stack(nstack+1)=bottom
				stack(nstack+2)=p+1
				bottom=p-1
				nstack=nstack+2
				if(nstack.gt.100) 
     $	stop ' Stack overflow in Qsort '
			end if
		end do
	end do

	return
	end

c ======================================================================
	function partition(i,j,n,d,v)
	implicit none
	integer i,j,n,partition, d(n)
	integer v(n)
c
	integer upper, lower, save

	upper=i
	lower=j
	save=d(i)

	do while (upper.ne.lower)
		do while((upper.lt.lower).and.(v(save).le.v(d(lower))))
			lower=lower-1
		end do
		if(upper.ne.lower) d(upper)=d(lower)
		do while((upper.lt.lower).and.(v(save).ge.v(d(upper))))
			upper=upper+1
		end do
		if(upper.ne.lower) d(lower)=d(upper)
	end do
	d(upper)=save
	partition=upper

	return
	end
c	
c     **********************************************************
c
	subroutine getpairpointers(npair, ipair, jpair, i_first, ni, 
     $		j_index, nresrot, ires_to_irot, nirot, rottores)
	implicit none
	include 'sizes.for'
	integer npair, j_index(npair)
	integer nresrot, ipair(npair), jpair(npair), 
     $		ni(nresrot), ires_to_irot(maxres), nirot(maxres),
     $		rottores(maxres)
	integer i_first(nresrot)
c
	integer i,j, ires
	integer irot

	write(*,*) npair
	do i=1,npair
		j_index(i)=i
	end do
	i=1
	j=npair
	call j_index_qsort(i, j, npair, j_index, jpair)
c	test
c	write(*,150) (jpair(j_index(i)),i=1,npair)
c150	format(20i4)

	do i=1,npair
		irot=ipair(i)
		ni(irot)=ni(irot)+1
		if(i_first(irot).eq.0) i_first(irot)=i
c		jrot=j_index(i)
	end do
c
c	pointer from ires --> i_first index (first irot of ires)
c	
	do i=1,nresrot
		ires=rottores(i)
		nirot(ires)=nirot(ires)+1
		if(ires_to_irot(ires).eq.0) ires_to_irot(ires)=i
	end do
	
	return
	end
c	
c     **********************************************************
c
	subroutine geteself(estatic, nres, nresrot,rottores,firstrotatm,
     $		resrot, firstatm, atmres, atmnam, resnam, xyz, mcmene,atmtype)
c
c	calculate interaction energy of all rotamers with their own residue
c	and with static atoms belonging to nearest neighbour residues
c	i.e. ires-1, ires, ires+1
c	this is the only place where nonbonded-ness needs to be checked !
c
	implicit none
	include 'sizes.for'
	integer nres, nresrot, rottores(maxres),
     $		firstrotatm(maxres),resrot(maxatm),firstatm(maxres),
     $		atmres(maxatm),atmtype(maxatm), fatm
	character*3 atmnam(maxatm), resnam(maxres)
	real estatic(maxres), efun, xyz(3,maxatm), mcmene(8,8,0:3600)
	logical nonbonded
c
	integer irot, ires, jres, iatm, jatm

	do irot=1, nresrot
		ires=rottores(irot)
c		wind iatm to start of irot
		iatm=firstrotatm(ires)
		if(iatm.eq.0) goto 99
10		if((resrot(iatm).ne.irot).and.(atmres(iatm).eq.ires)) then
			iatm=iatm+1
			goto 10
		end if
		fatm=iatm
		do while (resrot(iatm).eq.irot)
c	write(*,*) irot, ires, iatm, atmnam(iatm), resnam(atmres(iatm))
c 			for testing: backbones of 1..nres
c			do jres=1,nres
c			backbones of ires-1,ires,ires+1
			do jres=max(1,ires-1),min(nres,ires+1)
				jatm=firstatm(jres)
				do while(atmres(jatm).eq.jres)
c	write(*,*) ' backbone ', jatm, jres, atmnam(jatm), atmres(jatm),
c     $		resnam(atmres(jatm))
					if 
     $	(nonbonded(jatm,iatm,atmres,atmnam,resnam))
     $	estatic(irot)=estatic(irot)+efun(iatm,atmtype(iatm),
     $		jatm, atmtype(jatm), xyz, mcmene)
					jatm=jatm+1
				end do
			end do
c			sidechain atoms of irot
c			do jatm=fatm,iatm-1
c	write(*,*) ' sidechain ', jatm, atmnam(jatm), resnam(atmres(jatm))
c					if 
c    $	
c     $	estatic(irot)=estatic(irot)+efun(iatm,atmtype(iatm),
c     $		jatm, atmtype(jatm), xyz, mcmene)
c			end do
			iatm=iatm+1
		end do
99	end do

	return
	end
c	
c     **********************************************************
c
	function nonbonded(iatm, jatm,atmres, atmnam, resnam)
c
c	!!! ONLY FOR USE WITH GETESTATIC !!!
c	iatm is static, jatm is varying
c
c	for varying atoms (*G* - ...)
c	exclude:	
c		with the same ires
c			all - CB
c			*G* - all except O
c			*D* - CA
c		the following static atoms:
c		CD(ires, PRO) - CA(ires-1), C(ires-1), O(ires-1)
c		CG(ires, PRO) - C(ires-1)
c
	implicit none
	include 'sizes.for'
	logical nonbonded
	integer iatm, jatm, atmres(maxatm)
	character*3 atmnam(maxatm), resnam(maxres)
c
	integer ires, jres

	ires=atmres(iatm)
	jres=atmres(jatm)
	nonbonded=.true.
	if (jres.eq.ires) then
	  if(atmnam(iatm).eq.'CB ') nonbonded=.false.
	  if((atmnam(jatm)(2:2).eq.'G').and.(atmnam(iatm).ne.'O  '))
     $		 nonbonded=.false.
	  if((atmnam(jatm)(2:2).eq.'D').and.(atmnam(iatm).eq.'CA '))
     $		 nonbonded=.false.
	else if(ires.eq.jres-1) then
	  if(resnam(jres).eq.'PRO') then
	    if ((atmnam(jatm).eq.'CD ').and.((atmnam(iatm).eq.'C  ').or.
     $	      (atmnam(iatm).eq.'CB ').or.
     $	      (atmnam(iatm).eq.'O  ').or.(atmnam(iatm).eq.'CA ')))
     $		nonbonded=.false.
	    if ((atmnam(jatm).eq.'CG ').and.(atmnam(iatm).eq.'C  '))
     $		nonbonded=.false.
	  end if
	end if
	
	return
	end
c	
c     **********************************************************
c
	subroutine reject(estatic, activeatm, activerot, nres, nresrot,
     $		nstaticatm, natom, rottores, resrot, rejected)
	implicit none
	include 'sizes.for'
	real	estatic(maxres)
	integer nres, nresrot, nstaticatm, natom, rottores(maxres),
     $		resrot(maxatm)
	logical activeatm(maxatm), activerot(maxres), rejected(maxres)
c
	integer i
	real minestatic(maxres)

c	find best estatic(irot) for each residue
	do i=1,nres
		minestatic(i)=999999.999
	end do
	do i=1,nresrot
		if(estatic(i).lt.minestatic(rottores(i))) 
     $			minestatic(rottores(i))=estatic(i)
	end do
c	reject all rotamers whose estatic > max(emin+10.0, 0.0)
c	flag irot, iatms as .false.
	do i=1,nres
		if(minestatic(i).le.0.0) then
			minestatic(i)=10.0
		else
			minestatic(i)=minestatic(i)+10.0
		end if
	end do
	do i=1,nresrot
		if(estatic(i).gt.minestatic(rottores(i))) then
			activerot(i)=.false.
			rejected(i)=.true.
		else
			activerot(i)=.true.
			rejected(i)=.false.
		end if
	end do	
c	init activeatm flags: no static atoms; all varying atoms
	do i=1,nstaticatm
		activeatm(i)=.false.
	end do
	do i=nstaticatm+1,natom
		activeatm(i)=activerot(resrot(i))
	end do
	
	return
	end
c	
c     **********************************************************
c
	subroutine init_search(search, n)
c
c	get relative coordinates of 2 A grid cells for 6 A cutoff radius
c	n = number of cells to search
c	principle: inaccuracy of coordinates is +-1A
c		   condition: min((x+-1)-(+-1))**2+min(..y..)+min(..z..)<=R**2
c				<=> dx**2+dy**2+dz**2-2(dx+dy+dz)<=R**2-12
c
	implicit none
	integer search(3, 1000), n
c
	integer dx,dy,dz,d2,i
	logical x(0:1000)

	do i=0,1000
		x(i)=.false.
	end do
	do dx=-3,3
		do dy=-3,3
			do dz=-2,2
			  d2=dx*dx+dy*dy+dz*dz-2*(abs(dx)+abs(dy)+abs(dz))
			  if (d2.le.6) then
				x((dx+4)+(dy+4)*10+(dz+4)*100)=.true.
				x((dx+4)+(dz+4)*10+(dy+4)*100)=.true.
				x((dy+4)+(dx+4)*10+(dz+4)*100)=.true.
				x((dy+4)+(dz+4)*10+(dx+4)*100)=.true.
				x((dz+4)+(dy+4)*10+(dx+4)*100)=.true.
				x((dz+4)+(dx+4)*10+(dy+4)*100)=.true.
			  end if
			end do
		end do
	end do
	n=0
	do i=0,1000
		if(x(i)) then
			dx=mod(i,10)
			dy=mod(i,100)-dx
			dz=i-dy-dx
			dx=dx-4
			dy=dy/10-4
			dz=dz/100-4
			n=n+1
			search(1,n)=dx
			search(2,n)=dy
			search(3,n)=dz
		end if
	end do

	write(*,*) ' gridcells to be searched per atom: ',n
	return
	end
c	
c     **********************************************************
c
	subroutine getrotamers(firstatm, atmnam, natom, xyz, xyzrot, found,
     $		rotatmnam, rotresnam, nresrot, resrot, staticside, nrotlib,
     $		nres, resnam, atmres, rottores, rottolib, wala, firstrotatm)
c
c	wala is the standardized backbone of all rotamers, ready for u3b
c	1-c 2-n 3-ca 4-cb
c
	implicit none
	include 'sizes.for'
	character*3 resnam(maxres), atmnam(maxatm), rotresnam(maxrot)
	integer natom, atmres(maxatm), firstatm(maxres), nresrot, 
     $		  resrot(maxatm), nrotlib, nres, rottores(maxres),
     $		  rottolib(maxres), firstrotatm(maxres)
	real xyz(3, maxatm), xyzrot(3,maxrot, 14)
        real*8 wala(3,4)
	logical found(maxrot, 14), staticside(maxres)
	character*3 rotatmnam(maxrot, 14)
c
	integer i, j, n, ires
	real*8 wy(3,maxatm)
        real n0(3), ca(3), c0(3), cbeta(3)
	integer ier, k, irot
	real*8 w(4), u(3,3), t(3), ssq
	logical cbfound

c	weight array c-n-ca-cb for u3b
	do i=1,4
		w(i)=1
	end do
	do ires=1,nres
	  if((staticside(ires)).or.(resnam(ires).eq.'GLY')
     $		.or.(resnam(ires).eq.'UNK')) goto 99
c	  generate CBETA from n-ca-c
c		if it is missing only !!!
c	  find n-ca-c of static residue ires (wy)
c	  and CBETA !!!
	  n=0
	  cbfound=.false.
	  do i=firstatm(ires), natom
		if (atmres(i).eq.ires) then
			if(atmnam(i).eq.'CA ') then
				do j=1,3
					wy(j,3)=xyz(j,i)
				end do
				n=n+1
			end if
			if(atmnam(i).eq.'N  ') then
				do j=1,3
					wy(j,2)=xyz(j,i)
				end do
				n=n+1
			end if
			if(atmnam(i).eq.'C  ') then
				do j=1,3
					wy(j,1)=xyz(j,i)
				end do
				n=n+1
			end if
			if(atmnam(i).eq.'CB ') then
				do j=1,3
					wy(j,4)=xyz(j,i)
				end do
				n=n+1
				cbfound=.true.
			end if
			if (n.eq.4) goto 100
		end if
	  end do
c         incomplete backbone
	  write(*,*) n,' backbone atoms for ',ires
100	  continue
c	  generate CBETA if it is missing
	  if(.not.cbfound) then
		write(*,*) ' generating CB coordinates for residue ',ires
		do j=1,3
			n0(j)=wy(j,2)
			ca(j)=wy(j,3)
			c0(j)=wy(j,1)
		end do
		call generate_cb(n0,ca,c0,cbeta)
		do j=1,3
			wy(j,4)=cbeta(j)
		end do
	  end if
c	  get superposition matrices u,t for ca-cbeta-n-c of residue and 
c		rotamers
c	  for static residue ires they are in (wy)
c	  the rotamer ca-cb-n-c coordinates are in wala
c
	  call u3b(w,wala,wy,4,1,ssq,u,t,ier)
	write(*,*) ' u3b rms ',ssq/4, ires, resnam(ires)
c
c
c	  load rotamers of type resnam to residue ires
c 	  find all rotamers with same residue name
	  do irot=1,nrotlib
		  if (rotresnam(irot).eq.resnam(ires)) then
		    nresrot=nresrot+1
		    rottores(nresrot)=ires
		    rottolib(nresrot)=irot
c		    transrotate the whole rotamer sidechain cb-...
c		    ... and add to x,y,z,natom
c
c		    !!! CBETA is not copied !!!
c
	 	    do i=6,14
	  	      if(found(irot,i)) then
			natom=natom+1
			do j=1,3
			  xyz(j, natom)=t(j)
			  do k=1,3
			    xyz(j,natom)=xyz(j,natom)
     $					 +u(j,k)*xyzrot(k,irot,i)
			  end do
			end do
			atmnam(natom)=rotatmnam(irot,i)
			atmres(natom)=ires
			resrot(natom)=nresrot
			if(firstrotatm(ires).eq.0) firstrotatm(ires)=natom
		      end if
		    end do
		  end if
		end do
99	end do
c	check for too many atoms
c	if (natom.gt.maxatm)
c     $	stop ' getrotamers -- Too many atoms; increase MAXATM'
c	check for too many rotamers
c	if (nresrot+nres.gt.maxres)
c     $	stop ' getrotamers -- Too many rotamers; increase MAXRES'
	return
	end

c	
c     **********************************************************
c
	subroutine standardizelibrary(nrotlib, rotresnam, xyzrot, found, wy)
c
c	superpose all non-glycine rotamers on alanine
c	
	implicit none
	include 'sizes.for'
	integer nrotlib
	real xyzrot(3, maxrot, 14)
        real*8 wy(3,4)
	logical found(maxrot, 14)
	character*3 rotresnam(maxrot)
c
	integer iala, irot, i, j, k
	integer ier
	real*8 ssq, w(4), wx(3,4), u(3,3), t(3), rnew(3)
	
c	find alanine
	do irot=1,nrotlib
		if(rotresnam(irot).eq.'ALA') iala=irot
	end do
c	make iala the static vector wy for u3b; c(1)-n(3)-ca(4)-cb(5)
	do i=1,3
		wy(i, 1)=xyzrot(i, iala, 1)
		wy(i, 2)=xyzrot(i, iala, 3)
		wy(i, 3)=xyzrot(i, iala, 4)
		wy(i, 4)=xyzrot(i, iala, 5)
	end do
	do i=1,4
		w(i)=1
	end do
c	all non-glycines and non-alanine
	do irot=1,nrotlib
	  if((rotresnam(irot).ne.'GLY').and.(rotresnam(irot).ne.'ALA')) then
c		get moving vector wx for u3b
		do i=1,3
			wx(i, 1)=xyzrot(i, irot, 1)
			wx(i, 2)=xyzrot(i, irot, 3)
			wx(i, 3)=xyzrot(i, irot, 4)
			wx(i, 4)=xyzrot(i, irot, 5)
		end do
c		call u3b and transrotate all coordinates
		call u3b(w,wx,wy,4,1,ssq,u,t,ier)
		do i=1,14
		  if (found(irot,i)) then
		    do j=1,3
		      rnew(j)=t(j)
		      do k=1,3
			rnew(j)=rnew(j)+u(j,k)*xyzrot(k,irot,i)
		      end do
		    end do
		    do j=1,3
		      xyzrot(j,irot,i)=rnew(j)
		    end do
		  end if
		end do		
	  end if
	end do

	return 
	end

c
c     **********************************************************
c
	function efun(iatm, at1, jatm, at2, xyz, mcmene)
c
c	uses atom types at1, at2
c
	implicit none
	include 'sizes.for'
	real xyz(3, maxatm), a(3), b(3), d(3), efun, dx2,dy2,dz2,dr
	real mcmene(8,8,0:3600)
	integer iatm, jatm, i,d2, at1, at2

	do i=1,3
		a(i)=xyz(i, iatm)
		b(i)=xyz(i, jatm)
	end do
	call diff(a, b, d)
	dx2=d(1)*d(1)
	dy2=d(2)*d(2)
	dz2=d(3)*d(3)
	dr=dx2+dy2+dz2
	efun=0
	if((at1.eq.0).or.(at2.eq.0)) goto 10
	if (dr.lt.36.0) then
		d2=dr*100
		if (d2.eq.0) d2=1
		if (d2.lt.3600) efun=mcmene(at1, at2, d2)
	end if

10	return
	end
c
c     **********************************************************
c

	subroutine atomtype(iout, natom, atmres, atmnam,resnam,atmtype,nres)
c
c	1=all others; 2=polar N/O; 3=charged O; 4=charged N; 5=sulphur
c
      implicit none
      include 'sizes.for'
      integer i,natom, atmres(natom), iout, nres
      integer atmtype(natom), n(5)
      character*3 atmnam(natom),resnam(nres)

c
c	1=all others; 2=polar N/O; 3=charged O; 4=charged N; 5=sulphur
c
	do i=1,5
		n(i)=0
	end do
	do i=1,natom
		atmtype(i)=1
		if (atmnam(i).eq.'N  ') atmtype(i)=2
		if (atmnam(i).eq.'O  ') atmtype(i)=2
		if (atmnam(i)(1:2).eq.'OG') atmtype(i)=2
		if (atmnam(i).eq.'OH ') atmtype(i)=2
		if (atmnam(i)(1:2).eq.'ND') atmtype(i)=2
		if (atmnam(i)(1:2).eq.'NE') atmtype(i)=2
		if ((atmnam(i).eq.'OD1').and.
     $			(resnam(atmres(i)).eq.'ASN')) atmtype(i)=2
		if ((atmnam(i).eq.'OE1').and.(resnam(atmres(i)).eq.'GLN')) 
     $			atmtype(i)=2
		if ((atmnam(i)(1:2).eq.'OD').and.(resnam(atmres(i)).eq.'ASP')) 
     $			atmtype(i)=3
		if ((atmnam(i)(1:2).eq.'OE').and.(resnam(atmres(i)).eq.'GLU')) 
     $			atmtype(i)=3
		if (atmnam(i)(1:2).eq.'NH') atmtype(i)=4
		if (atmnam(i).eq.'NZ ') atmtype(i)=4
		if (atmnam(i)(1:1).eq.'S') atmtype(i)=5
		n(atmtype(i))=n(atmtype(i))+1
	end do

	write(iout,120) (n(i),i=1,5)
  120   format( 
     $	/'1=all others; 2=polar N/O; 3=charged O; 4=charged N; 5=sulphur',
     $		/' Summary of Atom Types:',/
     $          /' Number of atoms of each type :  ', 5i6)

	return
	end
c
c     **********************************************************
c
	subroutine init_mcmenergy(mcmene)
c	(* 6-9 potentials: e = A r(-9) - B r(-6)
c
c	1=all others; 2=polar N/O; 3=charged O; 4=charged N; 5=sulphur
c
c	atom types    	minimum	truncate   	    A	   B	comment
c	cys S-S (5,5)  	2.0    	1.25       	 131.1	 24.6	(and Met)
c	N-O (2-4) 	3.0    	1.75       	 5039	 280
c	all others     	4.0    	2.25       	67000	1573
c
c	all distances in 0.1 A
c	mcmenergy index is (r/0.1 A)**2

	implicit none
	real mcmene(8,8,0:3600)
c
	integer i, at1, at2
	real a, b, e, r, ebase

	ebase=0.0
c	(* all others *)
	a = log(67000.0)
	b = log(1573.0)
	r = log(2.25)
	e = ebase+exp(a-9*r)-exp(b-6*r)
	do i = 0,506
		do at1 = 1,5
			do at2 = at1,5
				mcmene(at1, at2, i)= e
				mcmene(at2, at1, i)= e
			end do
		end do
	end do
	do i = 507,3600 
		r = 0.5*log(float(i))-log(10.0)
		e = ebase+exp(a-9*r)-exp(b-6*r)
		do at1 = 1,5
			do at2 = at1,5
				mcmene(at1, at2, i)= e
				mcmene(at2, at1, i)= e
			end do
		end do
	end do
c	(* S-S *)
	a = log(131.1)
	b = log(24.6)
	r = log(1.25)
	e = ebase+exp(a-9*r)-exp(b-6*r)
	do i = 0,156 
		mcmene(5, 5, i) = e
	end do
	do i = 157, 3600 
		r = 0.5*log(float(i))-log(10.0)
		e = ebase+exp(a-9*r)-exp(b-6*r)
		mcmene(5, 5, i) = e
	end do
c	(* N-O *)
	a = log(5039.0)
	b = log(280.0)
	r = log(1.75)
	e = ebase+exp(a-9*r)-exp(b-6*r)
	do i =0, 306 
		do at1 = 2,4 
			do at2 = at1,4 
			  if(.not.((at1.gt.2).and.(at1.eq.at2))) then
				mcmene(at1, at2, i)= e
				mcmene(at2, at1, i)= e
			  end if
			end do
		end do
	end do
	do i = 307, 3600 
		r = 0.5*log(float(i))-log(10.0)
		e = ebase+exp(a-9*r)-exp(b-6*r)
		do at1 = 2,4
			do at2 = at1,4 
			  if(.not.((at1.gt.2).and.(at1.eq.at2))) then
				mcmene(at1, at2, i)= e
				mcmene(at2, at1, i)= e
			  end if
			end do
		end do
	end do
c	(* print energy profiles *)
	write(*,*) ' mcmenergy atom types'
     	write(*,*) 
     $	' 1=all others; 2=polar N/O; 3=charged O; 4=charged N; 5=sulphur'
c	write(*,*) ' mcmenergy:'
c	write(*,*) ' distance', '     S-S', '     N-O', '     C-C'
c	do i=1,3600,100
c		d=sqrt(float(i))/10
c		write(*,*) d, mcmene(5, 5, i), mcmene(2, 2, i),
c     $			mcmene(1, 1, i)
c	end do

	return
	end

c
c     **********************************************************
c
      subroutine fillgrid(iout, grid, overlay, xyz, natom, 
     $			  nocc, occ, activeatm)
c
c	center molecule at (0,0,0) by shifting coordinates x(),y(),z()
c	put atom numbers to respective grid cells or overlay array
c	clear grid of old atoms
c	insert only active atoms (activeatm(iatm)=.true.)
c
      implicit none
      include 'sizes.for'
      integer natom, iout
      real xyz(3, natom)
      integer grid(-maxgrid:maxgrid, -maxgrid:maxgrid, 
     $			   -maxgrid:maxgrid)
      integer overlay(natom), nocc, occ(3, maxatm)
	logical activeatm(natom)
c
      real xyz0(3)
	integer fung, gx,gy,gz
      integer hist(0:maxatm), i,j,m,n,k, s, nx

c
c	center molecule at 0,0,0
c
	do j=1,3
		xyz0(j)=0
	end do
	write(*,*) natom
	if (natom.eq.0) stop 'GRID -- No atoms !'
	do i=1,natom
	  do j=1,3
		xyz0(j)=xyz0(j)+xyz(j,i)
	  end do
	end do
	do j=1,3
		xyz0(j)=xyz0(j)/natom
	end do
c
c	shift all atoms by -x0,-y0,-z0 
c
	do i=1,natom
	  do j=1,3
		xyz(j,i)=xyz(j,i)-xyz0(j)
	  end do
	end do
c
c	zero grid and overlay -- all nocc atoms in occ()
c
	do i=1,nocc
		grid(occ(1,i), occ(2,i), occ(3,i))=0
	end do
	hist(0)=0
	do i=1,natom
		overlay(i)=0
		hist(i)=0
	end do
c
c	put atoms in grid -- only actives
c	(or pointers to overlay)
c	!!! warn if atoms are outside grid !!!
c
	nx=0
	m=0
	nocc=0
	do i=1,natom
	  if(activeatm(i)) then
		gx=fung(xyz(1,i))
		gy=fung(xyz(2,i))
		gz=fung(xyz(3,i))
c		check inside
		if ((abs(gx).gt.maxgrid).or.(abs(gy).gt.maxgrid).or.
     $		    (abs(gz).gt.maxgrid)) then
			nx=nx+1
			goto 19
		end if
		j=grid(gx,gy,gz)
c		write(*,*) i,j
		if (j.eq.0) then
			grid(gx,gy,gz)=i
			nocc=nocc+1
			occ(1,nocc)=gx
			occ(2,nocc)=gy
			occ(3,nocc)=gz
			hist(0)=hist(0)+1
		else
			n=0
10			k=overlay(j)
			n=n+1
			if (n.gt.m) m=n	
			if (k.eq.0) then
				overlay(j)=i
				hist(n)=hist(n)+1
			else
				j=k
				go to 10
			end if
		end if	
19	  end if		
	end do

	write(iout,40) m, nocc, natom, nx
40 	format (/' maximum overlap of atoms in one grid cell: ', i6,
     $	        /' number of occupied grid cells: ',i6,
     $		/' natom: ',i6, /' atoms outside grid: ', i6)
	s=0
	do i=0,m
		hist(i)=hist(i)-hist(i+1)
		s=s+hist(i)*(i+1)
	end do
	write(iout,*) ' sum: ',s
c	histogram of cell occupancy
c	write(iout,42) (hist(i),i=0,m)
c42	format(/' cell density histogram: ',/10(i8))
      
      return
      end

	function fung(x)
	real x
	integer fung

	fung=int(x/2)

	return
	end
c	
c     **********************************************************
c
	subroutine generate_cb(n,ca,c0,cb)
c
c 	a = ca to n; b = ca to c; c = right-vector; d = up-vector; 
c	cb = ca + 1.32 * right-vector + 0.765 * up-vector 
c	implicit none
	real n(3), ca(3), c0(3), cb(3)
	integer i
	real lc, ld, a(3), b(3), c(3), d(3)

	call diff(n,ca,a)
	call diff(c0,ca,b)
	call cross(a,b,c)
	lc=0
	ld=0
	do i=1,3
		lc=lc+c(i)*c(i)
		d(i)=a(i)+b(i)
		ld=ld+d(i)*d(i)
	end do
	lc=sqrt(lc)
	ld=sqrt(ld)
	if(lc.eq.0) lc=1
	if(ld.eq.0) ld=1
	do i=1,3
		c(i)=c(i)/lc
		d(i)=-d(i)/ld
	end do
	do i=1,3
		cb(i)=ca(i)+1.32*c(i)+0.765*d(i)
	end do

	return
	end

c	
c     **********************************************************
c
	subroutine getlibrary(iout, irotlib, nrotlib, rotkhi1, rotkhi2,
     $		rotatmnam, rotresnam, xyzrot, found, ierr, dfltfile)
c
c	reads rotamer library
c
	implicit none
	include 'sizes.for'

	character*20 dfltfile
	integer irotlib, nrotlib, ierr, iout
	integer rotkhi(maxrot, 5)
	character*3 rotresnam(maxrot)
	real xyzrot(3, maxrot, 14), rotkhi1(maxrot), rotkhi2(maxrot)
	logical found(maxrot, 14)
	character*3 rotatmnam(maxrot, 14)
c
	character*3 resname, atmnam(maxatm), resnam(maxres)
	integer nrot, nkhi, i, j, k, natm, rotnres, rotierr, irot, 
     $		  atmres(maxatm), ix1, ix2, ix3, ix4
	character start
	real xyz(3,maxatm),bval(maxatm)
	character*80 infile
	integer pdbresno(maxres)
c
c     open the data base file
c
	nrotlib=0
	ierr=0
	open(irotlib, err=99, status='old')
	goto 100
99	stop ' FATAL ERROR: cannot open rotamer library '	
c
c     read the data file
c
100	continue
	read(irotlib, 90, err=100, end=110) start, resname, nrot, nkhi
	if (start.eq.'>') then
		read(irotlib, 92, err=100, end=110) infile
c		fetch prototype rotamer coordinates
		call getcoor(99, iout, infile, xyz, natm, atmnam,
     $			atmres, rotnres,resnam,.false.,rotierr,bval,
     $ 			dfltfile,pdbresno)
		if (resname.eq.'GLY') then
			nrotlib=nrotlib+1
			rotresnam(nrotlib)=resname
		end if
		do irot=1,nrot
		  nrotlib=nrotlib+1
		  rotresnam(nrotlib)=resname
		  read(irotlib, 94) (rotkhi(irot,j),j=1,nkhi)
		  rotkhi1(nrotlib)=rotkhi(irot,1)
		  rotkhi2(nrotlib)=rotkhi(irot,2)
c		  put prototype coordinates to tables xyzrot, 
c			found, rotatmnam
		  do j=1,14
			found(nrotlib, j)=.false.
		  end do
		  do i=1,natm
			j=0
			if(atmnam(i).eq.'C  ') j=1
			if(atmnam(i).eq.'O  ') j=2
			if(atmnam(i).eq.'N  ') j=3
			if(atmnam(i).eq.'CA ') j=4
			if(atmnam(i).eq.'CB ') j=5
			if((atmnam(i).eq.'CG ').or.(atmnam(i).eq.'SG ').or.
     $			   (atmnam(i).eq.'OG ').or.
     $			   (atmnam(i).eq.'CG1').or.(atmnam(i).eq.'OG1')) j=6
			if((atmnam(i).eq.'CD ').or.(atmnam(i).eq.'OD1').or.
     $			   (atmnam(i).eq.'ND1').or.(atmnam(i).eq.'CD1').or.
     $			   (atmnam(i).eq.'SD ')) j=7
			if(j.eq.0) then
			  do j=8,14
				if (.not.found(nrotlib, j)) goto 80
			  end do
80			end if
			found(nrotlib, j)=.true.
			do k=1,3
				xyzrot(k, nrotlib, j)=xyz(k, i)
			end do
			rotatmnam(nrotlib, j)=atmnam(i)
		  end do
c		  set khi angles to those specified in rotamer library
c		  khi1=n-ca-cb-... khi2=ca-cb-cg-... khi3=cb-cg-cd-...
c		  table ix: n-3 ca-4 cb-5 'cg'-6 'cd'-7 '...'- up to 14
		  do j=1,nkhi
c	write(*,*) ' rotate ',nrotlib, rotresnam(nrotlib), irot, j,
c     $		   rotkhi(irot,j)
c
c		    dirty ILE patch: don't move CG2 in khi2 -- part 1
	  	    if ((rotresnam(nrotlib).eq.'ILE').and.(j.eq.2)) then
			found(nrotlib, 8)=.false.
		    end if
c		    end part 1 of dirty patch
		    ix1=2+j
		    ix2=3+j
		    ix3=4+j
		    ix4=5+j
		    call rotate(xyzrot,found,
     $			nrotlib,ix1, ix2, ix3, ix4 ,rotkhi(irot,j))
c		    dirty ILE patch: don't move CG2 in khi2 -- part2
	  	    if ((rotresnam(nrotlib).eq.'ILE').and.(j.eq.2)) then
			found(nrotlib, 8)=.true.
		    end if
c		    end part 2 of dirty patch
		  end do
		end do
	end if
	go to 100
110	close(irotlib)
	write(*,*) nrotlib, ' rotamers read from library'

c	check for too many rotamers
	  if (nrotlib.gt.maxrot) 
     $ 	  stop ' getlibrary -- Too many rotamers; increase MAXROT '

90	format (a1, a3, 2(i5))
92	format (a80)
94 	format (100i4)

999	return
      	end

c	
c     **********************************************************
c
	subroutine rotate(xyzrot, found, irot,
     $		ix1, ix2, ix3, ix4, targetangle)
	implicit none
	include 'sizes.for'
	real xyzrot(3, maxrot,14)
	integer targetangle, ix1,ix2,ix3, ix4, irot
	logical found(maxrot,14)
c
	real dot, dihedralangle
	real xyz1(3), xyz2(3), xyz3(3), xyz4(3), u0, angle, x(3),y(3),z(3),
     $ 		cacb(3), nca(3), a, b, c, lab, n, cosu, sinu, r, xyz(3)
	integer i,j

c	get initial angle, difference to target angle
	do i=1,3
		xyz1(i)=xyzrot(i,irot,ix1)
		xyz2(i)=xyzrot(i,irot,ix2)
		xyz3(i)=xyzrot(i,irot,ix3)
		xyz4(i)=xyzrot(i,irot,ix4)
	end do
	u0=dihedralangle(xyz1, xyz2, xyz3, xyz4)
	angle=(u0-targetangle)/180*3.14
c	write(*,*) ' rotate angle ',angle
c	coordinates around 'ca'-'cb'
	call diff(xyz2, xyz3, cacb)
	call diff(xyz1, xyz2, nca)
	lab=sqrt(dot(cacb, cacb))
	if (lab.eq.0) write(*,*) ' lab is zero'
	do i=1,3
		x(i)=cacb(i)/lab
	end do
	n=dot(nca, cacb)/lab
	do i=1,3
		y(i)=x(i)*n-nca(i)
	end do
	n=sqrt(dot(y,y))
	if (n.eq.0) write(*,*) ' n is zero'
	do i=1,3
		y(i)=y(i)/n
	end do
	call cross(x,y,z)
c	transform all atoms after ix4 inclusive
	do i=ix4,14
	  if (found(irot,i)) then
		do j=1,3
			xyz(j)=xyzrot(j,irot,i)-xyz2(j)
		end do
		a=dot(xyz, x)
		b=dot(xyz, y)
		c=dot(xyz, z)
		r=sqrt(b*b+c*c)
		u0=atan2(c,b)
		cosu=cos(angle+u0)
		sinu=sin(angle+u0)
		b=cosu*r
		c=sinu*r
		do j=1,3
			xyzrot(j,irot,i)=xyz2(j)+a*x(j)+b*y(j)+c*z(j)
		end do
	  end if
	end do
	return
	end
c	
c     **********************************************************
c
	subroutine putcoor(iout, xyz, natom, atmnam, atmres, resnam, 
     $		property)
c
c	write atom list out in Brookhaven format
c	first, sort atoms by residue number using j_index_Qsort
c	property is output in temp-factor column
c
	implicit none
	include 'sizes.for'
	real xyz(3, maxatm), property(maxatm)
	integer natom, atmres(maxatm), iout
	character*3 atmnam(maxatm), resnam(maxres)
c
	integer i, v(maxatm), code
	integer n, d(maxatm)

c	sort
	do i=1,natom
		d(i)=i
		v(i)=atmres(i)*14+code(atmnam(i))
	end do
	n=natom
	call j_index_Qsort(1, n, n, d, v)

c
c	write MaxSprout header !
c
	write(iout, 90) 
     $	'REMARK SIDECHAIN COORDINATES OPTIMIZED BY MAXSPROUT'
	write(iout, 90) 
     $	'REMARK REFERENCE: L.HOLM,C.SANDER (1991) J.MOL.BIOL.218:183-194'
	do i=1, natom
	  write(iout, 90) 'ATOM  ', i, atmnam(d(i)), resnam(atmres(d(i))),
     $		atmres(d(i)), xyz(1, d(i)), xyz(2,d(i)), xyz(3,d(i)),
     $		property(d(i))
	end do
   90   format (a6,i5,2x,a3,1x,a3,2x,i4,4x,3f8.3,6x,f6.2)

	return
	end

	function code(atmname)
	implicit none
	integer code
	character*3 atmname

	if(atmname.eq.'N  ') then
		code=1
	else if(atmname.eq.'CA ') then
		code=2
	else if(atmname.eq.'C  ') then
		code=3
	else if(atmname.eq.'O  ') then
		code=4
	else if(atmname.eq.'CB ') then
		code=5
	else if(atmname(2:2).eq.'G') then
		code=6
	else if(atmname(2:2).eq.'D') then
		code=7
	else if(atmname(2:2).eq.'E') then
		code=8
	else if(atmname(2:2).eq.'Z') then
		code=9
	else if(atmname(2:2).eq.'H') then
		code=10
	else
		code=11
	end if

	return
	end
c	
c     **********************************************************
c
	subroutine getcoor(ibhvn, iout, brkfile, xyz, natom, atmnam, 
     $		atmres, nres, resnam, verbose, ierr, bvalue, dfltfile,
     $		pdbresno)
c
c	read Brookhaven file
c	input: 	brkfile -- file name
c		ibhvn, iout -- i/o units
c		dfltfile -- open(..., defaultfile=dfltfile)
c	output: xyz -- coordinates
c		natom -- total # atoms
c		atmnam -- 'CA '
c		atmres -- residue index of atom
c		nres -- total # residues
c		resnam -- 'ALA'
c		ierr -- error code: 0 success; 1 file not opened
c		bvalue -- real
c		pdbresno -- original residue number
c
c	reads only standard amino acid names; N must be first atom in residue;
c	residues must be in one block; residues are numbered as read in;
c	hydrogens are ignored
c
      	implicit none
      	include 'sizes.for'
	character*80 brkfile
	character*20 dfltfile
	real xyz(3, maxatm), bvalue(maxatm)
	integer ibhvn
	integer natom, atmres(maxatm), nres, iout, ierr
	character*3 atmnam(maxatm), resnam(maxres)
	logical verbose
	integer pdbresno(maxres)
c	locals:
	character*60 record
	character*6 remark
	integer i, nig
	character*3 atmname, resname
	integer number
	real xx, yy, zz, bval
c
c     open the data base file
c
	ierr=1
	open (unit=ibhvn,file=brkfile,status='old',err=999)
	ierr=0
c
c     read the title, compound, source, etc. from data file
c
	if (verbose) then
	  do i=1,50
		read (ibhvn,70,err=60,end=65)  remark,record
      		if (remark.eq.'HEADER' .or. remark.eq.'COMPND' .or.
     $      	  remark.eq.'SOURCE' .or. remark.eq.'AUTHOR') then
            	  write (iout,80)  record
	        end if
60	  end do
	end if
65	rewind (unit=ibhvn)
70 	format (a6,4x,a60)
80	format (1x,a60)
c
c     initialize various values used below
c
	natom = 0
     	nres = 0
	nig=0
c
c     now, read individual atoms from the data file
c
      do i = 1, 30000
         read (ibhvn,90,err=100,end=110)  remark,atmname,resname,
     $                                    number,xx,yy,zz,bval
   90    format (a6,7x,a3,1x,a3,2x,i4,4x,3f8.3,6x,f6.2)
	 if ((remark.eq.'ATOM  ').and.(resname.ne.'GLY').and.
     $(resname.ne.'ALA').and.(resname.ne.'VAL').and.(resname.ne.'LEU')
     $.and.(resname.ne.'ILE').and.(resname.ne.'SER').and.
     $(resname.ne.'THR').and.(resname.ne.'CYS').and.
     $(resname.ne.'PRO').and.(resname.ne.'PHE').and.
     $(resname.ne.'TYR').and.(resname.ne.'TRP').and.
     $(resname.ne.'HIS').and.(resname.ne.'ASP').and.
     $(resname.ne.'ASN').and.(resname.ne.'GLU').and.
     $(resname.ne.'GLN').and.(resname.ne.'MET').and.
     $(resname.ne.'LYS').and.(resname.ne.'ARG').and.
     $(resname.ne.'UNK')) then
		if (nig.le.10) write(iout, 95) resname
		nig=nig+1
   95 		format(/' unknown residue type ', a3,' -- ignored ')
		go to 100
	end if
         if (remark.eq.'ATOM  ') then
            if (atmname .eq. 'N  ') then
               nres = nres + 1
               natom = natom + 1
               atmnam(natom) = atmname
               atmres(natom) = nres
               resnam(nres) = resname
	       pdbresno(nres) = number
               xyz(1,natom) = xx
               xyz(2,natom) = yy
               xyz(3,natom) = zz
	       bvalue(natom) = bval
            else if (atmname(1:1) .ne. 'H') then
               natom = natom + 1
               atmnam(natom) = atmname
               atmres(natom) = nres
               xyz(1,natom) = xx
               xyz(2,natom) = yy
               xyz(3,natom) = zz
	       bvalue(natom) = bval
	    else if (atmname(1:1).eq.'H') then
		nig=nig+1
            end if
         end if
  100    continue
      end do
  110 continue
      close (unit=ibhvn)
c
c     check for too many atoms or residues
c
c      if (natom .gt. maxatm)
c     $   stop  ' getcoor --  Too many atoms; increase MAXATM'
c      if (nres .gt. maxres)
c     $   stop  ' getcoor --  Too many residues; increase MAXRES'
c
c     print summary info on atoms read from the data file
c
      if (verbose) write (iout,120)  nres,natom, nig
  120 format (/' Summary of Atom Information read from Data File :',/
     $        /' Number of standard peptide units :  ',i6,
     $        /' Total number of atoms read in :     ',i6,
     $	      /' Number of ignored ATOM records :    ',i6)
c
c	error report
c
999	if (ierr.eq.1) write(iout, 130) brkfile
130	format (/' File not found : ', a80)

	return
      	end
c	
c     **********************************************************
c
	subroutine getstatic(natom, atmnam, atmres, staticside, staticmain,
     $			 nstaticatm, staticatm)
	implicit none
	include 'sizes.for'
	integer natom, nstaticatm, staticatm(maxatm), atmres(maxatm)
	character*3 atmnam(maxatm)
	logical staticmain(maxres), staticside(maxres)
	integer i

	nstaticatm=0
	do i=1,natom
		if ((atmnam(i).eq.'N  ').or.(atmnam(i).eq.'CA ').or.
     $		    (atmnam(i).eq.'CB ').or.
     $		    (atmnam(i).eq.'C  ').or.(atmnam(i).eq.'O  ')) then
			if (staticmain(atmres(i))) then
				nstaticatm=nstaticatm+1
				staticatm(nstaticatm)=i
			end if
		else if (staticside(atmres(i))) then
				nstaticatm=nstaticatm+1
				staticatm(nstaticatm)=i
		end if
	end do
	
	return
	end

c	
c     **********************************************************
c
	subroutine compress_static(xyz,natom, nstaticatm, staticatm,
     $					atmnam, atmres, firstatm)
	implicit none
	include 'sizes.for'
	integer i,j, natom, nstaticatm, staticatm(maxatm), atmres(maxatm),
     $		  firstatm(maxres)
	real xyz(3,maxatm)
	character*3 atmnam(maxatm)
c
c	leaves only static atoms in the x,y,z arrays; natom=nstaticatm
c
	write(*,*) nstaticatm, ' nstaticatm'
	do i=1,nstaticatm
		do j=1,3
			xyz(j, i)=xyz(j, staticatm(i))
		end do
		atmnam(i)=atmnam(staticatm(i))
		atmres(i)=atmres(staticatm(i))
		if(firstatm(atmres(i)).eq.0) firstatm(atmres(i))=i
		staticatm(i)=i
	end do
	natom=nstaticatm

	return
	end

c	
c     **********************************************************
c
	subroutine diff(a, b, c)
	implicit none
	real a(3), b(3), c(3)
	integer i

	do i=1,3
		c(i)=a(i)-b(i)
	end do	
	
	return
	end
c =========================================================================
	function dot(a, b)
	implicit none
	real dot, a(3), b(3)

	dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

	return
	end

c =========================================================================
	subroutine cross(a, b, c)
	implicit none
	real a(3), b(3), c(3)
	
	c(1)=a(2)*b(3)-b(2)*a(3)
	c(2)=a(3)*b(1)-a(1)*b(3)
	c(3)=a(1)*b(2)-a(2)*b(1)

	return
	end

c =========================================================================
	function atan2(y,x)
	implicit none
	real y, x, z, atan2

	if (x.ne.0.0) then
		z=atan(y/x)
	else if (y.gt.0.0) then
		z=1.570796
	else if (y.lt.0.0) then
		z=-1.570796
	else
		z=6.283185
	end if
	if (x.lt.0.0) then
	  if (y.gt.0.0) then
		z=z+3.141593
	  else
		z=z-3.141593
	  end if
	end if
	atan2=z

	return
	end
c =========================================================================
	function dihedralangle(v1, v2, v3, v4)
	implicit none
	real 	u,v, v12(3), v43(3), x(3), y(3), z(3), p(3),
     $		v1(3), v2(3), v3(3), v4(3), dihedralangle, dot

	call diff(v1, v2, v12)
	call diff(v4, v3, v43)
	call diff(v2, v3, z)
	call cross(z, v12, p)
	call cross(z, v43, x)
	call cross(z, x, y)
	u=dot(x,x)
	v=dot(y,y)
	dihedralangle=360.0
	if ((u.gt.0.0).and.(v.gt.0.0)) then
		u=dot(p,x)/sqrt(u)
		v=dot(p,y)/sqrt(v)
		if ((u.ne.0.0).or.(v.ne.0.0)) then
			dihedralangle=atan2(v,u)*57.29578
		end if
	end if

	return
	end

c	
c     **********************************************************
c
c	
c     **********************************************************
c
	subroutine singlet_freq(nresrot,freq,rejected)
	implicit none
	integer nresrot, freq(nresrot)
	logical rejected(nresrot)
c
	integer i

	do i=1,nresrot
		freq(i)=0
		if (.not.rejected(i)) freq(i)=1
	end do
	write(*,*) ' singlet bootstrap frequencies'
	write(*,*) (freq(i),i=1,nresrot)

	return
	end
c	
c     **********************************************************
c
	subroutine triplet_freq(nresrot,freq,rejected,npair,ipair,jpair,epair,
     $		i_first, w0, w1)
	implicit none
	integer nresrot, freq(nresrot)
	logical rejected(nresrot)
	integer npair, i_first(nresrot)
	integer ipair(npair), jpair(npair), w0
	real epair(npair), w1
c
	integer i,a,b,c
	integer k,l,m
	real e, w
c
c	ipair A A ... B
c	jpair B C ... C
c	epair * * ... *
c
c	if B>A then if e(A,B)<0 then
c		if C>B then if e(A,C)<0 then
c			if (B,C) found then if e(B,C)<0 then
c				e=e(A,B)+e(A,C)+e(B,C)
c				w=exp(-alfa*e)
c				increment weights of A,B,C by w
c
	do i=1,nresrot
		freq(i)=0
		if (.not.rejected(i)) freq(i)=w0
	end do
	do k=1,npair
	  a=ipair(k)
	  b=jpair(k)
	  if((b.gt.a).and.(epair(k).lt.0.0)) then
	    l=k+1
	    do while((ipair(l).eq.a).and.(l.le.npair))
	      c=jpair(l)
	      if((c.gt.b).and.(epair(l).lt.0.0)) then
c		find m, ipair(m)=b, jpair(m)=c
		m=i_first(b)
		do while((ipair(m).eq.b).and.(m.le.npair))
		  if(jpair(m).eq.c) then
		    if(epair(m).lt.0.0) then
		      e=(epair(k)+epair(l)+epair(m))/3
		      w=exp(-w1*e)
		      if(w.gt.100) w=100
		      freq(a)=freq(a)+w
		      freq(b)=freq(b)+w
		      freq(c)=freq(c)+w
		    end if
		    goto 111
		  end if
		  m=m+1
		end do
	      end if
111	      l=l+1
	    end do
	  end if
	end do
	write(*,*) ' triplet bootstrap frequencies'
	write(*,*) (freq(i),i=1,nresrot)

	return
	end
c	
c     **********************************************************
c
	subroutine doublet_freq(nresrot,freq,rejected,npair,ipair,jpair,epair,
     $		w0, w1)
	implicit none
	include 'sizes.for'
	integer nresrot, freq(nresrot)
	logical rejected(nresrot)
	integer npair
	integer ipair(npair), jpair(npair), w0
	real epair(npair), w1
c
	integer i, n(maxres)
	integer k
	real we(maxres)

	do i=1,nresrot
		freq(i)=0
		if (.not.rejected(i)) freq(i)=w0
		we(i)=0.0
		n(i)=0
	end do
	do k=1,npair
	  if(epair(k).lt.0.0) then
		we(ipair(k))=we(ipair(k))+epair(k)
		we(jpair(k))=we(jpair(k))+epair(k)
		n(ipair(k))=n(ipair(k))+1
		n(jpair(k))=n(jpair(k))+1
	  end if
	end do
	do i=1,nresrot
	  if(n(i).gt.0) freq(i)=freq(i)+int(exp(-we(i)/3/n(i)))
	end do

	write(*,*) ' doublet bootstrap frequencies '
	write(*,*) (freq(i), i=1,nresrot)

	return
	end
c	
c     **********************************************************
c
	subroutine bootstrap_freq(freq, wrot, wrottot, nres, nirot, 
     $		ires_to_irot, uniq, frozen, nresrot, rejected)
	implicit none
	include 'sizes.for'
	integer wrottot
	integer nres, nresrot, freq(nresrot), wrot(maxres*100),
     $		nirot(nres), ires_to_irot(nres)
	logical uniq(nres), frozen, rejected(nresrot)
c
	integer i,j,jrot,n
	integer k, s
	real norm 

	wrottot=0
	frozen=.true.
c	write(*,*) ' rotamer frequencies '
c	write(*,*) (freq(i),i=1,nresrot)
c
c	normalize freq(irot) per nresrot*100
c
	s=0
	do i=1,nresrot
	  s=s+freq(i)
 	end do
	norm=100.0/float(s)*nresrot
	do i=1,nres
          if(uniq(i)) goto 10
	  n=0
	  do j=1,nirot(i)
	    jrot=ires_to_irot(i)+j-1
	    if(freq(jrot).gt.0) n=n+1
	  end do
	  if(n.le.1) then
	    uniq(i)=.true.
	    goto 10
	  else
	    frozen=.false.
	    do j=1,nirot(i)
	      jrot=ires_to_irot(i)+j-1
	      freq(jrot)=int(norm*float(freq(jrot)))
	      do k=1,freq(jrot)
	   	wrottot=wrottot+1
		wrot(wrottot)=jrot
	      end do
	    end do
	  end if
10	end do

c	zero freq !
	do i=1,nresrot
		freq(i)=0
		if(.not.rejected(i)) freq(i)=1
	end do

	return
	end	
c	
c     **********************************************************
c
	subroutine update_freq(freq, wrot, wrottot, nres, nirot, ires_to_irot,
     $		uniq, frozen, nresrot, eres, rejected)
	implicit none
	include 'sizes.for'
	integer wrottot
	integer nres, nresrot, freq(nresrot), wrot(maxres*100),
     $		nirot(nres), ires_to_irot(nres)
	logical uniq(nres), frozen
	logical rejected(nresrot)
	real eres(nres)
c
	integer i,j,jrot,n
	integer k
	real norm, s
c
c	weight by eres(ires) !!!
c
	wrottot=0
	frozen=.true.
	do i=1,nres
	  norm=100.0
c	  if(eres(i).gt.10.0) then
c		norm=500
c	  else
c		norm=50
c	  end if
          if(uniq(i)) goto 10
	  n=0
	  s=0
	  do j=1,nirot(i)
	    jrot=ires_to_irot(i)+j-1
	    s=s+freq(jrot)
	    if(freq(jrot).gt.0) n=n+1
	  end do
	  if(n.le.1) then
	    uniq(i)=.true.
	    goto 10
	  else
	    frozen=.false.
	    norm=norm/s
	    do j=1,nirot(i)
	      jrot=ires_to_irot(i)+j-1
	      freq(jrot)=int(float(freq(jrot))*norm)
	      do k=1,freq(jrot)
	   	wrottot=wrottot+1
		wrot(wrottot)=jrot
	      end do
	    end do
	  end if
10	end do

c	write(*,*) (i, uniq(i),i=1,nres)

c	zero freq !
	do i=1,nresrot
		if(rejected(i)) then
			freq(i)=0
		else
			freq(i)=1
		end if
	end do

	return
	end	
c	
c     **********************************************************
c
	function khipercent(nres, nresrot, rotamer, khi1, rotkhi1, 
     $		rotkhi2, rottolib)
	implicit none
	include 'sizes.for'
	integer nres, nresrot, rotamer(nres), rottolib(nresrot)
	real khi1(nres), rotkhi1(maxrot), rotkhi2(maxrot), khipercent
c
	integer ires, irot, nok, n
	real dkhi

	nok=0
	n=0
	do ires=1,nres
		irot=rotamer(ires)
c		don't want alanines or glycines 
c			nor prolines
		if((irot.gt.0).and.(khi1(ires).ne.0.0)
     $			.and.(rotkhi1(rottolib(irot)).ne.0.0)) then
			n=n+1
			dkhi=abs(rotkhi1(rottolib(irot))-khi1(ires))
			dkhi=min(abs(dkhi), abs(dkhi+360.0), abs(dkhi-360.0))
c	if(dkhi.gt.30.0) write(*,*) ires, irot, rottolib(irot), 
c     $		rotkhi1(rottolib(irot)), khi1(ires), dkhi
			if(dkhi.le.30.0) nok=nok+1
		end if
	end do
	write(*,*) nok, n, nres, nresrot
	if(n.gt.0) khipercent=100*nok/n
	if(n.eq.0) khipercent=-100.0

	return
	end
c	
c     **********************************************************
c
	subroutine protomcmsolv(npair, j_index, nres, nresrot, rotamer,
     $		rejected, i_first, ni, rottores, ipair, jpair, nirot, bestemin,
     $		eres, etot, epair, estatic, jactive, ires_to_irot,
     $		activeatm, activerot, nstaticatm, natom, firstocc, occstatic,
     $		niocc, occjrot, occs, occ_index, maxocc, minocc, solvatomtype,
     $		gsol, gtrf, atmres, resrot, firstrotatm, occiatm,
     $		occjrot_index, njrot, firstjrot, noccs, solvate,
     $		khi1, khi2, rotkhi1, rotkhi2, rottolib)
c
c	prototype
c
c	Monte Carlo with native sequence
c	no rotamer preferences
c	no residue preferences
c
c	include solvation (logical solvate)
c
	implicit none
	include 'sizes.for'
	integer npair, j_index(npair), noccs
	integer natom, nres, nresrot, rotamer(nres), ni(nresrot), 
     $		rottores(nresrot), ipair(npair), jpair(npair), nirot(nres),
     $		ires_to_irot(nres), atmres(natom), resrot(nresrot),
     $		firstrotatm(natom), nstaticatm
	integer i_first(nresrot)
	real bestemin, eres(nres), etot, epair(npair), estatic(maxres)
	logical jactive(npair), rejected(maxres)
c
	integer occ_index(noccs), niocc(natom), firstocc(natom),
     $		occjrot_index(noccs),  firstjrot(nresrot), njrot(nresrot)
	integer occjrot(noccs),
     $		solvatomtype(natom), occiatm(noccs)
	real occstatic(natom), occs(noccs), maxocc(1000), minocc(1000),
     $		esoltot, etrftot, gsol(1000), gtrf(1000)
	logical activeatm(natom), activerot(nresrot)
c
	real khi1(nres), khi2(nres), rotkhi1(maxrot), rotkhi2(maxrot),
     $		khipercent
	integer rottolib(nresrot)
c
	integer i, j, rot1, rot2, ires, nprint, irun, nrun,
     $		best(maxres), ndivide, improve, iatmlist(maxatm),niatm,
     $		curbest(maxres), freq(maxres), updatefreq
	real emin, p, alfa, alfamax, x, x1, delta, 
     $		esol(maxatm), etrf(maxatm), cs(maxatm),
     $		csxx(maxatm), dsol, delta0,alfamax0,alfa0
	integer seed, k
	logical solvate
	integer tslp(maxatm), nrestore
c
	integer wrot(maxres*100), w0, bootstrap, felsenstein,
     $		start0, start1
	integer wrottot, niter, steplimit, naccept, nup, nchecksim
	logical uniq(maxres), frozen
	real d0, w1, cutsim, solscale
c
	real ran
c       vsi: real for khipercent call
        real khiresult
c
c	initialize global variables from for081
c
	call initialize_global(solvate, natom, nres, seed, w0, w1,
     $		dsol, frozen, tslp, uniq, nrun, ndivide, nrestore, delta0,
     $		alfa0, alfamax0, nprint, updatefreq, steplimit, bestemin, esol,
     $		bootstrap, felsenstein, start0, start1, cutsim, nchecksim,
     $		solscale)
        write (*,*) 'initialize_global finished'
c
c	initialize rotamer weights
c
c	initialize with equal rotamer weights
	  if(bootstrap.eq.1) call singlet_freq(nresrot, freq, rejected)
c	get initial rotamer frequencies from 1+exp(epair)
	  if(bootstrap.eq.2) call doublet_freq(nresrot,freq,rejected,npair,
     $		ipair,jpair,epair,w0,w1)
	  if(bootstrap.eq.3) call triplet_freq(nresrot,freq,rejected,npair,
     $		ipair,jpair,epair,i_first,w0,w1)
	  call bootstrap_freq(freq, wrot, wrottot, nres, 
     $		nirot, ires_to_irot, uniq, frozen, nresrot, rejected)
c
c	do nrun runs
c
	irun=0
	do while ((irun.lt.nrun).and.(.not.frozen))
	  irun=irun+1
c
c	  initialize variables for the run
c
	  call initialize_local(irun, nrestore, delta0, alfamax0, ndivide,
     $		nstaticatm, nres, nirot, ires_to_irot, rejected, npair, 
     $		i_first, ni, j_index, natom, seed, alfa0, nresrot, 
     $		delta, alfamax, alfa, niter, naccept, improve, activeatm,
     $		rotamer, curbest, jactive, nup, start0, start1)
          write (*,*) 'initialize_local finished'
	  d0=delta
c	  get initial energy
	  call initialenergy(etot, esoltot, emin, eres,  etrftot,esol,etrf,
     $		solvate, nres, rotamer, estatic, ni, i_first, jactive, epair,
     $		nstaticatm, natom, resrot, atmres, activeatm, activerot, npair,
     $		firstocc, occstatic, niocc, occjrot, occs, occ_index, 
     $		maxocc, minocc, solvatomtype, nresrot, gsol, gtrf, cs, noccs,
     $		solscale)
          write (*,*) 'initialenergy finished'
c
c	  iterate until frozen or alfa equals alfamax
c
	  do while ((alfa.le.alfamax).and.(niter-improve.lt.steplimit))
	    alfa=alfa+delta
	    niter=niter+1
c
c	    every nchecksim steps check that rotamer() is further than
c	      cutsim from best() -- if irun.gt.1
c
	    if((irun.gt.1).and.(mod(niter, nchecksim).eq.0)) then
		j=0
		do i=1,nres
			if(rotamer(i).eq.best(i)) j=j+1
		end do
		p=float(j)/float(nres)
		write(*,*) p, j, niter, irun, alfa
		if(p.gt.cutsim) goto 199
	    end if
c
c	    get test rotamer rot2, corresponding old rotamer rot1
c	      loop here until rot2 is really new and not rejected !!!
c	      max 50 trials
c
	    j=0
10	    continue
	    rot2=wrot(int(ran(seed)*wrottot)+1)
	    j=j+1
	    if(j.gt.50) goto 110
	    if(rejected(rot2)) goto 10
	    ires=rottores(rot2)
	    rot1=rotamer(ires)
	    if(rot2.eq.rot1) goto 10
c	    get energy x for test rotamer rot2
	    x=estatic(rot2)
	    do j=0, ni(rot2)-1
	      k=i_first(rot2)+j
	      if(jactive(k)) x=x+epair(k)
	    end do
c	    if(solvate) then
c	      dsol is new-old
c	      call testsolv(rot1, rot2, dsol, dtrf, iatmlist, niatm, cs, csxx,
c     $		natom, occjrot_index, njrot, occiatm, firstjrot,
c     $		nresrot, noccs, occstatic, firstocc, occ_index, niocc,
c     $		activerot, solvatomtype, maxocc, minocc, gsol, gtrf, occs,
c     $		occjrot, activeatm, resrot, atmres, ires, firstrotatm, tslp)
c	      write(*,*) ' dsol ', dsol
c	    end if
c
c	    Metropolis criterion: p=1/exp(x1)
c	    Felsenstein criterion: p=1/(1+exp(x1))
c		
	    x1=-alfa*(eres(ires)-x-dsol)
	    if(x1.gt.50) then
	      p=0.0
	    else
	      p=1/(felsenstein+exp(x1))
	    end if
c
c	    accept change
c	    p is probability of acceptance
c
	    if(ran(seed).le.p) then
	      call changerotamer(rot1, rot2, eres, etot, jactive, i_first, ni,
     $		j_index, rottores, nresrot, epair, ipair, jpair, estatic, 
     $		nres, npair)
c		write(*,*) rot1,rot2,' changed'
	      call accept(rotamer, ires, rot2, eres, x, esoltot, dsol, emin,
     $		nres, curbest, improve, niter, bestemin, best, solvate,
     $		firstrotatm, atmres, natom, resrot, rot1, activeatm, activerot,
     $		niatm, iatmlist, cs, csxx, naccept, x1, nup, nresrot, etot,
     $		solscale)
c	      log only accepted steps !
c	       if(mod(naccept, nprint).eq.0) write(13,131) 
c     $		niter,char(9),alfa,char(9),etot,char(9),emin
c     $		, char(9), khipercent(nres, nresrot, rotamer, khi1, rotkhi1, 
c     $			rotkhi2, rottolib)
c     $		,char(9),esoltot,char(9),etot+esoltot
	      if(mod(naccept, nprint).eq.0) then
                khiresult = khipercent(nres, nresrot, rotamer, khi1,
     $			rotkhi1, rotkhi2, rottolib)
                write(13,131) 
     $		niter,char(9),alfa,char(9),etot,char(9),emin
     $		, char(9), khiresult, char(9),esoltot,char(9),etot+esoltot
              end if
131	      format(i10,6(a1,f7.2))
	    end if
110	  end do
c	  end of iteration
c	  update rotamer freq from curbest !!! during an optimization run !!!
	  do i=1,nres
		j=curbest(i)
		if (j.gt.0) freq(j)=freq(j)+1
	  end do
c	  write final energy to for013 and acceptance ratio to for014.dat
c	  write final rotamer() to for084 (binary file)
c         vsi: separated khipercent call
          khiresult = khipercent(nres, nresrot, curbest, khi1,
     $			rotkhi1, rotkhi2, rottolib)
199	  write(14, 144) irun, char(9), d0, char(9), naccept, char(9),
     $		nup, char(9), niter, char(9), emin, char(9)
     $		, khiresult
c199	  write(14, 144) irun, char(9), d0, char(9), naccept, char(9),
c     $		nup, char(9), niter, char(9), emin, char(9)
c     $		, khipercent(nres, nresrot, curbest, khi1, rotkhi1, rotkhi2, 
c     $			rottolib)
144	  format(i5,a1,f20.10,2(a1,i6),a1,i10,2(a1,f10.3))
	  write(84) emin, nres, (curbest(i),i=1,nres)
c	  now and then update wrot from freq; zero freq array
	  if(mod(irun,updatefreq).eq.0) 
     $	    call update_freq(freq, wrot, wrottot, nres, nirot, 
     $		ires_to_irot, uniq, frozen, nresrot, eres, rejected)	
	end do
c	end of optimization run

c	return best configuration in rotamer()
	do i=1,nres
	  rotamer(i)=best(i)
	end do

	return
	end
c	
c     **********************************************************
c
	subroutine initialize_global(solvate, natom, nres, seed, w0, w1,
     $		dsol, frozen, tslp, uniq, nrun, ndivide, nrestore, delta0,
     $		alfa0,alfamax0, nprint, updatefreq, steplimit, bestemin, esol,
     $		bootstrap, felsenstein, start0, start1, cutsim, nchecksim,
     $		solscale)
	implicit none
	integer natom, nres, tslp(natom), nrun, ndivide, nrestore, nprint,
     $		updatefreq, w0, bootstrap, felsenstein,
     $		start0, start1
	integer seed, steplimit, nchecksim
	real dsol, delta0, alfa0, alfamax0, bestemin,esol(natom),w1,cutsim,
     $		solscale
	logical solvate, frozen, uniq(nres)
c
	integer i, ierr
	character c
	
	dsol=0.0
	frozen=.false.
	do i=1,natom
		tslp(i)=0
	end do
	do i=1,nres
		uniq(i)=.false.
	end do
c	
c	read parameters from for081
c
	ierr=-1
	open(81, status='old', err=185)
	ierr=0
	c=' '
	do while (c.ne.'#')
		read(81, 181) c
	end do
181	format(a1)
	read(81, 182) seed
	read(81, 182) nrun
	read(81, 182) ndivide
	read(81, 182) nrestore
	read(81, 183) delta0
	read(81, 183) alfa0
	read(81, 183) alfamax0
	read(81, 182) nprint
	read(81, 182) updatefreq
	read(81, 182) steplimit
	read(81, 182) felsenstein
	read(81, 182) bootstrap
	read(81, 182) w0
	read(81, 183) w1
	read(81, 182) start0
	read(81, 182) start1
	read(81, 183) cutsim
	read(81, 182) nchecksim
	read(81, 183) solscale
182	format(i20)
183	format(f20.1)
185	close(81)
	if (ierr.eq.-1) 
     $	  write(*,*) ' unable to open Monte Carlo parameters file '
	bestemin=999999999.9999
	do i=1,natom
		esol(i)=0.0
	end do

	return
	end
c	
c     **********************************************************
c
	subroutine initialize_local(irun, nrestore, delta0, alfamax0, ndivide,
     $		nstaticatm, nres, nirot, ires_to_irot, rejected, npair, 
     $		i_first, ni, j_index, natom, seed, alfa0,nresrot,
     $		delta, alfamax, alfa, niter, naccept, improve, activeatm,
     $		rotamer, curbest, jactive, nup, start0, start1)
	implicit none
	integer irun, nrestore, ndivide, nstaticatm, nres, nirot(nres),
     $		nresrot, ires_to_irot(nres), ni(nresrot),natom,start0,start1
	integer niter, npair, i_first(nresrot), j_index(npair), seed
	logical rejected(nresrot)
	real delta0, alfamax0, alfa0
c	
	real delta, alfamax, alfa
	integer improve, rotamer(nres), curbest(nres)
	logical activeatm(natom), jactive(npair)
	integer nup, naccept
c
	integer i, j, irot 
	integer k
c
	real ran
        integer n, list(100), random

c	  initialize variables for the run
	  if(mod(irun,nrestore).eq.1) then
		delta=delta0
		alfamax=alfamax0
	  end if
	  if(mod(irun,ndivide).eq.0) then
	    alfamax=alfamax/10
	    delta=delta/10
	  end if
	  alfa=alfa0
	  niter=0
	  naccept=0
	  nup=0
	  improve=0
	  do i=1,nstaticatm
		activeatm(i)=.true.
	  end do
c         vsi: initialize list
          do i=1,100
                list(i)=0
          end do

          write(*,*) ' initialize_local nirot=', nirot

c	  write(*,*) ' run ',irun,' of ',nrun,' with delta, alfamax', delta,
c     $			alfamax
c	  select random starting configuration and initialize jactive
c	  check that rotamers are not rejected
c	  curbest is best combination for each optimization run
	  if(((irun.eq.1).and.(start0.eq.0)).or.
     $		((irun.gt.1).and.(start1.eq.0))) then
            do i=1,nres
              rotamer(i)=1
              if(nirot(i).gt.0) then
c               shortlist acceptable rotamers
                n=0
                do j=1,nirot(i)
                  if(.not.rejected(rotamer(i))) then
                    n=n+1
                    list(n)=j
                  end if
                end do
c                write(*,*) ' initialize_local list=', list
                random = int(ran(seed)*n)+1
c                write (*,*) 'initialize_local n=', n, ' random=', random
                if(n.gt.0) then
                  rotamer(i)=list(random)+ires_to_irot(i)
c                  write (*,*) 'i=', i
c                  write (*,*) 'list(random)=', list(random)
c                  write (*,*) 'ires_to_irot(i)=', ires_to_irot(i)
c                  write (*,*) 'rotamer(i)=', rotamer(i)
                end if
              end if
              curbest(i)=rotamer(i)
            end do
          else
c          vsi: old version of above if, infinite loop in Linux
c	   if(((irun.eq.1).and.(start0.eq.0)).or.
c     $		((irun.gt.1).and.(start1.eq.0))) then
c	     do i=1,nres
c	       if(nirot(i).eq.0) then
c		 rotamer(i)=1
c	       else
cc		! this could be a never-ending loop !
c12	        rotamer(i)=int(ran(seed)*nirot(i))+ires_to_irot(i)
c		 if(rejected(rotamer(i))) goto 12
c	       end if
c	       curbest(i)=rotamer(i)
c	     end do
c	   else
	    do i=1,nres
	      rotamer(i)=curbest(i)
	    end do
	  end if
	  do k=1,npair
	    jactive(k)=.false.
	  end do
	  do i=1,nres
	    irot=rotamer(i)
	    if(irot.gt.0) then
	      do k=i_first(irot),i_first(irot)+ni(irot)-1
		jactive(j_index(k))=.true.
	      end do
	    end if
	  end do
	
	return
	end
c	
c     **********************************************************
c
	subroutine initialenergy(etot, esoltot, emin, eres, etrftot,esol,etrf,
     $		solvate, nres, rotamer, estatic, ni, i_first, jactive, epair,
     $		nstaticatm, natom, resrot, atmres, activeatm, activerot, npair,
     $		firstocc, occstatic, niocc, occjrot, occs, occ_index, 
     $		maxocc, minocc, solvatomtype, nresrot, gsol, gtrf, cs, noccs,
     $		solscale)
	implicit none
	include 'sizes.for'
	integer noccs
	integer nres, rotamer(nres), nresrot, ni(nresrot),nstaticatm, natom,
     $		resrot(nresrot), atmres(natom), occjrot(noccs),
     $		solvatomtype(natom)
	integer i_first(nresrot), npair, firstocc(natom), niocc(natom),
     $		occ_index(noccs)
	real etot, esoltot, emin, eres(nres), estatic(maxres), epair(npair),
     $		occstatic(natom), occs(noccs), cs(natom), etrftot, esol(natom),
     $		maxocc(1000), minocc(1000), gsol(1000), gtrf(1000),etrf(natom),
     $		solscale
	logical solvate, jactive(npair), activeatm(natom), activerot(nresrot)
c
	integer i, irot, j
	integer k

c	  get initial energy
	  etot=0.0
	  do i=1,nres
	    irot=rotamer(i)
	    eres(i)=estatic(irot)
	    do j=0, ni(irot)-1
	      k=i_first(irot)+j
	      if(jactive(k)) eres(i)=eres(i)+epair(k)
	    end do
	    etot=etot+eres(i)
	  end do	
c	  if(solvate) then
c	    get initial solvation energy and atomic occupancies cs(iatm)
c	    do i=nstaticatm+1,natom
c	      if(resrot(i).eq.rotamer(atmres(i))) then
c		activeatm(i)=.true.
c	      else
c		activeatm(i)=.false.
c	      end if
c	    end do
c	    do i=1,nresrot
c		activerot(i)=.false.
c	    end do
c	    do i=1,nres
c		if(rotamer(i).gt.0) activerot(rotamer(i))=.true.
c	    end do
c	    call getesolv(natom, activeatm, activerot, firstocc, occstatic,
c    $		niocc, occjrot, occs, occ_index, maxocc, minocc, solvatomtype,
c     $		nresrot, esoltot, etrftot, gsol, gtrf, esol, etrf, cs, noccs,
c     $		solscale)
c	  else
		esoltot=0.0
c	  end if
	  emin=etot+esoltot
	  write(*,*) ' initial energy ',etot, esoltot, etrftot

	return
	end
c	
c     **********************************************************
c
	subroutine accept(rotamer, ires, rot2, eres, x, esoltot, dsol, emin,
     $		nres, curbest, improve, niter, bestemin, best, solvate,
     $		firstrotatm, atmres, natom, resrot, rot1, activeatm, activerot,
     $		niatm, iatmlist, cs, csxx, naccept, x1, nup, nresrot, etot,
     $		solscale)
	implicit none
	integer nres, rotamer(nres), ires, rot2, curbest(nres), improve,
     $		natom, niter, best(nres), firstrotatm(natom), atmres(natom),
     $		nresrot, resrot(nresrot), rot1, niatm, iatmlist(natom)
	integer naccept, nup
	logical solvate, activeatm(natom), activerot(nresrot)
	real eres(nres), x, esoltot, dsol, emin, bestemin, cs(natom),
     $		csxx(natom), etot, x1, solscale
c
	integer i
c
c	      accept change -- remember to update rotamer(), esoltot, etrftot
c
	      rotamer(ires)=rot2
	      eres(ires)=x
	      if(x1.gt.0.0) nup=nup+1
	      esoltot=esoltot+solscale*dsol
c	      save curbest
	      if(etot+esoltot.lt.emin) then
		emin=etot+esoltot
		do i=1,nres
			curbest(i)=rotamer(i)
		end do
		improve=niter
	      end if
	      if(emin.lt.bestemin) then
		bestemin=emin
		do i=1,nres
		  best(i)=rotamer(i)
		end do
	      end if
c	      if(solvate) then
c		i=firstrotatm(ires)
c		do while((atmres(i).eq.ires).and.(i.le.natom))
c			if(resrot(i).eq.rot1) activeatm(i)=.false.
c			if(resrot(i).eq.rot2) activeatm(i)=.true.
c			i=i+1
c		end do
c		activerot(rot1)=.false.
c		activerot(rot2)=.true.
c		replace cs with csx for atoms in iatmlist
c		do i=1,niatm
c			iatm=iatmlist(i)
c			cs(iatm)=csxx(iatm)
c		end do
c	      end if
	      naccept=naccept+1

	return
	end
c	
c     **********************************************************
c
	subroutine changerotamer(irot1, irot2, eres, etot, jactive,
     $		i_first, ni, j_index, rottores, nresrot,
     $		epair, ipair, jpair, estatic, nres, npair)
c
c	old rotamer is irot1, new rotamer is irot2
c	subtract all pair interactions of irot1 -- and its estatic
c	add all pair interactions of irot2 -- and its estatic
c	
c	eres(ires) is per-residue energy
c	etot is total energy
c	jactive says which rotamers are active and has the same index as jpair
c	i_first is a pointer irot --> ipair index
c	ni goes with i_first
c	j_index contains all jpair addresses ordered by irot
c
	implicit none
	include 'sizes.for'
	integer npair, j_index(npair)
	integer irot1, irot2, nres, nresrot,
     $		ni(nresrot), rottores(nresrot),
     $		ipair(npair), jpair(npair)
	integer i_first(nresrot)
	real eres(nres), etot, epair(npair), estatic(maxres)
	logical jactive(npair)
c
	integer k
	integer i, ires, jres
	real de

c
c	subtract irot1
c
	do i=0,ni(irot1)-1
		k=i_first(irot1)+i
		if(jactive(k)) then
			de=epair(k)
			etot=etot-de-de
			ires=rottores(ipair(k))
			jres=rottores(jpair(k))
			eres(ires)=eres(ires)-de
			eres(jres)=eres(jres)-de
		end if
	end do	
	de=estatic(irot1)
	etot=etot-de
	ires=rottores(irot1)
	eres(ires)=eres(ires)-de
c
c	add irot2
c
	do i=0,ni(irot2)-1
		k=i_first(irot2)+i
		if(jactive(k)) then
			de=epair(k)
			etot=etot+de+de
			ires=rottores(ipair(k))
			jres=rottores(jpair(k))
			eres(ires)=eres(ires)+de
			eres(jres)=eres(jres)+de
		end if
	end do	
	de=estatic(irot2)
	etot=etot+de
	ires=rottores(irot2)
	eres(ires)=eres(ires)+de
c
c	deactivate irot1 
c
	do k=i_first(irot1), i_first(irot1)+ni(irot1)-1
		jactive(j_index(k))=.false.
	end do
c
c	activate irot2
c	
	do k=i_first(irot2), i_first(irot2)+ni(irot2)-1
		jactive(j_index(k))=.true.
	end do

	return
	end
c
c==========================================================================
c
	function atomid(atmnam, resnam)
	implicit none
	character*3 atmnam, resnam
c
	integer restype, getatomtype, atomid

	atomid=getatomtype(atmnam)+(restype(resnam)-1)*20

	return
	end

c==========================================================================
	function getatomtype(atmnam)
	implicit none
	character*3 atmnam
	integer getatomtype, j

		j=0
		if(atmnam.eq.'C  ') j=1
		if((atmnam.eq.'O  ').or.(atmnam.eq.'OXT')) j=2
		if((atmnam.eq.'N  ').or.(atmnam.eq.'NT ')) j=3
		if(atmnam.eq.'CA ') j=4
		if(atmnam.eq.'CB ') j=5
		if((atmnam.eq.'CG ').or.(atmnam.eq.'SG ').or.
     $		   (atmnam.eq.'OG ').or.
     $		   (atmnam.eq.'CG1').or.(atmnam.eq.'OG1')) j=6
		if((atmnam.eq.'CD ').or.(atmnam.eq.'CG2').or.
     $		   (atmnam.eq.'CZ3').or.(atmnam.eq.'SD ')) j=7
		if((atmnam.eq.'CD1').or.(atmnam.eq.'OD1').or.
     $		   (atmnam.eq.'ND1')) j=8
		if((atmnam.eq.'CD2').or.(atmnam.eq.'OD2').or.
     $		   (atmnam.eq.'ND2')) j=9
		if((atmnam.eq.'OE1').or.(atmnam.eq.'CE ').or.
     $		   (atmnam.eq.'CE1').or.(atmnam.eq.'NE1')) j=10
		if((atmnam.eq.'OE2').or.(atmnam.eq.'CE2').or.
     $		   (atmnam.eq.'NZ ').or.(atmnam.eq.'NH1').or.
     $		   (atmnam.eq.'NE2')) j=11
		if((atmnam.eq.'CE3').or.(atmnam.eq.'CZ ')) j=12
		if((atmnam.eq.'CZ2').or.(atmnam.eq.'OH ').or.
     $		   (atmnam.eq.'NH2')) j=13
		if((atmnam.eq.'CH2').or.(atmnam.eq.'NE ')) j=14

	getatomtype=j

	return
	end
c==========================================================================
	function reference(aa)
	implicit none
	character aa
	real reference

	reference=-1.0
	if(aa.eq.'A') then 
		reference=106.3
	else if(aa.eq.'P') then
		reference=135.9
	else if(aa.eq.'S') then
		reference=123.1
	else if(aa.eq.'C') then
		reference=138.5
	else if(aa.eq.'T') then
		reference=141.7
	else if(aa.eq.'V') then
		reference=148.4
	else if(aa.eq.'I') then
		reference=171.5
	else if(aa.eq.'L') then
		reference=163.6
	else if(aa.eq.'D') then
		reference=149.5
	else if(aa.eq.'N') then
		reference=149.8
	else if(aa.eq.'H') then
		reference=182.2
	else if(aa.eq.'F') then
		reference=200.9
	else if(aa.eq.'Y') then
		reference=212.4
	else if(aa.eq.'W') then
		reference=245.4
	else if(aa.eq.'M') then
		reference=193.7
	else if(aa.eq.'E') then
		reference=182.8
	else if(aa.eq.'Q') then
		reference=186.6
	else if(aa.eq.'K') then
		reference=200.8
	else if(aa.eq.'R') then
		reference=239.5
	else if(aa.eq.'G') then
		reference=83.6
	end if

	return
	end

c==========================================================================
	function restype(aa)
	implicit none
	character*3 aa
	integer restype

	if(aa.eq.'ALA') then 
		restype=1
	else if(aa.eq.'PRO') then
		restype=2
	else if(aa.eq.'SER') then
		restype=3
	else if(aa.eq.'CYS') then
		restype=4
	else if(aa.eq.'THR') then
		restype=5
	else if(aa.eq.'VAL') then
		restype=6
	else if(aa.eq.'ILE') then
		restype=7
	else if(aa.eq.'LEU') then
		restype=8
	else if(aa.eq.'ASP') then
		restype=9
	else if(aa.eq.'ASN') then
		restype=10
	else if(aa.eq.'HIS') then
		restype=11
	else if(aa.eq.'PHE') then
		restype=12
	else if(aa.eq.'TYR') then
		restype=13
	else if(aa.eq.'TRP') then
		restype=14
	else if(aa.eq.'MET') then
		restype=15
	else if(aa.eq.'GLU') then
		restype=16
	else if(aa.eq.'GLN') then
		restype=17
	else if(aa.eq.'LYS') then
		restype=18
	else if(aa.eq.'ARG') then
		restype=19
	else if(aa.eq.'GLY') then
		restype=20
	else
		restype=0
	end if

	return
	end

c==========================================================================
	function restype_1(aa)
	implicit none
	character aa
	integer restype_1

	if(aa.eq.'A') then 
		restype_1=1
	else if(aa.eq.'P') then
		restype_1=2
	else if(aa.eq.'S') then
		restype_1=3
	else if(aa.eq.'C') then
		restype_1=4
	else if(aa.eq.'T') then
		restype_1=5
	else if(aa.eq.'V') then
		restype_1=6
	else if(aa.eq.'I') then
		restype_1=7
	else if(aa.eq.'L') then
		restype_1=8
	else if(aa.eq.'D') then
		restype_1=9
	else if(aa.eq.'N') then
		restype_1=10
	else if(aa.eq.'H') then
		restype_1=11
	else if(aa.eq.'F') then
		restype_1=12
	else if(aa.eq.'Y') then
		restype_1=13
	else if(aa.eq.'W') then
		restype_1=14
	else if(aa.eq.'M') then
		restype_1=15
	else if(aa.eq.'E') then
		restype_1=16
	else if(aa.eq.'Q') then
		restype_1=17
	else if(aa.eq.'K') then
		restype_1=18
	else if(aa.eq.'R') then
		restype_1=19
	else if(aa.eq.'G') then
		restype_1=20
	else
		restype_1=0
	end if

	return
	end

c==========================================================================
	function distance(a1,a2,a3,b1,b2,b3)
	implicit none
	real a1,a2,a3,b1,b2,b3,distance

	distance=sqrt((a1-b1)*(a1-b1)+(a2-b2)*(a2-b2)+(a3-b3)*(a3-b3))

	return
	end

c==========================================================================
