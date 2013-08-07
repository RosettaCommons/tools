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
c**************************************************************************
c	sizes.for
c**************************************************************************
c
c	include file for TORSO
c
c       parameter           maximum number of
c
c     	maxatm          atoms in the protein
c     	maxres          amino acid residues + rotamers
c     	maxgrid         grid points in each dimension (-maxgrid,+maxgrid)
c     	maxrot          amino acid rotamers in library
c	maxpair		rotamer-rotamer interactions
c	maxnoccs	partial atomic occupancies 
c
	integer*2 maxatm,maxres,maxgrid,maxrot
	integer*4 maxpair, maxnoccs

	parameter (maxatm=30000)
	parameter (maxres=4500)
	parameter (maxgrid=40)
	parameter (maxrot=400)
	parameter (maxpair=240000)
	parameter (maxnoccs=2)
