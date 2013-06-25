c This module/program is part of DaliLite (c) L. Holm 1999
c

	function ran(seed) 
	implicit none
	real ran 
	integer*8 seed,x,y,a,c
	data a,c,y/155,138,65536/
c        write (*,*) 'ran: seed=', seed
	x=mod (a*seed+c, y)
	seed=x
	ran=float(x)/float(y)
c        write (*,*) 'ran: x=', x, ' y=', y, ' ran=', ran
	return
	end
