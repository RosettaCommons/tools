


hide lines
show cartoon
cartoon loop
select bb, (not resn PRO and (name n+ca+c+o )) or (resn PRO and (name ca+c+o))
select bbn , (not resn PRO and (name n+c+o )) or (resn PRO and (name c+o)) or elem h*
bg_color white
select water, resname WAT
select light, chain L
select heavy, chain H
select ab, light or heavy
select antigen, not ab
select doc, resname DOC
select l1, light and resi 24-34
select l2, light and resi 50-56
select l3, light and resi 89-97
select h1, heavy and resi 26-35
select h2, heavy and resi 50-56
select h2new, heavy and resi 50-56
select h3, heavy and resi 95-102
select kinkplus, heavy and resi 94+100A-103
select kinkhbond, heavy and resi 94+101
select cdrl, l1 or l2 or l3
select cdrh, h1 or h2 or h3
select cdr, cdrl or cdrh
select chargedcdr, cdr and resname ASP+GLU+LYS+ARG
#select lightframework, light and not cdr
#select heavyframework, heavy and not cdr

color yellow, light
color blue, heavy
color brightorange, cdrl
color cyan, cdrh
color green, antigen
color greencyan, h3
color orange, l3

#select lightframework, light and (resi 4-23 or resi 35-49 or resi 57-88 or resi 98-106)
select lightframework, light and (resi 4-6 or resi 10-23 or resi 35-38 or resi 45-49 or resi 57-66 or resi 71-88 or resi 98-105)

#select heavyframework, heavy and (resi 3-25 or resi 36-49 or resi 66-94 or resi 103-110)
select heavyframework, heavy and (resi 4-6 or resi 10-25 or resi 36-39 or resi 46-49 or resi 66-94 or resi 103-110)



select framework, lightframework or heavyframework
disable framework

#select lbeta, light and (resi 35-38 or resi 45-49 or resi 84-88 or resi 98-105)
#select hbeta, heavy and (resi 36-39 or resi 46-49 or resi 88-94 or resi 103-110)
#select beta, lbeta or hbeta
#select svrl, light and (not(framework or cdr))
#select svrh, heavy and (not(framework or cdr))
#select svr, svrl or svrh
#select inner, heavy and name ca and (resi 36-39 or resi 46-49 or resi 85-94 or resi 103-110)
#select outer, heavy and name ca and (resi 4-6 or resi 10-25 or resi 66-84)
#select vl, light and resi 38+98+44+95+36+43+46
#select vh, heavy and resi 39+91+45+46+103+47+37+101
select l1s, light and (resi 21-23 or resi 35-37)
select l2s, light and (resi 47-49 or resi 57-59)
select l3s, light and (resi 86-88 or resi 98-100)
select h1s, heavy and (resi 23-25 or resi 36-38)
select h2s, heavy and (resi 47-49 or resi 66-68)
select h3s, heavy and (resi 92-94 or resi 103-105)

select ls, l1s or l2s or l3s
select hs, h1s or h2s or h3s
select stem, ls or hs
disable stem
#select hotspot_l, light and resi 50+53+92+96
#select hotspot_h, heavy and resi 33+50+52+100+99
#select hotspot_ab, hotspot_l or hotspot_h
#select hotspot_anti, resi 682+684+685+687+688+654
#select hotspot, hotspot_anti or hotspot_ab
#disable hotspot
#color red, hotspot
#select heavystem, heavy and (resi 23-25 or resi 47-49 or resi 92-94 or resi 36-38 or resi 66-68 or resi 103-105)
#disable heavystem
zoom all
hide everything, not ab
zoom ab
