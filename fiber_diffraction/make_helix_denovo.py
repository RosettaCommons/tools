#!/usr/bin/env python
"""
Helix denovo symmetry definition file generator
Usage example: make_helix_denovo.py -p 3.0 -n 40 -v 5 -u 27 
where p - a helical pitch, n - number of subunits to be included in symm def file
u - number of helical units, v - number of helical turns 
"""
__author__ = "Wojtek Potrzebowski"
__credits__ = ["Wojtek Potrzebowski and Ingemar Andre"]
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

#Optparse and os import
import optparse
import os
#Numpy imports
from numpy import array
from numpy import matrix
from numpy import dot
from numpy import cross
from numpy import pi
from numpy import cos
from numpy import sin
from numpy.linalg import norm
from math import radians

class Virtual(object):
    """
    Utility functions for virtual resiudes
    """
    def __init__(self, vrtX, vrtY, vrtOrig, name=""):
        self.name = name
        self.vrtX = vrtX
        self.vrtY = vrtY
        self.vrtOrig = vrtOrig
    
    def display(self):
        return "xyz %s %6f,%6f,%6f %6f,%6f,%6f %6f,%6f,%6f\n"%(self.name,self.vrtX[0],self.vrtX[1],self.vrtX[2], self.vrtY[0],self.vrtY[1],self.vrtY[2],self.vrtOrig[0],self.vrtOrig[1],self.vrtOrig[2])
    
    def print_pdbline(self,resnum,atomnum):
        pdblines = []
        pdblines.append('ATOM    %3d  C     X Z %3d     %7.3f %7.3f %7.3f  1.00  0.00\n'%(resnum, atomnum, self.vrtX[0]+self.vrtOrig[0],self.vrtX[1]+self.vrtOrig[1], self.vrtX[2]+self.vrtOrig[2] ))
        pdblines.append('ATOM    %3d  C     Y Z %3d     %7.3f %7.3f %7.3f  1.00  0.00\n'%(resnum, atomnum+1, self.vrtY[0]+self.vrtOrig[0],self.vrtY[1]+self.vrtOrig[1], self.vrtY[2]+self.vrtOrig[2] ))
        pdblines.append('ATOM    %3d  C   ORI Z %3d     %7.3f %7.3f %7.3f  1.00  0.00\n'%(resnum, atomnum+2, self.vrtOrig[0],self.vrtOrig[1], self.vrtOrig[2] ))
        return pdblines
    
    def set_name(self, name):
        self.name = name
    
    def normalize(self):
        normval = norm( self.vrtX )
        if normval != float(0.0):
            self.vrtX = self.vrtX/normval
        normval = norm( self.vrtY )
        if normval != float(0.0):
            self.vrtY = self.vrtY/normval

class Helixer(object):
    """
    Main function to generate helical symmetry
    """
    def __init__(self, pitch, nsub, vturns, unit, name="helix_denovo", virtuals_name = "virtuals"):
        self.pitch = pitch
        self.nsub = nsub
        self.vturns = vturns
        self.unit = unit
        self.name = name
        self.virtuals_name = virtuals_name
        self.vrt_names = {}
        self.jump_names = {}
        self.jump_sub_names = {}
        for i in range(-self.nsub/2, self.nsub/2+1):
            if i<0:
                vrt_name = "VRT_0_n"+str(-i)+"_0_base"
                jump_name = "JUMP_0_n"+str(-i)+"_0"
                jump_sub_name = "JUMP_0_n"+str(-i)+"_0_to_subunit"
            else:
                vrt_name = "VRT_0_"+str(i)+"_0_base"
                jump_name = "JUMP_0_"+str(i)+"_0"
                jump_sub_name = "JUMP_0_"+str(i)+"_0_to_subunit"
            self.vrt_names[i] = vrt_name
            self.jump_names[i] =jump_name
            self.jump_sub_names[i] =jump_sub_name
        #List that stores all output lines
        self.output = []
        self.virtuals = []
        self.output.append("symmetry_name "+name+"\n")

    def __rotateMatrix(self, axis, angle ):
        unit = axis/norm(axis)
        cos_theta = cos( radians(angle) )
        one_minus_cos_theta = 1.0 - cos_theta
        sin_theta = sin( radians(angle) )
        xx = cos_theta + unit[0]*unit[0]*one_minus_cos_theta
        xy = unit[0]*unit[1]*one_minus_cos_theta - unit[2]*sin_theta
        xz = unit[0]*unit[2]*one_minus_cos_theta + unit[1]*sin_theta
        yx = unit[1]*unit[0]*one_minus_cos_theta + unit[2]*sin_theta
        yy = cos_theta + unit[1]*unit[1]*one_minus_cos_theta
        yz = unit[1]*unit[2]*one_minus_cos_theta - unit[0]*sin_theta
        zx = unit[2]*unit[0]*one_minus_cos_theta - unit[1]*sin_theta
        zy = unit[2]*unit[1]*one_minus_cos_theta + unit[0]*sin_theta
        zz = cos_theta + unit[2]*unit[2]*one_minus_cos_theta
        Rot = array( [ [xx, xy, xz], [yx, yy,yz ], [zx, zy, zz] ] )
        return Rot
    
    def __rotateVrt(self, VRT, R ):
        x = dot(VRT.vrtX,R)
        y = dot(VRT.vrtY,R)
        o = dot(VRT.vrtOrig,R)
        newVRT = Virtual( x,y,o )
        return newVRT

    def __setupEnergyLines(self):
        root_name = "VRT_0_0_0_base"
        E_line = "E = 1*"+root_name
        for i in range(self.nsub/2+1):
            if i!=0:
                E_line+=" + 1*("+root_name+":"+self.vrt_names[i]+")"
        self.output.append(E_line+"\n")
    
    def __setupVirtualCoordsLines(self):
        self.output.append("anchor_residue COM\n")
        self.output.append("recenter\n")
        self.output.append("virtual_coordinates_start\n")
        # First root virtual
        x = array( [1,0,0] )
        y = array( [0,1,0] )
        o = array( [0,0,0] )
        root =  Virtual(x,y,o,"VRT_0_0_0_base")
        #And here will be for loop from 1 to n iterating over negative and positive virtuals
        atomnum = 1
        resnum = 1
        for i in range(-self.nsub/2, self.nsub/2+1):
            if i!=0:
                zaxis = array([0,0,1])
                angle = -i*360*float(self.vturns)/float(self.unit)
                Rz = self.__rotateMatrix( zaxis, angle )
                vrt_obj = self.__rotateVrt( root, Rz )
                vrt_obj.set_name(self.vrt_names[i])
                vrt_obj.vrtOrig = root.vrtOrig+[0,0,i*self.pitch]
                # Normalize
                vrt_obj.normalize()
                x = vrt_obj.vrtX
                y = vrt_obj.vrtY
                z = cross( x, y )
                #TODO: Check it!
                y = cross (z, x )
                vrt_obj.vrtY = y
                self.virtuals.extend(vrt_obj.print_pdbline(resnum,atomnum))
                self.output.append(vrt_obj.display())
            else:
                self.virtuals.extend(root.print_pdbline(resnum,atomnum))
                self.output.append(root.display())
            atomnum+=3
            resnum+=1
        self.output.append("virtual_coordinates_stop\n")

    def __setupConnectLines(self):
        #inter-vrts connects
        for i in range(-self.nsub/2, self.nsub/2):
            self.output.append("connect_virtual "+self.jump_names[i]+" "+self.vrt_names[i]+" "+self.vrt_names[i+1]+"\n")
        #vrts-subunits connects
        for i in range(-self.nsub/2, self.nsub/2+1):
            self.output.append("connect_virtual "+self.jump_sub_names[i]+" "+self.vrt_names[i]+" SUBUNIT\n")

    def __setupDofs(self):
        #TODO: x value?
        self.output.append("set_dof "+self.jump_sub_names[0]+ " x angle_x angle_y angle_z\n")
    
    def __setupJumpGroups(self):
        jump_line1 = "set_jump_group JUMPGROUP1"
        jump_line3 = "set_jump_group JUMPGROUP3"
        for i in range(self.nsub/2):
            jump_line1+="  "+self.jump_names[i]+":"+str(i+1)
            jump_line3+="  "+self.jump_sub_names[i]
            if i>0:
                jump_line1+="  "+self.jump_names[-i]
                jump_line3+="  "+self.jump_sub_names[-i]
        #Ading last negative subuint - isn't a bit inconsitent?
        jump_line1 +="  "+self.jump_names[-self.nsub/2]
        jump_line3 +="  "+self.jump_sub_names[-self.nsub/2]+"  "+self.jump_sub_names[self.nsub/2]
        self.output.append(jump_line1+"\n")
        self.output.append(jump_line3+"\n")
    
    def writePDBlines(self):
        out_file = open(self.virtuals_name+".pdb", "w")
        out_file.writelines(self.virtuals)
        out_file.close()

    def __str__(self):
        return "".join(self.output)
    
    def write(self):
        out_file = open(self.name+".sdef", "w")
        out_file.writelines(self.output)
        out_file.close()
    
    def execute(self):
        self.__setupEnergyLines()
        self.__setupVirtualCoordsLines()
        self.__setupConnectLines()
        self.__setupDofs()
        self.__setupJumpGroups()

if '__main__' == __name__:
    doc = """
    Generates helices denovo based on parameters that can be obtained from fiber diffraction experiment
    usage: python make_helix_denovo.py --help
    """
    print doc
    usage = "usage: %prog [options] args"
    option_parser_class = optparse.OptionParser
    parser = option_parser_class( usage = usage, version='0.1' )

    parser.add_option("-p", "--rise", dest="rise", default = None,
                      type = 'float',
                      help="Helical axial rise [OBLIGATORY]")
    parser.add_option("-u", "--unit", dest="unit", default =None,
                      type = 'float',
                      help="Unit rise, subunit axial transition")
    parser.add_option("-v", "--vturns", dest="vturns", default =None,
                      type = 'int',
                      help="Number of turns in one repeat")
    parser.add_option("-n", "--nsub", dest="nsub", default = None,
                      type = 'int',
                      help="number of subuints in final oligomeric model [OBLIGATORY]")
    parser.add_option("-o", "--output", dest="name", default = "helix_denovo",
                      help="symmetry output file name")
    parser.add_option("-r", "--virtuals", dest="virtuals_name", default = "virtuals",
                      help="virtuals residues in PDB format")
    options, args = parser.parse_args()
    helixer = Helixer(options.rise, options.nsub, options.vturns, options.unit, options.name, options.virtuals_name)
    #Executing pipeline
    helixer.execute()
    helixer.write()
    helixer.writePDBlines()
        

