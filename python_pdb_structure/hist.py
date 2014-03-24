from optparse import OptionParser
import os
import sys
import math

def initialize_options_parser():
   parser = OptionParser()
   parser.add_option( "-d", "--data-file", dest="data_file",
      help="The data file, which should consist of a single column of numbers" )
   parser.add_option( "-s", "--half-bin", action="store_true", dest="half_bin", default=False,
      help="Represent each bin as the range +/- half a bin width from the bin representative;" +
      "if absent then the range is from the bin representative to the bin representative + the bin width." )
   parser.add_option( "-w", "--bin-width", dest="bin_width",
      help="The bin width to use", default=1.0 )
   parser.add_option( "-l", "--lower-boundary", dest="lower_boundary", type="float",
      help="The lower boundary for the histogram.  Points lower than the boundary are not counted.", default=None )
   parser.add_option( "-u", "--upper-boundary", dest="upper_boundary", type="float",
      help="The upper boundary for the histogram.  Points above the uper boundary are not counted.", default=None )
   parser.add_option( "-n", "--normalize", action="store_true", dest="normalize", default=False,
      help="report the frequency of occurances per bin instead of the total number of occurances" )
   parser.add_option( "-m", "--normalize-over-range", action="store_true", dest="normalize_over_range",
      help="report the frequency of occurances per bin instead of the total number; the frequency is defined as the percentage " +
      "of data points that fall in the range specified by the --upper-boundary or the --lower-boundary flags.  Behaves " +
      "just like --normalize if niether --upper-boundary nor --lower-boundary are specified." )
   parser.add_option( "-q", "--quiet", action="store_true", dest="quiet",
      help="print out the histogram values without the bin center column" )
   return parser

class Histogramer() :
   def __init__( self ) :
      self.data = {}
      self.n_values = 0
      self.n_values_considered = 0
      self.lower_boundary = None
      self.upper_boundary = None
      self.normalize = False
      self.normalize_over_range = False
      self.bin_width = 1.0
      self.half_bin = False
      self.quiet = False

   def setup_from_options( self, opts ) :
      self.bin_width = opts.bin_width
      self.bin_width_over_2 = self.bin_width / 2.0
      self.normalize = opts.normalize
      self.normalize_over_range = opts.normalize_over_range
      self.half_bin = opts.half_bin
      self.quiet = opts.quiet
      if opts.lower_boundary != None :
         self.lower_boundary = opts.lower_boundary
         self.lower_remainder = self.lower_boundary - self.bin_width * math.floor( self.lower_boundary * self.bin_width )
      if opts.upper_boundary != None :
         self.upper_boundary = opts.upper_boundary
      if self.upper_boundary != None and self.lower_boundary != None :
         rep = self.lower_boundary
         while rep <= self.upper_boundary :
            self.data[ rep ] = []
            rep += self.bin_width

   def representative_for_value( self, v ) :
      if self.lower_boundary == None :
         if self.half_bin :
             v2 = v + self.bin_width_over_2
         else :
             v2 = v
         rep = math.floor( v2 * self.bin_width ) / self.bin_width
         return rep
      else :
         if self.half_bin :
            v2 = v - self.bin_width_over_2
         else :
            v2 = v
         rep = math.floor( v2 * self.bin_width ) / self.bin_width + self.lower_remainder 
         return rep

   def insert_value( self, v ) :
      rep = self.representative_for_value( v )
      self.add_value( rep, v )

   def add_value( self, rep, v ) :
      self.n_values_considered += 1.0
      if self.lower_boundary != None :
         if v < self.lower_boundary :
            return
      if self.upper_boundary != None :
         if v > self.upper_boundary :
            return
      if rep not in self.data :
         self.data[ rep ] = []
      self.data[ rep ].append( v )
      self.n_values += 1.0

   def print_histogram( self  ) :
      keys = self.data.keys()
      keys.sort()
      for key in keys :
         print self.first_column_for_key( key ) + self.printval_for_key( key )

   def printval_for_key( self, key ) :
      len_for_key = len(self.data[key])
      if self.normalize :
         return str( float( len_for_key ) / self.n_values_considered )
      elif self.normalize_over_range :
         return str( float( len_for_key ) / self.n_values )
      else :
         return str( float ( len_for_key ) )
   def first_column_for_key( self, key ) :
      if self.quiet :
         return ""
      return str( key ) + "\t"

if __name__ == "__main__" :
   parser = initialize_options_parser()
   (options, args ) = parser.parse_args()
   hist = Histogramer()
   hist.setup_from_options( options )
   lines = open( options.data_file ).readlines()
   for line in lines :
      hist.insert_value( float( line.strip() ) )
   hist.print_histogram()
