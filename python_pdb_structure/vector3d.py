import math

class vector3d :
    def __init__( self, x=0, y=0, z=0 ) :
        self.x_ = x
        self.y_ = y
        self.z_ = z
    def __add__( self, other ) :
        v = vector3d()
        v.x_ = self.x_ + other.x_
        v.y_ = self.y_ + other.y_
        v.z_ = self.z_ + other.z_
        return v
    def __sub__( self, other ) :
        v = vector3d()
        v.x_ = self.x_ - other.x_
        v.y_ = self.y_ - other.y_
        v.z_ = self.z_ - other.z_
        return v
    def __mul__( self, other ) :
        dotprod = 0
        dotprod += self.x_ * other.x_
        dotprod += self.y_ * other.y_
        dotprod += self.z_ * other.z_
        return dotprod
    def __str__( self ) :
        me = "("
        me += str( self.x_ ) + ", "
        me += str( self.y_ ) + ", "
        me += str( self.z_ ) + ")"
        return me
    def x(self) :
        return self.x_
    def y(self) :
        return self.y_
    def z(self) :
        return self.z_
    def copy(self, other ) :
        self.x_ = other.x_;
        self.y_ = other.y_;
        self.z_ = other.z_;
    def set( self, x, y, z ):
        self.x_ = x
        self.y_ = y
        self.z_ = z
    def set_x( self, val ) :
        self.x_ = val
    def set_y( self, val ) :
        self.y_ = val
    def set_z( self, val ) :
        self.z_ = val
    def distance_squared( self, other ) :
        return ( self.x_ - other.x_ ) * (self.x_ - other.x_ ) + \
               ( self.y_ - other.y_ ) * (self.y_ - other.y_ ) + \
               ( self.z_ - other.z_ ) * (self.z_ - other.z_ )
    def distance( self, other ) :
        return math.sqrt( self.distance_squared(other) )
    def normalize( self ) :
        len = self.length()
        if len == 0.0 :
            return
        invlength = 1 / len
        self.scale( invlength )
        return self
    def scale( self, stretch ) :
        self.x_ *= stretch
        self.y_ *= stretch
        self.z_ *= stretch
        return self
    def length( self ) :
        return math.sqrt( self * self )
    def max( self, other ) :
        #print "max start: ", self, other,
        if self.x_ < other.x_ : self.x_ = other.x_
        if self.y_ < other.y_ : self.y_ = other.y_
        if self.z_ < other.z_ : self.z_ = other.z_
        #print "max end: ", self
    def min( self, other ) :
        #print "min start: ", self, other,
        if self.x_ > other.x_ : self.x_ = other.x_
        if self.y_ > other.y_ : self.y_ = other.y_
        if self.z_ > other.z_ : self.z_ = other.z_
        #print "min end: ", self


if __name__ == "__main__" :
    p1 = vector3d()
    p1.set( 1, 2, 3 )
    p2 = vector3d()
    p2.set( 4, 5, 6 )
    print(p1 + p2)
    print(p1 * p2)
    p3 = p1 - p2
    p3.normalize()
    print(p3)
