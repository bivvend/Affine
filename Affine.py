# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 22:30:55 2017

@author: Simon Henley
"""
## Modified from version written by Jarno Elonen
import math
import sys
from ast import literal_eval # used to convert string to list of tuples

def Affine_Fit( from_pts, to_pts ):
    """Fit an affine transformation to given point sets.
      More precisely: solve (least squares fit) matrix 'A'and 't' from
      'p ~= A*q+t', given vectors 'p' and 'q'.
      Works with arbitrary dimensional vectors (2d, 3d, 4d...).

      Written by Jarno Elonen <elonen@iki.fi> in 2007.
      Placed in Public Domain.

      Based on paper "Fitting affine and orthogonal transformations
      between two sets of points, by Helmuth Späth (2003)."""

    q = from_pts
    p = to_pts
    if len(q) != len(p) or len(q)<1:
        print("from_pts and to_pts must be of same size.")
        return False

    dim = len(q[0]) # num of dimensions
    if len(q) < dim:
        print("Too few points => under-determined system.")
        return False

    # Make an empty (dim) x (dim+1) matrix and fill it
    c = [[0.0 for a in range(dim)] for i in range(dim+1)]
    for j in range(dim):
        for k in range(dim+1):
            for i in range(len(q)):
                qt = list(q[i]) + [1]
                c[k][j] += qt[k] * p[i][j]

    # Make an empty (dim+1) x (dim+1) matrix and fill it
    Q = [[0.0 for a in range(dim)] + [0] for i in range(dim+1)]
    for qi in q:
        qt = list(qi) + [1]
        for i in range(dim+1):
            for j in range(dim+1):
                Q[i][j] += qt[i] * qt[j]

    # Ultra simple linear system solver. Replace this if you need speed.
    def gauss_jordan(m, eps= 1.0/(10**10)):
        #Puts given matrix (2D array) into the Reduced Row Echelon Form.
        #Returns True if successful, False if 'm' is singular.
        #NOTE: make sure all the matrix items support fractions! Int matrix will NOT work!
        #Written by Jarno Elonen in April 2005, released into Public Domain
        (h, w) = (len(m), len(m[0]))
        for y in range(0,h):
            maxrow = y
            for y2 in range(y+1, h):    # Find max pivot
                if abs(m[y2][y]) > abs(m[maxrow][y]):
                    maxrow = y2
            (m[y], m[maxrow]) = (m[maxrow], m[y])
            if abs(m[y][y]) <= eps:     # Singular?
                return False
            for y2 in range(y+1, h):    # Eliminate column y
                c = m[y2][y] / m[y][y]
                for x in range(y, w):
                    m[y2][x] -= m[y][x] * c
        for y in range(h-1, 0-1, -1): # Backsubstitute
            c  = m[y][y]
            for y2 in range(0,y):
                for x in range(w-1, y-1, -1):
                    m[y2][x] -=  m[y][x] * m[y2][y] / c
            m[y][y] /= c
            for x in range(h, w):       # Normalize row y
                m[y][x] /= c
        return True

    # Augement Q with c and solve Q * a' = c by Gauss-Jordan
    M = [ Q[i] + c[i] for i in range(dim+1)]
    if not gauss_jordan(M):
        print("Error: singular matrix. Points are probably coplanar.")
        return False

    # Make a result object
    class Transformation:
        """Result object that represents the transformation
           from affine fitter."""

        def To_Str(self):
            res = ""
            for j in range(dim):
                str = "x%d' = " % j
                for i in range(dim):
                    str +="x%d * %f + " % (i, M[i][j+dim+1])
                str += "%f" % M[dim][j+dim+1]
                res += str + "\n"
            return res
        
        def Get_Dim(self):
            return dim
            
        def Get_Element(self,i,j):
            return M[i][j]
            
        def Transform(self, pt):
            res = [0.0 for a in range(dim)]
            for j in range(dim):
                for i in range(dim):
                    res[j] += pt[i] * M[i][j+dim+1]
                res[j] += M[dim][j+dim+1]
            return res
    return Transformation()

if __name__ == '__main__':
   
    print("Running script " + sys.argv[0] + "....")

    if len(sys.argv)>2:
        print("Using EXTERNAL data")
        from_pt =literal_eval(sys.argv[1])
        to_pt = literal_eval(sys.argv[2])
    else:
        print("Using INTERNAL data")
        #from_pt = ((1,1),(1,2),(2,2),(2,1))   
        #to_pt =  ((4,4),(6,6),(8,4),(6,2))

        #from_pt = ((-3, 0), (0, 3), (3, 0))
        #to_pt =  ((2, 3), (3, 2), (4, 3)) 

        #Scale X1.1 Y1.3 rotate 80 degs translate X0.1 Y-80   -  OK
        #from_pt = ((-30,30),(30,-30),(-30,-30))
        #to_pt =  ((-44.0379,-105.726),(44.23789,-54.2736),(32.77711,-119.271))   

        #Scale X1.1 Y1.3 rotate 175 degs translate X0.1 Y-80   -  GETS SCALE NEGATIVE
        #from_pt = ((-30,30),(30,-30),(-30,-30))
        #to_pt =  ((29.57535,-121.728),(-29.3754,-38.2723),(36.3735,-44.0245))   

        #Scale X1.1 Y1.3 rotate 175 degs translate X0.1 Y-80   -  OK
        from_pt = ((-30,10),(30,-10),(30,10))
        to_pt =  ((31.8414,-95.8267),(-31.6414,-64.1733),(-33.9074,-90.0744))  

    trn = Affine_Fit(from_pt, to_pt)
    dim = trn.Get_Dim()

    print("Dimension: %d" %dim)
    print("Transformation is:")
    print(trn.To_Str())
    
    print(" ---------- Full Matrix   ---------- ")
    for j in range(dim+dim+1):
        for i in range(dim+1):
            print( str(i) +"," + str(j)  + " " + str(trn.Get_Element(i,j))  )       
    print(" ------------------------------------ ")
    print("\n")
    
    a = trn.Get_Element(0,3)
    b = trn.Get_Element(1,3)
    c = trn.Get_Element(0,4)
    d = trn.Get_Element(1,4)
    
    print("a = " + str(a))
    print("b = " + str(b))
    print("c = " + str(c))
    print("d = " + str(d))
    

    #ROTATE SCALE TRANSLATE SOLUTION
    print(" -------------------------------------- ")
    print(" --- ROTATE --- SCALE --- TRANSLATE --- ")
    
    offset_x = trn.Get_Element(2, 3)
    offset_y = trn.Get_Element(2, 4)
    print("Offset x:" + str(offset_x))
    print("Offset y:" + str(offset_y))
    
    #EXCLUDE NEGATIVE SCALES
    signa=1.0
    signd=1.0
    """
    if a<0:
        signa = -1.0
    else:
        signa = 1.0
    if d<0:
        signd = -1.0
    else:
        signd = 1.0
    """  
       
    sx = signa*math.sqrt(a*a + b*b)
    sy = signd*math.sqrt(c*c + d*d)
    angle1 = math.atan2(-1.0*b,a)
    angle2 = math.atan2(c,d)

    angle = (angle1+angle2)/2.0
    
    print("Scale x = " + str(sx) )
    print("Scale y = " + str(sy) )
    print("Angle (-b,a) = " + str(angle1) )
    print("Angle (c,d) = " + str(angle2) )
    print("Average Angle = " + str(angle))

    #SCALE ROTATE TRANSLATE SOLUTION
    print(" -------------------------------------- ")
    print(" --- SCALE --- ROTATE --- TRANSLATE --- ")

    offset_x = trn.Get_Element(2, 3)
    offset_y = trn.Get_Element(2, 4)
    print("Offset x:" + str(offset_x))
    print("Offset y:" + str(offset_y))
    
    #EXCLUDE NEGATIVE SCALES
    signa=1.0
    signd=1.0
    """
    if a<0:
        signa = -1.0
    else:
        signa = 1.0
    if d<0:
        signd = -1.0
    else:
        signd = 1.0  
    """ 
       
    sx = signa*math.sqrt(a*a + c*c)
    sy = signd*math.sqrt(b*b + d*d)

    angle1 = math.atan2(-1.0*b,d)
    angle2 = math.atan2(c,a)
    angle = (angle1+angle2)/2.0

    print("Scale x = " + str(sx) )
    print("Scale y = " + str(sy) )
    print("Angle (-b,d) = " + str(angle1) )
    print("Angle (c,a) = " + str(angle2) )
    print("Average Angle = " + str(angle))


    print(" -------------------------------------- ")

    err = 0.0
    for i in range(len(from_pt)):
        fp = from_pt[i]
        tp = to_pt[i]
        t = trn.Transform(fp)
        print ("%s => %s ~= %s" % (fp, tuple(t), tp))
        err += ((tp[0] - t[0])**2 + (tp[1] - t[1])**2)**0.5
    
    print("Fitting error = %f" % err)
