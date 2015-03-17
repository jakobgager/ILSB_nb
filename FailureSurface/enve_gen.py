# -*- coding: utf-8 -*-

import sympy as sp

# FPF Surfaces

# Tsai-Wu Failure Envelope

def Tsai_Wu():
    # Tsai-Wu (taken from http://de.wikipedia.org/wiki/Tsai-Wu-Kriterium)

    # Input stress
    s11, s22, s12 = sp.symbols('S_{11} S_{22} S_{12}')
    # Strength quantities
    Rlt, Rlc, Rqt, Rqc, Rql = sp.symbols('R_{lt},R_{lc},R_{qt},R_{qc}, R_{ql}')
    # Tsai-Wu factor
    F12s = sp.symbols('F_{12}^*')

    #####
    a11, a22 = sp.symbols('a11 a22')
    b1111, b2222, b1122, b1212 = sp.symbols('b1111 b2222 b1122 b1212')
    #
    #
    a11 = 1/Rlt - 1/Rlc
    a22 = 1/Rqt - 1/Rqc
    b1111 = 1/(Rlt*Rlc)
    b2222 = 1/(Rqt*Rqc)
    b1212 = 1/(4*Rql**2)
    b1122 = F12s/sp.sqrt(Rlt*Rlc*Rqt*Rqc)

    F_tw = a11*s11 + a22*s22 + b1111*s11**2 + b2222*s22**2 + 2*b1122*s11*s22 + 4*b1212*s12**2

    return F_tw



'''Create a unitsphere recursively by subdividing all triangles in an octahedron recursivly.

A unitsphere has a radius of 1, which also means that all points in this sphere
have an absolute value of 1. Another feature of an unitsphere is that the normals 
of this sphere are exactly the same as the vertices.

This recursive method will avoid the common problem of the polar singularity, 
produced by 2d parameterization methods.

If you wish a sphere with another radius than that of 1, simply multiply every single
value in the vertex array with this new radius 
(although this will break the "vertex array equal to normal array" property)
'''
import numpy

octahedron_vertices = numpy.array( [ 
    [ 1.0, 0.0, 0.0], # 0 
    [-1.0, 0.0, 0.0], # 1
    [ 0.0, 1.0, 0.0], # 2 
    [ 0.0,-1.0, 0.0], # 3
    [ 0.0, 0.0, 1.0], # 4 
    [ 0.0, 0.0,-1.0]  # 5                                
    ] )
octahedron_triangles = numpy.array( [ 
    [ 0, 4, 2 ],
    [ 2, 4, 1 ],
    [ 1, 4, 3 ],
    [ 3, 4, 0 ],
    [ 0, 2, 5 ],
    [ 2, 1, 5 ],
    [ 1, 3, 5 ],
    [ 3, 0, 5 ]] )

def normalize_v3(arr):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    lens = numpy.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr

def divide_all( vertices, triangles ):    
    #new_triangles = []
    new_triangle_count = len( triangles ) * 4
    # Subdivide each triangle in the old approximation and normalize
    #  the new points thus generated to lie on the surface of the unit
    #  sphere.
    # Each input triangle with vertices labelled [0,1,2] as shown
    #  below will be turned into four new triangles:
    #
    #            Make new points
    #                 a = (0+2)/2
    #                 b = (0+1)/2
    #                 c = (1+2)/2
    #        1
    #       /\        Normalize a, b, c
    #      /  \
    #    b/____\ c    Construct new triangles
    #    /\    /\       t1 [0,b,a]
    #   /  \  /  \      t2 [b,1,c]
    #  /____\/____\     t3 [a,b,c]
    # 0      a     2    t4 [a,c,2]    
    v0 = vertices[ triangles[:,0] ]
    v1 = vertices[ triangles[:,1] ]
    v2 = vertices[ triangles[:,2] ]
    a = ( v0+v2 ) * 0.5
    b = ( v0+v1 ) * 0.5
    c = ( v1+v2 ) * 0.5  
    normalize_v3( a )
    normalize_v3( b )
    normalize_v3( c )
    
    #Stack the triangles together.
    # modified the following line (add indices)
    vertices = numpy.vstack( [(v0[i],b[i],a[i],  b[i],v1[i],c[i],  a[i],b[i],c[i], a[i],c[i],v2[i]) for i in range(len(v0))] )
    #Now our vertices are duplicated, and thus our triangle structure are unnecesarry.    
    return vertices, numpy.arange( len(vertices) ).reshape( (-1,3) )

def create_unit_sphere( recursion_level=2 ):
    vertex_array, index_array = octahedron_vertices, octahedron_triangles
    for i in range( recursion_level - 1 ):
        vertex_array, index_array  = divide_all(vertex_array, index_array)
    return vertex_array, index_array


def vertex_array_only_unit_sphere( recursion_level=2 ):
    vertex_array, index_array = create_unit_sphere(recursion_level)
    if recursion_level > 1:    
        return vertex_array.reshape( (-1) )
    else:
        return vertex_array[index_array].reshape( (-1) )



def readvtk(filename, points=[], cells=[], celltype='POLYGONS'):
    offset = len(points)
    with open(filename) as fil:

        lin = 1
        while lin:
            lin = fil.readline()
            if lin.startswith('POINTS'):
                numpt = int(lin.split()[1])
                for i in range(numpt):
                    points.append(np.array(fil.readline().split(),dtype=float))
            if lin.startswith(celltype):
                numcells = int(lin.split()[1])
                for i in range(numcells):
                    cells.append(np.array(fil.readline().split()[1:],dtype=int)+offset)

    return points, cells


# VTK-Write Routine

def writevtk(filename,nodes, translist, nr_sets, nr_levels, sym_idx=0, data=[]):
    nonr=len(nodes)
    elnr=len(translist)
    f=open(filename,'w')
    f.writelines('# vtk DataFile Version 3.0\n')
    f.write('{0} {1} {2} {3} {4}\nASCII\nDATASET POLYDATA\n'.format(nonr/nr_sets, elnr/nr_sets, nr_sets, nr_levels, sym_idx))
    f.write('POINTS %d float\n'%nonr)
    for i in range(nonr):
        f.write('%f %f %f\n'%tuple(nodes[i]))
    f.write('POLYGONS %d %d\n'%(elnr,elnr*4))
    for n in range(elnr):
        f.write('3 %d %d %d \n'%tuple(translist[n]))
    # write data tags
    if len(data) > 0:
        f.write('POINT_DATA {0}\n'.format(nonr))
        for dset in data:
            f.write('SCALARS {0} int\nLOOKUP_TABLE default\n'.format(dset[0]))
            for i in dset[1]:
                f.write('{0}\n'.format(i))
            
    f.close()

