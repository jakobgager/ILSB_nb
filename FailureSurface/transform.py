#! /usr/bin/env python
"""
Created on Tue Feb  7 12:14:41 2012

@author: gager

Transform input file using an affine transformation
"""
import os.path
import numpy as np

def trans_inp(filename, matrix):
    """Method to transform nodes of an input file

    Reads the nodes of an Abaqus input file (filename.inp) and transforms 
    these using a transformation matrix. The transformed nodes are written 
    into a new file (filename_trans.inp). All other lines of the original 
    file are copied (but not modfied) as well.

    The transformation matrix has the following form:

    .. note::
        Currently a 3x3 matrix is implemented, thus no shifting is possible

    .. math::

      \left( \begin{array}{cccc}
      r1 & r2 & r3 & x1  \\
      r4 & r5 & r6 & y1 \\
      r7 & r8 & r9 & z1 \\
      0  &  0 & 0 &  1
      \end{array} \right)

    with:
        (x1,y1,z1) ...translation
        (r1 .. r9) ...rotation/scaleing/shearing matrix
        
    **Args:**
      - **filename** *(string)* - filename of the input file
      - **matrix** *(numpy matrix)* - transformation matrix
    """
    #Node switch (True if Node)
    ns = False

    #input file
    ifile = open(filename,'r')
    #output file
    filename = os.path.splitext(filename)[0]
    ofile = open(filename + '_trans.inp','w')

    #read file
    for line in ifile:
        # if no node active
        if not ns:
            ofile.write(line)
        else:
        #check if still node
            if line.startswith('*'):
                if line.startswith('**'):
                    #comment line
                    continue
                else:
                    #new command
                    ns = False
                ofile.write(line)
            else:
            #transform
                node = line.strip().split(',')
                #print node
                node = map(float,node)
                nodenr = int(node[0])
                node = np.array(node[1:])
                onode = np.dot(matrix,node)
                #print nodenr,onode
                #ofile.write('war: {0}, {1[0]}, {1[1]}, {1[2]} \n'.format(nodenr, node))
                ofile.write(' {0}, {1[0]}, {1[1]}, {1[2]} \n'.format(nodenr, onode))

        if line.lower().startswith('*node'): 
            ns = True

    ifile.close()
    ofile.close()

def trans_vtk(filename, matrix):
    """Method to transform nodes of an vtk file

    Reads the nodes of an VTK file (filename.vtk) and transforms 
    these using a transformation matrix. The transformed nodes are written 
    into a new file (filename_trans.vtk). All other lines of the original 
    file are copied (but not modfied) as well.

    The transformation matrix has the following form:

    .. note::
        Currently a 3x3 matrix is implemented, thus no shifting is possible

    .. math::

      \left( \begin{array}{cccc}
      r1 & r2 & r3 & x1  \\
      r4 & r5 & r6 & y1 \\
      r7 & r8 & r9 & z1 \\
      0  &  0 & 0 &  1
      \end{array} \right)

    with:
        (x1,y1,z1) ...translation
        (r1 .. r9) ...rotation/scaleing/shearing matrix
        
    **Args:**
      - **filename** *(string)* - filename of the input file
      - **matrix** *(numpy matrix)* - transformation matrix
    """
    #Node switch (True if Node)
    ns = False

    #input file
    ifile = open(filename,'r')
    #output file
    filename = os.path.splitext(filename)[0]
    ofile = open(filename + '_trans.vtk','w')

    #read file
    for line in ifile:
        # if no node active
        if not ns:
            ofile.write(line)
        else:
        #check if still node
            if not line[0].isdigit() and line[0]!='-':
                ns = False
                ofile.write(line)
            else:
            #transform
                node = line.strip().split(' ')
                #print node
                node = map(float,node)
                node = np.array(node)
                onode = np.dot(matrix,node)
                #print nodenr,onode
                #ofile.write('war: {0}, {1[0]}, {1[1]}, {1[2]} \n'.format(nodenr, node))
                ofile.write('{0[0]} {0[1]} {0[2]}\n'.format(onode))

        if line.lower().startswith('points'): 
            ns = True

    ifile.close()
    ofile.close()

# some transformation matrices
def scalmat(sf):
    """Transformation matrix isotropic scaling

    ! Needs Numpy as np

    Matrix:

    .. math::

      \left( \begin{array}{cccc}
      sf & 0 & 0  \\
      0 & sf & 0  \\
      0 & 0 & sf  
      \end{array} \right)

    **Args:**
      - **sf** *(float)* - scale factor

    **Returns:**
      - Transformation matrix as numpy array
    """
    return np.eye(3)*sf
##############################
def rotmat(phi):
    """Transformation matrix rotation

    ! Needs Numpy as np

    Matrix:

    .. math::

      \left( \begin{array}{cccc}
      \cos \phi & \sin \phi & 0  \\
      -\sin \phi & \cos \phi & 0  \\
      0 & 0 & 1  
      \end{array} \right)

    **Args:**
      - **phi** *(float)* - rotation angle in deg

    **Returns:**
      - Transformation matrix as numpy array
    """
    phirad = phi/180. * np.pi
    return np.array([[np.cos(phirad), np.sin(phirad), 0.0],
                     [-np.sin(phirad),np.cos(phirad), 0.0],
                     [0.0,         0.0,         1.0, ]])
##############################
def shearmat(gamma):
    """Transformation matrix pure shear

    ! Needs Numpy as np

    Matrix:

    .. math::

      \left( \begin{array}{cccc}
      1 & \tan\frac{\gamma}{2} & 0  \\
      \tan\frac{\gamma}{2} & 1 & 0  \\
      0 & 0 & 1  
      \end{array} \right)

    **Args:**
      - **gamma** *(float)* - shear angle

    **Returns:**
      - Transformation matrix as numpy array
    """
    gamma = gamma/180.*np.pi
    a = np.tan(gamma/2.)
    return np.array([[1, a, 0.0],
                     [a, 1, 0.0],
                     [0.0, 0.0, 1.0, ]])

##################################################
def main():
# the main code goes here

        filename = 'tows.vtk'

        idmat = np.eye(3)*2
        rot = rotmat(45)
        she = shearmat(30)
        rot2 = rotmat(-15)
        she2 = np.array(np.matrix(rot2)*np.matrix(she))

        if os.path.splitext(filename)[1] == '.inp':
            trans_inp(filename, she2)
        elif os.path.splitext(filename)[1] == '.vtk':
            trans_vtk(filename, she2)
        else:
            print 'Error: Filetype not defined.'

if __name__ == "__main__":
    main()
