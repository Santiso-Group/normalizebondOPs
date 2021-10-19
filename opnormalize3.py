import numpy as np
import mdtraj as md
import sys
np.set_printoptions(threshold=sys.maxsize)
import argparse
from scipy.spatial import cKDTree

def getargs():
    '''
    Description: This pulls command line arguments and parses them using argparser
    arguments
    the arguments passed in the command line
    returns
    the arguments
    '''
    parser = argparse.ArgumentParser() #define parser
    #get command line arguments.  includes default values and help text
    parser.add_argument("--frame", help="name of the frame", default='', type=str)
    parser.add_argument("--origin", help="origin atom", default=0, type=int)
    parser.add_argument("--atomspermol", help="number of atoms per molecule", default=1, type=int)
    parser.add_argument("--nmols", help="number of molecules in the frame", default=1, type=int)
    parser.add_argument("--cutoff", help="cutoff distance in angstrom,", default=1, type=float)
    parser.add_argument("--cop", help="name of the cop file", default='', type=str)
    parser.add_argument("--copnorm", help="name of the normalized cop file", default='', type=str)
    parser.add_argument("--resname", help="name of the residue", default='', type=str)
    args = parser.parse_args() #put command line arguments into args variable
    frame = args.frame
    origin = args.origin
    atomspermol = args.atomspermol
    nmols = args.nmols
    cutoff = args.cutoff
    cop = args.cop
    copnorm = args.copnorm
    resname = args.resname
    return frame,origin,atomspermol,nmols,cutoff,cop,copnorm,resname #return variables

def distance(a,b,bounds):
    a = np.asarray((a))
    b = np.asarray((b))
    min_dists = np.min(np.dstack(((a - b) % bounds, (b - a) % bounds)), axis = 2)
    dist = np.sqrt(np.sum(min_dists ** 2, axis = 1))
    return dist

def normalize(frame,origin,atomspermol,nmols,cutoff,cop,copnorm,resname):
    t = md.load(frame)
    top = md.load(frame).topology
    atom_indices = np.zeros(nmols, dtype=int)
    res_string = 'resname ' + resname
    #table, bonds = top.to_dataframe()
    #print(table)
    resindices = top.select(res_string)
    t = t.atom_slice(resindices)
    for i in range(len(atom_indices)):
        atom_indices[i] = int(i*atomspermol+origin)
    t = t.atom_slice(atom_indices)
        #get PBCs bounds.  This only works for orthohombic PBCs for the time being.  Might implement more general solution later.
    pbc_vectors = np.diag(np.asarray((t.unitcell_vectors*10)[0]))
    coords = t.xyz*10 #defines coordinates.  multiply by 10 to get angstroms, since mdtraj converts to nm.
    coords = coords[0]
    xm = min(coords[:,0])
    ym = min(coords[:,1])
    zm = min(coords[:,2])
    #pbc_vectors = [max(coords[:,0]),max(coords[:,1]),max(coords[:,2])]
    for i in range(len(coords)):
        coords[i][0] = coords[i][0]-xm
        if coords[i][0] > pbc_vectors[0]:
            coords[i][0] = coords[i][0]-pbc_vectors[0]
        coords[i][1] = coords[i][1]-ym
        if coords[i][1] > pbc_vectors[1]:
            coords[i][1] = coords[i][1] - pbc_vectors[1]
        coords[i][2] = coords[i][2]-zm
        if coords[i][2] > pbc_vectors[2]:
            coords[i][2] = coords[i][2] - pbc_vectors[2]
    try:
        tree = cKDTree(coords,boxsize = pbc_vectors)
    except:
        pbc_vectors = [max(coords[:,0]),max(coords[:,1]),max(coords[:,2])]
        tree = cKDTree(coords,boxsize = pbc_vectors)
        
    neighbors = tree.query_ball_tree(tree,cutoff)
    nneighbors = []
    for molecule in neighbors:
        nneighbors.append(len(molecule))
    f = open(cop,'r')
    g = open(copnorm,'w')
    orderparameters = []
    for x in f:
        orderparameters.append(x)
    f.close()
    bondops = np.array(orderparameters[(10+nmols):(10+2*nmols)],dtype=float)
    normalops = bondops/(nneighbors)
    np.savetxt(copnorm,normalops)
    np.savetxt(copneighbors,nneighbors)

def main():
    frame,origin,atomspermol,nmols,cutoff,cop,copnorm,resname = getargs()
    normalize(frame,origin,atomspermol,nmols,cutoff,cop,copnorm,resname)

if __name__ == '__main__':
  main() #run main
