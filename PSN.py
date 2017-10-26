#!/Users/erol/Code/2.7/bin/python
#Change this before running

import sys, getopt, glob, csv
from Bio import Struct
from Bio.Struct.Geometry import center_of_mass
from Bio.PDB import *
import numpy as np

def string_to_int(val):
  try:
    return int(val)
  except ValueError as e:
    raise


# Returns the number of atoms per residue pair that have a
# euclidean distance below a threshold
# Inputs:
# model: Structure object
# threshold: Threshold
# excluded: Excluded atoms
# Returns: |R| x |R| matrix |R| = # of residues, Rij = # of connected atom bet Ri and Rj
def compute_atom_contact(model, threshold, excluded=[]):
  res_count = 0
  # Get residue counts on all chain
  for chain in model:
    res_count += len(chain)
  # generate mod contact matrix
  mod_contact_matrix = np.zeros((res_count,res_count))
  last_i = 1
  for chain in model:
    i = last_i
    # Assume for now na single chain lang
    for res_i in range(len(chain)):
      index_i = res_i + i # add i to consider multiple chains
      last_i = index_i # we note the last value of i for multiple chain
      for res_j in range(len(chain)-res_i):
        index_j = res_j + index_i
        if index_i < index_j: # Wag na magrecompute
          # for ptm
          try:
              for at_i in chain[index_i]:
                if at_i not in excluded:
                  for at_j in chain[index_j]:
                    if at_j not in excluded:
                      euclidean = np.linalg.norm(at_i.get_coord()-at_j.get_coord())
                      if euclidean <= threshold:
                        # add one interacting node to both residues
                        mod_contact_matrix[index_i-1][index_j-1] += 1
                        mod_contact_matrix[index_j-1][index_i-1] += 1
          except KeyError as e:
            break
          
  return mod_contact_matrix


# Returns the ave number of atoms per residue pair that have a
# euclidean distance below a threshold
# Inputs:
# structure: Structure object
# threshold: Threshold
# excluded: Excluded atoms
# Returns: |R| x |R| matrix |R| = # of residues
def compute_ave_atom_contact(structure, threshold, excluded=[]):
  res_count = 0
  # Get residue counts on all chain
  for chain in structure[0]:
    res_count += len(chain)
  struct_contact_matrix = np.zeros((res_count,res_count))
  for model in structure:
    mod_contact_matrix = compute_atom_contact(model, threshold, excluded)
    struct_contact_matrix += mod_contact_matrix
  ave_struct_contact_matrix = struct_contact_matrix / len(structure)

  return ave_struct_contact_matrix


def compute_max_residue(structure, ave_struct_contact_matrix):
  res_dict = {}
  i = 0
  for chain in structure[0]:
    for res in chain:
      contact_atoms = np.sum(ave_struct_contact_matrix[i])
      i+=1
      if res.get_resname() not in res_dict:
        res_dict[res.get_resname()] = contact_atoms
      elif res_dict[res.get_resname()] < contact_atoms:
        res_dict[res.get_resname()] = contact_atoms
  return res_dict


def compute_normalization_factor(pdb_dir, threshold, output):
  parser = PDBParser()
  tot_files = 0.0
  res_dir_dict = {}
  for filename in glob.iglob(pdb_dir + "/*.pdb"):
    tot_files += 1
    label = filename.split('/')[-1].split('.')[0]
    print label
    structure = parser.get_structure(label,filename)
    ave_struct_contact_matrix = compute_ave_atom_contact(structure, 2)
    res_dict = compute_max_residue(structure, ave_struct_contact_matrix)
    for key in res_dict:
      if key not in res_dir_dict:
        res_dir_dict[key] = res_dict[key]
      else:
        res_dir_dict[key] += res_dict[key]
  for key in res_dir_dict:
    res_dir_dict[key] = res_dir_dict[key]/tot_files
  with open(output, "w") as outfile:
    csvfile = csv.writer(outfile)
    for key in res_dir_dict:
      csvfile.writerow([key,res_dir_dict[key]])
  print res_dir_dict




def main(argv):
  pdb_dir = ''
  norm_fact = ''
  output = ''
  norm = 1
  interaction_thresh = 4
  contact_thresh = 4.5
  normalize = False

  bb = ['N','CA','C','O']

  try:
    opts, args = getopt.getopt(argv,"hnp:f:c:i:o:",["pdb_dir=",
                                                    "norm_fact=",
                                                    "interaction_thresh=",
                                                    "contact_thresh=",
                                                    "normalize",
                                                    "output="])
  except getopt.GetoptError:
    print 'blast_parser.py -i <input fasta> -x <blast xml> -o <output>'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'Computes a network given a threshold or computes normalization factor'
      print '''PSN.py 
                      -p --pdb_dir\t\tpdb directory 
                      -f --norm_fact\t\tnormalization factor file
                      -i --interaction_thresh\tcomma separated interaction strengths
                      -c --contact_thresh\tcontact threshold
                      -n --normalize\t\tComputes normalization factor of input pdbs
                      -o --output\t\tName of output'''
      sys.exit()
    elif opt in ("-p", "--pdb_dir"):
      pdb_dir = arg
    elif opt in ("-f", "--norm_fact"):
      norm_fact = arg
    elif opt in ("-c", "--contact_thresh"):
      
      contact_thresh = string_to_int(arg)
    elif opt in ("-i", "--interaction_thresh"):
      interaction_thresh = []
      for i in arg.split(","):
        interaction_thresh.append(string_to_int(i))
    elif opt in ("-o", "--output"):
      output = arg
    elif opt in ("-n","--normalize"):
      normalize = True

  # parser = PDBParser()
  # structure = parser.get_structure('1akg','../native/1akg.pdb')
  # a = compute_ave_atom_contact(structure, contact_thresh)
  compute_normalization_factor(pdb_dir, contact_thresh, output)



if __name__ == "__main__":
  main(sys.argv[1:])
