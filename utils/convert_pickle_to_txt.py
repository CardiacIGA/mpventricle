## To  be written
# Might be useful to have a human-readable file unlike the current .pickle files

import os, pickle
import numpy as np

def load_pickle(filename : str = 'filename'):
    '''load_pickle
    Load an existing 'filename.pickle' file and output its results.
    '''
    with open(filename + '.pickle', 'rb') as f: 
      return pickle.load(f) 
    
def main(file : str, name : str, filewrite : str):
    cps, w, patchverts, patches, nelems, knotval, knotmult, boundaries = load_pickle(file)
    Title   = f"{name} multipatch geometry data"
    Data    = (cps, w, patchverts, patches, nelems, knotval, knotmult, boundaries)
    Headers = "Control points", "Weights", "Patch vertices", "Patch connectivity", \
              "Number of elements per boundary", "Knot values per boundary", \
              "Knot multiplicity per boundary", "Boundary names"

    with open(filewrite+".txt", 'w') as f:
       f.write(Title+"\n\n")
       for head, data in zip(Headers,Data):
           f.write(head+"\n")
           if type(data) == np.ndarray:
               lines = "\n".join( [ str(row.tolist()) for row in data ] )
           else:
              lines = ""
           f.write(lines) 
 
           if type(data) == dict:
              for key, value in data.items(): 
                  f.write('%s:%s\n' % (key, value))
           f.write("\n\n")        
    return

if __name__ == "__main__":
   geomVar = True
   if geomVar:
    files  = "BV_GEOMETRY_DATA_PAPER_REF", "BV_GEOMETRY_DATA_PAPER_LONG", "BV_GEOMETRY_DATA_PAPER_THK"
    names  = "Bi-ventricle (reference variant)", "Bi-ventricle (elongated variant)", "Bi-ventricle (thickness variant)"
    folder_read  = "examples/output (geom variations)/pickle" 
    folder_write = "examples/output (geom variations)/txt" 
   else:
    files  = "BV_GEOMETRY_DATA", "LV_GEOMETRY_DATA"
    names  = "Bi-ventricle", "Left ventricle" 
    folder_read  = "examples/output/pickle"
    folder_write = "examples/output/txt"
   for file, name in zip(files,names):
    main(os.path.join(folder_read,file),name,os.path.join(folder_write,file))                