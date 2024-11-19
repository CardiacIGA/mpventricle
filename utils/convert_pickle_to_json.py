## Converts (existing) pickle files to json files. Script was used when building the repo
# Might be useful to have a human-readable file unlike the current .pickle files

import os, pickle, json
import numpy as np

def load_pickle(filename : str = 'filename'):
    '''load_pickle
    Load an existing 'filename.pickle' file and output its results.
    '''
    with open(filename + '.pickle', 'rb') as f: 
      return pickle.load(f) 
    
def remap_keys(mapping, convListfloat=False, convValfloat=False):
    if convValfloat:
        return [{'key':[int(i) for i in k ], 'value': v} for k, v in mapping.items()]    
    else:
       return [{'key':k, 'value': v} for k, v in mapping.items()]
    
def main(file : str, name : str, filewrite : str):
    cps, w, patchverts, patches, nelems, knotval, knotmult, boundaries = load_pickle(file)

    # Convert np.arrays to lists
    cps = cps.tolist()
    w   = w.tolist()
    patches    = patches.tolist()
    patchverts = patchverts.tolist()
    
    # Assign everything to one dictionary
    Data = {}
    Data["cps"] = cps
    Data["w"]   = w
    Data["patches"]    = patches
    Data["patchverts"] = patchverts
    Data["nelems"]  = remap_keys(nelems, convValfloat=True)
    Data["knotval"] = remap_keys(knotval, convValfloat=True)
    Data["knotmult"]   = remap_keys(knotmult, convValfloat=True)
    Data["boundaries"] = boundaries
    json_object = json.dumps(Data, indent=4)

    with open(filewrite + ".json", "w") as outfile:
        outfile.write(json_object)
  
    return

if __name__ == "__main__":
   geomVar = True
   if geomVar:
    files  = "BV_GEOMETRY_DATA_PAPER_REF", "BV_GEOMETRY_DATA_PAPER_LONG", "BV_GEOMETRY_DATA_PAPER_THK"
    names  = "Bi-ventricle (reference variant)", "Bi-ventricle (elongated variant)", "Bi-ventricle (thickness variant)"
    folder_read  = "examples/output (geom variations)/pickle" 
    folder_write = "examples/output (geom variations)/json" 
   else:
    files  = "BV_GEOMETRY_DATA", "LV_GEOMETRY_DATA"
    names  = "Bi-ventricle", "Left ventricle" 
    folder_read  = "examples/output/pickle"
    folder_write = "examples/output/json"
   for file, name in zip(files,names):
    main(os.path.join(folder_read,file),name,os.path.join(folder_write,file))   