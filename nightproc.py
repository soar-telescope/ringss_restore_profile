# coding: utf-8


import os
import sys
# import cube2
import getzen
import getweight5
import profrest5
import json
import codecs

# Process one night with cubes in indir

#indir =  '/home/andrei/ASICAP/CapObj/20210321/'
#outdir = indir + 'outpython/'
#   par = profrest5.read_json('par-mar21.json')

def Nightproc(indir, outdir, parfile):

    # list of files to process
    #files =  [f for f in os.listdir(indir) if f.endswith('_cube.fits')]
    #n = len(files)
    #print("Files to process: ",n)

    par = profrest5.read_json(parfile)

    # getweight5.computeweight(par) # saves in weight.json

    # Crunch all data cubes
    #datafiles = []
    #for cubefile in files:
    #    datafile = cube2.cubeproc(cubefile,indir,outdir) # Creates the data dictionary in outdir
    #    print("Created "+datafile)

    # list of files to process
    files =  [f for f in os.listdir(outdir) if f.endswith('.json')]
    n = len(files)
    print("Files to process: ",n)

    starcat = profrest5.read_json(par["profrest"]["starcat"])  # Catalog
    weight = profrest5.read_json(par["profrest"]["weightfile"])
    zmat = profrest5.read_json(par["profrest"]["zmat"])

    # Compute zenith distance
    for f in files:
        data = profrest5.read_json(outdir+f)
        hr = data["cubepar"]["star"] 
        star = getzen.getstar(hr,starcat)
        if star[0] > 0:
            starpar = getzen.getstarpar(par,data,star,hr)
        # print(starpar)
        else:
            starpar='none'
            data["starpar"] = starpar
        profile = profrest5.Restore(par, data, weight, zmat) # returns profile or None
        data["profile"] = profile # add to the data dictionary
    # Save the data with starpar and profile    
    json.dump(data, codecs.open(outdir+f, 'w', encoding='utf-8'), separators=(',', ':'))
    print("Finished night processing!")

    
# ### Main module. usage: > python <par> <data>
# ---------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python nightproc.py <indir> <outdir> <par-file>")
        sys.exit()

    indir = sys.argv[1]
    outdir = sys.argv[2]
    parfile = sys.argv[3]

    Nightproc(indir,outdir,parfile)
    

