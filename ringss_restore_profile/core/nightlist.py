#!/usr/bin/env python
# coding: utf-8

import sys
import os

from ringss_restore_profile.restore_profile.core import read_json


# ### Main module. usage: > python <par> <data>
# ---------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python nightlist.py <directory>")
        sys.exit()
        
    dir = sys.argv[1]
    print(dir)
    # list of files to process
    files =  [f for f in os.listdir(dir) if f.endswith('.json')]
    n = len(files)
    print("Files to process: ",n)
    if (n == 0):
        sys.exit()
    # Make the list!
    
    list = open(dir+"list.txt", 'w')
    for f in files:
        data = read_json(dir+f)
        profile = data["profile"]
        sp = ""  # print profile in readable form, in 10^(-13) units
        for p in profile["prof"]:
            sp += " {:.2f}".format(p)
    
        s0 = (f[0:17] + " {} {:.2f}").format(data["starpar"]["HR"],data["starpar"]["zen"])
        eladu = 3.60*pow(10, -data["cubepar"]["gain"]/200 )
        flux = eladu*data["image"]["impar"][1] # flux in electrons
        s0 += " {:.3e} {:.2f}".format(flux,data["image"]["impar"][3] ) # flux, ring radius
        s0 += "  {:.2f}  {:.2f}  {:.2f}".format(profile["see2"],profile["see"],profile["fsee"]) # seeing values
        s0 += "  {:.3f}  {:.3f}  {:.2f}".format(profile["totvar"],profile["erms"],profile["wind"]) # scint,rms, wind
        print(s0, sp)  
        list.write(s0+sp+'\n')
    list.close()
    print(" Created the listing: "+dir+"list.txt")
