# RINGSS Restore Profile
[![RINGSS Restore Profile](https://github.com/soar-telescope/ringss_restore_profile/actions/workflows/python-package.yml/badge.svg)](https://github.com/soar-telescope/ringss_restore_profile/actions/workflows/python-package.yml)
[![Upload to PYPI](https://github.com/soar-telescope/ringss_restore_profile/actions/workflows/python-publish.yml/badge.svg)](https://github.com/soar-telescope/ringss_restore_profile/actions/workflows/python-publish.yml)

## Overview

Restores atmospheric turbulence profile


## About RINGSS
Short description +links


# Usage

## First commit
Starting from code5.tar.gz from Andrei Tokovinin email

### Project requirements:
1) par_tololo.json: added from Andrei's email at May 24, 202 (par-mar21.json has incorrect pixel scale)
2) readstm.py: added from Andrei's email at May 20, 2022
3) profrest5.py: updated from Andrei's email at May 20, 2022
4) starcat.json: added from Andrei's email May 28, 2021
5) zmat.json: matrix constant

Notes:
1) **test9b.py**: it **IS** actually asi290mm_imageacq project
2) **cube2.py**: it is **IN** asi290mm_imageacq

### Python Requirements:
1) zernike package, done.
2) h5py required by zernike, done.

### Running profile restoration
#### Parameters file par-*.json
par-tololo.json is the initial parameters file provided by Andrei 
It needs to be modified to use site parameters,  
from this: "site":{"name":"LaSerena","lon":-71.2425,"lat":-29.91744},
to this:   "site":{"name":"CTIO","lon":-70.8065,"lat":-30.1684},

#### Weights File weights.json
The file is created by getweight5.py

These are the weights used by profrest5.py 

*[ebustos@localhost pythonProject]$ python3 getweight5.py par-tololo.json*

#### zmat.json
This is a "constant". It is an empirical model that depends on the optics used (telescope)

#### Calculating profiles *.prof
[ebustos@localhost ringss_restore_profile]$ python3 readstm.py par-tololo.json 2022-05-18

