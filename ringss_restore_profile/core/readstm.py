# read .stm file, return data structure of the first line. 
# Returns the data dictionary with moments, impar, noisepar, starpar
# see cube2.py for the data structure:
# moments = {'var': power.tolist(), 'cov': cov.tolist(),'rnoise': drnoise[0],'rvar': meanrvar,'mcoef':mcoef.tolist()}
# data = {'image': {'impar': impar, 'noisepar': noisepar}, 'moments': moment
# cubepar = {'nx': nx,'nz':nz,'texp':texp,'gain':gain,'date':DATE,'star':star}

import sys
from os.path import splitext
from .getzen import get_star, get_starpar
from .profrest5 import restore
from ringss_restore_profile.restore_profile.core import read_json


# Parse O-line, return cubepar dictionary
def parseO(line: list):
    nz = int(line[3])
    texp = float(line[4]) * 1e-3
    gain = float(line[5])
    nx = int(line[7])
    cubepar = {'nx': nx, 'nz': nz, 'texp': texp, 'gain': gain, 'date': line[1], 'star': line[2]}
    #    print(cubepar)
    return cubepar


# Parse M-line of .stm file, return the data dictionary
def parseM(line: list, cubepar):
    mcoef = 20  # fixed for now
    date = line[1][1:20]  # remove first space, insert 'T'
    date1 = date[0:10] + "T" + date[11:19]
    cubepar['date'] = date1
    cubepar['star'] = line[2]
    k0 = 4
    impar = {'backgr': line[k0], 'flux': line[k0 + 1], 'fluxvar': line[k0 + 2], 'meanrad': line[k0 + 3],
             'rwidth': line[k0 + 4],
             'xc': line[k0 + 5], 'yc': line[k0 + 6], 'xcvar': line[k0 + 7], 'ycvar': line[k0 + 8],
             'coma': line[k0 + 9], 'angle': line[k0 + 10], 'contrast': line[k0 + 11]}
    noisepar = line[k0 + 12:k0 + 16]
    m = mcoef + 1
    k1 = k0 + 16
    var = line[k1:k1 + m]
    cov = line[k1 + m:k1 + 2 * m]
    k2 = k1 + 2 * m
    rnoise = line[k2]
    rrms = line[k2 + 1]
    meancoef = line[k2 + 2:k2 + 16]
    moments = {'var': var, 'cov': cov, 'rnoise': rnoise, 'rvar': rrms, 'mcoef': meancoef}
    data = {'image': {'impar': impar, 'noisepar': noisepar}, 'moments': moments, 'cubepar': cubepar}
    return data


# read and process .stm file
def read_stm(par: dict, stm: str):
    # Get starcat, weight, zmat
    starcat = read_json(par["profrest"]["starcat"])  # Catalog
    weights = read_json(par["profrest"]["weightsfile"])
    zmat = read_json(par["profrest"]["zmat"])

    try:
        f = open(stm, 'r')  # input .stm file
    except FileNotFoundError as err:
        print(err)
        return None

    filename, file_extension = splitext(stm)
    try:
        fout = open(filename + '.prof', 'w')  # output profile file
    except FileNotFoundError as err:
        print(err)
        return None

    for line in f:
        alist = line.split(",")
        tag = alist[0]
        if tag == 'O':
            cubepar = parseO(alist)
        if tag == 'M':
            data = parseM(alist, cubepar)
            hr = data["cubepar"]["star"]
            star = get_star(hr, starcat)
            #  Compute zenith distance
            if star[0] > 0:
                star_par = get_starpar(par, data, star, hr)
            else:
                star_par = 'none'
            data["starpar"] = star_par

            profile = restore(par, data, weights, zmat)  # returns profile or None
            if profile:
                # data["profile"] = profile # add to the data dictionary
                # Output
                s = data['cubepar']['date'] + ',' + data['cubepar']['star'] + ',{:.2f}'.format(data['starpar']['zen'])
                s = s + ',' + data["image"]["impar"]["flux"]
                s = s + ',  {:.3f}'.format(profile['see2']) + ',{:.3f}'.format(profile['see']) + ',{:.3f}'.format(
                    profile['fsee'])
                s = s + ',{:.2f}'.format(profile['wind']) + ',{:.3f}'.format(profile['tau0']) + ',{:.3f}'.format(
                    profile['theta0'])
                s = s + ',  {:.3f}'.format(profile['totvar']) + ',{:.3f}'.format(profile['erms'])
                for x in profile['prof']:
                    s += ",{:.2f}".format(x)
                print(data['cubepar']['date'])
                fout.write(s + '\n')
    f.close()
    fout.close()
    print('Processing finished!')


# ### Main module. usage: > python <par> <data>
# ---------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python readstm.py <par-file.json> <file.stm>")
        sys.exit()

    parameters_filename = sys.argv[1]
    stm_filename = sys.argv[2]
    print('STM file: ' + stm_filename)

    # Ingest all information needed
    parameters = read_json(parameters_filename)
    if parameters is None:
        sys.exit()

    read_stm(parameters, stm_filename)
