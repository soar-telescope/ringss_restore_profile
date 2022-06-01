# read .stm file, return data structure of the first line. 
# Returns the data dictionary with moments, impar, noisepar, starpar
# see cube2.py for the data structure:
# moments = {'var': power.tolist(), 'cov': cov.tolist(),'rnoise': drnoise[0],'rvar': meanrvar,'mcoef':mcoef.tolist()}
# data = {'image': {'impar': impar, 'noisepar': noisepar}, 'moments': moment
# cubepar = {'nx': nx,'nz':nz,'texp':texp,'gain':gain,'date':DATE,'star':star}

import sys
from os.path import splitext
import numpy as np
from scipy import optimize
from ringss_restore_profile.restore_profile.star import Star
from ringss_restore_profile.restore_profile.core import read_json


class ReadMoments(object):
    def __init__(self, ringss_parameters_filename: str):
        self.parameters = read_json(ringss_parameters_filename)
        self.star_catalog = read_json(self.parameters["profrest"]["starcat"])
        self.weights = read_json(self.parameters["profrest"]["weightsfile"])
        self.zmat = read_json(self.parameters["profrest"]["zmat"])
        self.star = Star(ringss_parameters_filename)

    def parse_line(self, moments_line: str) -> dict:
        moments_list = moments_line.split(",")

        nz = int(moments_list[5])
        texp = float(moments_list[2]) * 1e-3
        gain = float(moments_list[4])
        nx = int(moments_list[6])
        cubepar = {'nx': nx, 'nz': nz, 'texp': texp, 'gain': gain, 'date': moments_list[0],
                   'star': moments_list[1].strip()}

        mcoef = 20  # fixed for now

        k0 = 9
        impar = {'backgr': moments_list[k0], 'flux': moments_list[k0 + 1], 'fluxvar': moments_list[k0 + 2],
                 'meanrad': moments_list[k0 + 3], 'rwidth': moments_list[k0 + 4],
                 'xc': moments_list[k0 + 5], 'yc': moments_list[k0 + 6],
                 'xcvar': moments_list[k0 + 7], 'ycvar': moments_list[k0 + 8],
                 'coma': moments_list[k0 + 9], 'angle': moments_list[k0 + 10], 'contrast': moments_list[k0 + 11]}
        noisepar = moments_list[k0 + 12:k0 + 16]

        m = mcoef + 1
        k1 = k0 + 16
        var = moments_list[k1:k1 + m]
        cov = moments_list[k1 + m:k1 + 2 * m]
        k2 = k1 + 2 * m
        rnoise = moments_list[k2]
        rrms = moments_list[k2 + 1]
        meancoef = moments_list[k2 + 2:k2 + 16]
        momentos = {'var': var, 'cov': cov, 'rnoise': rnoise, 'rvar': rrms, 'mcoef': meancoef}
        data_dict = {'image': {'impar': impar, 'noisepar': noisepar}, 'moments': momentos, 'cubepar': cubepar}

        return data_dict

    # read & process file with averaged moments only
    # based in read_stm()
    # data_filename, starcat, weights and zmatrix are inputs
    def read_from_file(self, data_filename: str):
        try:
            fh_in = open(data_filename, 'r')  # input .stm file
        except FileNotFoundError as err:
            print(err)
            return None

        filename, file_extension = splitext(data_filename)
        try:
            fh_out = open(filename + '.prof', 'w')  # output profile file
        except FileNotFoundError as err:
            print(err)
            return None

        for line in fh_in:
            data = self.parse_line(line)

            #print(data)
            #print(type(data["cubepar"]["star"]))

            # Star
            star_par = self.star.get_starpar(data["cubepar"]["date"], data["cubepar"]["star"])

            data["starpar"] = star_par

            profile = self.restore(data)  # returns profile or None

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
                fh_out.write(s + '\n')
        fh_in.close()
        fh_out.close()
        print('Processing finished!')

    # Profile restoration. Inputs: dictionaries of parameters, data, weights, and Z-matrix
    # Output: dictionary of profile parameters
    # Before calling restore(), run getzen.py to define the zenith distance and star color in data
    #
    def restore(self, data: dict) -> dict:
        #    print(data["image"]["impar"])
        var = data["moments"]["var"]  # variance of a-coefficients
        cov = data["moments"]["cov"]  # covariance of a-coefficients
        # impar = data["image"]["impar"]  # ring parameters
        # noisepar = data["image"]["noisepar"]  # noise parameters
        starpar = data["starpar"]  # star parameters
        zen = starpar["zen"]
        bv = starpar["BV"]
        z0 = self.parameters["profrest"]["zgrid"]  # grid of heights representing the turbulence profile

        z = np.array(self.weights["z"])  # nominal distance grid of weights
        nz = len(self.weights["z"])  # number of layers in weights
        nm = len(self.weights["wt0"][0])  # number of coefficients per layer
        wt = np.reshape(np.array(self.weights["wt0"]), (nz, nm))  # wt.shape = (16,21)
        wt += bv * np.reshape(np.array(self.weights["wtslope"]), (nz, nm))  # must be >0!

        # interpolate weight to the Z0 grid, trim to mmax
        cosz = np.cos(zen / 180 * np.pi)  # cos(z)
        mmax = 15  # hard-coded number of terms to use
        nz0 = len(z0)
        z00 = np.array(z0) / cosz  # stretched array of nominal heights
        wt1 = np.zeros((nz0, mmax))  # interpolated and trimmed weights
        for m in range(0, mmax):
            wt1[:, m] = np.interp(z00, z, wt[:, 1 + m])
        wtsect = np.interp(z00, z, wt[:, 0])  # sector-motion weight on z0 grid

        # Check that moments and weights have the same number of coefficients
        mcoef = len(var)  # 21 for mmax=20
        if mcoef != nm:
            print("Error! Mismatching number of coefficients: ", mcoef, nm)
            return

        # noise bias, see allcubes5.pro => noisecubes
        gain = data["cubepar"]["gain"]  # electrons per ADU
        eladu = 3.60 * pow(10, -gain / 200)
        noisepar = data["image"]["noisepar"]  # list of 4 numbers
        fluxadu = float(data["image"]["impar"]["flux"])
        #    print(fluxadu)
        flux = eladu * fluxadu  # flux in electrons

        anoise = float(noisepar[0]) / flux + float(noisepar[1]) * pow(self.parameters["telescope"]["ron"] / flux,
                                                                      2)  # noise variance of a-coef
        rnoise = 2. * (float(noisepar[2]) / flux + float(noisepar[3]) * pow(self.parameters["telescope"]["ron"] / flux,
                                                                            2) / 8)  # noise on radius, pix^2

        #    anoise = noisepar[0]/flux + noisepar[1]*pow(par["telescope"]["ron"]/flux,2) # noise variance of a-coef
        #    rnoise = 2.*(noisepar[2]/flux + noisepar[3]*pow(par["telescope"]["ron"]/flux,2)/8) # noise on radius, pix^2

        # Select max frequency used in the restoration
        var1 = np.array(var[1:mmax + 1], float) - anoise
        var1 *= (var1 > 0)  # prevent negative
        cov1 = np.array(cov[1:mmax + 1], float)
        rho = cov1 / var1
        varcorr = var1 / (0.8 + 0.2 * rho)  # correction for finite exposure time
        totvar = np.sum(np.array(var, float))  # full scintillation power, incl. noise

        # Z-matrix correction. Avoid negative values!
        z1 = np.array(self.zmat)
        ncol = z1.shape[1]
        var2 = np.array(var[1:ncol + 1], float)  # indices used for correction
        varz = varcorr / (1. + np.dot(z1[0:mmax, :], var2))  # correct for saturation

        # weighted nnls. Weight is proportional to 1/var (before noise subtraction, non-negative)
        varwt = np.power(np.array(var[1:mmax + 1], float), -1)
        # varwt = np.power(np.array(var[1:mmax+1],float), -0.5)
        a2 = np.transpose(wt1.copy())  # matrix of (mmax-1,nz) dimension
        for i in range(0, mmax):  # weighted system matrix
            a2[i, :] *= varwt[i]
        varz2 = varz * varwt  # weighted right-hand vector
        prof, resvar = optimize.nnls(a2, varz2)  # turbulence integrals in  [m^1/3] in z0 layers

        varmod = np.dot(a2, prof) / varwt
        erms = np.std(1. - varmod / varz)
        print("RMS residual: {:.3f}".format(erms))

        prof1 = prof * cosz  # zenith-corrected integrals
        jtot = np.sum(prof1)  # turbulence integral
        seeconst = 6.83e-13  # integral for 1" seeing at 500nm
        see = pow(jtot / seeconst, 0.6)  # seeing in arcseconds
        jfree = np.sum(prof1[2:nz0 - 1])  # start at 0.5km
        fsee = pow(jfree / seeconst, 0.6)  # FA seeing in arcseconds
        prof1 = prof * cosz * 1e13  # scaled

        # Compute seeing from radius var.
        d = self.parameters["telescope"]["D"]
        pixel = self.parameters["telescope"]["pixel"]
        lameff = self.weights["lameff"]
        lam0 = lameff[0] + bv * (lameff[1] - lameff[0])  # effective wavelength for the star color
        rvar = float(data["moments"]["rvar"]) - rnoise  # radius variance in pix^2, noise-subtracted
        lamd = lam0 / d * 206265.  # lambda/D in arcseconds
        rvarnorm = rvar * pow(pixel / lamd, 2)  # variance in (lam/D^2) units

        wcoef = np.sum(wtsect * prof) / np.sum(prof)  # profile-adjusted weight of sector variance
        jtot2 = rvarnorm / wcoef / 4  # turbulence integral, m^1/3. Explain factor 4!
        see2 = pow(jtot2 * cosz / seeconst, 0.6)  # seeing at zenith, arcsecs
        see2 *= 1. / (1. - 0.4 * totvar)  # saturation correction
        print("Seeing (sect,tot,FA): {:.3f} {:.3f} {:.3f}".format(see2, see, fsee))

        # Wind measurement
        texp = 1e-3  # hard-coded exposure time
        delta = (var1 - cov1) * (texp ** (-2))
        delta = delta * (delta > 0)
        ucoef = np.array(self.weights["ucoef0"]) + bv * np.array(self.weights["ucoefslope"])
        umm = np.array(self.weights["umm"], int) - 1
        delta = np.array(delta[umm])
        # print(delta)  # debugging
        v2mom = np.sum(ucoef * delta)
        jwind = np.sum(prof[2:nz0])  # exclude ground layer, the speed is not zen-corrected
        #    jwind = np.sum(prof[1:nz0])  # exclude ground layer, the speed is not zen-corrected
        v2 = pow(v2mom / jwind, 0.5)  # use profile uncorrected for zenith
        print("Wind speed [m/s]: {:.3f}".format(v2))

        # tau_0 at 500nm
        r0 = pow(6.680e13 * jwind, -0.6)
        tau0 = 310. * r0 / v2  # in [ms]
        # isoplanatic angle at 500nm at zenith
        r0 = 0.101 / see
        tmp = pow(np.array(z0), 1.6667)
        heff = pow(np.sum(tmp * prof1) / np.sum(prof1), 0.6)
        theta0 = 205265. * 0.31 * r0 / heff
        # Output dictionary
        profile_dict = {"z0": z0, "prof": prof1.tolist(), "see": see, "fsee": fsee, "see2": see2, "wind": v2, "erms": erms,
                        "totvar": totvar, "tau0": tau0, "theta0": theta0}

        return profile_dict


# ### Main module. usage: > python <par> <data>
# ---------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python read_moments.py <par-file.json> <file.stm>")
        sys.exit()

    parameters_filename = sys.argv[1]
    moments = ReadMoments(parameters_filename)

    stm_filename = sys.argv[2]
    print('STM file: ' + stm_filename)

    # Ingest all information needed

    moments.read_from_file(stm_filename)
