# coding: utf-8

# based on getzen.py
# ### Zenith distance and JD calculation. May 27, 2021, AT

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from ringss_restore_profile.restore_profile.core import read_json


class Star(object):
    def __init__(self, ringss_parameters_filename):
        self.parameters = read_json(ringss_parameters_filename)
        self.star_catalog = read_json(self.parameters["profrest"]["starcat"])

    # for a given date and HR compute JD and zenith distance and returns the starpar dictionary
    def get_starpar(self, date: str, hr_number: str) -> dict:
        star_by_hr = self.star_catalog[hr_number]
        time = Time(date, format='isot')
        jd = time.jd
        location = EarthLocation(lat=self.parameters["site"]["lat"] * u.deg, lon=self.parameters["site"]["lon"] * u.deg)
        coo = SkyCoord(ra=float(star_by_hr[0]) * u.degree, dec=float(star_by_hr[1]) * u.degree)
        alt_az = coo.transform_to(AltAz(obstime=time, location=location))
        zen = 90. - alt_az.alt.degree
        az = alt_az.az.degree
        star_parameters = {"JD": jd, "HR": int(hr_number), "Vmag": float(star_by_hr[2]), "BV": float(star_by_hr[3]), "zen": zen,
                           "az": az}  # starpar dictionary
        return star_parameters

