# code from Mike that interpolates over tables of photometric models
# such as AMES-COND, BT-SETTI, SONORA, and the like, allowing us to
# estimate planetary masses from a given contrast

import re

import numpy as np
from astropy.io import fits
from astropy.table import Table
import pandas as pd
from scipy.interpolate import interp2d

__all__ = ["MassModel", "SONORAModel"]


def Mj_to_solar(mass_MJ):
    return mass_MJ / 1047.9


def solar_to_Mj(mass_ms):
    return mass_ms * 1047.9


def distance_modulus(dist_pc):
    return 5.0 * np.log10(dist_pc) - 5


def extract_first_number_from_string(stringy):
    allnums = re.findall(r"-?\d+\.?\d*", stringy)
    firstnum = float(allnums[0])
    return firstnum


def get_headers(headerline):
    # return the line column labels
    # unfucking the M/MsTeff thing
    return list(filter(None, headerline.replace("Teff", " Teff").split(" ")))


def split_row(row):
    bb = filter(None, row.split(" "))
    return bb


def parse_allard_file(fname):
    with open(fname) as f:
        content = f.readlines()
        content = [x.strip() for x in content]

    output = []
    age = None
    for line in content:
        if "Gyr" in line:
            age = extract_first_number_from_string(line)
            continue
        if "M/Ms" in line:
            headers = get_headers(line)
        else:
            # split the row by spaces
            row = split_row(line)
            try:
                numrow = [float(x) for x in row]
                # if not empty
                if numrow:
                    output.append([1000 * age] + numrow)
            except:
                pass
    tab = Table(np.array(output), names=["Age(Myr)"] + headers)
    return tab.to_pandas()


class MassModel:
    def __init__(self, filename):
        self.table = parse_allard_file(filename)

    def get_below(self, item, amt, subtable=None):
        # returns a subtable
        # eg get_below('Age(Myr)', 1000) will return
        # a subtable of the closest age at or below 1000 Myrs
        if subtable is None:
            subtable = self.table
        ret = subtable.loc[subtable[item] <= amt][item].max()
        if np.isnan(ret):
            ret = subtable[item].min()
        return subtable.loc[subtable[item] == ret]

    def get_above(self, item, amt, subtable=None):
        # returns a subtable
        # eg get_above('Age(Myr)', 1000) will return
        # a subtable of the closest above 1000 Myrs
        if subtable is None:
            subtable = self.table
        ret = subtable.loc[subtable[item] > amt][item].min()
        if np.isnan(ret):
            ret = subtable[item].max()
        return subtable.loc[subtable[item] == ret]

    def interp_between(
        self, prop_0=None, val_0=None, prop_1=None, val_1=None, prop_2=None
    ):
        """interpolate between any three properties.
        prop0,val_0 -> dependent variable 0, with value val_0
        prop1,val_1 -> dependent variable 1, with value val_1
        prop2,      -> output variable.  for example:
        m = MassModel('/Users/mbottom/Desktop/Coronagraphy/MassModels/model.BT-Settl.M-0.0.MKO.AB.txt')
        m.interp_between(prop_0 = 'Age(Myr)', val_0 = 1000, prop_1 = "L'", val_1 = 7.091, prop_2 = 'M/Ms')
        array([0.8])
        says that at 1000 Myrs, for an L' value of 7.091, you would expect the mass to be 0.8 M*
        """

        prop_0_below_table = self.get_below(prop_0, val_0)
        x0 = prop_0_below_table[prop_0].iloc[0]
        x1 = x0.copy()

        sub0 = self.get_below(prop_1, val_1, subtable=prop_0_below_table)
        y0 = float(sub0[prop_1])
        z0 = float(sub0[prop_2])

        sub1 = self.get_above(prop_1, val_1, subtable=prop_0_below_table)
        y1 = float(sub1[prop_1])
        z1 = float(sub1[prop_2])

        prop_0_above_table = self.get_above(prop_0, val_0)
        x2 = prop_0_above_table[prop_0].iloc[0]
        x3 = x2.copy()

        sub2 = self.get_below(prop_1, val_1, subtable=prop_0_above_table)
        y2 = float(sub2[prop_1])
        z2 = float(sub2[prop_2])

        sub3 = self.get_above(prop_1, val_1, subtable=prop_0_above_table)
        y3 = float(sub3[prop_1])
        z3 = float(sub3[prop_2])

        f = interp2d([x0, x1, x2, x3], [y0, y1, y2, y3], [z0, z1, z2, z3])
        return float(f(val_0, val_1))

    def band_abbrevs(self, band):
        if band in ["Lp", "L"]:
            theband = "L'"
        elif band in ["M", "Mp", "Ms"]:
            theband = "M'"
        elif band in ["Kp"]:
            theband = "Ks"
        else:
            theband = band
        return theband

    def mass_to_absolute_mag(self, mass_MJ=None, band=None, age_Myr=None):
        theband = self.band_abbrevs(band)
        mass_ms = Mj_to_solar(mass_MJ)
        return self.interp_between(
            prop_0="Age(Myr)",
            val_0=age_Myr,
            prop_1="M/Ms",
            val_1=mass_ms,
            prop_2=theband,
        )

    def mass_to_apparent_mag(self, mass_MJ=None, band=None, age_Myr=None, dist_pc=None):
        absmag = self.mass_to_absolute_mag(mass_MJ=mass_MJ, band=band, age_Myr=age_Myr)
        dist_modulus = distance_modulus(dist_pc)
        return absmag + dist_modulus

    def mass_to_contrast(
        self,
        mass_MJ=None,
        band=None,
        age_Myr=None,
        stellar_apparent_mag=None,
        dist_pc=None,
    ):
        planet_appmag = self.mass_to_apparent_mag(
            mass_MJ=mass_MJ, band=band, age_Myr=age_Myr, dist_pc=dist_pc
        )
        return planet_appmag - stellar_apparent_mag

    def absolute_mag_to_mass(self, absmag=None, age_Myr=None, band=None):
        theband = self.band_abbrevs(band)
        mass_ms = self.interp_between(
            prop_0="Age(Myr)",
            val_0=age_Myr,
            prop_1=theband,
            val_1=absmag,
            prop_2="M/Ms",
        )
        return solar_to_Mj(mass_ms)

    def apparent_mag_to_mass(self, appmag=None, age_Myr=None, band=None, dist_pc=None):
        abs_mag = appmag - distance_modulus(dist_pc)
        return self.absolute_mag_to_mass(absmag=abs_mag, age_Myr=age_Myr, band=band)

    def contrast_to_mass(
        self,
        deltaMag=None,
        age_Myr=None,
        band=None,
        dist_pc=None,
        stellar_apparent_mag=None,
    ):
        planet_app_mag = deltaMag + stellar_apparent_mag
        return self.apparent_mag_to_mass(
            appmag=planet_app_mag, age_Myr=age_Myr, band=band, dist_pc=dist_pc
        )


class SONORAModel(MassModel):
    def __init__(self, filename, condgrid):
        super().__init__(condgrid)
        self.sonoragrid = pd.read_csv(
            filename, header=0, delim_whitespace=True, comment="#"
        )

    def get_below_sonora(self, item, amt, subtable=None):
        if subtable is None:
            subtable = self.sonoragrid
        ret = subtable.loc[subtable[item] <= amt][item].max()
        if np.isnan(ret):
            ret = subtable[item].min()
        return subtable.loc[subtable[item] == ret]

    def get_above_sonora(self, item, amt, subtable=None):
        if subtable is None:
            subtable = self.sonoragrid
        ret = subtable.loc[subtable[item] > amt][item].min()
        if np.isnan(ret):
            ret = subtable[item].max()
        return subtable.loc[subtable[item] == ret]

    def interp_sonora(self, prop_0, val_0, prop_1, val_1, prop_2):
        prop_0_below_table = self.get_below_sonora(prop_0, val_0)
        x0 = prop_0_below_table[prop_0].iloc[0]
        x1 = x0.copy()

        sub0 = self.get_below_sonora(prop_1, val_1, subtable=prop_0_below_table)
        y0 = float(sub0[prop_1])
        z0 = float(sub0[prop_2])

        sub1 = self.get_above_sonora(prop_1, val_1, subtable=prop_0_below_table)
        y1 = float(sub1[prop_1])
        z1 = float(sub1[prop_2])

        prop_0_above_table = self.get_above_sonora(prop_0, val_0)
        x2 = prop_0_above_table[prop_0].iloc[0]
        x3 = x2.copy()

        sub2 = self.get_below_sonora(prop_1, val_1, subtable=prop_0_above_table)
        y2 = float(sub2[prop_1])
        z2 = float(sub2[prop_2])

        sub3 = self.get_above_sonora(prop_1, val_1, subtable=prop_0_above_table)
        y3 = float(sub3[prop_1])
        z3 = float(sub3[prop_2])

        f = interp2d([x0, x1, x2, x3], [y0, y1, y2, y3], [z0, z1, z2, z3])
        return float(f(val_0, val_1))

    def absolute_mag_to_mass(self, absmag=None, age_Myr=None, band=None):
        theband = self.band_abbrevs(band)
        teff = self.interp_between(
            prop_0="Age(Myr)",
            val_0=age_Myr,
            prop_1=theband,
            val_1=absmag,
            prop_2="Teff(K)",
        )
        mass_ms = self.interp_sonora(
            prop_0="age(Gyr)",
            val_0=age_Myr / 1000,
            prop_1="Teff(K)",
            val_1=teff,
            prop_2="M/Msun",
        )
        return solar_to_Mj(mass_ms)
