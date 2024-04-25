"""
Created on 2024/04/18
@author: Timo Brenninkmeyer
Description:
"""
from astropy.io import fits
from pathlib import Path
import numpy as np
from dataclasses import dataclass

_ = "Not a script, do not run"



class Hdu:
    def __init__(self, path: Path, index: int, header: dict, data: np.ndarray):
        self.path: Path = path
        self.index: int = index
        self.header: dict = header
        self.data: np.ndarray = data
    def __repr__(self):
        return f"<{self.path}> {self.header['TIME-OBS']}{self.index}"


class Observations:
    """
        WORKS WITH FILTERLESS DATA (SPECTOGRAPH), DOESNT DESTINGUISH BETWEEN FILTERS
    """
    def __init__(self, hdus: list):
        self.hdus = hdus
        d = {}
        print(hdus)
        for hdu in hdus:
            type = hdu.header.get("IMAGETYP", "NO_IMTYPE")
            assert not [hdu] is None
            if not d.get(type, None):
                d[type] = [hdu]
            print(d[type] + [hdu])
            d[type] = d[type] + [hdu]
        print(d)
        self.hdu_sep = d

    def include_type(self, type: str, start: int, stop: int):
        self.hdu_sep[type] = self.hdu_sep[type][start:stop]

    def create_masters(self):
        self.masterbias = self.get_masterbias()
        self.masterdark = self.get_masterdark()
        self.masterflat = self.get_masterflat()
        self.masterlight = self.get_masterlight()

    def get_masterbias(self):
        biasstack = []
        for bias in self.hdu_sep["Bias Frame"]:
            biasstack.append(bias.data)
        stack = np.stack(biasstack)
        return np.median(stack, axis=0)

    def get_masterdark(self):
        darkstack = []
        for dark in self.hdu_sep["Dark Frame"]:
            dark_adj = (dark.data - self.masterbias)/dark.header["EXPTIME"]
            darkstack.append(dark_adj)
        stack = np.stack(darkstack)
        return np.median(stack, axis=0)

    def get_masterflat(self):
        flatstack = []
        for flat in self.hdu_sep["Light Frame"]:
            flat_adj = (flat.data - self.masterbias - self.masterdark*flat.header["EXPTIME"])
            flat_median = np.median(flat_adj)
            flat_norm = flat_adj/flat_median
            flatstack.append(flat_norm)
        stack = np.stack(flatstack)
        return np.median(stack, axis=0)/np.median(stack)

    def get_masterlight(self):
        lightstack = []
        for light in self.hdu_sep["Light Frame"]:
            light_adj = ((light.data - self.masterbias - self.masterdark*light.header["EXPTIME"])
                         /np.where(self.masterflat!=0, self.masterflat, 0.001)) #getting /0 errors
            lightstack.append(light_adj)
        stack = np.stack(lightstack)
        return np.median(stack, axis=0)

    def display_final(self):
        pass

class FitsHandler:
    filetypes = [".FIT", ".fits"]
    def __init__(self,
                 folder_directory: str,
                 max_files: int=-1):
        self.directory = folder_directory
        self.files = list(Path(self.directory).iterdir())
        self.max_files = max_files


    def get_hdus(self):
        files = self.files.copy()
        hdus = []

        for file in files:
            if not file.suffix in self.filetypes:
                print(f"Skipping {file}")
                continue

            with fits.open(file) as hdulist:
                for i, hdu in enumerate(hdulist):
                    if len(hdus) > self.max_files and self.max_files != -1:
                        return hdus

                    new_hdu = Hdu(file, i, hdu.header, hdu.data)
                    hdus.append(new_hdu)
        return Observations(hdus)