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


class Observations:
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