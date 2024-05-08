"""


"""


import numpy as np
from numpy.typing import NDArray

class Spectrum:
    def __init__(self, array: NDArray[float], type_: str) -> None:
        self.array = array
        self.type = self.set_type(type_)
        self.flat_array = None


    def magic_numbers(self):
        """Forgive me father for I have fudged the dice"""
        testx = np.linspace(177, 1000, self.array.size)
        testx *= 1e-9
        testy = self.blackbody_curve(testx)
        offset = 350
        y = testy*1.1e-11+offset
        return y

    def get_peaks(self, *args, **kwargs) -> tuple[list[float], list[float]]:
        cut_array = self.peak_cut(*args, **kwargs)
        seperated_peaks = self.separate_peaks(cut_array)
        max_peaks = self.peak_maxima(seperated_peaks)
        self.peaks = max_peaks
        xcoords = []
        ycoords = []
        for val in max_peaks.values():
            xcoords.append(val[0])
            ycoords.append(val[1])
        return xcoords, ycoords

    def peak_cut_function(self):
        raise NotImplementedError("Abstract function, only cut peaks from a subclass of Spectrum")

    def peak_cut(self, *args, **kwargs) -> NDArray:
        data = self.flat_array
        cut_array = np.where(data>self.peak_cut_function(*args, **kwargs), data, np.nan)
        return cut_array

    def separate_peaks(self, array: NDArray) -> dict[int, list[tuple[float, int]]]:
        sections = {}
        was_nan = True
        i = 0
        for index, element in enumerate(array):
            if not np.isnan(element):
                if was_nan:
                    sections[i] = []
                    was_nan = False
                sections[i].append((element, index))
            else:
                if not was_nan:
                    i += 1
                was_nan = True
        return sections

    def peak_maxima(self, separated_array: dict[int, list[tuple[float, int]]]) -> dict[int, tuple[float, float]]:
        max_vals = {}
        for key in separated_array:
            val, indx = zip(*separated_array[key])
            val = np.array(val)
            indx = np.array(indx)
            m: float = np.max(val)
            i: float = indx[val.argmax()]
            max_vals[key] = (i, m)
        return max_vals

    def set_type(self, type: str) -> None:
        if not (type == 'E' or type == 'A'):
            raise ValueError(f"invalid input: {type}\nInput must be 'E' for emission or 'A' for absorption")
        self.type = type


class EmissionSpectrum(Spectrum):
    def __init__(self, array: NDArray[float]):
        super().__init__(array, 'E')
        self.flat_array: NDArray[float] = array

    def peak_cut_function(self, slope=0.07) -> NDArray:
        data = self.array
        x_max = data.argmax()
        shape = np.linspace(x_max, data.size, data.size-x_max)
        base = np.median(np.abs(data))
        start = 0.65 * data.max()
        right_curve = (start * np.exp(-(shape-x_max)/(slope*data.size)) + base)
        left_curve = (start * np.exp(-(shape-x_max)/(slope*data.size)) + base)[x_max-1::-1]
        extract = np.array(list(left_curve)+list(right_curve))
        return extract

class AbsorptionSpectrum(Spectrum):
    def __init__(self, array: NDArray[float]):
        super().__init__(array, 'A')
        self.blackbody = self.magic_numbers()
        self.flat_array: NDArray[float] = self.magic_numbers()-self.array

    def get_peaks(self, *args, **kwargs) -> tuple[list[float], list[float]]:

        cut_array = self.peak_cut(*args, **kwargs)
        seperated_peaks = self.separate_peaks(cut_array)
        self.send_help = seperated_peaks
        max_peaks = self.peak_maxima(seperated_peaks)
        self.peaks = max_peaks
        y = self.blackbody
        xcoords = []
        ycoords = []
        for val in max_peaks.values():
            xcoords.append(val[0])
            ycoords.append(-val[1]+y[val[0]])
        return xcoords, ycoords

    def blackbody_curve(self, wavelength: NDArray, temperature: float = 11000, base: float = 0):
        c = 3*10**8
        h = 6.26*10**-34
        kB = 1.38*10**-23
        c1, c2 = 2*h*c**2, h*c/kB
        wl = np.array(wavelength)
        curve = c1 / wl ** 5 * (1 / (np.exp(c2 / (wl * temperature)) - 1)) + base
        return curve

    def peak_cut_function(self, left_slope: float = 0.05, right_slope: float = 0.2) -> NDArray:
        data = self.blackbody - self.array
        x_max = data.argmax()
        right_shape = np.linspace(x_max, data.size, data.size-x_max)
        base = np.median(np.abs(data))
        right_curve = 0.65*data.max() * np.exp(-(right_shape-x_max)/(right_slope*data.size)) + base
        left_curve = (0.65*data.max() * np.exp(-(right_shape-x_max)/(left_slope*data.size)) + base)[x_max-1::-1]
        extract = np.array(list(left_curve)+list(right_curve))
        return extract