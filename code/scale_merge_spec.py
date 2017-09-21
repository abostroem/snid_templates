import os
import numpy as np
from astropy.io import fits, ascii

from matplotlib import pyplot as plt

import spectroscopy

class fits_spec():
    def __init__(self, directory, filename, ext=0):
        self.filename=filename
        if os.path.splitext(self.filename)[1]!='.fits':
            raise IOError('Ascii file passed to ascii class')
        spec_file = fits.open(os.path.join(directory, self.filename))
        self.flux = spec_file[ext].data
        if len(self.flux.shape) > 1:
            self.flux = self.flux[0,0,:]
        self.hdr = spec_file[ext].header
        pix = np.arange(1, len(self.flux)+1)
        self.wavelength = spectroscopy.calc_wavelength(self.hdr, pix)
        
        
class ascii_spec():
    def __init__(self, directory, filename):
        self.filename = filename
        if os.path.splitext(self.filename)[1]=='.fits':
            raise IOError('Fits file passed to ascii class')
        spec_table = ascii.read(os.path.join(directory, self.filename), data_start=0)
        if 'col0' in spec_table.colnames:
            self.flux = spec_table['col1']
            self.wavelength = spec_table['col0']
        elif 'col2' in spec_table.colnames:
            self.flux = spec_table['col2']
            self.wavelength = spec_table['col1']
        else:
            print('WARNING, unsure how to read {}'.format(os.path.join(directory,filename)))
            sys.exit()


def read_spec(directory, filename):
    if os.path.splitext(filename)[1] == '.fits':
        if 'wifes' in filename.lower():
            spec = fits_spec(directory, filename, ext=1)
        else:
            spec = fits_spec(directory, filename)
    else:
        spec = ascii_spec(directory, filename)
    return spec

# Trim wavelength range
def trim_wavelengths(spec):
    '''
    Remove the first and last few points, these are most often bad
    Also remove any points at wavelengths greater than 10000
    '''
    spec.wavelength = spec.wavelength[10:-10]
    spec.flux = spec.flux[10:-10]
    if (spec.wavelength>10000).any():
        indx = spec.wavelength<10000
        spec.wavelength = spec.wavelength[indx]
        spec.flux = spec.flux[indx]
    return spec


# Find overlap range
def find_overlap_values(spec1, spec2):
    '''
    Find the overlaping end points between 2 wavelength arrays
    '''
    wl1_min = np.min(spec1.wavelength)
    wl1_max = np.max(spec1.wavelength)
    wl2_min = np.min(spec2.wavelength)
    wl2_max = np.max(spec2.wavelength)
    
    #No overlap
    if (wl1_min > wl2_max) or (wl2_min > wl1_max):
        overlap=None
    #Overlap
    else:
        overlap = [max(wl1_min, wl2_min), min(wl1_max, wl2_max)]
    return overlap


def get_overlap_data(spec, overlap):
    '''
    Given a range of overlapping wavelengths, determine 
    the flux and wavelength in the overlap region
    '''
    overlap_index = [(spec.wavelength <= overlap[1])&
                     (spec.wavelength>=overlap[0])]
    overlap_wave = spec.wavelength[overlap_index]
    overlap_flux = spec.flux[overlap_index]
    return overlap_wave, overlap_flux


# Find scale factor
def integrate_flux(wavelength, flux):
    '''
    Integrate the flux of a spectrum
    '''
    spacing = wavelength[1:] - wavelength[0:-1]
    spacing = np.append(spacing, spacing[-1])
    integral = np.sum(flux*spacing)
    return integral
    
def find_scalefactor(wavelength1, flux1, wavelength2, flux2):
    '''
    Find the scale factor needed to scale the smaller flux
    to the larger flux and the scale array - which array s
    should be scaled
    '''
    sum1 = integrate_flux(wavelength1, flux1)
    sum2 = integrate_flux(wavelength2, flux2)
    #Scale everything to larger flux value (even net counts)
    if (sum1 > sum2) or ((flux1 < 10**-10).all() and (flux2 > 10**-10).all()):
        scalefactor = sum1/sum2
        scale_array = 2
    else:
        scalefactor=sum2/sum1
        scale_array = 1
    return scalefactor, scale_array


# Scale data
def scale_data(spec, scalefactor):
    '''
    Scale flux based on scale factor
    '''
    spec.flux = spec.flux*scalefactor
    return spec


# Put it all together and visually test it works
def scale_spectrum(spec1, spec2, spectrum_to_scale= None):
    '''
    Scale one spectrum to another. The default is to scale the lower
    value to the higher value, however, spectrum to scale allows the user
    to set which spectrum is scaled
    '''
    spec1 = trim_wavelengths(spec1)
    spec2 = trim_wavelengths(spec2)
    #Scale only when trimming doesn't cut whole wavelength range (e.g. IR data)
    if (len(spec1.wavelength)>0) & (len(spec2.wavelength)>0):
        overlap = find_overlap_values(spec1, spec2)
        if overlap is not None:
            wl_1, flux_1 = get_overlap_data(spec1, overlap)
            wl_2, flux_2 = get_overlap_data(spec2, overlap)

            scale_factor, scale_array = find_scalefactor(wl_1, flux_1, wl_2, flux_2)

            if spectrum_to_scale is None:
                if scale_array == 2:
                    spec2.flux = spec2.flux*scale_factor
                elif scale_array == 1:
                    spec1.flux = spec1.flux*scale_factor
            elif spectrum_to_scale == 2:
                if scale_array == 2:
                    spec2.flux = spec2.flux*scale_factor
                elif scale_array == 1:
                    spec2.flux = spec2.flux/scale_factor
            elif spectrum_to_scale == 1:
                if scale_array == 2:
                    spec1.flux = spec1.flux/scale_factor
                elif scale_array == 1:
                    spec1.flux = spec1.flux*scale_factor
    return spec1, spec2

def median_combine_data(spec_list):
    wl_range = []
    disp = []
    for spec in spec_list:
        wl_range.append(np.min(spec.wavelength))
        wl_range.append(np.max(spec.wavelength))
        disp.append(np.median(spec.wavelength[1:]-spec.wavelength[0:-1]))
    min_wl = min(wl_range)
    max_wl = max(wl_range)
    disp = np.array(disp)
    avg_disp = np.mean(disp)
    if ((disp - avg_disp) > .1*avg_disp).any():
        print('WARNING, Disperions differ by more than 10\% from average {}'.format(avg_disp))
        print(disp)
        print(spec_list[0].filename)
    wl_array = np.arange(min_wl, max_wl, avg_disp)
    stack_spec = np.interp(wl_array, spec_list[0].wavelength, spec_list[0].flux, left=np.nan, right=np.nan)
    for spec in spec_list[1:]:
        stack_spec = np.vstack((stack_spec, np.interp(wl_array, spec.wavelength, spec.flux, left=np.nan, right=np.nan)))
    median_spec = np.nanmedian(stack_spec, axis=0)
    return wl_array, median_spec



def test_basic_functionality():
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    spec1 = ascii_spec('../2013ej', '2013ej_20130811_wifes_BR.dat')
    spec2 = fits_spec('../2013ej', '2013ej_20130811_fts.fits')
    ax1.plot(spec1.wavelength, spec1.flux, label='spec1-before')
    ax1.plot(spec2.wavelength, spec2.flux, label='spec2-before')

    scale_spectrum(spec1, spec2)

    ax2.plot(spec1.wavelength, spec1.flux, label='spec1-after')
    ax2.plot(spec2.wavelength, spec2.flux, label='spec2-after')
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax1.set_title('Before')
    ax2.set_title('After')


# Test net to flux scaling works
def test_flux_net_scaling():
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    spec1 = fits_spec('../2013ej', '2013ej_20130812_fts.fits')
    spec2 = ascii_spec('../2013ej', 'sn2013ej-20130812.508-ui.flm')
    ax1.plot(spec1.wavelength, spec1.flux, label='spec1-before')
    ax1.plot(spec2.wavelength, spec2.flux, label='spec2-before')

    scale_spectrum(spec1, spec2)

    ax2.plot(spec1.wavelength, spec1.flux, label='spec1-after')
    ax2.plot(spec2.wavelength, spec2.flux, label='spec2-after')
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    ax1.set_title('Before')
    ax2.set_title('After')





