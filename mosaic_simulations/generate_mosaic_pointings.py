import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import astropy.units as u
from astropy.io import fits
from matplotlib.colors import SymLogNorm

def generate_rectangle(x_fov,y_fov,RA_centre,Dec_centre,theta):
    fov = np.array([x_fov/2.,y_fov/2.])
    xy = np.array([fov,fov*np.array([-1,1]),fov*-1,fov*np.array([1,-1])]).T
    return translate(rotate(xy,theta*(np.pi/180.)),[RA_centre,Dec_centre]).T

def mosaic_pointings_square(centre_ra, centre_dec, centre_freq, ra_fov, dec_fov, theta, pointing_file_path, pb_fwhm):

    # Initialise some things...

    ra_mosaic_buffer = []   # array of RA coordinates for plotting   
    dec_mosaic_buffer = []  # array of Dec coordinates for plotting
    pointing_counter = 0    # the total number of pointings

    margin = 0.0           # determines edge-of-field margin, consistent with simdata

    # Everything in degrees:

    if ra_fov < 0.0:
        ra_fov = ra_fov * -1.0

    coords = generate_rectangle(ra_fov,dec_fov,centre_ra,centre_dec,theta)
    
    ra_fov = 1.3*np.abs((np.max(coords.T[0])-np.min(coords.T[0])))
    dec_fov = 1.3*np.abs((np.max(coords.T[1])-np.min(coords.T[1])))

    # Airy disk of a 12-m dish in degrees

    #pb_fwhm = 1.2 * 3.0e8/1.4e11/12.0*180.0/np.pi
    half_pb  = pb_fwhm * 0.6

    ra_spacing = half_pb 
    dec_spacing = half_pb * 0.866025404  # cos 60 for hex

    n_rows = 1 + int(np.floor((dec_fov / dec_spacing) - 2.0 * margin / 0.866025404))

    float_cols = 1 + (ra_fov / ra_spacing) - (2.0 * margin)
    n_cols = int(np.floor(float_cols))

    if float_cols - n_cols >= 0.5 and n_rows > 1:
        even_cols = n_cols
        n_cols_min = 0.5 * (n_cols - 0.5)
    else:
        even_cols = n_cols -1
        n_cols_min = 0.5 * (n_cols - 1)

    current_dec = centre_dec + (0.5 * (n_rows-1) * dec_spacing)

    for i in range(0,n_rows):
        ra_spacing = half_pb / np.cos(current_dec*(np.pi/180.))

        if i % 2:
            ra_min = (centre_ra - (n_cols_min * ra_spacing))
            stop_col = n_cols
        else:
            ra_min = (centre_ra - ((n_cols_min - 0.5) * ra_spacing))
            stop_col = even_cols

        for j in range(0,stop_col):
            current_ra = ra_min + j * ra_spacing
            ra_mosaic_buffer.append(current_ra)
            dec_mosaic_buffer.append(current_dec)
            pointing_counter += 1
        current_dec = current_dec - dec_spacing
        
    temp_ra = []
    temp_dec =[]
    for i in range(len(ra_mosaic_buffer)):
        truth = determine_in_out(coords,[ra_mosaic_buffer[i],dec_mosaic_buffer[i]])
        if truth == True:
            temp_ra.append(ra_mosaic_buffer[i])
            temp_dec.append(dec_mosaic_buffer[i])
            
    ra_mosaic_buffer = temp_ra
    dec_mosaic_buffer = temp_dec
    # Write out a pointings file and generate a list of beams
    os.system('rm %s'%pointing_file_path)
    ptgfile = open(pointing_file_path,'w')

    print('#Epoch     RA          DEC      RANGE',file=ptgfile)

    if ra_mosaic_buffer:
        for index in range(0, len(ra_mosaic_buffer)):
            tmp_ra = str(ra_mosaic_buffer[index])
            tmp_dec = str(dec_mosaic_buffer[index])
            ptgstring = 'J2000 '+str(tmp_ra)+' '+str(tmp_dec)+' '+str(pb_fwhm)
            print(ptgstring, file=ptgfile)
    else:
        ptgstring = 'J2000 '+str(centre_ra)+' '+str(centre_dec)+' '+str(pb_fwhm)
        print(ptgstring,file=ptgfile)
    ptgfile.close()
    return coords
    #return ra_mosaic_buffer, dec_mosaic_buffer, half_pb, pointing_counter, coords

def hyp(co_a,co_b):
    return np.sqrt(np.abs(co_a[0]-co_b[0])**2. + np.abs(co_a[1]-co_b[1])**2)

def tri_area(co_a,co_b,co_c):
    ##Get sizes of each side - assume co_d=[x_d,y_d]
    a = hyp(co_a,co_b)
    b = hyp(co_b,co_c)
    c = hyp(co_c,co_a)
    s = (a+b+c)/2
    return (s*(s-a)*(s-b)*(s-c)) ** 0.5

def rect_area(co_a,co_b,co_c):
    a = hyp(co_a,co_b)
    b = hyp(co_b,co_c)
    return a*b

def rotate(xy, theta):
    # https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
    theta=theta*-1
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)

    return np.array([
        xy[0] * cos_theta - xy[1] * sin_theta,
        xy[0] * sin_theta + xy[1] * cos_theta])


def translate(xy, offset):
    return np.array([xy[0] + offset[0], xy[1] + offset[1]])


def determine_in_out(coords,point):
    triangle_areas = []
    for i in range(len(coords)):
        if i == (len(coords)-1):
            triangle_areas.append(tri_area(coords[i],coords[0],point))
        else:
            triangle_areas.append(tri_area(coords[i],coords[i+1],point))
    tri_sum = np.sum(triangle_areas)
    rect_sum = rect_area(coords[0],coords[1],coords[2])
    return np.isclose(tri_sum,rect_sum)