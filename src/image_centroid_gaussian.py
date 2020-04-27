import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import sys
#from ds9 import ds9

def twoD_Gaussian(x, y, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """Returns a two-d gaussian function
    coordinate variables:
    x - the x coordinate.  may be an array
    y - the y coordinate.  may be an array
    
    function variables
    amplitude - the peak amplitude of the gaussian
    xo - the x centroid of the gaussian
    yo - the y centroid of the gaussian
    sigma_x - the std dev of the gaussian in x
    sigma_y - the std dev of the gaussian in y
    theta - the rotation angle (in radians) - Note--degenerate with sigmaxy when off by 180
    offset - the background level

    Returns
    a scalar or vector or array of z-values at each x, y
    """
    #x, y = xytuple
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( -(a*((x-xo)**2) 
                                    + 2*b*(x-xo)*(y-yo) 
                                    + c*((y-yo)**2)))
    return g

def twoD_Gaussian_fitfunc(xytuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    #passes a single parameter for x, y
    #ravels hte output
    return twoD_Gaussian(xytuple[0], xytuple[1], amplitude, xo, yo, sigma_x, sigma_y, theta, offset).ravel()

def fitgaussian(image, x=None, y=None):
    """Returns a fit for a gaussian function given an image input"""
    if (x is None) or (y is None):
        x, y = np.meshgrid( np.arange(image.shape[1]),
                            np.arange(image.shape[0]))
    xindguess= np.argwhere(np.nansum(image, axis=0)==np.max(np.nansum(image, axis=0)))[0][0]
    yindguess= np.argwhere(np.nansum(image, axis=1)==np.max(np.nansum(image, axis=1)))[0][0] 
    
    xoguess = x[yindguess, xindguess]
    yoguess = y[yindguess, xindguess]
    offsetguess = np.nanmean(image)
    amplitudeguess = np.nanmax(image)-offsetguess
    initguess = (amplitudeguess,
                 xoguess,
                 yoguess,
                 1.0,
                 1.0,
                 0.0,
                 offsetguess)
    inds = np.isfinite(image)
    input_image = image[inds]
    input_x = x[inds]
    input_y = y[inds]
    popt, pcov = opt.curve_fit(twoD_Gaussian_fitfunc, (input_x, input_y), 
                   input_image.ravel(),
                   p0 = initguess, maxfev = 100000000)
    return popt

def image_centroid_gaussian(image, x=None, y=None):
    """Returns the coordinates of the center of a gaussian blob in the image"""
    popt = fitgaussian(image, x=x, y=y)
    return popt[1], popt[2]

if __name__ == "__main__":
    print("Fitting a noisy image with known x and y coordinates")
    xcoords = np.arange(30, 60)
    ycoords = np.arange(20, 80)
    x, y = np.meshgrid(xcoords, ycoords)
    print("labels".ljust(15), [x.ljust(10) for x in "amp, x0, y0, sigx, sigy, th, off".split(',')])
    true_params = 100, 55, 70, 2.0, 1.0, np.pi/180*60, 100
    print("True params: ".ljust(15), [('%.2f'%x).ljust(10) for x in true_params])
    testimage = twoD_Gaussian(x, y, *true_params)
    noisy_testimage = np.random.poisson(np.abs(testimage))*1.0
    noisy_testimage[np.random.randint(np.shape(noisy_testimage)[0], size=10),
                     np.random.randint(np.shape(noisy_testimage)[1], size=10)]=np.nan 
    popt = fitgaussian(noisy_testimage, x=x, y=y)
    print("Fit params: ".ljust(15), [('%.2f'%x).ljust(10) for x in popt])
    print("\n")
    
    print("Fitting the same index with no x and y coordinates")
    popt2 = fitgaussian(noisy_testimage)
    print("Fit params: ".ljust(15), [('%.2f'%x).ljust(10) for x in popt2])
    fit = twoD_Gaussian(x,y, *popt)
    
    fig, ax = plt.subplots(1,3)
    imargs = {'origin':'lower', 'cmap':'Greys',
              'vmin':np.min(testimage)/2, 'vmax':np.max(testimage)}
    titles = ['Original', 'Noisy', 'Fit to Noisy']
    data = [testimage, noisy_testimage, fit]
    for i in range(0, 3):
        ax[i].imshow(data[i], **imargs)
        ax[i].set_title(titles[i])
        ax[i].set_xticks(range(len(x[0,:])))
        ax[i].set_yticks(range(len(y[:,0])))
        ax[i].set_xticklabels(x[0,:])
        ax[i].set_yticklabels(y[:,0])
        ax[i].set_xlabel('x')
        ax[i].set_ylabel('y')
    plt.show()
    #ds9([testimage, noisy_testimage,fit])
