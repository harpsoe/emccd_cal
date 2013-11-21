'''
Created on Nov 21, 2013

@author: Kennet Harpsoee

This project contains function to callibrate the parameters of an EMCCD gain curve .

'''


from scipy.stats import scoreatpercentile, expon
from scipy.interpolate import interp1d
import numpy as np

def gaussian(x,mu,sigma):
    '''Gaussian PDF'''
    return (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(x-mu)**2/(2*sigma**2))

def expon2(x,gamma):
    '''Exponential PDF'''
    return expon.pdf(x,loc=0,scale=gamma)

def cal_line(x,mu,sigma,gamma,p):
    '''Calibrated a line according to the expectation maximisation
       procedure described in [KH2012]_

       .. [KH2012] High Frame-rate Imaging Based Photometry
          http://arxiv.org/abs/1202.3814'''
    def M_step(x,mu,sigma,gamma,p):
        return p*gaussian(x,mu,sigma)/((1-p)*expon2(x,gamma)+p*gaussian(x,mu,sigma))
    def E_step(x,lamb):
        return np.mean(lamb)
    lamb = M_step(x,mu,sigma,gamma,p)
    k = 0
    chi_old = np.sum((1-lamb)*expon2(x,gamma)+lamb*gaussian(x,mu,sigma))
    for i in range(50):
        lamb = M_step(x,mu,sigma,gamma,p)
        p = E_step(x,lamb)
        #mu = np.sum(lamb*(x-mu)) / np.sum(lamb)
        sigma = np.sqrt(np.sum(lamb*(x-mu)**2) / np.sum(lamb))
        gamma = np.sum(1-lamb)/np.sum((1-lamb)*(x-mu))
        chi_temp = np.sum((1-lamb)*expon2(x,gamma)+lamb*gaussian(x,mu,sigma))
        delta_chi = abs(chi_old-chi_temp)
        chi_old = chi_temp
        k+=1
        print p, mu, sigma, delta_chi
    return mu,sigma,gamma,1-p
    
def callibration(dataCube):
    '''Calibrates a data cube of EMCCD bias frames''' 
    shape = dataCube.shape
    biasArr = np.zeros((shape[1],shape[2]))
    ampArr = np.zeros((shape[1],shape[2]))
    ronArr = np.zeros((shape[1],shape[2]))
    spArr = np.zeros((shape[1],shape[2]))
    for i in range(shape[1]):
        for j in range(shape[2]):
            l = dataCube[:,i,j]
            mu,sigma,gamma,p = cal_line(l,0.1,4.2,21,0.6)
            biasArr[i,j] = mu
            ronArr[i,j] = sigma
            ampArr[i,j] = gamma
            spArr[i,j] = p
        print 'Line',i
    return biasArr,ronArr,ampArr,spArr

def tmean(a,axis=None):
    '''Returns the truncated mean from above of an array a, should not be used'''
    percentile = scoreatpercentile(a.reshape(-1),85)
    masked = np.where(a<percentile,a,0)
    number = len(np.flatnonzero(masked))
    return np.sum(masked) / float(number)

def neighbors(arr,x,y,n=3):
    ''' Given a 2D-array, returns an n x n array whose "center" element is arr[x,y]'''
    arr=np.roll(np.roll(arr,shift=-x+1,axis=0),shift=-y+1,axis=1)
    return arr[:n,:n]
