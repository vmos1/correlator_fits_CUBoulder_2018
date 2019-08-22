
import lsqfit
import gvar as gv
import matplotlib.pyplot as plt
import numpy as np


# Functions for performing Fit and plot.
   
def f_make_p0(par):
    p0=par.copy()
    return p0

 

def f_fit_plot(x,y,all_dta,fit,error_band=False,semilog=False,full_data=True):
    '''
    Function for plotting data with the fit lines and error bands.
    For correlators, using a semi-log plot.
    full_data=True, plots the entire data and the best-fit in the fit region
    '''

    plt.figure()
    # Plots data points
    if full_data: # Plot all the correlators even those not used in the fit.
        plt.errorbar(x=all_dta['x'],y=gv.mean(all_dta['y']),yerr=gv.sdev(all_dta['y']),color='black',linestyle='None',marker='o')
    else:         # Plot the data points used in the fit
        plt.errorbar(x=x['x'],y=gv.mean(y['y']),yerr=gv.sdev(y['y']),linestyle='None',color='red',marker='s')
    
    # Plot the best fit line
    # # Using a finer grid to get a continuous curve.
    curvex=dict(x)
    curvex['x']=np.linspace(min(fit.x['x']),max(fit.x['x']),500)
    curvey=gv.mean(fit.fcn(curvex,fit.p))
    plt.plot(curvex['x'],curvey,color='blue')

    if error_band:
        obs_fit=gv.mean(fit.fcn(curvex,fit.p))
        err_fit=gv.sdev(fit.fcn(curvex,fit.p))  
        sigma=2.0
        plt.fill_between(curvex['x'],obs_fit-sigma*err_fit,obs_fit+sigma*err_fit,color='yellow') # providing an error band.

    if semilog: plt.semilogy()
   
    plt.title("Plot")

 

def f_perform_fit(all_dta,x,y,par,f_func,verbose=True,concise=False,plot=True,full_data=False,error_band=False,semilog=False):
    
    '''
    Function wrapper to perform correlator fit for mesons.
    Reads dictionary dta with x and y data.
    Performs both cosh and sinh fits.
    '''
    
    # Performing meson fit
    
    p0=f_make_p0(par)
    fit = lsqfit.nonlinear_fit(data=(x, y['y']), fcn=f_func, p0=p0,extend=True,svdcut=1e-8)
    
    # Print the fit results
    
    if concise : 
        pass
    else:
        print "*****************************",'\n'
        if verbose:
            print f_func.__doc__.strip('\n').strip('  ') # Prints the functional form of fit function.
            print(fit.format(maxline=True)),"\n\n"

        else :
            print fit,"\n\n"
    
    # Plot fit
    if plot: f_fit_plot(x,y,all_dta,fit,full_data=full_data,error_band=error_band,semilog=semilog)

    return fit

# Code to find outlier
def f_outlier_list(fit):
    ''' Code to find outliers in a fit
    Reads fit and returns list of indices of outliers.
    '''
    y_th=fit.fcn(fit.x,fit.p)
    y_exp=fit.y
    outlier_list=[]
    for i,(j,k) in enumerate(zip(y_exp,y_th)):
        # Find points with chi-square more than 3.
        if gv.chi2(j,k)>3.0:
            outlier_list.append(i)

    return outlier_list









