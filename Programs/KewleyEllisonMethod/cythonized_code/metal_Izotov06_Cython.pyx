'''Calculates metallicity based on method outlined in Izotov et al. (2006)'''


################################################################################
#
#                                IMPORT LIBRARIES
#
################################################################################


from scipy.optimize import bisect, brentq
from math import log10, log, exp
import numpy as np
import warnings

import sys
sys.path.insert(1, '/Users/kellydouglass/Documents/Drexel/Research/Programs/KewleyEllisonMethod/')
from Op_approximation import Op_ratio_approx # Op_approx
from ICF import N_ICF_I06, Ne_ICF_I06


################################################################################
#
#                        SOLVE ITERATIVELY FOR TEMPERATURE
#
################################################################################


def equations(p, Ne, ratio, bad_temp):
    t_OIII,Ct,x = p
    
    if t_OIII > 6.48:
        # temperature has reached / exceeded the asymptote
        t_OIII = 6.48
        bad_temp[0] = 1
        
    if t_OIII < 0.:
        # temperature has gone negative (for no reason)
        t_OIII = 0.1
        bad_temp[0] = 1
        
    A = 1e-4*Ne*t_OIII**-0.5 - x
    B = (8.44 - 1.09*t_OIII + 0.5*t_OIII**2 - 0.08*t_OIII**3)*((1. + 0.0004*x)/(1. + 0.044*x)) - Ct
    C = 1.432/(ratio - log10(Ct)) - t_OIII

    return (A, B, C)
    
    
cpdef double equation_nolog(double t_OIII, double Ne, double ratio) except *:
    
    cdef double A
    
    A = (10**(1.432/t_OIII))*(8.44 - 1.09*t_OIII + 0.5*t_OIII**2 - 0.08*t_OIII**3)*((1 + 0.0004e-4*Ne*t_OIII**-0.5)/(1 + 0.044e-4*Ne*t_OIII**-0.5)) - ratio
    
    return A
    
    
cpdef double equation(double t_OIII, double Ne, double ratio) except *:
        
    #print t_OIII
    cdef double A
        
    A = (1.432/t_OIII) + log10(8.44 - 1.09*t_OIII + 0.5*t_OIII**2 - 0.08*t_OIII**3) + log10(1 + 0.0004e-4*Ne*t_OIII**-0.5) - log10(1 + 0.044e-4*Ne*t_OIII**-0.5) - ratio
    
    return A
    

################################################################################
#
#                                CALCULATE METALLICITY
#
################################################################################


cpdef Izotov06(galaxy):
    '''galaxy is a row from an astropy table containing all required emission line 
    fluxes
    
    flag_3727 is boolean - 0 if 3727 cannot be used, 1 if it can
    flag_NH is boolean - 0 if N abundance cannot be calculated, 1 if it can
    flag_NeH is boolean - 0 if Ne abundance cannot be calculated, 1 if it can'''
    
    #print '---------', galaxy['index'],'---------'
    
    ############################################################################
    #    Define global variables
    ############################################################################
    
    cdef double ratio, ratio2, Ne, t_OIII, logOppHp12, OppHp, t_OII, a, b, x, c
    cdef double logOpHp12, OpHp, OH, Z, v, w, logNH, logNpHp12, logNeH, logNeppHp12
    
    
    ratio = (galaxy['OIII_4959_FLUX'] + galaxy['OIII_5007_FLUX'])/galaxy['OIII_4363_FLUX']
    ratio2 = (galaxy['OIII_4959_FLUX'] + galaxy['OIII_5007_FLUX'])/galaxy['H_BETA_FLUX']
    
    ############################################################################
    # Calculate electron number density in [OIII]
    ############################################################################
    
    Ne = 100.    # cm^-3

    ############################################################################    
    # Calculate electron temperature in [OIII]
    ############################################################################
    
    t_OIII = bisect(equation_nolog, 0.01, 7, args=(Ne, ratio))
    '''
    warnings.filterwarnings('error')
    
    try:
        #t_OIII,Ct,x = fsolve(equations, (1.,10.,0.01), args=(Ne, log10(ratio), bad_temp))
        #t_OIII = fsolve(equation_nolog, 1.)
        #t_OIII = root(equation, 0.01, args=(Ne, ratio, bad_temp), method='linearmixing')['x']
        #t_OIII = bisect(equation, 0.01, 6.60531, args=(Ne, log10(ratio)))
        t_OIII = bisect(equation_nolog, 0.01, 7, args=(Ne, ratio))
        #t_OIII = brentq(equation_nolog, 0.01, 7, args=(Ne, ratio))
    except ValueError as myError:
        print '---------', galaxy['index'],'---------'
        print myError
        print galaxy
        #t_OIII = float(input('t3 = '))
        raise
        exit()
    '''
    ############################################################################
    # Calculate metallicity of O++
    ############################################################################
    
    logOppHp12 = log10(ratio2) + 6.200 + (1.251/t_OIII) - 0.55*log10(t_OIII) - 0.014*t_OIII
    
    ############################################################################
    # Calculate oxygen abundance
    ############################################################################
    
    OppHp = 10**(logOppHp12 - 12)
    
    ########################################################################
    # Calculate electron temperature in [OII]
    ########################################################################
    
    '''if logOppHp12 <= 7.2:
        t_OII = -0.577 + t_OIII*(2.065 - 0.498*t_OIII)
    elif logOppHp12 <= 7.6:
        t_OII = -0.744 + t_OIII*(2.338 - 0.610*t_OIII)
    else:
        t_OII = 2.967 + t_OIII*(-4.797 + 2.827*t_OIII)'''
        
    t_OII = 0.7*t_OIII + 0.3
    
    

    ############################################################################
    #    Calculate ion abundance of O+
    ############################################################################

    if galaxy['flag3727'] == 1:
        
        ########################################################################
        # Calculate metallicity of O+
        ########################################################################
        
        a = log10(galaxy['OII_3727_FLUX']/galaxy['H_BETA_FLUX'])
        b = log10(t_OII)
        x = 1e-4*Ne*t_OII**-0.5
        c = log10(1. + 1.35*x)
        
        logOpHp12 = a + 5.961 + (1.676/t_OII) - 0.40*b - 0.034*t_OII + c
        logOpHp12_alt = np.NaN
        
        ########################################################################
        # Calculate total oxygen abundance
        ########################################################################
        
        OpHp = 10**(logOpHp12 - 12.)
        
        ########################################################################
        # Calculate metallicity
        ########################################################################
    
        OH = OppHp + OpHp
        Z = 12. + log10(OH)
        
    # [OII] 3727 is not available
    elif galaxy['flagOppSp'] == 1 and (galaxy['BPTclass'] in [1,3,4]):
    
        ########################################################################
        # Uses derived relationship between logOppHp12 and Z (with BPT separation and O++/S+ ratio)
        ########################################################################

        OppSp_ratio = galaxy['OIII_5007_FLUX']/(galaxy['SII_6717_FLUX'] + galaxy['SII_6731_FLUX'])
        
        #Z_arr,logOpHp12_arr,OpHp_arr = Op_approx(np.array([galaxy['BPTclass']]), np.array([galaxy['rabsmag']]), np.array([logOppHp12]), np.array([OppHp]))
        Z_arr,logOpHp12_arr,OpHp_arr = Op_ratio_approx(np.array([galaxy['BPTclass']]), np.array([galaxy['rabsmag']]), np.array([logOppHp12]), np.array([OppHp]), OppSp_ratio)
        Z = Z_arr[0]
        logOpHp12 = logOpHp12_arr[0]
        OpHp = OpHp_arr[0]
        
        '''
    else:
        # Assumes existence of OII_7325 flux
        
        a = log10(galaxy['OII_7325_FLUX']/galaxy['H_BETA_FLUX'])
        b = log10(t_OII)
        x = 1e-4*Ne*t_OII**-0.5
        c = log10(1. - 3.48*x)
        
        logOpHp12 = a + 6.901 + (2.487/t_OII) - 0.483*b - 0.013*t_OII + c
        logOpHp12_alt = np.NaN
        
        OpHp = 10**(logOpHp12 - 12.)
        
        OH = OppHp + OpHp
        Z = 12. + log10(OH)
        '''
    else:
        Z = np.NaN
        OpHp = np.NaN
        logOpHp12 = np.NaN
    
    
    ############################################################################
    # Calculate ion abundance of N+
    ############################################################################

    if (galaxy['flagNH'] == 1) and (np.isfinite(Z)):

        ########################################################################
        # Calculate N+ abundance
        ########################################################################
    
        a = log10((galaxy['NII_6548_FLUX'] + galaxy['NII_6584_FLUX'])/galaxy['H_BETA_FLUX'])
        b = log10(t_OII)
        x = 1e-4*Ne*t_OII**-0.5
        c = log10(1. + 0.116*x)
    
        logNpHp12 = a + 6.234 + (0.950/t_OII) - 0.42*b - 0.027*t_OII + c
        NpHp = 10**(logNpHp12 - 12.)
    
        ########################################################################
        # ICF
        ########################################################################
    
        v = OpHp/(OpHp + OppHp)
    
        ICF = N_ICF_I06(np.array([Z]),np.array([v]))
        
        ####################################################################
        # Calculate N abundance
        ####################################################################
    
        NH = ICF[0]*NpHp
        logNH = 12. + np.log10(NH)
        '''
        try:
            logNH = 12. + np.log10(NH)
        except RuntimeWarning as myWarning:
            print '---------', galaxy['index'],'---------'
            print myWarning
            print 'N+/H+:', NpHp
            print 'O++/H+:', OppHp
            print 'O+/H+:', OpHp
            print 'v:', v
            print '12 + log(N+/H+):', logNpHp12
            print 'ICF:', ICF[0]
            print 'N/H: ', NH
            raise
        '''
    else:
    
        logNH = np.NaN
        logNpHp12 = np.NaN
        
        
    ############################################################################
    # Calculate ion abundance of Ne++
    ############################################################################

    if (galaxy['flagNeH'] == 1) and (np.isfinite(Z)):

        ########################################################################
        # Calculate Ne++ abundance
        ########################################################################
    
        a = log10(galaxy['NeIII_3869_FLUX']/galaxy['H_BETA_FLUX'])
        b = log10(t_OII)
    
        logNeppHp12 = a + 6.444 + (1.606/t_OIII) - 0.42*b - 0.009*t_OIII
        NeppHp = 10**(logNeppHp12 - 12.)
    
        ########################################################################
        # ICF
        ########################################################################
    
        w = OppHp/(OpHp + OppHp)
    
        ICF = Ne_ICF_I06(np.array([Z]),np.array([w]))
        
        ####################################################################
        # Calculate N abundance
        ####################################################################
    
        NeH = ICF[0]*NeppHp
        logNeH = 12. + np.log10(NeH)
        '''
        try:
            logNH = 12. + np.log10(NH)
        except RuntimeWarning as myWarning:
            print '---------', galaxy['index'],'---------'
            print myWarning
            print 'N+/H+:', NpHp
            print 'O++/H+:', OppHp
            print 'O+/H+:', OpHp
            print 'v:', v
            print '12 + log(N+/H+):', logNpHp12
            print 'ICF:', ICF[0]
            print 'N/H: ', NH
            raise
        '''
    else:
    
        logNeH = np.NaN
        logNeppHp12 = np.NaN

    
    '''
    #print 'log(OIII/4363) =', ratio
    print '10^-4 Te =', t_OIII
    #print 'log(OII/Hbeta) =', a
    #print 'log(OIII/Hbeta) =', math.log10(ratio2)
    #print 'log(OIII/Hbeta) =', math.log(ratio2)
    print 'log(O+/H+) =', logOpHp12
    print 'log(O++/H+) =', logOppHp12
    print '12+log(O/H) =', Z
    print '12+log(N/H) =', logNH
    print 'log(N/O) =', logNH - Z
    '''
    
    ############################################################################
    #    Return abundance estimates
    ############################################################################

    return Z,logOppHp12,logOpHp12,t_OIII,logNH,logNpHp12,logNeH,logNeppHp12
