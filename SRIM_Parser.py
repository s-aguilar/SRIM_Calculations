import re
from scipy.interpolate import interp1d
import numpy as np

class SRIM_Table(object):
    def __init__(self,fname):

        # TO DO: 1. Incorporate some way of determining the stopping units used (EASY/ INCOMLPETE)
        #        2. convert lists to numpy arrays (NOT STARTED)

        sci = '(\d+\.\d+E[+-]*\d+)' # check sign of power
        flo = '(\d+\.\d+|\d+)'      # '|' is OR. Check whether it is a float or an int
        unit = '(\w+)'
        fmt = ('\s*'+flo+'\s'+unit+'\s+'+sci+'\s+'+sci
                +'\s+'+flo+'\s+'+unit
                +'\s+'+flo+'\s+'+unit
                +'\s+'+flo+'\s+'+unit
                +'.*')

        data = []

        # eV / A
        units_eVa = {'eV':1, 'keV':1e3, 'MeV':1e6, \
                        'm':1e10, 'mm':1e7, 'um':1e4, 'A':1}

        # keV / um
        units_keVum = {'eV':1e-3, 'keV':1, 'MeV':1e3, \
                        'm':1e-6, 'mm':1e-3, 'um':1, 'A':1e-4}

        # MeV / mm
        units_MeVmm = {'eV':1e-6, 'keV':1e-3, 'MeV':1, \
                        'm':1e3, 'mm':1, 'um':1e-3, 'A':1e-7}

        # keV / (ug/cm2)
        units_keVugcm2 = {'eV':1e-3, 'keV':1, 'MeV':1e3}

        # MeV / (mg/cm2)
        units_MeVmgcm2 = {'eV':1e-6, 'keV':1e-3, 'MeV':1}

        # keV / (mg/cm2)
        units_keVmgcm2 = {'eV':1e-3, 'keV':1, 'MeV':1e3}

        # eV / (1E15 atoms/cm^2)
        units_eVatomscm2 = {'eV':1, 'keV':1e3, 'MeV':1e6}


        # dictionary of units
        dictUnits = {'eV / A':units_eVa}

        with open(fname) as f:
            for line in f:
                res = re.match(fmt,line)
                stopUnit = re.search(' Stopping Units = ',line)
                # if stopUnit: print(line[19:-1])
                if stopUnit: units = str(line[19:-1])
                if res is None: continue

                # The '-u' suffix denotes unit and L is length/range
                E,Eu, dE_dX1, dE_dX2, L, Lu = res.groups()[:6]

                # print(E,Eu,dE_dX1, dE_dX2, L, Lu)
                Total_dEdX = float(dE_dX1) + float(dE_dX2)
                data.append({'E':float(E)*units_keVum[Eu],'Range':float(L)*units_keVum[Lu],\
                'dXdE':1./float(Total_dEdX),'dEdX':float(Total_dEdX)})

        if not data:
            raise Exception("Input a valid table")

        self.MaxE0 = max(_['E'] for _ in data)
        self.MaxR = max(_['Range'] for _ in data)
        self.dXdE = interp1d(*zip(*[(_['E'],_['dXdE'],) for _ in data]),
                fill_value='extrapolate' )
        self.dEdX = interp1d(*zip(*[(_['E'],_['dEdX'],) for _ in data]),
                fill_value='extrapolate' )
        self.R2E = interp1d(*zip(*[(_['Range'],_['E'],) for _ in data]),
                fill_value='extrapolate' )
        self.E2R = interp1d(*zip(*[(_['E'],_['Range'],) for _ in data]),
                fill_value='extrapolate' )


    def GetE(self,E0,X):
        if np.logical_or.reduce(np.array(E0)>self.MaxE0):
            raise Exception("exceed the maximum energy range of the table")

        R = self.E2R(E0)
        R_ = R-X
        E = self.R2E(R_)
        return np.clip(E,0,float('inf'))

    def GetE0(self,x):
        if np.logical_or.reduce(np.array(x)>self.MaxR):
            raise Exception("exceed the maximum distance range of the table")

        return self.R2E(x)
