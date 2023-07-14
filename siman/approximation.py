#!/usr/bin/env python3
""" 
Author: Kartamyshev A.I. (Darth Feiwante)
"""

# Fitting of the E(a,c) dependence for the equilibrium c/a searching

import xalglib

class ALGLIB:

    def build_2d_bicubic_spline(self, x, m, y, n, z, d): self.bicubicv2d = xalglib.spline2dbuildbicubicv(x, m, y, n, z, d)

    def calc(self, x, y, ind): 
        l = xalglib.spline2ddiff(self.bicubicv2d,x,y)
        if ind==0: return l[0]    # z
        elif ind==1: return l[1]  # dz/dx
        elif ind==2: return l[2]  # dz/dy
        elif ind==3: return l[3]  # d2z/dxdy
        else: raise RuntimeError ('Unknown ind = '+str(ind))


class Approximation:

    def aprx_lsq(self, fun_aprx, num_coef, xx, yy):

        from scipy.optimize import leastsq

        coefs = ''
        for i in range(num_coef): coefs += 'a'+str(i)+', '
        coefs+=' = par'
        
        def f_aprx(par, yy1, xx1):
            exec(coefs)
            x1 = []
            for x in xx1: x1.append(eval(fun_aprx))
            return [yy1[i] - x1[i] for i in range(len(x1))]
            
        plsq = leastsq(f_aprx, [1 for i in range(num_coef)], args=(yy, xx))
        self.coefs = list(plsq[0])