import os, subprocess
import math
import numpy as np
import copy
import datetime
import shutil
import traceback
import glob
from random import randint, random

from siman.header import printlog, kB






def metropolis(E1, E2, T = 1):
    """
    Metropolis algorithm
    """
    decrease = False # energy reduction
    
    # kB = 1.3806488*10**-23  / 1.6 * 10 **19
    dE = E2 - E1
    
    printlog("metropolis(): dE is ", dE)
    r = random()
    e = math.exp(-dE/kB/T)
    printlog("random number is ", r)
    printlog("exponent is  ", e)
    if dE < -0.000001:
        printlog( "dE is ", dE, "Accept!")
        decrease = True
    elif  1 > e > r:
        printlog ("Accepted due to the temperature; exponent is ", e )
        decrease = True
    else:
        printlog('Not accepted')

    return decrease


