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
    if dE < -0.000001:
        printlog( "dE is ", dE, "Accept!")
        decrease = True
    elif  1 > math.exp(-dE/kB/T) > random():
        printlog ("Accepted due to the temperature; exponent is ", math.exp(-dE/kB/T) )
        decrease = True
    else:
        printlog('Not accepted')

    return decrease


