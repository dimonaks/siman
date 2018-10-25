#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import 
import os, math, re, sys
import numpy as np
from collections import Iterable
import shutil, gzip
import traceback
from contextlib import contextmanager

try:
    from six import string_types
except:
    print('Warning! six module was not found, I use only str and basestring as string_types; errors can be expected')
    try:
        string_types = basestring #for python 2.7
    except NameError:
        string_types = str # for python 3


from siman.header import printlog




class TracePrints(object):
  def __init__(self):    
    self.stdout = sys.stdout
  def write(self, s):
    self.stdout.write("Writing %r\n" % s)
    traceback.print_stack(file=self.stdout)


def angle(v1, v2):
  return math.acos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))) / math.pi * 180



def red_prec(value, precision = 100.):
    a = value * precision
    return round(a)/1./precision



def is_list_like(obj): 
    
	return not isinstance(obj, string_types) and isinstance(obj, Iterable)

def is_string_like(s):
    return isinstance(s, string_types)

def list2string(ilist):
    #' '.join(['{:}']*len(lis)).format(*lis)
    return ' '.join(np.array(ilist).astype(str))


def merge_dics(dic1, dic2):
	"""
	return dic
	"""
	dic_new = dic1.copy()
	dic_new.update(dic2)
	return dic_new


def cat_files(files, output_file):
    #files (list) - file names
    #

    with open(output_file,'wb') as wfd:
        for f in files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd, 1024*1024*10)
    return


def grep_file(string, file, reverse = False):
    #reverse(bool) - read file in reversed order - require upload to memory
    out = ''
    with open(file, 'r') as f:
        
        if reverse:
            f = reversed(f.readlines())

        for line in f:
            if string in line:
                out = line
                

    return str(out.strip() )

def gunzip_file(filename):
    printlog('unzipping file', filename)
    with open(filename.replace('.gz', ''), 'wb') as f_out:
        with gzip.open(filename, 'rb') as f_in:
            shutil.copyfileobj(f_in, f_out)

    return



def makedir(path):
    """
    *path* - path to some file 
    Make dirname(path) directory if it does not exist
    """
    dirname = os.path.dirname(path)

    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)
        printlog("Directory", dirname, " was created", imp = 'y')
    return


def latex_chem(formula):
    """ """
    # print(formula)
    if '$' not in formula:
        formula =re.sub("([0-9]{1,3})", "$_{\\1}$", formula)
    return formula

def latex_spg(spg):

    # print(spg)
    # spg = spg.replace('_1', '$_1$')
    spg = spg.replace('p', 'P')
    if '-' in spg:
        pos = spg.find('-')
        dig = spg[pos+1]
        spg = spg.replace('-'+dig, '\\bar{'+dig+'}')
    spg = '$'+spg+'$'
    return spg


def bash_chk_file_cmd(file):
    #bash returns empty string if file exist
    return " [ -e "+   file     +""" ] || echo "NO"     ; """

# def find_transition_atom(elements):
#     #return list of unique transition elements
#     for elements

def get_common_chemical_base(st1,st2):
    from difflib import SequenceMatcher
    s1 = st1.get_reduced_formula()
    s2 = st2.get_reduced_formula()
    match = SequenceMatcher(None, s1, s2).find_longest_match(0, len(s1), 0, len(s2))
    base  = s1[match.a: match.a + match.size]
    return latex_chem(base)


@contextmanager
def cwd(path):
    oldpwd=os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)
