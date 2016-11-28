
from six import string_types
from collections import Iterable
import shutil, gzip

from header import printlog


def is_list_like(obj): 
	return not isinstance(obj, string_types) and isinstance(obj, Iterable)

def is_string_like(s):
	return isinstance(s, string_types)



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