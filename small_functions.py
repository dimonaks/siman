from six import string_types

from collections import Iterable
def is_list_like(obj): 
	return not isinstance(obj, str) and isinstance(obj, Iterable)

def is_string_like(s):
	return isinstance(s, string_types)




