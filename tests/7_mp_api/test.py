"""
Test mp-api

Author: Dmitry Aksyonov
"""


from siman.header import _update_configuration
from siman.calc_manage import get_structure_from_matproj_new


_update_configuration('~/simanrc.py') # read configuration, required to run job
try:
    st = get_structure_from_matproj_new(mat_proj_id = 'mp-100')
    # st.printme()
    exit('success')
except:
    exit('failure')
