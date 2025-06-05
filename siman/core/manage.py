# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import io, json, sys
from siman.classes import CalculationVasp
def setup_shiftk_average(cl, params, list_to_copy):
    """
    Prepare configuration file and pickle file for shiftk_average calculation regime
    """
    print('ngkpt serialized', cl.set.ngkpt)
    # sys.exit()
    cl.serialize(cl.dir+'/shiftk')

    # cl_test = CalculationVasp().deserialize(cl.dir+'/shiftk.pickle')
    # print('ngkpt deserialized', cl_test.set.ngkpt)


    pm = params['shiftk_average']
    pm['vasp_com'] = cl.cluster['vasp_com']
    pm['v'] = cl.id[2]
    pm['chgcar'] = params['chgcar']

    with io.open('shiftk_conf.json', 'w', newline = '') as fp:
        json.dump(pm, fp)


    list_to_copy.append(cl.dir+'/shiftk.pickle')
    list_to_copy.append('shiftk_conf.json')

    return

    