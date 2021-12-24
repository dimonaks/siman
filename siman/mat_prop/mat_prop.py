# -*- coding: utf-8 -*- 
# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU General Public License v2.0.

"module contains tabulatated information about Materials propertis"

import json, os


def ionic_radius(el, state, coordination, radius_type = 'r_ionic'):
    """
    https://github.com/prtkm/ionic-radii
    R. D. Shannon, Revised Effective Ionic Radii and Systematic Studies of Interatomic Distances in Halides and Chalcogenides,
    Acta Crystallographica A32 (1976) 751-767.

    INPUT:
        el (str) - name of element
        state (int) - oxidation state, -1,1,2,3 ...
        coordination (str) - roman numbers, VI, IV, and spin state for transition metals, e.g. VIHS or VILS
        spin (str) 
            'HS' - high spin
            'LS' - low spin 
        radius_type (str) - radius type and two additional fields
            'r_ionic'
            'r_crystal'
            'spin' - check possible spin
            'remark'



    """
    # print(os.path.dirname(__file__))
    with open(os.path.dirname(__file__)+"/shannon-radii.json") as f:
        out = f.read()

    d = json.loads(out)

    # Enter Element, Charge, Coordination and one of - r_crystal, r_ionic, spin, remark

    return d[el][str(state)][coordination][radius_type]


def mendeleev():
    from mendeleev import element

    Fe = element(26)

    # for ir in Fe.ionic_radii:
        # print(ir)
    print(Fe.ionic_radii)
    

    from mendeleev.fetch import fetch_ionic_radii
    # irs = fetch_ionic_radii(radius="ionic_radius")
    # print(irs.head(10))