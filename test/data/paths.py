# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

UNIT_TEST_PATHS = {
    'BOMD': {
        'Gaussian': {
            'basicGaussian09': ('dvb_bomd.out', ),
        },
        'NWChem': {
            'basicNWChem6.6': ('dvb_bomd_ks.out',),
        },
        'QChem': {
            'basicQChem4.2': ('dvb_bomd.out',),
            'basicQChem5.1': ('dvb_bomd.out',),
        },
    },
}
