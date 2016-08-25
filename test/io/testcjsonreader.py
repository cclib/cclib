# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Unit tests for reader cjsonreader module."""

import os
import unittest
import numpy as np
import tempfile

import cclib


__filedir__ = os.path.dirname(__file__)
__filepath__ = os.path.realpath(__filedir__)
__datadir__ = os.path.join(__filepath__, "..", "..")


class CJSONReaderTest(unittest.TestCase):

    def setUp(self):
        self.CJSON = cclib.io.CJSONReader

    def test_cjson_read(self):
        """File->ccData->CJSON->ccData, the two instances of ccData should be same"""
        fpath = os.path.join(__datadir__, "data/ADF/basicADF2007.01/dvb_gopt.adfout")
        data = cclib.io.ccopen(fpath).parse()

        fp = tempfile.NamedTemporaryFile(mode='w')
        cjson_data = cclib.io.cjsonwriter.CJSON(data, terse=True).generate_repr()
        fp.write(cjson_data)
        fp.flush()

        cjson_Reader = cclib.io.cjsonreader.CJSON(fp.name)
        read_cjson_data = cjson_Reader.read_cjson()

        # The attributes in the object returned by cjson reader should be the same as the
        # attribute in the ccData instance
        ccdata_dict = data.getattributes()

        is_equal = True

        # The attributes inserted into the CJSON is a subset of the total attributes present
        # in the ccData instance
        for key in read_cjson_data:
            ccdata_value = ccdata_dict[key]
            cjson_value = read_cjson_data[key]

            # The values in the ccData object might be of numpy types whereas the values
            # obtained by the cjsonReader are of the inbuilt python types.
            # Conversion of numpy types into python types happens here:
            if isinstance(ccdata_value, np.ndarray) or isinstance(ccdata_value,list):
                ccdata_value = (np.asarray(ccdata_value)).tolist()
            if isinstance(ccdata_value,dict):
                temp_dict = {}
                for ccdata_key in ccdata_value:
                    dict_value = ccdata_value[ccdata_key]
                    if isinstance(dict_value,np.ndarray):
                        dict_value = (np.asarray(dict_value)).tolist()
                        temp_dict[ccdata_key] = dict_value
                ccdata_value = temp_dict

            # The 'moments' attribute present in the CJSON is a post processed value obtained from the
            # moments key within ccData.
            if ccdata_value != cjson_value and key != 'moments':
                is_equal = False
                break

        self.assertEqual(is_equal, True)

if __name__ == "__main__":
    unittest.main()