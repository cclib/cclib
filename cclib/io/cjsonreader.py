# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import json

from cclib.parser.data import ccData


class CJSON:
    """ CJSON log file"""

    def __init__(self, source, *args, **kwargs):

        # Set the filename to source if it is a string.
        # To Do: Add functionality to accept multiple cjson files and streams
        if isinstance(source, str):
            self.filename = source
        else:
            raise ValueError

        self.datatype = {}

    def read_cjson(self):
        inputfile = self.filename

        json_data = json.loads(open(inputfile).read())

        # Actual update of attribute dictionary happens here
        self.construct(json_data)

        return self.datatype

    def construct(self,json_data):
        for Key, Value in ccData._attributes.items():
            jsonKey = Value.jsonKey
            attributePath = Value.attributePath.split(":")

            if attributePath[0] == 'N/A':
                continue

            levels = len(attributePath)
            if attributePath[0] in json_data:
                l1_data_object = json_data[attributePath[0]]
                if levels == 1:
                    if jsonKey in l1_data_object:
                        self.datatype[Key] = l1_data_object[jsonKey]

                elif levels >= 2:
                    if attributePath[1] in l1_data_object:
                        l2_data_object = l1_data_object[attributePath[1]]
                        if jsonKey in l2_data_object:
                            self.datatype[Key] = l2_data_object[jsonKey]

                        if levels == 3 and attributePath[2] in l2_data_object:
                            l3_data_object = l2_data_object[attributePath[2]]
                            if jsonKey in l3_data_object:
                                self.datatype[Key] = l3_data_object[jsonKey]
