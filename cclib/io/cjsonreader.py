# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""A reader for chemical JSON (CJSON) files."""

import json

from cclib.io import filereader
from cclib.parser.data import ccData


class CJSON(filereader.Reader):
    """A reader for chemical JSON (CJSON) log files."""

    def __init__(self, source, *args, **kwargs):
        super().__init__(source, *args, **kwargs)

        self.representation = dict()

    def parse(self):
        super().parse()

        json_data = json.loads(self.filecontents)

        self.generate_repr(json_data)

        return self.representation

    def generate_repr(self, json_data):
        for k, v in ccData._attributes.items():
            json_key = v.json_key
            attribute_path = v.attribute_path.split(":")

            if attribute_path[0] == 'N/A':
                continue

            levels = len(attribute_path)
            if attribute_path[0] in json_data:
                l1_data_object = json_data[attribute_path[0]]
                if levels == 1:
                    if json_key in l1_data_object:
                        self.representation[k] = l1_data_object[json_key]

                elif levels >= 2:
                    if attribute_path[1] in l1_data_object:
                        l2_data_object = l1_data_object[attribute_path[1]]
                        if json_key in l2_data_object:
                            self.representation[k] = l2_data_object[json_key]

                        if levels == 3 and attribute_path[2] in l2_data_object:
                            l3_data_object = l2_data_object[attribute_path[2]]
                            if json_key in l3_data_object:
                                self.representation[k] = l3_data_object[json_key]
