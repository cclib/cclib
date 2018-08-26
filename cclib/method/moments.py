# -*- coding: utf-8 -*-
#
# Copyright (c) 2018, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculation of electric multipole moments based on data parsed by cclib."""

import numpy

from cclib.parser.utils import convertor
from cclib.method.calculationmethod import Method


class Moments(Method):
    def __init__(self, data):
        self.required_attrs = ('atomcoords', 'atomcharges')

        super(Moments, self).__init__(data)

    def __str__(self):
        """Returns a string representation of the object."""
        return "Multipole moments of %s" % (self.data)

    def __repr__(self):
        """Returns a representation of the object."""
        return 'Moments("%s")' % (self.data)

    def calculate(self, gauge='nuccharge', population='mulliken', **kwargs):
        """Calculate electric dipole moment using parsed partial atomic
        charges.
        
        Args:
            gauge - a choice of gauge origin. Can be either a string or
                a three-element iterable. If iterable, then it
                explicitly defines the origin. If string, then the
                value can be any one of the following and it describes
                what is used as the origin:
                    * 'nuccharge' -- center of positive nuclear charge
                    * 'mass' -- center of mass

        Keyword Args:
            masses - if None, then use default atomic masses. Otherwise,
                the user-provided will be used.
            decimals - a number of decimal places that atomic charges are
                to be truncated to without rounding.

        Returns:
            A list where the first element is the origin of coordinates,
            while other elements are multipole moments expressed in Debye
            units.
        """
        coords = self.data.atomcoords[-1]
        try:
            charges = self.data.atomcharges[population]
        except KeyError:
            raise ValueError

        if hasattr(gauge, '__iter__') and not isinstance(gauge, str):
            origin = numpy.asarray(gauge)
        elif gauge == 'nuccharge':
            origin = numpy.average(coords, weights=self.data.atomnos, axis=0)
        elif gauge == 'mass':
            if 'masses' in kwargs:
                atommasses = numpy.asarray(kwargs.get('masses'))
            else:
                atommasses = self.data.atommasses
            origin = numpy.average(coords, weights=atommasses, axis=0)
        else:
            raise ValueError

        one_debye = convertor(1, 'Angstrom', 'bohr') * convertor(1, 'ebohr', 'Debye')
        if self.data.charge == 0:
            dipole = numpy.dot(charges, coords) * one_debye
        else:
            dipole = numpy.dot(charges, coords - origin) * one_debye

        moments = [origin, dipole]
        return moments
