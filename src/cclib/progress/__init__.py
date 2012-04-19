# This file is part of cclib (http://cclib.sf.net), a library for parsing
# and interpreting the results of computational chemistry packages.
#
# Copyright (C) 2006, the cclib development team
#
# The library is free software, distributed under the terms of
# the GNU Lesser General Public version 2.1 or later. You should have
# received a copy of the license along with cclib. You can also access
# the full license online at http://www.gnu.org/copyleft/lgpl.html.


__revision__ = "$Revision$"


import sys

if 'qt' in sys.modules.keys():
    from qtprogress import QtProgress
if 'PyQt4' in sys.modules.keys():
    from qt4progress import Qt4Progress

from textprogress import TextProgress
