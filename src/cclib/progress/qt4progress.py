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

from PyQt4 import QtGui, QtCore


class Qt4Progress(QtGui.QProgressDialog):

    def __init__(self, title, parent=None):

        QtGui.QProgressDialog.__init__(self, title, QtCore.QString(), 0, 100, parent)

        self.nstep = 0
        self.text = None
        self.loop = QtCore.QEventLoop(self)
        self.setWindowTitle(title)
        self.setRange(0, 100)
        self.setLabelText("Parsing...")

    def update(self, step, text=None):

        if text:
            self.setLabelText(text)
        self.setValue(step)
        self.loop.processEvents(QtCore.QEventLoop.AllEvents)

if __name__ == "__main__":
    """ An example of the Qt4 progress use. """

    import sys, logging
    from cclib.parser import ccopen

    app = QtGui.QApplication(sys.argv)
    dialog = Qt4Progress("Parsing %s..." % sys.argv[1])
    parser = ccopen(sys.argv[1], dialog.update, loglevel=logging.CRITICAL)
    data = parser.parse()

    sys.exit(app.exec_())

