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


from qt import QProgressDialog


class QtProgress(QProgressDialog):

    def __init__(self, parent):

        QProgressDialog.__init__(self, parent, "progress", True)

        self.nstep = 0
        self.text = None
        self.oldprogress = 0
        self.progress = 0
        self.calls = 0

        self.setCaption("Progress...")

    def initialize(self, nstep, text=None):

        self.nstep = nstep
        self.text = text
        self.setTotalSteps(nstep)
        if text:
            self.setLabelText(text)
        self.setProgress(1)
        #sys.stdout.write("\n")

    def update(self, step, text=None):

        self.setLabelText(text)
        self.setProgress(step)

        return
