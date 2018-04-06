# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

from PyQt4 import QtGui, QtCore


class Qt4Progress(QtGui.QProgressDialog):

    def __init__(self, title, parent=None):

        QtGui.QProgressDialog.__init__(self, parent)

        self.nstep = 0
        self.text = None
        self.oldprogress = 0
        self.progress = 0
        self.calls = 0
        self.loop=QtCore.QEventLoop(self)
        self.setWindowTitle(title)

    def initialize(self, nstep, text=None):

        self.nstep = nstep
        self.text = text
        self.setRange(0,nstep)
        if text:
            self.setLabelText(text)
        self.setValue(1)
        #sys.stdout.write("\n")

    def update(self, step, text=None):

        if text:
            self.setLabelText(text)
        self.setValue(step)
        self.loop.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)

