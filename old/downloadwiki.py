# -*- coding: utf-8 -*-
#
# Copyright (c) 2016, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import urllib

from BeautifulSoup import BeautifulSoup


wiki = "http://sourceforge.net/apps/mediawiki/cclib"
# http://openbabel.sourceforge.net/w/index.php?title=Special:Export&action=submit&pages=FAQ
allpages = urllib.urlopen(wiki + "/index.php?title=Special:AllPages").read()

soup = BeautifulSoup(allpages)

hrefs = []
for anchor in soup('table')[2]('a'):
    hrefs.append(anchor['href'].split("/")[-1])
print 'Pages to download:\n', ', '.join(hrefs)

params = urllib.urlencode({"action":"submit",
                           "pages":"\n".join(hrefs)})

query = urllib.urlopen(wiki + "/index.php/Special:Export", params)
outputfile = open("backup.xml","w")
print >> outputfile, query.read()
outputfile.close()
