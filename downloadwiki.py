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
