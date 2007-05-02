import urllib

from BeautifulSoup import BeautifulSoup

wiki = "http://cclib.sourceforge.net/wiki"
# http://openbabel.sourceforge.net/w/index.php?title=Special:Export&action=submit&pages=FAQ
allpages = urllib.urlopen(wiki + "/index.php/Special:Allpages").read()

soup = BeautifulSoup(allpages)

hrefs = []
for anchor in soup('table')[2]('a'):
    hrefs.append(anchor['href'].split("/")[-1])
print hrefs

params = urllib.urlencode({"action":"submit",
                           "pages":"\n".join(hrefs)})

query = urllib.urlopen(wiki + "/index.php/Special:Export", params)
outputfile = open("backup.xml","w")
print >> outputfile, query.read()
outputfile.close()
