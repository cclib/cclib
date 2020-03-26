import scholarly

# Retrieve the author's data and fill-in
search_query = scholarly.search_author('Karol M. Langner')
author = next(search_query).fill()

# Take a closer look at the cclib paper
for i in range(len(author.publications)):
    if author.publications[i] == "Cclib: a library for package‐independent computational chemistry algorithms":
        pub = author.publications[i].fill()

# Which papers cited that publication?
print([citation.bib['title'] for citation in pub.get_citedby()])
