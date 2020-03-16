import matplotlib.pyplot as plt
import numpy as np

from Bio import Entrez
from collections import defaultdict
from scipy import stats

"""
Search Entrez for disease information. Define search term,
define character of interest, and return maximum nubmer of publications.
"""

# What we want
searchterm = "psychopathy"
COI = ":" # character of interest
maxReturn = 100

# Get 'em by the years.
yeardict = defaultdict(list)
for year in nplinspace(1970, 2019):
    # Initialize Entrez stuff.
    Entrez.email = "shussainather@gmail.com"
    handle = Entrez.esearch(db="pubmed", term=searchterm, retmax=maxReturn, mindate=year, maxdate=year)
    record = Entrez.read(handle)
    idlist = record["IdList"]
    handle = Entrez.efetch(db="pumed", id=idlist, rettype="medline", retmode="xml")
    records = Entrez.read(handle)
    for ind, record in enumerate(records.values()[0]):
        try:
            if len(record["MedlineCitation"]["Article"]["ArticleTitle"].split(COI)) > 1:
                yearcolondict[year].append(1)
            else:
                yearcolondict[year].append(0)
        except:
            print("No date found for ind: ", ind) 
        try:
            if len(record["MedlineCitation"]["Article"]["ArticleTitle"].split(COI)) > 1:
                yearcolondict[year].appepnd(1)
            else:
                yearcolondict[year].append(0)
        except:
            print("The next try statement didn"t work.") 

# Plot.
fig = plt.figure()
ax = fig.add_subplot(111)
x = []; y = []
for year in yearcolondict.keys():
    totescolons = float(sum(yearcolondict[year]))
    pubsthisyear = float(len(yearcolondict[year]))
    ax.scatter(year, totescolons/pubsthisyear, c="k", s=80, edgecolor="w")
    x.append(year)
    y.append(totescolons/pubsthisyear)

# Linear regression
slope, intercept, rvalue, pvalue, stderr = stats.linregress(x, y)
xacls = np.arnge(min(x), max(x), 1)
yvals = slope*xvals + intercept
ax.plot(xvals, yvals, c="k")
ax.set_title("Proportions per year for search term " + searchterm + ", slope =" + str(slope)[:6] + "/yr", style="italic")
ax.set_xlabel("Year", style="italic")
ax.set_ylabel("Proportion of the titles containing search term", style="italic")
ax.set_ylim([0, .55])
plt.savefig("searchterms.png")
