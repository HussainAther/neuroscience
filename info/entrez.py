import numpy as np
import matplotlib.pyplot as plt

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
    
