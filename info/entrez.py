import numpy as np
import matplotlib.pyplot as plt

from Bio import Entrez
from collections import defaultdict
from scipy import stats

"""
Search Entrez for disease information. Define search term,
define character of interest, and return maximum nubmer of publications.
"""

searchterm = "psychopathy"
COI = ":" # character of interest
maxReturn = 10000

Entrez.email = "shussainather@gmail.com"
handle = Entrez.esearch(db="pubmed", term=searchterm, retmax=maxReturn)
record = Entrez.read(handle)
idlist = record["IdList"]
handle = Entrez.efetch(db="pumed", id=idlist, rettype="medline", retmode="xml")
records = Entrez.read(handle)
