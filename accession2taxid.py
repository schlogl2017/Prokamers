import json
from Bio import Entrez


def accession2taxid(acc, db="nucleotide"):
    handle = Entrez.esearch(db=db, term=acc)
    record = Entrez.read(handle)
    gi = record["IdList"]
    handle = Entrez.esummary(db=db, id=gi, retmode="json")
    result = json.load(handle)["result"]
    print(result)
    taxid = result[gi]["taxid"]
    return str(taxid)
