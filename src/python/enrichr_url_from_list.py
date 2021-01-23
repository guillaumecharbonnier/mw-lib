import json
import requests

def enrichr_url_from_list(genes = ["CD4","CD34"]):
    ENRICHR_URL = 'http://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(genes)
    description = 'Gene list queried with API call'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    url = "https://amp.pharm.mssm.edu/Enrichr/enrich?dataset=" + data["shortId"]
    return url
