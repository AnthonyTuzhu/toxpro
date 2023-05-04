"""
interface to interact with PubChem's PUGRest API
"""
import requests
import pandas as pd
import numpy as np
from io import StringIO
from typing import List

def get_raw_bioactivity_data(aid: int):
    """ uses PubChem PUGRest to revieve Bioassay information for a given
     pubchem assay identifier (AID) """
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{}/concise/csv".format(aid)

    r = requests.get(url)

    if r.status_code == 200:
        s = StringIO(r.text)
        df = pd.read_csv(s)
        return df[['CID', 'Activity Outcome']], None
    else:
        return r.status_code, r.reason


def clean_bioactivity_frame(df):
    """ converts activity to 1, 0 """
    df['Activity Outcome'] = df['Activity Outcome'] \
                            .replace('Inactive', 0) \
                            .replace('Active', 1) \
                            .replace('Probe', 1) \
                            .replace('Inconclusive', np.nan) \
                            .replace('Unspecified', np.nan) \
                            .dropna()
    # if there are multiple responses for a CID
    # take the conservative, "Active" response
    df = df.sort_values('Activity Outcome')
    df = df.drop_duplicates('CID', keep='last')
    return df


def get_inchi_from_cids(cids: List):
    """ just as the name says, gets inchis from a list of cids """

    identifier_list = list(map(str, cids))

    # make the base URL for the PubChem POST Request
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/inchi/csv'
    #url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/MolecularFormula,MolecularWeight/CSV'

    headers = {'Content-Type': 'Content-Type: application/x-www-form-urlencoded'}
    data = {'cid': ','.join(identifier_list)}

    r= requests.post(url, data=data)

    if r.status_code == 200:
        s = StringIO(r.text)
        df = pd.read_csv(s)
        return df, None
    else:
        return r.status_code, r.reason



def get_assay_name(aid: int):
    """ uses PubChem PUGRest to revieve Bioassay name """
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{}/description/json".format(aid)

    r = requests.get(url)

    if r.status_code == 200:

        json_data = r.json()
        name = json_data["PC_AssayContainer"][0]['assay']['descr']['name']
        return name, None
    else:
        return r.status_code, r.reason

def import_pubchem_aid(aid: int):
    """ given a pubchem identifier import all the bioassay data and
    structure information in a format ready for importing to the sqlite database """

    result, reason = get_raw_bioactivity_data(aid)
    if reason != None:
        return 'Failed', reason

    bioactivity = clean_bioactivity_frame(result)
    result, reason = get_inchi_from_cids(bioactivity.CID.values.tolist())
    if reason != None:
        return 'Failed', reason

    result = result.merge(bioactivity)

    result = result.rename({
        'CID': 'compound_id',
        'InChI': 'inchi',
        'Activity Outcome': 'activity'
    },
    axis=1
    )

    return result, None

print(get_assay_name(1000))