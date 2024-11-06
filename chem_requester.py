import os
import requests
import re
from time import sleep

from datetime import datetime
from libchebipy import ChebiEntity
import numpy as np
import pandas as pd
from rdkit import Chem

from enums import AllowedRequestMethods, AllowedChEBIRelations
from utils import rhea_query

def expand_search_chebi(chebi_id):
    ce = ChebiEntity(chebi_id)
    outgoings = ce.get_outgoings()
    outgoings_of_interest = []
    for outgoing in outgoings:
        if outgoing.get_type() in [relation.value for relation in AllowedChEBIRelations]:
            outgoings_of_interest.append(outgoing.get_target_chebi_id())
    return outgoings_of_interest

def execute_request(url, handle_response, method="GET", params=None, back_off_time=0.2):
    # Setting minimum back_off_time to 0.2s to ensure no more than 5 requests per second
    # Read more here: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
    if back_off_time < 0.2:
        raise ValueError("Back off time must be greater than or equal to 0.2")
    sleep(back_off_time)
    back_off_time *= 10
    response = requests.request(method=method, url=url, params=params)
    if response.status_code != 200 and back_off_time > 20:
        print(f"Execution stopped due to exceeding back off time on url: {response.url}")
        return {"CID": [], "ChEBI": []}
    elif response.status_code != 200:
        print("Request failed, trying again with back off time = {} seconds".format(back_off_time))
        return execute_request(url, handle_response, method=method, params=params, back_off_time=back_off_time)
    else:
        return handle_response(response)

def process_pubchem_synonyms(response):
    cids = []
    chebi_ids = set()
    cid_to_synonyms = response.json().get('InformationList').get('Information')
    for items in cid_to_synonyms:
        if cid := items.get('CID'):
            cids.append(cid)
        if terms := items.get('Synonym') or items.get('RegistryID'):
            for term in terms:
                if re.match("CHEBI:-?\\d+", term):
                    chebi_ids.add(term)
    return {"CID": cids, "ChEBI": list(chebi_ids)}

def expand_chebi(data):
    chebi_ids = set(data["ChEBI"])
    # If expand_chebi or expand_all methods are specified then expand results on ChEBI's side
    for chebi_id in data["ChEBI"]:
        for expanded_chebi_id in expand_search_chebi(chebi_id):
            chebi_ids.add(expanded_chebi_id)
    data["ChEBI"] = list(chebi_ids)
    return data

def process_rhea_IDs(response):
    df = pd.json_normalize(response.json()['results']['bindings'])
    # No results
    if df.empty:
        return df
    df.drop(columns=[column for column in df.columns if "type" in column], inplace=True)
    df.columns = df.columns.str.replace(".value", "")
    expected_columns = response.json()['head']['vars']
    # Force the dataframe to have defined columns
    df = df.reindex(columns=expected_columns, fill_value=np.nan)
    # Handle together all cases
    df = df.groupby(["rhea", "ecNumber"], dropna=False)["protein"].apply(lambda x: [] if pd.isna(x).any() else list(x), include_groups=False).reset_index()
    return df


def map_smiles_to_proteins(smiles, search_method, to_csv=True):
    compound_data = {}
    reaction_data = []
    similar_reaction_data = []
    if type(smiles) == str:
        smiles = [smiles]
    if type(smiles) != list:
        raise ValueError("smiles must be str or list")
    for _smiles in smiles:
        if not Chem.MolFromSmiles(_smiles):
            raise ValueError(f"smiles: {_smiles} is not valid")

    if to_csv:
        folder_name = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
        path = os.getcwd() + "/" + folder_name
        os.mkdir(path)
        print(f"Saving files on the following path {path}")

    for parent_smiles in smiles:
        if parent_smiles in compound_data.keys():
            # Skipping unnecessary requests
            continue
        # Using exact match SMILES method
        compound_data[parent_smiles] = execute_request(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/synonyms/JSON", process_pubchem_synonyms,
            params={"smiles": parent_smiles})
        if search_method == AllowedRequestMethods.STRICT_MATCH.value:
            pass
        elif search_method == AllowedRequestMethods.EXPAND_PUBCHEM.value:
            # Expand search according to CID match
            compound_data[parent_smiles]["related_results"] = execute_request(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/synonyms/JSON",
                process_pubchem_synonyms, params={"smiles": parent_smiles})
        elif search_method == AllowedRequestMethods.EXPAND_CHEBI.value:
            compound_data[parent_smiles] = expand_chebi(compound_data[parent_smiles])
        elif search_method == AllowedRequestMethods.EXPAND_ALL.value:
            # Expand search according to CID match
            compound_data[parent_smiles]["related_results"] = execute_request(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/synonyms/JSON",
                process_pubchem_synonyms, params={"smiles": parent_smiles})
            # First update chebi IDs and rhea IDs of similar structures
            compound_data[parent_smiles]["related_results"] = expand_chebi(
                compound_data[parent_smiles]["related_results"])
            compound_data[parent_smiles] = expand_chebi(compound_data[parent_smiles])
        else:
            raise ValueError(
                f"{search_method} is not a valid AllowedRequestMethods, allowed methods are {[meth.value for meth in AllowedRequestMethods]}")

        if chebi_IDs := compound_data[parent_smiles].get("ChEBI"):
            reaction_data_df = execute_request("https://sparql.rhea-db.org/sparql/", process_rhea_IDs,
                                               params={"format": "json",
                                                       "query": rhea_query(chebi_IDs).replace("\n", " ")})
            if not reaction_data_df.empty:
                reaction_data_df.insert(0, "smiles", parent_smiles)
                reaction_data += reaction_data_df.to_dict("records")

        if search_method in [AllowedRequestMethods.EXPAND_ALL.value, AllowedRequestMethods.EXPAND_PUBCHEM.value]:
            if chebi_IDs := compound_data[parent_smiles]["related_results"].get("ChEBI"):
                similar_reaction_data_df = execute_request("https://sparql.uniprot.org/sparql/", process_rhea_IDs,
                                                           params={"format": "json",
                                                                   "query": rhea_query(chebi_IDs).replace("\n", " ")})
                if not similar_reaction_data_df.empty:
                    similar_reaction_data_df.insert(0, "smiles", parent_smiles)
                    similar_reaction_data += similar_reaction_data_df.to_dict("records")

    smiles = list(compound_data.keys())
    df = pd.json_normalize([compound_data[key] for key in smiles])
    df.index = pd.Series(smiles, name="smiles")

    reaction_data_df = pd.json_normalize(reaction_data)
    similar_reaction_data_df = pd.json_normalize(similar_reaction_data)
    if to_csv:
        df.to_csv(f"{path}/compounds_data.tsv", sep="\t")
        if not reaction_data_df.empty:
            reaction_data_df.to_csv(f"{path}/reaction_data.tsv", sep="\t", index=False)
        if not similar_reaction_data_df.empty:
            similar_reaction_data_df.to_csv(f"{path}/similar_reaction_data.tsv", sep="\t", index=False)

    return df, reaction_data_df, similar_reaction_data_df

if __name__ == "__main__":
    smiles = ["C1=C(C(=CC(=C1Cl)Cl)OCC(=O)[O-])Cl", "c1ccccc1"]
    search_method = "expand_all"
    df, rd, srd = map_smiles_to_proteins(smiles, search_method)