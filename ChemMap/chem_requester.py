import os
import requests
import re
from time import sleep

from datetime import datetime
from libchebipy import ChebiEntity
import numpy as np
import pandas as pd

from ChemMap.utils import is_valid_search_method, is_valid_smiles, expand_search_chebi
from enums import AllowedRequestMethods, AllowedChEBIRelations
from utils import uniprot_query

class ChemMap:
    PUBCHEM_REST_SMILES_INPUT = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/synonyms/JSON"
    PUBCHEM_REST_SIMILARITY_OPERATION = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/synonyms/JSON"
    UNIPROT_ENDPOINT = "https://sparql.uniprot.org/sparql/"
    def __init__(self):
        self.compound_data = {}
        self.reaction_data = []
        self.similar_reaction_data = []

    def __reset(self):
        self.compound_data = {}
        self.reaction_data = []
        self.similar_reaction_data = []

    def __execute_request(self, url, handle_response, method="GET", params=None, back_off_time=0.2):
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
            return self.__execute_request(url, handle_response, method=method, params=params, back_off_time=back_off_time)
        else:
            return handle_response(response)

    def process_pubchem_synonyms(self, response):
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

    def expand_chebi(self, data):
        chebi_ids = set(data["ChEBI"])
        # If expand_chebi or expand_all methods are specified then expand results on ChEBI's side
        for chebi_id in data["ChEBI"]:
            for expanded_chebi_id in expand_search_chebi(chebi_id):
                chebi_ids.add(expanded_chebi_id)
        data["ChEBI"] = list(chebi_ids)
        return data

    def __process_rhea_IDs(self, response):
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


    def request_pubchem_and_chebi(self, smiles, search_method):
        # Using exact match SMILES method
        self.compound_data[smiles] = self.__execute_request(self.PUBCHEM_REST_SMILES_INPUT, self.process_pubchem_synonyms,
                                                            params={"smiles": smiles})
        if search_method == AllowedRequestMethods.EXPAND_PUBCHEM.value:
            # Expand search according to CID match
            self.compound_data[smiles]["related_results"] = self.__execute_request(self.PUBCHEM_REST_SIMILARITY_OPERATION,
                                                                                   self.process_pubchem_synonyms,
                                                                                   params={"smiles": smiles})
        elif search_method == AllowedRequestMethods.EXPAND_CHEBI.value:
            self.compound_data[smiles] = self.expand_chebi(self.compound_data[smiles])
        elif search_method == AllowedRequestMethods.EXPAND_ALL.value:
            # Expand search according to CID match
            self.compound_data[smiles]["related_results"] = self.__execute_request(self.PUBCHEM_REST_SMILES_INPUT,
                                                                                   self.process_pubchem_synonyms,
                                                                                   params={"smiles": smiles})
            # First update chebi IDs and rhea IDs of similar structures
            self.compound_data[smiles]["related_results"] = self.expand_chebi(self.compound_data[smiles]["related_results"])
            self.compound_data[smiles] = self.expand_chebi(self.compound_data[smiles])
        return self.compound_data[smiles]

    def request_to_uniprot(self, smiles, chebi_IDs):
        reaction_data_df = self.__execute_request(self.UNIPROT_ENDPOINT, self.__process_rhea_IDs,
                                                  params={"format": "json",
                                                          "query": uniprot_query(chebi_IDs).replace("\n", " ")})
        if not reaction_data_df.empty:
            reaction_data_df.insert(0, "smiles", smiles)
            self.reaction_data += reaction_data_df.to_dict("records")
        return reaction_data_df

    def map_smiles_to_proteins(self, smiles, search_method, to_tsv=True):
        self.__reset()
        is_valid_smiles(smiles)
        is_valid_search_method(search_method)

        for parent_smiles in smiles:
            if smiles in self.compound_data.keys():
                # Skipping unnecessary requests
                continue

            self.request_pubchem_and_chebi(parent_smiles, search_method)

            if chebi_IDs := self.compound_data[parent_smiles].get("ChEBI"):
                reaction_data_df = self.__execute_request(self.UNIPROT_ENDPOINT, self.__process_rhea_IDs,
                                                          params={"format": "json",
                                                                  "query": uniprot_query(chebi_IDs).replace("\n", " ")})
                if not reaction_data_df.empty:
                    reaction_data_df.insert(0, "smiles", parent_smiles)
                    self.reaction_data += reaction_data_df.to_dict("records")

            if search_method in [AllowedRequestMethods.EXPAND_ALL.value, AllowedRequestMethods.EXPAND_PUBCHEM.value]:
                if chebi_IDs := self.compound_data[parent_smiles]["related_results"].get("ChEBI"):
                    similar_reaction_data_df = self.__execute_request("https://sparql.uniprot.org/sparql/", self.__process_rhea_IDs,
                                                               params={"format": "json",
                                                                       "query": uniprot_query(chebi_IDs).replace("\n", " ")})
                    if not similar_reaction_data_df.empty:
                        similar_reaction_data_df.insert(0, "smiles", parent_smiles)
                        self.similar_reaction_data += similar_reaction_data_df.to_dict("records")

        smiles_list = list(self.compound_data.keys())
        compound_data_df = pd.json_normalize([self.compound_data[key] for key in smiles_list])
        compound_data_df.index = pd.Series(smiles_list, name="smiles")

        reaction_data_df = pd.json_normalize(self.reaction_data)
        similar_reaction_data_df = pd.json_normalize(self.similar_reaction_data)
        if to_tsv:
            folder_name = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
            path = os.getcwd() + "/" + folder_name
            os.mkdir(path)
            print(f"Saving files on the following path {path}")
            compound_data_df.to_csv(f"{path}/compounds_data.tsv", sep="\t")
            if not reaction_data_df.empty:
                reaction_data_df.to_csv(f"{path}/reaction_data.tsv", sep="\t", index=False)
            if not similar_reaction_data_df.empty:
                similar_reaction_data_df.to_csv(f"{path}/related_reaction_data.tsv", sep="\t", index=False)

        return compound_data_df, reaction_data_df, similar_reaction_data_df

if __name__ == "__main__":
    smiles = ["C1=C(C(=CC(=C1Cl)Cl)OCC(=O)[O-])Cl", "c1ccccc1"]
    search_method = "expand_all"
    cm = ChemMap()
    print(cm.request_to_uniprot("C1=NC(=O)NC(=C1F)N", chebi_IDs=["CHEBI:58413"]))
    # df, rd, srd = cm.map_smiles_to_proteins(smiles, search_method)

