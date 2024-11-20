import requests
import re
from time import sleep

import numpy as np
import pandas as pd

from ChemMap.utils import expand_search_chebi, add_or_append_values_to_dict
from enums import AllowedRequestMethods
from utils import uniprot_query

class ChemRequester:
    PUBCHEM_REST_SMILES_INPUT = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/synonyms/JSON"
    PUBCHEM_REST_SIMILARITY_OPERATION = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/synonyms/JSON"
    PUBCHEM_REST_REGISTRY_INPUT = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/xrefs/RegistryID/JSON"
    UNIPROT_ENDPOINT = "https://sparql.uniprot.org/sparql/"

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

    def __process_pubchem_synonyms(self, response):
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

    def expand_chebi(self, ChEBI_IDs, excluded_ChEBI_IDs=None):
        unique_chebi_ids = set(ChEBI_IDs)
        # If expand_chebi or expand_all methods are specified then expand results on ChEBI's side
        for chebi_id in ChEBI_IDs:
            for expanded_chebi_id in expand_search_chebi(chebi_id):
                unique_chebi_ids.add(expanded_chebi_id)
        if excluded_ChEBI_IDs:
            unique_chebi_ids = unique_chebi_ids - set(excluded_ChEBI_IDs)
        return list(unique_chebi_ids)

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
        compound_data = self.__execute_request(self.PUBCHEM_REST_SMILES_INPUT, self.__process_pubchem_synonyms,
                                              params={"smiles": smiles})
        compound_data_temp = self.__execute_request(self.PUBCHEM_REST_REGISTRY_INPUT, self.__process_pubchem_synonyms,
                                                    params={"smiles": smiles})
        add_or_append_values_to_dict(new_dictionary=compound_data_temp, reference_dictionary=compound_data)
        if (search_method == AllowedRequestMethods.EXPAND_PUBCHEM.value or
                search_method == AllowedRequestMethods.EXPAND_ALL.value):
            # Expand search according to CID match
            temp_related_results = self.__execute_request(self.PUBCHEM_REST_SIMILARITY_OPERATION,
                                                          self.__process_pubchem_synonyms, params={"smiles": smiles})
            compound_data["related_results"] = add_or_append_values_to_dict(new_dictionary=temp_related_results,
                                                                            reference_dictionary=compound_data,
                                                                            empty_dictionary={"CID": [], "ChEBI": []})
        if (search_method == AllowedRequestMethods.EXPAND_CHEBI.value or
                search_method == AllowedRequestMethods.EXPAND_ALL.value):
            compound_data["ChEBI"] = self.expand_chebi(compound_data["ChEBI"])

        if search_method == AllowedRequestMethods.EXPAND_ALL.value:
            # First update chebi IDs and rhea IDs of similar structures
            compound_data["related_results"]["ChEBI"] = self.expand_chebi(compound_data["related_results"]["ChEBI"],
                                                                          excluded_ChEBI_IDs=compound_data["ChEBI"])
        return compound_data

    def request_to_uniprot(self, smiles, chebi_IDs, old_reaction_data):
        reaction_query = uniprot_query(chebi_IDs).replace("\n", " ")
        reaction_data_df = self.__execute_request(self.UNIPROT_ENDPOINT, self.__process_rhea_IDs,
                                                  params={"format": "json",
                                                          "query": reaction_query})
        if not reaction_data_df.empty:
            reaction_data_df.insert(0, "smiles", smiles)
            old_reaction_data += reaction_data_df.to_dict("records")
        return reaction_data_df