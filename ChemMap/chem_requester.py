from typing import Union

import requests
import re
from time import sleep

import numpy as np
import pandas as pd

from ChemMap.utils import expand_search_chebi, add_or_append_values_to_dict
from ChemMap.enums import AllowedRequestMethods
from ChemMap.utils import uniprot_query

class ChemRequester:
    """A class representing the requester side of ChemMap"""
    PUBCHEM_COMPOUND_DOMAIN = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    PUBCHEM_SMILES_NAMESPACE = "/smiles"
    PUBCHEM_2D_SIMILARITY_NAMESPACE = "/fastsimilarity_2d"
    PUBCHEM_XREFS_NAMESPACE = "/xrefs"

    PUBCHEM_REGISTRY_OPERATION = "/RegistryID"
    PUBCHEM_REST_PROPERTY_OPERATION = "/property"
    PUBCHEM_SYNONYMS_OPERATION = "/synonyms"

    PUBCHEM_JSON_OUTPUT = "/JSON"
    UNIPROT_SPARQL_ENDPOINT = "https://sparql.uniprot.org/sparql/"

    def request_pubchem_and_chebi(self, smiles: str, search_method: str):
        """
        Given a SMILES and search_method strings, this method delegates to the proper data workflow.

        :param smiles: a SMILES string
        :param search_method: a string indicating a valid search method to use, according to AllowedRequestMethods
        :return: A dictionary containing all the extracted
         data for this compound
        """
        url = (self.PUBCHEM_COMPOUND_DOMAIN + self.PUBCHEM_SMILES_NAMESPACE + self.PUBCHEM_SYNONYMS_OPERATION +
               self.PUBCHEM_JSON_OUTPUT)
        # Using exact match SMILES method
        compound_data = self.__execute_request(url, self.__process_pubchem_synonyms, params={"smiles": smiles})
        url = (self.PUBCHEM_COMPOUND_DOMAIN + self.PUBCHEM_SMILES_NAMESPACE + self.PUBCHEM_XREFS_NAMESPACE +
               self.PUBCHEM_REGISTRY_OPERATION + self.PUBCHEM_JSON_OUTPUT)
        compound_data_temp = self.__execute_request(url, self.__process_pubchem_synonyms, params={"smiles": smiles})
        add_or_append_values_to_dict(new_dictionary=compound_data_temp, reference_dictionary=compound_data)
        if (search_method == AllowedRequestMethods.EXPAND_PUBCHEM.value or
                search_method == AllowedRequestMethods.EXPAND_ALL.value):
            # Expand search according to CID match
            url = (self.PUBCHEM_COMPOUND_DOMAIN + self.PUBCHEM_2D_SIMILARITY_NAMESPACE + self.PUBCHEM_SMILES_NAMESPACE +
                   self.PUBCHEM_SYNONYMS_OPERATION + self.PUBCHEM_JSON_OUTPUT)
            temp_related_results = self.__execute_request(url, self.__process_pubchem_synonyms, params={"smiles": smiles})
            compound_data["related_results"] = add_or_append_values_to_dict(new_dictionary=temp_related_results,
                                                                            reference_dictionary=compound_data,
                                                                            empty_dictionary={"CID": [], "ChEBI": []})
        if (search_method == AllowedRequestMethods.EXPAND_CHEBI.value or
                search_method == AllowedRequestMethods.EXPAND_ALL.value):
            compound_data["ChEBI"] = self.expand_chebi(compound_data["ChEBI"])

        if search_method == AllowedRequestMethods.EXPAND_ALL.value:
            # First update chebi IDs and rhea IDs of similar structures
            compound_data["related_results"]["ChEBI"] = self.expand_chebi(compound_data["related_results"]["ChEBI"],
                                                                          excluded_chebi_ids=compound_data["ChEBI"])
        return compound_data

    def request_to_uniprot(self, smiles: str, chebi_ids: list[str], old_reaction_data: list[dict],
                           reference_reaction_data: Union[pd.DataFrame | None]=None):
        """
        Method to perform a request to uniProt given a list of ChEBI IDs

        :param smiles: the SMILES for which the request is being performed
        :param chebi_ids: the ChEBI IDs for that given SMILES
        :param old_reaction_data: a dataframe containing all the stored reaction data until now
        :param reference_reaction_data: None or a dataframe, in which case only the data not present in
            reference_reaction_data will be added to the output
        :return: None or a non-empty reaction dataframe
        """
        # TODO: Enable chunk based queries, this should improved request times
        reaction_query = uniprot_query(chebi_ids).replace("\n", " ")
        reaction_data_df = self.__execute_request(self.UNIPROT_SPARQL_ENDPOINT, self.__process_uniprot_IDs,
                                                  params={"format": "json",
                                                          "query": reaction_query})
        if not reaction_data_df.empty:
            if reference_reaction_data is not None:
                reaction_data_df = reaction_data_df[~reaction_data_df.rhea.isin(reference_reaction_data.rhea)]
            reaction_data_df.insert(0, "smiles", smiles)
            old_reaction_data += reaction_data_df.to_dict("records")
            return reaction_data_df

    def request_pubchem_properties(self, smiles: str, properties: Union[str | list[str]]):
        """
        Given a SMILES and search_method strings, this method delegates to the proper data workflow.

        :param smiles: a SMILES string or a list of SMILES strings
        :param properties: a valid property of PubChem compounds or a list of them (read more
            here: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=Compound-Property-Tables)
        :return: A dictionary containing all the extracted
         data for this compound
        """
        if isinstance(properties, str):
            properties = [properties]
        url = (self.PUBCHEM_COMPOUND_DOMAIN + self.PUBCHEM_SMILES_NAMESPACE + self.PUBCHEM_REST_PROPERTY_OPERATION +
               "/" + ",".join(properties) + self.PUBCHEM_JSON_OUTPUT)
        # Using exact match SMILES method
        compound_data = self.__execute_request(url, self.__process_pubchem_properties,
                                               params={"smiles": smiles})
        return compound_data

    @staticmethod
    def expand_chebi(chebi_ids: list[str], excluded_chebi_ids: Union[list[str] | None]=None):
        """
        Expands the ChEBI IDs provided in a list according to AllowedChEBIRelations

        :param chebi_ids: a list of ChEBI IDs
        :param excluded_chebi_ids: a list of ChEBI IDs to exclude
        :return: A list of non-repeated ChEBI IDs
        """
        unique_chebi_ids = set(chebi_ids)
        for chebi_id in chebi_ids:
            for expanded_chebi_id in expand_search_chebi(chebi_id):
                unique_chebi_ids.add(expanded_chebi_id)
        if excluded_chebi_ids:
            unique_chebi_ids = unique_chebi_ids - set(excluded_chebi_ids)
        return list(unique_chebi_ids)

    def __execute_request(self, url: str, handle_response: callable, method: str="GET", params: Union[dict[str, str], None]=None,
                          back_off_time:int=0.2):
        """
        The general method to generate requests internally. The method expects a URL and a method to handle the response
        if the request has a status_code==200.

        In order to fulfill PubChem's requirements on maximum number of requests per second, I set up a minimum
        back_off_time to 0.2s. Read more here: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

        :param url: A string indicating the endpoint for the request
        :param handle_response: a function that handles the response when the request has status code 200
        :param method: the method of the request, since we are retrieving data this is always "GET"
        :param params: a dictionary of parameters for the request
        :param back_off_time: the number of seconds to back off after each failed request.
        :return: the preprocessed response
        """
        if back_off_time < 0.2:
            raise ValueError("Back off time must be greater than or equal to 0.2")
        sleep(back_off_time)
        back_off_time *= 10
        response = requests.request(method=method, url=url, params=params)
        if response.status_code != 200 and back_off_time > 20:
            print(f"Execution stopped due to exceeding back off time on url: {response.url}")
            return {"CID": [], "ChEBI": []}
        elif response.status_code != 200:
            print(f"Request with status code {response.status_code} failed with error message {response.reason}")
            if 500 <= response.status_code < 600:
                print("trying again with back off time = {} seconds".format(back_off_time))
                return self.__execute_request(url, handle_response, method=method, params=params, back_off_time=back_off_time)
            else:
                response.raise_for_status()
        else:
            return handle_response(response)

    @staticmethod
    def __process_pubchem_synonyms(response: requests.Response):
        """
        Handler for PubChem request

        :param response: the response from the request
        :return: a dictionary containing the CIDs and the ChEBI IDs found on the PubChem endpoint for the given request
        """
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

    @staticmethod
    def __process_pubchem_properties(response: requests.Response):
        """
        Handler for property PubChem request

        :param response: the response from the request
        :return: a list containing dictionaries with the properties found on the PubChem endpoint for the given request
        """
        output = []
        cid_to_synonyms = response.json().get('PropertyTable').get('Properties')
        for items in cid_to_synonyms:
            output.append(items)
        return output

    @staticmethod
    def __process_uniprot_IDs(response: requests.Response):
        """
        Handler for uniProt request

        :param response: the response from the request
        :return: a dataframe containing the rhea, ecNumber and proteins found on the UniProt endpoint for the given request
        """
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
