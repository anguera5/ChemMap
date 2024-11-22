import os

from datetime import datetime
from typing import Union

import pandas as pd
from tqdm import tqdm

from ChemMap.chem_requester import ChemRequester
from ChemMap.utils import is_valid_search_method, is_valid_smiles
from ChemMap.enums import AllowedRequestMethods

class ChemMap:
    """
    A class that represents the main functionality of ChemMap.
    """
    def __init__(self):
        self.requester = ChemRequester()
        self.compound_data = {}
        self.reaction_data = []
        self.similar_reaction_data = []

    def __reset(self):
        """Resets the state of the instance"""
        self.__init__()

    def map_smiles_to_proteins(self, smiles: Union[str | list[str]], search_method: str, to_tsv=True):
        """
        A method that accepts SMILES and a valid search_method as an input and search them in PubChem, ChEBI and UniProt.

        :param smiles: a string or list of strings representing valid SMILES
        :param search_method: one of the values in AllowedRequestMethods
        :param to_tsv: whether to save the output locally as .tsv files
        :return: Three dataframes containing the identifiers from compound databases, the reaction data and the
            reaction data for similar structures (in case the search_method is one of "expand_all" or "expand_pubchem"),
            respectively.
        """
        is_valid_smiles(smiles)
        is_valid_search_method(search_method)

        for parent_smiles in tqdm(smiles):
            if parent_smiles in self.compound_data.keys():
                # Skipping unnecessary requests
                continue

            self.compound_data[parent_smiles] = self.requester.request_pubchem_and_chebi(parent_smiles, search_method)

            chebi_ids = self.compound_data[parent_smiles]["ChEBI"]
            current_reaction_data = self.requester.request_to_uniprot(parent_smiles, chebi_ids, self.reaction_data)

            if search_method in [AllowedRequestMethods.EXPAND_ALL.value, AllowedRequestMethods.EXPAND_PUBCHEM.value]:
                chebi_ids = self.compound_data[parent_smiles]["related_results"]["ChEBI"]
                self.requester.request_to_uniprot(parent_smiles, chebi_ids, self.similar_reaction_data,
                                                  reference_reaction_data=current_reaction_data)

        smiles_list = list(self.compound_data.keys())
        compound_data_df = pd.json_normalize([self.compound_data[key] for key in smiles_list])
        compound_data_df.index = pd.Series(smiles_list, name="smiles")

        reaction_data_df = pd.json_normalize(self.reaction_data)
        similar_reaction_data_df = pd.json_normalize(self.similar_reaction_data)
        self.__reset()
        if to_tsv:
            folder_name = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
            current_path = os.getcwd()
            if not "data" in os.listdir(current_path):
                os.mkdir(current_path + "/data/")
            path = current_path + "/data/" + folder_name
            os.mkdir(path)
            print(f"Saving files on the following path {path}")
            compound_data_df.to_csv(f"{path}/compounds_data.tsv", sep="\t")
            if not reaction_data_df.empty:
                reaction_data_df.to_csv(f"{path}/reaction_data.tsv", sep="\t", index=False)
            if not similar_reaction_data_df.empty:
                similar_reaction_data_df.to_csv(f"{path}/related_reaction_data.tsv", sep="\t", index=False)

        return compound_data_df, reaction_data_df, similar_reaction_data_df
