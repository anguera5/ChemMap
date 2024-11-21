import os

from datetime import datetime
import pandas as pd
from tqdm import tqdm

from ChemMap.chem_requester import ChemRequester
from ChemMap.utils import is_valid_search_method, is_valid_smiles
from ChemMap.enums import AllowedRequestMethods

class ChemMap:
    def __init__(self):
        self.requester = ChemRequester()
        self.compound_data = {}
        self.reaction_data = []
        self.similar_reaction_data = []

    def __reset(self):
        self.__init__()

    def map_smiles_to_proteins(self, smiles, search_method, to_tsv=True):
        is_valid_smiles(smiles)
        is_valid_search_method(search_method)

        for parent_smiles in tqdm(smiles):
            if parent_smiles in self.compound_data.keys():
                # Skipping unnecessary requests
                continue

            self.compound_data[parent_smiles] = self.requester.request_pubchem_and_chebi(parent_smiles, search_method)

            chebi_IDs = self.compound_data[parent_smiles].get("ChEBI")
            current_reaction_data = self.requester.request_to_uniprot(parent_smiles, chebi_IDs, self.reaction_data)

            if search_method in [AllowedRequestMethods.EXPAND_ALL.value, AllowedRequestMethods.EXPAND_PUBCHEM.value]:
                chebi_IDs = self.compound_data[parent_smiles].get("ChEBI")
                self.requester.request_to_uniprot(parent_smiles, chebi_IDs, self.similar_reaction_data,
                                                  reference_reaction_data=current_reaction_data)

        smiles_list = list(self.compound_data.keys())
        compound_data_df = pd.json_normalize([self.compound_data[key] for key in smiles_list])
        compound_data_df.index = pd.Series(smiles_list, name="smiles")

        reaction_data_df = pd.json_normalize(self.reaction_data)
        similar_reaction_data_df = pd.json_normalize(self.similar_reaction_data)
        self.__reset()
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
    smiles = ["CN(C)C(=O)NC1=CC=C(C(=C1)Cl)Cl", "CC(C)C1=CC=C(C=C1)NC(=O)N(C)C"]
    search_method = "expand_all"
    cm = ChemMap()
    print(cm.map_smiles_to_proteins(smiles, search_method="expand_all"))
