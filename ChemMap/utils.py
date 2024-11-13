from libchebipy import ChebiEntity
from rdkit import Chem
from enums import AllowedRequestMethods, AllowedChEBIRelations


def is_valid_smiles(smiles):
        if type(smiles) == str:
                smiles = [smiles]
        if type(smiles) != list:
                raise ValueError("smiles must be str or list")
        for _smiles in smiles:
                if not Chem.MolFromSmiles(_smiles):
                        raise ValueError(f"smiles: {_smiles} is not valid")


def is_valid_search_method(search_method):
        try:
                AllowedRequestMethods(search_method)
        except Exception:
                raise ValueError(f"{search_method} is not a valid AllowedRequestMethods, allowed methods "
                                 f"are {[meth.value for meth in AllowedRequestMethods]}")

def expand_search_chebi(chebi_id):
        ce = ChebiEntity(chebi_id)
        outgoings = ce.get_outgoings()
        outgoings_of_interest = []
        for outgoing in outgoings:
            if outgoing.get_type() in [relation.value for relation in AllowedChEBIRelations]:
                outgoings_of_interest.append(outgoing.get_target_chebi_id())
        return outgoings_of_interest

def uniprot_query(ChEBiIDs):
        return ("PREFIX rh: <http://rdf.rhea-db.org/>\nPREFIX CHEBI: <http://purl.obolibrary.org/obo/CHEBI_>\n "
                "PREFIX up: <http://purl.uniprot.org/core/> \n"
            " SELECT Distinct ?rhea ?ecNumber ?protein\n"
            "WHERE { \n ?rhea rh:side ?reactionSide1 . \n ?reactionSide1  rh:contains / rh:compound / rh:chebi "
            "?chebi .\n ?reactionSide1 rh:transformableTo ?reactionSide2 .\n "
            "OPTIONAL{?rhea rh:ec ?ecNumber . \n"
            "?protein ( up:enzyme | up:domain/up:enzyme | up:component/up:enzyme ) ?ecNumber . \n"
            "?protein up:sequence ?isoform .}\n"
            "VALUES (?chebi) {" + " ".join(["(" + ChEBiID + ")" for ChEBiID in ChEBiIDs]) + "}\n}")
