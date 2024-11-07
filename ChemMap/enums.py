from enum import Enum


class AllowedChEBIRelations(Enum):
    IS_CONJUGATE_ACID_OF = "is_conjugate_base_of"
    IS_CONJUGATE_BASE_OF = "is_conjugate_acid_of"
    IS_A = "is_a"
    IS_TAUTOMER_OF = "is_tautomer_of"
    IS_ENANTIOMER_OF = "is_enantiomer_of"


class AllowedRequestMethods(Enum):
    EXPAND_ALL = "expand_all"
    EXPAND_CHEBI = "expand_chebi"
    EXPAND_PUBCHEM = "expand_pubchem"
    STRICT_MATCH = "strict_match"

