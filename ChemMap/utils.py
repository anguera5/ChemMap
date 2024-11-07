def rhea_query(ChEBiIDs):
    return ("PREFIX rh: <http://rdf.rhea-db.org/>\nPREFIX CHEBI: <http://purl.obolibrary.org/obo/CHEBI_>\n "
            "PREFIX up: <http://purl.uniprot.org/core/> \n"
            " SELECT Distinct ?rhea ?ecNumber ?protein\n"
            "WHERE { \n ?rhea rh:side ?reactionSide1 . \n ?reactionSide1  rh:contains / rh:compound / rh:chebi "
            "?chebi .\n ?reactionSide1 rh:transformableTo ?reactionSide2 .\n "
            "OPTIONAL{?rhea rh:ec ?ecNumber . \n"
            "?protein ( up:enzyme | up:domain/up:enzyme | up:component/up:enzyme ) ?ecNumber . \n"
            "?protein up:sequence ?isoform .}\n"
            "VALUES (?chebi) {" + " ".join(["(" + ChEBiID + ")" for ChEBiID in ChEBiIDs]) + "}\n}")
