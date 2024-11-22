# ChemMap
ChemMap is a Python library that tries to bridge the gap from metabolomics to proteomics using existing databases.
A sketch of the main method of ChemMap can be found on the following diagram.

| ![app_schema.png](assets/app_schema.png) | 
|:----------------------------------------:| 
| *Schema showing the workflow of ChemMap* |

The main functionality of ChemMap, the function `map_smiles_to_proteins`, accepts a 
[SMILES](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html) or a list of them and on the first phase 
tries to extract PubChem's and ChEBI's chemical identifiers of this molecule using the 
[PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial). Should you select `"expand_all"`
or `"expand_pubchem"` as parameters of the `search_method`, ChemMap would then find molecules that are structurally 
similar using PUG REST API `fastsimilarity_2d` endpoint, which uses Tanimoto similarity scores. It is noteworthy that
in order to extract ChEBI's identifiers at this stage we are relying on them being reported on PubChem, which might not
be the case for newly reported ChEBI substances.

On the second phase, if either `"expand_all"` or `"expand_chebi"` where selected as input for the parameter 
`search_method`. The workflow will use [libChEBIpy](https://github.com/libChEBI/libChEBIpy) to find substances that are 
related to the ones found by one of the following relationships `is_conjugate_base_of` `is_conjugate_acid_of`, `is_a`, 
`is_tautomer_of` or `is_enantiomer_of`.

On the last step, the ChEBI identifiers are used to search for the presence of the compound on a [Rhea](https://www.rhea-db.org/) 
reaction as a substrate. If we found one, we retrieve the [EC Number](https://en.wikipedia.org/wiki/Enzyme_Commission_number) and 
UniProt protein identifier, if available. On the background we are using the [UniProt SPARQL Endpoint](https://sparql.uniprot.org/) and
the fact that Rhea and UniProt are synchronized on every UniProt release (more 
[here](https://www.uniprot.org/help/synchronization)).

## How to Download

This library can be downloaded through pip 

```bash
pip install chemmap
```

or by direct clone using

```bash
git clone git@github.com:anguera5/ChemMap.git
```

