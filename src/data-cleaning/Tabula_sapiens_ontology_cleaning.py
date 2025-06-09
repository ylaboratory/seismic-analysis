import obonet
import pandas as pd


cl_class = pd.read_table("data/expr/Tabula_sapiens/all_cell_ontology_class.csv")
# This can be downloaded from: https://obofoundry.org/ontology/cl.html
cl_graph = obonet.read_obo("data/ref/cell_ontology/cl.obo")

# fix name for cell ontology
cl_class["official_name"] = cl_class["cell_ontology_class"]

# Apply general regex replacements first
cl_class["official_name"] = cl_class["official_name"].replace(
    to_replace=["cd", "nk", " t cell", "^t cell", " b cell", "^b cell", "type i ", "type ii "],
    value=["CD", "NK", " T cell", "T cell", " B cell", "B cell", "type I ", "type II "],
    regex=True
)

# Define all specific replacements in a dictionary
name_map = {
    "cd8-positive alpha-beta t cell": "CD8-positive, alpha-beta T cell",
    "cd4-positive alpha-beta t cell": "CD4-positive, alpha-beta T cell",
    "muller cell": "Mueller cell",
    "nkt cell": "natural killer cell",
    "langerhans cell": "Langerhans cell",
    "erythroid progenitor": "erythroid progenitor cell",
    "myeloid progenitor": "myeloid lineage restricted progenitor cell",
    "hillock-club cell of prostate epithelium": "hillock cell of prostatic urethral epithelium",
    "t follicular helper cell": "T follicular helper cell",
    "dn4 thymocyte": "DN4 thymocyte",
    "bronchial vessel endothelial cell": "bronchial epithelial cell",
    "respiratory mucous cell": "respiratory goblet cell",
    "club cell of prostate epithelium": "club-like cell of the urethral epithelium",
    "hillock cell of prostate epithelium": "hillock cell of prostatic urethral epithelium",
    "dn3 thymocyte": "DN3 thymocyte",
    "dn1 thymic pro-t cell": "DN1 thymic pro-T cell",
    "schwann cell": "Schwann cell",
    "artery endothelial cell": "endothelial cell of artery",
    "pancreatic pp cell": "pancreatic PP cell"
}

# Apply all specific replacements at once using the dictionary
cl_class["official_name"] = cl_class["official_name"].replace(name_map)


# find cell ontology that is poorly annotated
name_to_id = {data.get('name'): id for id, data in cl_graph.nodes(data=True)}
# synonym to id mapping
synonym_to_id = {}
for node_id, data in cl_graph.nodes(data=True):
    if 'synonym' in data:
        for synonym_entry in data['synonym']:
            # Extract the actual synonym text from between the quotes
            if '"' in synonym_entry:
                # Find the text between the first and second quote
                synonym_text = synonym_entry.split('"')[1]
                synonym_to_id[synonym_text] = node_id
# update
name_to_id.update(synonym_to_id)

cl_class["cl_term"] = cl_class["official_name"].map(name_to_id)
cl_class = cl_class.dropna(subset=["cl_term"])
cl_class = cl_class[["cl_term", "cell_ontology_class", "official_name"]]
cl_class = cl_class[cl_class["cl_term"].str.startswith("CL")]

cl_class.to_csv("data/expr/Tabula_sapiens/all_cell_ontology_class.annotated.new.csv", index=False)
