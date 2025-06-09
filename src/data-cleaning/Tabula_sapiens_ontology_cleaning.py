import obonet
import pandas as pd


cl_class = pd.read_table("data/expr/Tabula_sapiens/all_cell_ontology_class.csv")
# This can be downloaded from: https://obofoundry.org/ontology/cl.html
cl_graph = obonet.read_obo("data/ref/cell_ontology/cl.obo")

# fix name for cell ontology
cl_class["official_name"] = cl_class["cell_ontology_class"]
cl_class["official_name"] = cl_class["official_name"].replace(to_replace="cd", value="CD", regex=True)
cl_class["official_name"] = cl_class["official_name"].replace(to_replace="nk", value="NK", regex=True)
cl_class["official_name"] = cl_class["official_name"].replace(to_replace=" t cell", value=" T cell", regex=True)
cl_class["official_name"] = cl_class["official_name"].replace(to_replace="^t cell", value="T cell", regex=True)
cl_class["official_name"] = cl_class["official_name"].replace(to_replace=" b cell", value=" B cell", regex=True)
cl_class["official_name"] = cl_class["official_name"].replace(to_replace="^b cell", value="B cell", regex=True)
cl_class["official_name"] = cl_class["official_name"].replace(to_replace="type i ", value="type I ", regex=True)
cl_class["official_name"] = cl_class["official_name"].replace(to_replace="type ii ", value="type II ", regex=True)
cl_class.loc[cl_class.cell_ontology_class == "cd8-positive alpha-beta t cell", "official_name"] = "CD8-positive, alpha-beta T cell"
cl_class.loc[cl_class.cell_ontology_class == "cd4-positive alpha-beta t cell", "official_name"] = "CD4-positive, alpha-beta T cell"
cl_class.loc[cl_class.cell_ontology_class == "muller cell", "official_name"] = "Mueller cell"
cl_class.loc[cl_class.cell_ontology_class == "nkt cell", "official_name"] = "natural killer cell"
cl_class.loc[cl_class.cell_ontology_class == "langerhans cell", "official_name"] = "Langerhans cell"
cl_class.loc[cl_class.cell_ontology_class == "erythroid progenitor", "official_name"] = "erythroid progenitor cell"
cl_class.loc[cl_class.cell_ontology_class == "myeloid progenitor", "official_name"] = "myeloid lineage restricted progenitor cell"
cl_class.loc[
    cl_class.cell_ontology_class == "hillock-club cell of prostate epithelium", "official_name"
] = "hillock cell of prostatic urethral epithelium"
cl_class.loc[cl_class.cell_ontology_class == "t follicular helper cell", "official_name"] = "T follicular helper cell"
cl_class.loc[cl_class.cell_ontology_class == "dn4 thymocyte", "official_name"] = "DN4 thymocyte"
cl_class.loc[cl_class.cell_ontology_class == "bronchial vessel endothelial cell", "official_name"] = "bronchial epithelial cell"
cl_class.loc[cl_class.cell_ontology_class == "respiratory mucous cell", "official_name"] = "respiratory goblet cell"
cl_class.loc[
    cl_class.cell_ontology_class == "club cell of prostate epithelium", "official_name"
] = "club-like cell of the urethral epithelium"
cl_class.loc[
    cl_class.cell_ontology_class == "hillock cell of prostate epithelium", "official_name"
] = "hillock cell of prostatic urethral epithelium"
cl_class.loc[cl_class.cell_ontology_class == "dn3 thymocyte", "official_name"] = "DN3 thymocyte"
cl_class.loc[cl_class.cell_ontology_class == "dn1 thymic pro-t cell", "official_name"] = "DN1 thymic pro-T cell"
cl_class.loc[cl_class.cell_ontology_class == "schwann cell", "official_name"] = "Schwann cell"
cl_class.loc[cl_class.cell_ontology_class == "artery endothelial cell", "official_name"] = "endothelial cell of artery"
cl_class.loc[cl_class.cell_ontology_class == "pancreatic pp cell", "official_name"] = "pancreatic PP cell"


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
