import sys
import os

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from baseutil import *

import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
from modelseedpy import MSPackageManager, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem
from modelseedpy.core.annotationontology import convert_to_search_role, split_role
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.helpers import get_template

class AnnoOntUtil(BaseUtil):
    def __init__(self):
        BaseUtil.__init__(self)
        self.get_kbdevutil("AnnoOnt")

    def build_tigrfam_dictionary(self):
        data = pd.read_csv("TIGRFAMs_improvements_2018-10-11.tsv",sep="\t")
        output = {"term_hash":{},"date":"2018-10-11","ontology":"TIGR"}
        for [i,row] in data.iterrows():
            output["term_hash"]["TIGR:"+row["# accession"]] = {
                "id":row["# accession"],
                "name":row["product_name"],
                "type":row["familytype"],
                "synonyms":[row["TC1"],row["TC2"]],
                "length":row["length"],
            }
            if row["genesymbol"] != "NULL":
                output["term_hash"]["TIGR:"+row["# accession"]]["synonyms"].append(row["genesymbol"])
        with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/TIGR_dictionary.json", 'w') as f:
            json.dump(output,f,indent=4,skipkeys=True)
    
    def build_panther_dictionary(self):
        data = pd.read_csv("PANTHER.txt",sep="\t")
        output = {"term_hash":{},"date":"2018-10-11","ontology":"PANTHER"}
        for [i,row] in data.iterrows():
            output["term_hash"]["PANTHER:"+row["ID"]] = {
                "id":row["ID"],
                "name":row["GO"],
            }
        with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/PANTHER_dictionary.json", 'w') as f:
            json.dump(output,f,indent=4,skipkeys=True)

    def build_pfam_dictionary(self):
        data = pd.read_csv("PFAM_GO.txt",sep="\t")
        output = {"term_hash":{},"date":"2018-10-11","ontology":"PFAM"}
        for [i,row] in data.iterrows():
            output["term_hash"]["PFAM:"+row["ID"]] = {
                "id":row["ID"],
                "name":row["GO_name"]
            }
        with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/PFAM_dictionary.json", 'w') as f:
            json.dump(output,f,indent=4,skipkeys=True)

util = AnnoOntUtil()