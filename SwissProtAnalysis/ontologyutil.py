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

class OntologyUtil(BaseUtil):
    def __init__(self):
        BaseUtil.__init__(self)
        self.obs_ec = None
        self.get_kbdevutil("Ontology")

    #Code for translating obsolete EC numbers
    def trans_ec(self,ec):
        if not self.obs_ec:
            with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/obsolete_ec.json") as json_file:
                self.obs_ec = json.load(json_file)
        original_ec = ec
        count=0
        while ec in self.obs_ec:
            count += 1
            if count == 20:
                #print("Circular reference:",original_ec,"->",ec)
                return original_ec
            ec = self.obs_ec[ec]
        return ec

    def create_msrxn_data(self):
        biochem = ModelSEEDBiochem.get()

        filtered_reactions = pd.read_csv(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/FilteredReactions.csv",sep="\t")
        filtered_reaction_hash = {}
        for [i,row] in filtered_reactions.iterrows():
            filtered_reaction_hash[row["id"]] = row["reason"]

        msrxn_data = {}
        for rxn in biochem.reactions:
            msrxn_data[rxn.id] = {
                "id":rxn.id,
                "name":rxn.name,
                "equation":rxn.build_reaction_string(use_metabolite_names=True),
                "ec":[],
                "filtered":None
            }
            ecnums = rxn.ec_numbers
            
            for ec in ecnums:
                ec = self.trans_ec(ec)
                if ec not in msrxn_data[rxn.id]["ec"]:
                    msrxn_data[rxn.id]["ec"].append(ec)

            if rxn.id in filtered_reaction_hash:
                msrxn_data[rxn.id]["filtered"] = filtered_reaction_hash[rxn.id]

        reaction_ec = pd.read_csv(self.kbdevutil.config["data"]+"/ModelSEEDDatabase/Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt",sep="\t")
        for [i,row] in reaction_ec.iterrows():
            if row["ModelSEED ID"] in msrxn_data:
                ec = row["External ID"]
                ec = self.trans_ec(ec)
                if ec not in msrxn_data[row["ModelSEED ID"]]["ec"]:
                    msrxn_data[row["ModelSEED ID"]]["ec"].append(ec)

        self.kbdevutil.save("msrxn_data",msrxn_data)

    def create_rhea_data(self):
        rhea_data = {}
        msrxn_data = self.kbdevutil.load("msrxn_data")
        #Loading GO terms to get names for rhea IDs because GO has a single reaction resolution that corresponds to Rhea
        #A problem with this is that not every Rhea ID has a GO term; also, we are assuming the GO mappings are right, and they may not be
        with open(self.kbdevutil.codebase+'/cb_annotation_ontology_api/data/GO_dictionary.json') as json_file:
            go_dictionary = json.load(json_file)

        #Here I'm also reading in the EC mappings for Rhea so I can get good EC numbers for the Rhea IDs
        ec_trans = pd.read_csv(self.kbdevutil.config["data"]+"/TemplateFunctions/rhea2ec.tsv",sep="\t")
        for [i,row] in ec_trans.iterrows():
            row["RHEA_ID"] = str(row["RHEA_ID"])
            if row["RHEA_ID"] not in rhea_data:
                rhea_data[row["RHEA_ID"]] = {
                    "id":row["RHEA_ID"],
                    "ec":[],
                    "name":None,
                    "genes":[],
                    "msrxn":[]
                }
            ecnum = self.trans_ec(row["ID"])
            if ecnum not in rhea_data[row["RHEA_ID"]]["ec"]:
                rhea_data[row["RHEA_ID"]]["ec"].append(row["ID"])

        #Here I load alternative EC mappings for Rhea so I can get good EC numbers for the Rhea IDs
        ec_trans = pd.read_csv(self.kbdevutil.config["data"]+"/TemplateFunctions/rhea-ec-iubmb.tsv",sep="\t")
        for [i,row] in ec_trans.iterrows():
            row["RHEA_ID"] = str(row["RHEA_ID"])
            if row["RHEA_ID"] not in rhea_data:
                rhea_data[row["RHEA_ID"]] = {
                    "id":row["RHEA_ID"],
                    "ec":[],
                    "name":None,
                    "genes":[],
                    "msrxn":[]
                }
            ecnum = self.trans_ec(row["EC"])
            if ecnum not in rhea_data[row["RHEA_ID"]]["ec"]:
                rhea_data[row["RHEA_ID"]]["ec"].append(ecnum)

        #Here I use the GO mappings to assign names to the Rhea IDs                                      
        go_trans = pd.read_csv(self.kbdevutil.config["data"]+"/TemplateFunctions/rhea2go.tsv",sep="\t")
        for [i,row] in go_trans.iterrows():
            row["RHEA_ID"] = str(row["RHEA_ID"])
            if row["RHEA_ID"] not in rhea_data:
                rhea_data[row["RHEA_ID"]] = {
                    "id":row["RHEA_ID"],
                    "ec":[],
                    "name":None,
                    "genes":[],
                    "msrxn":[]
                }
            rhea_data[row["RHEA_ID"]]["name"] = row["ID"]
            if row["ID"] in go_dictionary["term_hash"]:
                rhea_data[row["RHEA_ID"]]["name"] += ":"+go_dictionary["term_hash"][row["ID"]]["name"]

        ms_aliases = pd.read_csv(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/ModelSEED_Reaction_Aliases.txt",sep="\t")
        for [i,row] in ms_aliases.iterrows():
            if row["Source"] == "Rhea":
                if row["External ID"] not in rhea_data:
                    rhea_data[row["External ID"]] = {
                        "id":row["External ID"],
                        "ec":[],
                        "name":None,
                        "genes":[],
                        "msrxn":[]
                    }
                msrxn = row["ModelSEED ID"]
                if msrxn not in rhea_data[row["External ID"]]["msrxn"]:
                    rhea_data[row["External ID"]]["msrxn"].append(msrxn)
                if msrxn in msrxn_data:
                    if msrxn_data[msrxn]["name"] and not rhea_data[row["External ID"]]["name"]:
                        rhea_data[row["External ID"]]["name"] = msrxn_data[msrxn]["name"]
                    for ec in msrxn_data[msrxn]["ec"]:
                        ec = self.trans_ec(ec)
                        if ec not in rhea_data[row["External ID"]]["ec"]:
                            rhea_data[row["External ID"]]["ec"].append(ec)

        #TODO: We need a new more robust naming procedure for Rhea reactions
        #Step one, if the Rhea is mapped to a ModelSEED reaction with a name, use that name
        #Step two, if there is no MS rxn or the MS rxn has no name or just and rxn ID or Rhea ID for a name, do the following:
        #Take the "activity" from the first number in the EC number (e.g. 1 = oxidoreductase, 2 = transferase, etc.)
        #Then print the reactant list and produce list with the activity of the first EC (e.g. Glucose-6-phosphate hydrolase)
        #Possible explore filtering out obvious cofactors from the lists above (e.g. O2, H+, H2O)

        self.kbdevutil.save("rhea_data",rhea_data)
    
    def create_sso_data(self):
        sso_data = {}
        msrxn_data = self.kbdevutil.load("msrxn_data")

        with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/SSO_dictionary.json") as json_file:
            sso = json.load(json_file)

        with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/SSO_reactions.json") as json_file:
            sso_rxns = json.load(json_file)

        for sso_id in sso["term_hash"]:
            sso_id = sso_id[4:]
            sso_data[sso_id] = {
                "id":sso_id,
                "name":sso["term_hash"]["SSO:"+sso_id]["name"],
                "ec":[],
                "genes":[],
                "msrxn":[],
                "class":None
            }
            if sso_id == "000009137":
                sso_data[sso_id]["class"] = "hypothetical"
            match = re.search(r'(\d+\.[\d-]+\.[\d-]+\.[\d-]+)',sso["term_hash"]["SSO:"+sso_id]["name"])
            if match:
                ec = match.group(0)
                ec = self.trans_ec(ec)
                if ec not in sso_data[sso_id]["ec"]:
                    sso_data[sso_id]["ec"].append(ec)
                
        for sso_id in sso_rxns:
            orig = sso_id
            sso_id = sso_id[4:]
            if sso_id in sso_data:
                for rxn in sso_rxns[orig]:
                    if rxn not in sso_data[sso_id]["msrxn"]:
                        sso_data[sso_id]["msrxn"].append(rxn)
                    if rxn in msrxn_data:
                        for ec in msrxn_data[rxn]["ec"]:
                            if ec not in sso_data[sso_id]["ec"]:
                                sso_data[sso_id]["ec"].append(ec)

        self.kbdevutil.save("sso_data",sso_data)

    def create_ko_data(self):
        ko_data = {}
        msrxn_data = self.kbdevutil.load("msrxn_data")

        with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/KO_dictionary.json") as json_file:
            kodict = json.load(json_file)

        for ko in kodict["term_hash"]:
            if ko not in ko_data:
                ko_data[kodict["term_hash"][ko]["id"]] = {
                    "id": kodict["term_hash"][ko]["id"],
                    "name": kodict["term_hash"][ko]["name"],
                    "ec":[],
                    "genes":[],
                    "msrxn":[]
                }
                match = re.search(r'(\d+\.[\d-]+\.[\d-]+\.[\d-]+)',kodict["term_hash"][ko]["name"])
                if match:
                    ec = match.group(0)
                    ec = self.trans_ec(ec)
                    if ec not in ko_data[kodict["term_hash"][ko]["id"]]["ec"]:
                        ko_data[kodict["term_hash"][ko]["id"]]["ec"].append(ec)

        with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/kegg_95_0_ko_seed.tsv") as f:
            korxn = pd.read_csv(f,sep="\t")

        for index, row in korxn.iterrows():
            if row["ko_id"] in ko_data:
                rxns = row["seed_ids"].split(";")
                for rxn in rxns:
                    if rxn not in ko_data[row["ko_id"]]["msrxn"]:
                        ko_data[row["ko_id"]]["msrxn"].append(rxn)
                    if rxn in msrxn_data:
                        for ec in msrxn_data[rxn]["ec"]:
                            if ec not in ko_data[row["ko_id"]]["ec"]:
                                ko_data[row["ko_id"]]["ec"].append(ec)

        self.kbdevutil.save("ko_data",ko_data)

    def create_ec_data(self):
        ec_data = {}

        with open(self.kbdevutil.codebase+"/cb_annotation_ontology_api/data/EC_dictionary.json") as json_file:
            ecdict = json.load(json_file)

        for ec in ecdict["term_hash"]:
            if ec not in ec_data:
                ec_data[ecdict["term_hash"][ec]["id"]] = {
                    "id": ecdict["term_hash"][ec]["id"],
                    "name": ecdict["term_hash"][ec]["name"],
                    "ec":[ec],
                    "dram_genes":[],
                    "prokka_genes":[],
                    "msrxn":[],
                    "rhea":[],
                    "sso":[],
                    "ko":[]
                }

        rhea_data = self.kbdevutil.load("rhea_data")
        sso_data = self.kbdevutil.load("sso_data")
        msrxn_data = self.kbdevutil.load("msrxn_data")
        ko_data = self.kbdevutil.load("ko_data")

        all_data = {
            "rhea":rhea_data,
            "sso":sso_data,
            "msrxn":msrxn_data,
            "ko":ko_data
        }

        for type in all_data:
            for id in all_data[type]:
                for ec in all_data[type][id]["ec"]:
                    if ec not in ec_data:
                        ec_data[ec] = {"id":ec,"name":None,"ec":[ec],"rhea":[],"sso":[],"msrxn":[],"ko":[],"dram_genes":[],"prokka_genes":[]}
                    ec_data[ec][type].append(id)
                    
        self.kbdevutil.save("ec_data",ec_data)

    def pull_swissprot_annotations(self):
        # Load annotation ontology events
        output = annoapi.get_annotation_ontology_events({
            "input_ref": "154984/SwissProtCuratedProteins.RAST.Prokka.DRAM.Rhea2.glm4ec"
        })
        kbdevutil.save("swiss_prot_anno", output)

        annotations_by_gene = {}
        for event in output["events"]:
            print(event["event_id"])
            for gene in event["ontology_terms"]:
                if gene not in annotations_by_gene:
                    annotations_by_gene[gene] = {}
                if event["event_id"].startswith("RAST"):
                    annotations_by_gene[gene]["sso"] = event["ontology_terms"][gene]
                elif event["event_id"].startswith("DRAM:KO"):
                    annotations_by_gene[gene]["ko"] = event["ontology_terms"][gene]
                elif event["event_id"].startswith("DRAM:EC"):
                    annotations_by_gene[gene]["dec"] = event["ontology_terms"][gene]
                elif "RHEA" in event["id"]:  # Adjusted to check "RHEA" in "id"
                    annotations_by_gene[gene]["rhea"] = event["ontology_terms"][gene]
                elif event["event_id"].startswith("Prokka"):
                    annotations_by_gene[gene]["pec"] = event["ontology_terms"][gene]
                elif event["event_id"].startswith("GLM4EC"):
        #            print("TEST")
                    annotations_by_gene[gene]["glm"] = event["ontology_terms"][gene]

        # Save updated data
        kbdevutil.save("annotations_by_gene", annotations_by_gene)

    def add_genes_to_annotation_data(self):
        annotations_by_gene = kbdevutil.load("annotations_by_gene")

        rhea_data = kbdevutil.load("rhea_data")
        sso_data = kbdevutil.load("sso_data")
        msrxn_data = kbdevutil.load("msrxn_data")
        ko_data = kbdevutil.load("ko_data")
        ec_data = kbdevutil.load("ec_data")
        all_data = {
            "rhea": rhea_data,
            "sso": sso_data,
            "msrxn": msrxn_data,
            "ko": ko_data,
            "ec": ec_data,
            "glm": ec_data  # Assuming this is intended to be 'ec_data' as placeholder for 'glm' data
        }

        # Update annotations in data files
        for gene in annotations_by_gene:
            for source in annotations_by_gene[gene]:
                if source != "swiss_prot_name" and source != "swiss_prot_EC": 
                    for item in annotations_by_gene[gene][source]:
                        term = item["term"].split(":")[1]
                        if source in all_data:
                            if term not in all_data[source]:
                                print(f"Term not found: {term} in {source}, for gene: {gene}")
                                continue  # Skip to next item
                            # Initialize 'genes' key if missing
                            if "genes" not in all_data[source][term]:
                                all_data[source][term]["genes"] = []
                            if gene not in all_data[source][term]["genes"]:
                                all_data[source][term]["genes"].append(gene)
                        elif source == "dec" or source == "pec":
                            if term not in all_data["ec"]:
                                print(f"Term not found: {term} in {source} (handled as 'ec'), for gene: {gene}")
                                continue
                            key = "dram_genes" if source == "dec" else "prokka_genes"
                            if key not in all_data["ec"][term]:
                                all_data["ec"][term][key] = []
                            if gene not in all_data["ec"][term][key]:
                                all_data["ec"][term][key].append(gene)
                        else:
                            print(f"Unknown source: {source}, for gene: {gene}")

        kbdevutil.save("rhea_data", rhea_data)
        kbdevutil.save("sso_data", sso_data)
        kbdevutil.save("msrxn_data", msrxn_data)
        kbdevutil.save("ko_data", ko_data)
        kbdevutil.save("ec_data", ec_data)

    def load_swissprot_names(self):
        # eventualy we can change to use Filipe's arango database, using the uniprot API with the code below it takes ~9 MINUTES
        import requests
        from concurrent.futures import ThreadPoolExecutor, as_completed
        import time

        def query_uniprot(accession_id):
            try:
                # Remove '.CDS' suffix if present
                cleaned_accession_id = accession_id.split('.CDS')[0]
                url = f'https://rest.uniprot.org/uniprotkb/{cleaned_accession_id}.json'
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    protein_name = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Name not found')
                    return accession_id, protein_name  # Return original accession_id to maintain consistency
                else:
                    return accession_id, 'Name not found'
            except Exception as e:
                return accession_id, 'Request failed'

        def process_accessions_batch(accession_ids, batch_size=5000, sleep_time=10):
            results = {}
            total_batches = (len(accession_ids) + batch_size - 1) // batch_size
            start_time = time.time()  # Start timing

            for batch_number in range(total_batches):
                start_index = batch_number * batch_size
                end_index = start_index + batch_size
                batch = list(accession_ids)[start_index:end_index]
                print(f"Processing batch {batch_number + 1} of {total_batches}...")

                with ThreadPoolExecutor(max_workers=20) as executor:
                    future_to_accession = {executor.submit(query_uniprot, accession_id): accession_id for accession_id in batch}
                    for future in as_completed(future_to_accession):
                        accession_id = future_to_accession[future]
                        try:
                            _, protein_name = future.result()
                            results[accession_id] = protein_name
                        except Exception as exc:
                            results[accession_id] = 'Error retrieving name'
                print(f"Batch {batch_number + 1} completed. Sleeping for {sleep_time} seconds...")
                time.sleep(sleep_time)

            end_time = time.time()  # End timing
            print(f"All batches processed. Total time: {end_time - start_time:.2f} seconds.")
            return results

        def enrich_annotations_with_swiss_prot_names(annotations_by_gene):
            # Extract and clean accession IDs, removing any '.CDS' suffix
            accession_ids = set(gene.split('.CDS')[0] for gene in annotations_by_gene.keys())
            swiss_prot_names = process_accessions_batch(accession_ids)

            # Update annotations_by_gene with fetched Swiss-Prot names, accounting for '.CDS' in original keys
            for gene, data in annotations_by_gene.items():
                cleaned_gene = gene.split('.CDS')[0]
                swiss_prot_name = swiss_prot_names.get(cleaned_gene, "Name not found")
                data["swiss_prot_name"] = swiss_prot_name

        # Assuming kbdevutil is your custom utility for loading and saving data
        # Load existing annotations
        annotations_by_gene = kbdevutil.load("annotations_by_gene")

        # Enrich annotations with Swiss-Prot names
        enrich_annotations_with_swiss_prot_names(annotations_by_gene)

        # Save updated annotations
        kbdevutil.save("annotations_by_gene", annotations_by_gene)

    def create_domain_gene_lists(self):
        domain_specific_lists = {
            "Fungi" : "154984/SwissProt_Rhea_Fungi",
            "Other" : "154984/SwissProt_Rhea_Other",
            "Viridiplantae" : "154984/SwissProt_Rhea_Viridiplantae",
            "Archaea" : "154984/SwissProt_Rhea_Archaea",
            "Bacteria" : "154984/SwissProt_Rhea_Bacteria",
            "Metazoa" : "154984/SwissProt_Rhea_Metazoa"
        }
        domain_proteins = {}
        for domain in domain_specific_lists:
            data = kbdevutil.get_object(domain_specific_lists[domain])
            for item in data["data"]["sequences"]:
                domain_proteins[item["id"]] = domain

        kbdevutil.save("domain_proteins",domain_proteins)

    def print_ontology_names_for_comparison(self):
        import pandas as pd
        import time

        # Load the annotations and domain_proteins data
        annotations_by_gene = kbdevutil.load("annotations_by_gene")
        domain_proteins = kbdevutil.load("domain_proteins")  # Load the domain information
        all_data = {
            "rhea": kbdevutil.load("rhea_data"),
            "sso": kbdevutil.load("sso_data"),
            "msrxn": kbdevutil.load("msrxn_data"),
            "ko": kbdevutil.load("ko_data"),
            "dec": kbdevutil.load("ec_data"),
            "pec": kbdevutil.load("ec_data"),
            "glm": kbdevutil.load("ec_data")
        }

        records = {"Term1": [], "Term2": [], "Name1": [], "Name2": [], "Source1": [], "Source2": [], "Gene": [], "Domain": [], "ReactionMatch": [], "ECMatch": []}

        start_time = time.time()

        for gene, gene_data in annotations_by_gene.items():
            ontology_sources = [source for source in gene_data if source != "swiss_prot_name"]
            swiss_prot_name = gene_data.get("swiss_prot_name", "Name not found") or "Name not found"
            gene_id = gene.split('.')[0]  # Adjust as necessary
            domain = domain_proteins.get(gene_id, "Unknown")

            for i, source1 in enumerate(ontology_sources):
                for j in range(i + 1, len(ontology_sources)):
                    source2 = ontology_sources[j]
                    for item in gene_data.get(source1, []):
                        for oitem in gene_data.get(source2, []):
                            term1 = item.get("term")
                            term2 = oitem.get("term")
                            term1_id = term1.split(":")[1] if term1 else ""
                            term2_id = term2.split(":")[1] if term2 else ""

                            # For rhea terms, use swiss_prot_name, otherwise fetch the name from all_data
                            name1 = swiss_prot_name if source1 == "rhea" else all_data.get(source1, {}).get(term1_id, {}).get("name", "Name not found") or "Name not found"
                            name2 = swiss_prot_name if source2 == "rhea" else all_data.get(source2, {}).get(term2_id, {}).get("name", "Name not found") or "Name not found"

                            # Append EC details if available and not rhea terms
                            if source1 != "rhea" and "ec" in all_data.get(source1, {}).get(term1_id, {}):
                                ec_details1 = ";".join(all_data[source1].get(term1_id, {}).get('ec', []))
                                if ec_details1:
                                    name1 += f" ({ec_details1})"
                            if source2 != "rhea" and "ec" in all_data.get(source2, {}).get(term2_id, {}):
                                ec_details2 = ";".join(all_data[source2].get(term2_id, {}).get('ec', []))
                                if ec_details2:
                                    name2 += f" ({ec_details2})"

                            reaction_match = "No"
                            ec_match = "No"

                            # Check for reaction and EC matches
                            msrxn1 = set(all_data.get(source1, {}).get(term1_id, {}).get("msrxn", []))
                            msrxn2 = set(all_data.get(source2, {}).get(term2_id, {}).get("msrxn", []))
                            reaction_match = "Yes" if msrxn1 & msrxn2 else "No"

                            ec1 = set(all_data.get(source1, {}).get(term1_id, {}).get("ec", []))
                            ec2 = set(all_data.get(source2, {}).get(term2_id, {}).get("ec", []))
                            ec_match = "Yes" if ec1 & ec2 else "No"

                            records["Term1"].append(term1)
                            records["Term2"].append(term2)
                            records["Name1"].append(name1)
                            records["Name2"].append(name2)
                            records["Source1"].append(source1)
                            records["Source2"].append(source2)
                            records["Gene"].append(gene)
                            records["Domain"].append(domain)
                            records["ReactionMatch"].append(reaction_match)
                            records["ECMatch"].append(ec_match)

        # Convert records to DataFrame and save
        df = pd.DataFrame.from_dict(records)

        # Specify the output directory and filename
        output_filename = "/annotation_pairs.csv"
        output_path = kbdevutil.output_dir + output_filename
        df.to_csv(output_path, index=False)

        end_time = time.time()
        print(f"Process completed in {end_time - start_time:.2f} seconds. Output saved to {output_path}")

    def print_condensed_table_of_unique_name_pairings(self):
        import pandas as pd

        # Assuming you have the path to the generated annotation_pairs.csv
        input_filename = "/annotation_pairs.csv"
        annotation_pairs_path = kbdevutil.output_dir + input_filename
        output_filename = "/annotation_pairs_condensed.csv"
        condensed_output_path = kbdevutil.output_dir + output_filename

        # Read the annotation_pairs.csv into a DataFrame
        annotation_pairs_df = pd.read_csv(annotation_pairs_path)

        # Condense the DataFrame by grouping on 'Name1' and 'Name2' and concatenating other columns
        condensed_df = annotation_pairs_df.groupby(['Name1', 'Name2'], as_index=False).agg({
            'Term1': lambda x: '; '.join(sorted(set(x))),
            'Term2': lambda x: '; '.join(sorted(set(x))),
            'Source1': lambda x: '; '.join(sorted(set(x))),
            'Source2': lambda x: '; '.join(sorted(set(x))),
            'Gene': lambda x: '; '.join(sorted(set(x))),
            'Domain': lambda x: '; '.join(sorted(set(x))),
            'ReactionMatch': lambda x: 'Yes' if any(x == 'Yes') else 'No',
            'ECMatch': lambda x: 'Yes' if any(x == 'Yes') else 'No'
        })

        # Save the condensed DataFrame to a new CSV file
        condensed_df.to_csv(condensed_output_path, index=False)

        print(f"Condensed annotation pairs saved to: {condensed_output_path}")

    def print_annotation_comparison_table(self):
        annotations_by_gene = kbdevutil.load("annotations_by_gene")

        all_data = {
            "rhea":kbdevutil.load("rhea_data"),
            "sso":kbdevutil.load("sso_data"),
            "msrxn":kbdevutil.load("msrxn_data"),
            "ko":kbdevutil.load("ko_data"),
            "dec":kbdevutil.load("ec_data"),
            "pec":kbdevutil.load("ec_data"),
            "glm": kbdevutil.load("ec_data"),
            "domain":kbdevutil.load("domain_proteins")
        }

        translation = {
            "rhea":"Rhea",
            "sso":"RAST",
            "ko":"DramKO",
            "dec":"DramEC",
            "pec":"Prokka",
            "glm":"GLM4EC"
        }
        records = {"Gene":[],"Domain":[],"Rhea":[],"Rhea rxn list":[],"RAST":[],"Prokka":[],"DramKO":[],"DramEC":[],"GLM4EC":[],"RAST rxn":[],"Prokka rxn":[],"DramKO rxn":[],"DramEC rxn":[],"GLM4EC rxn":[],"RAST ec":[],"Prokka ec":[],"DramKO ec":[],"DramEC ec":[],"GLM4EC ec":[],"RAST rxn list":[],"Prokka rxn list":[],"DramKO rxn list":[],"DramEC rxn list":[],"GLM4EC rxn list":[]}
        found = {}
        notfound = {}
        for gene in annotations_by_gene:
            records["Gene"].append(gene)
            if gene in all_data["domain"]:
                records["Domain"].append(all_data["domain"][gene])
            else:
                records["Domain"].append("None")
            all_rhea_rxn = {}
            all_rhea_ec = {}
            if "rhea" in annotations_by_gene[gene]:
                for item in annotations_by_gene[gene]["rhea"]:
                    term = item["term"].split(":").pop()
                    if term in all_data["rhea"]:
                        found[term] = 1
                        if "msrxn" in all_data["rhea"][term]:
                            for rxn in all_data["rhea"][term]["msrxn"]:
                                if rxn in all_data["msrxn"] and not all_data["msrxn"][rxn]["filtered"] == "yes":
                                    all_rhea_rxn[rxn] = 1
                        if "ec" in all_data["rhea"][term]:
                            for ec in all_data["rhea"][term]["ec"]:
                                all_rhea_ec[ec] = 1
                        if "swiss_prot_EC" in annotations_by_gene[gene]:
                            for ec in annotations_by_gene[gene]["swiss_prot_EC"]:
                                all_rhea_ec[ec] = 1
                            if "ec" not in all_data["rhea"][term] or len(all_data["rhea"][term]["ec"]) == 0:
                                all_data["rhea"][term]["ec"] = annotations_by_gene[gene]["swiss_prot_EC"]
                        if "swiss_prot_name" in annotations_by_gene[gene] and ("name" not in all_data["rhea"][term] or all_data["rhea"][term]["name"] == None):
                            all_data["rhea"][term]["name"] = annotations_by_gene[gene]["swiss_prot_name"]
                    else:
                        notfound[term] = 1
                        if "swiss_prot_name" in annotations_by_gene[gene]:
                            all_data["rhea"][term] = {
                                "id":term,
                                "ec":[],
                                "name":annotations_by_gene[gene]["swiss_prot_name"],
                                "genes":[gene],
                                "msrxn":[]
                            }
                            if "swiss_prot_EC" in annotations_by_gene[gene]:
                                all_data["rhea"][term]["ec"] = annotations_by_gene[gene]["swiss_prot_EC"]
                        elif "swiss_prot_EC" in annotations_by_gene[gene]:
                            all_data["rhea"][term] = {
                                "id":term,
                                "ec":annotations_by_gene[gene]["swiss_prot_EC"],
                                "name":None,
                                "genes":[gene],
                                "msrxn":[]
                            }
            for source in translation:
                if source in annotations_by_gene[gene]:
                    name = ""
                    other_ec = {}
                    other_rxn = {}
                    rxnmatch = "No"
                    ecmatch = "No"
                    for item in annotations_by_gene[gene][source]:
                        if len(name) > 0:
                            name += "\n"
                        term = item["term"].split(":").pop()
                        name += term
                        if term in all_data[source]:
                            if all_data[source][term]["name"]:
                                name += ":"+all_data[source][term]["name"]
                            if (len(all_data[source][term]["ec"]) > 0):
                                name += "["+";".join(all_data[source][term]["ec"])+"]"
                                for ec in all_data[source][term]["ec"]:
                                    other_ec[ec] = 1
                                    if ec in all_rhea_ec:
                                        ecmatch = "Yes"
                            for rxn in all_data[source][term]["msrxn"]:
                                other_rxn[rxn] = 1
                                if rxn in all_rhea_rxn:
                                    rxnmatch = "Yes"
                    records[translation[source]].append(name)
                    filtered_list = []
                    for rxn in other_rxn:
                        if rxn in all_data["msrxn"] and not all_data["msrxn"][rxn]["filtered"] == "yes":
                            filtered_list.append(rxn)
                    records[translation[source]+" rxn list"].append(";".join(filtered_list))
                    if source != "rhea":
                        if len(all_rhea_rxn) == 0:
                            rxnmatch = "No Rhea rxn"
                        if len(all_rhea_ec) == 0:
                            ecmatch = "No Rhea ec"
                        if len(other_rxn) == 0:
                            rxnmatch = "No other rxn"
                        if len(other_ec) == 0:
                            ecmatch = "No other ec"
                        records[translation[source]+" rxn"].append(rxnmatch)
                        records[translation[source]+" ec"].append(ecmatch)
                else:
                    records[translation[source]].append("None")
                    records[translation[source]+" rxn list"].append("")
                    if source == "rhea":
                        records[translation[source]+" rxn"].append("No Rhea")
                        records[translation[source]+" ec"].append("No Rhea")
                    else:
                        records[translation[source]+" rxn"].append("No function")
                        records[translation[source]+" ec"].append("No function")
        df = pd.DataFrame.from_dict(records)
        df.to_csv(kbdevutil.output_dir+"/annotation_comparison.csv",index=False)
        print(str(len(found))+" terms found")
        print(str(len(notfound))+" terms not found")

ontologyutil = OntologyUtil()