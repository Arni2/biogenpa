from flask import Flask, request, jsonify, send_file

import json
import re
import pymongo
from pymongo import MongoClient

#plot
import matplotlib.pyplot as plt
from collections import Counter
import io

import threading
import random
import string

from get_ref_fasta_value_on_position import get_sequence_around

client = MongoClient('mongodb://localhost:27017/')
db = client['vcf_database']
vcf_table = db['vcf_data']
dis_table = db['dis_data']

annotation_frequency_sessions = {}

app = Flask(__name__)


def validate_input(searched):
    if re.match(r"^[-ążźćęńłóa-zA-Z0-9\s]*$", searched):
        return searched
    raise ValueError("Invalid input")

def validate_input_integer(searchedRsPos):
    if re.match(r"^[0-9]*$", searchedRsPos):
        return searchedRsPos
    raise ValueError("Invalid input")

def mergeURLs(links):
    res = []
    for url in links:
        if url:
            if isinstance(url, str):
                res.append(url)
            elif url != []:
                res.append(', '.join(url))
    return res

def create_annotations(info):
    randomId = ''.join(random.choices(string.ascii_letters + string.digits, k=20))
    annotation_frequency_sessions[randomId] = info
    return randomId

def group_by_rs(session_id):
    data = annotation_frequency_sessions[session_id]
    print("Len: " + str(len(data)))
    output = {}
    counter_len = 0
    for item in data:
        if not output.get(item[1], None): #rsId
            counter_len += 1
            output[item[1]] = []
        output[item[1]].append(item)
    annotation_frequency_sessions[session_id] = output
    return output

def collect_allels(poses):
    try:
        print(len(poses))
        output = {}
        print(poses)
        with pymongo.timeout(100):
            query = {"POS": {"$in": poses}}
            projection = {
                "POS": 1,  # Include column1
                "frequency": 1,  # Include column2
                "ALT_VALUE": 1,  # Include column2
                "REF": 1,  # Include column2
                "CHROM": 1,  # Include column2
                "gene": 1,  # Include column2
                "RS_ID": 1,  # Include column2
                "oper_mutation": 1,  # Include column2
                "_id": 0       # Exclude the default _id field (optional)
            }
            x = dis_table.find(query, projection).limit(const_limit_data * 6)
            for document in x:
                if not output.get(document["POS"]):
                    output[document["POS"]] = []
                output[document["POS"]].append(document)
        return output
    except Exception as e:
        print(format(e))
        return "Error: {}".format(e)
    
def generate_report(searched, searchedRsPos, or_and):
    output = []
    try:
        with pymongo.timeout(100):
            searched = validate_input(searched.lower().replace("\n", " ")).split()
            searchedRsPos = validate_input_integer(searchedRsPos)
            rgxes = []
            for ser in searched:
                rgxes.append(re.compile('.*' + re.escape(ser) + '.*', re.IGNORECASE))
            fields_to_search = ["disease_ids", "disease_name", "gene", 
                        "mutation", "REF", "ALT_VALUE"]
            searchedBy = ""
            if searched:
                searchedBy = "szukano po tekscie"
                query = {
                "$or": [  # Ensure that at least one term is matched in each field
                        {
                            ("$or" if or_and == "or" else "$and"): [
                                {field: rgx}
                                for rgx in rgxes  # Apply OR within each field for the search terms
                            ]
                        }
                        for field in fields_to_search
                    ]
                }
            else:
                searchedBy = "szukano rs i pozycji"
                query = {
                    "$or": [
                        {"RS_ID": {"$regex": re.compile('.*' + searchedRsPos + '.*')}},
                        {"POS": {"$regex": re.compile('.*' + searchedRsPos + '.*')}}
                    ]
                }
            x = dis_table.find(query).limit(const_limit_data)
            for document in x:
                document.pop('_id', None)
                document["urls"] = mergeURLs([document.get("rs_link"), document.get("omim_links"), document.get("orphanet_links"),
                          document.get("medgen_links"), document.get("hpo_links"),
                          document.get("MONDO_links"), document.get("MeSH_links"), document.get("UniProt_links"), document.get("ClinVar_links"),
                          document.get("Ensembl_links"), document.get("LOVD_links"), document.get("GeneCards_links"), document.get("GnomAD_links")])
                document.pop('rs_link', None)
                document.pop('omim_links', None)
                document.pop('orphanet_links', None)
                document.pop('medgen_links', None)
                document.pop('hpo_links', None)
                document.pop('MONDO_links', None)
                document.pop('MeSH_links', None)
                document.pop('UniProt_links', None)
                document.pop('ClinVar_links', None)
                document.pop('Ensembl_links', None)
                document.pop('LOVD_links', None)
                document.pop('GeneCards_links', None)
                document.pop('GnomAD_links', None)
                output.append(document)
            session_id = create_annotations(list([value["frequency"], value["RS_ID"][0], 
                                             value["REF"], value["ALT_VALUE"],
                                             value["gene"], value["CHROM"], value["POS"]] for value in output))
            als = collect_allels(list(value["POS"] for value in output))
            output.append({'alleles': als})
            output.append({'session_id': session_id, "data": group_by_rs(session_id), "searchedBy": searchedBy})
    except Exception as e:
        print(format(e))
        return "Error: {}".format(e)
    return output

def process_request(searched, searchedRsPos, or_and):
    try:
        output = generate_report(searched, searchedRsPos, or_and)
        return output
    except Exception as e:
        print(format(e))
        return "Error: {}".format(e)

@app.route('/api/reference', methods=['POST'])
def allels_sequence():
    with lock:
        try:
            data = request.get_json()
            got = data.get('allels', None)
            if got:
                return jsonify({'success': True, 'report': { "sequence" : get_sequence_around(got[0], got[1]), "isVisible" : 1 } })
            else:
                return jsonify({'success': False, 'report': "-"})
        except Exception as e:
            return "Error: {}".format(e)

@app.route('/api/report', methods=['POST'])
def generate_report_for_fastq():
    with lock:
        try:
            data = request.get_json()
            searched = data.get('searched', '')
            searchedRsPos = data.get('rsOrPos', '')
            or_and = data.get('or_and', 'and')
            report = process_request(searched, searchedRsPos, or_and)
            return jsonify({'success': True, 'report': report})
        except Exception as e:
            return "Error: {}".format(e)

if __name__ == '__main__':
    try:
        print("Starting Flask server. Press Ctrl+C to stop and print results.")
        app.run(host='0.0.0.0', port=5000)
    except KeyboardInterrupt:
        print("\nCtrl+C detected! Here are the results:")
