import time
import pandas as pd
from Bio import Entrez

# Set your email and optional API key here
Entrez.email = "subhajitbn.maths@gmail.com"  # Replace with your email
Entrez.api_key = "33aaabf511eb2a854b32fd07aa18553a3b09"     # Optional but recommended

def build_query(gene, tumor_region, tumor_type_synonyms):
    # Always include cancer
    if tumor_type_synonyms:
        # term_query = " OR ".join(terms)
        term_query = " OR ".join(f'{t}[Title/Abstract]' for t in tumor_type_synonyms)
        return f'"{gene}"[Title/Abstract] AND cancer AND ({term_query})'
    else:
        return f'"{gene}"[Title/Abstract] AND cancer'


def search_pubmed_data(query, 
                       recent_pubmedids_cutoff, 
                       start_date,
                       end_date):
    """
    Search PubMed using the given query and return the total hit count
    and the article IDs of the most recent hits (up to max_hits).
    """
    try:
        handle = Entrez.esearch(db = "pubmed", 
                                term = query, 
                                retmax = recent_pubmedids_cutoff, 
                                mindate = start_date.strftime("%Y/%m/%d"),
                                maxdate = end_date.strftime("%Y/%m/%d"),
                                sort = "most_recent")
        record = Entrez.read(handle)
        count = int(record["Count"])
        ids = record.get("IdList", [])  # Get the list of PubMed IDs
        return count, ids
    except Exception as e:
        print(f"Error during search: {e}")
        return -1, []


def validate_searched_pubmed_data(genes_list, tumor_region, tumor_type_synonyms, start_date, end_date, novel_cutoff, recent_pubmedids_cutoff, sleep_time=0.2):
    """
    Validate genes in the CSV files under csv_dir using PubMed searches,
    and save the results to output_dir.
    """
    results = []
    for gene in genes_list:
        query = build_query(gene, tumor_region, tumor_type_synonyms)
        count, ids = search_pubmed_data(query, recent_pubmedids_cutoff, start_date, end_date) 
        status = "known" if count >= novel_cutoff else "novel"
        id_string = ", ".join(ids)  # Combine IDs into a single comma-separated string
        results.append({"Gene": gene, "PubMedHits": count, "Status": status, "RecentPubMedIDs": id_string})
        time.sleep(sleep_time)

    return pd.DataFrame(results)