import time
import pandas as pd
from Bio import Entrez


def build_query(gene, tumor_region, tumor_type_synonyms):
    """
    Constructs a PubMed query for a given gene, tumor region, and optional tumor type synonyms.

    Args:
        gene (str): The gene symbol to include in the query.
        tumor_region (str): The selected tumor region, not directly used in the query.
        tumor_type_synonyms (list of str): A list of synonyms for the tumor type to refine the query.

    Returns:
        str: A formatted PubMed query string that includes the gene, 'cancer', and optionally
             synonyms for the tumor type in the title or abstract fields.
    """

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
    Searches PubMed for articles based on the specified query and date range.

    Parameters:
    query (str): The search query to be used in PubMed.
    recent_pubmedids_cutoff (int): Maximum number of recent PubMed IDs to retrieve.
    start_date (datetime.date): The minimum date for the search cutoff.
    end_date (datetime.date): The maximum date for the search cutoff.

    Returns:
    tuple: A tuple containing the count of articles found and a list of PubMed IDs. 
           Returns (-1, []) in case of an exception.
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


def validate_searched_pubmed_data(genes_list, tumor_region, tumor_type_synonyms, start_date, end_date, novel_cutoff, recent_pubmedids_cutoff, entrez_email, entrez_api_key, sleep_time=0.2):
    """
    Validates the results of a PubMed search for a given list of genes, 
    tumor region, and optional tumor type synonyms. This function takes into
    account the number of PubMed hits as well as the date range for the search.

    Parameters:
    genes_list (list of str): A list of gene symbols to search for.
    tumor_region (str): The tumor region to search for.
    tumor_type_synonyms (list of str): A list of synonyms for the tumor type.
    start_date (datetime.date): The minimum date for the search cutoff.
    end_date (datetime.date): The maximum date for the search cutoff.
    novel_cutoff (int): The number of PubMed hits required for a gene to be considered "known".
    recent_pubmedids_cutoff (int): The maximum number of recent PubMed IDs to retrieve.
    entrez_email (str): Email address for Entrez API authentication.
    entrez_api_key (str): API key for Entrez API authentication.
    sleep_time (float, optional): The time in seconds to wait between searches. Defaults to 0.2.

    Returns:
    pandas.DataFrame: A DataFrame with columns for the gene symbol, number of PubMed hits, status (known/novel), and recent PubMed IDs.
    """
    if entrez_email and entrez_api_key:
        Entrez.email = entrez_email
        Entrez.api_key = entrez_api_key
        sleep_time = 0.2
    else:
        sleep_time = 0.5
    
    results = []
    for gene in genes_list:
        query = build_query(gene, tumor_region, tumor_type_synonyms)
        count, ids = search_pubmed_data(query, recent_pubmedids_cutoff, start_date, end_date) 
        status = "known" if count >= novel_cutoff else "novel"
        id_string = ", ".join(ids)  # Combine IDs into a single comma-separated string
        results.append({"Gene": gene, "PubMedHits": count, "Status": status, "RecentPubMedIDs": id_string})
        time.sleep(sleep_time)

    return pd.DataFrame(results)