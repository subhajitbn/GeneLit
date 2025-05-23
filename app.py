from search import validate_searched_pubmed_data
import datetime
import streamlit as st

st.title("GeneLit")
st.subheader("PubMed literature status for your list of cancer-associated genes")
with st.expander("ℹ️ How does the PubMed search work?"):
    st.markdown(
        """
        We build a query for searching PubMed, using your inputs and the following criteria:

        - ✅ We look for articles that mention the **gene symbol** and **at least one cancer-type synonym** in the **title or abstract**.
        - ✅ The article must also mention the word **"cancer"** **anywhere** in the content.
        - 📅 Only articles published within your **selected date range** are considered.

        Simple, fast, and tailored for your literature review needs.
        """
    )

st.header("Inputs")

# INPUT 1: Do you have an email and api key for the PubMed search?
have_email_and_api_key = st.checkbox("Do you have an email and api key for the PubMed search? This is optional, but provides faster results.")

if have_email_and_api_key:
    entrez_email = st.text_input("Enter your email", "")
    entrez_api_key = st.text_input("Enter your API key", "")
else:
    entrez_email = None
    entrez_api_key = None

# Create two columns: left narrow, right wider
left_col, right_col = st.columns([1, 3])

with st.form("Inputs"):  # Start the form

    with left_col:
        # INPUT 2: Gene symbols list
        genesymbols = st.text_area("Paste your gene symbols here", height=400, placeholder="One gene symbol per line")
        genes_list = [gene.strip() for gene in genesymbols.split("\n") if gene.strip()]
        
    with right_col:
        # INPUT 3: tumor region
        tumor_region_opts = ["lung", "breast", "colorectum", "prostate", "stomach", "liver", "cervix", "thyroid", "esophagus", "ovary", "pancreas", "bladder", "kidney", "lymph nodes", "bone marrow", "skin", "brain", "oral cavity", "uterus"]
        tumor_region_selected_opt = st.selectbox("Select tumor region", 
                                                 tumor_region_opts)
        
        # INPUT 4: tumor type synonyms
        tumor_region_synonyms_dict = {
            "lung": ["lung", "pulmonary", "pulmonic", "pulmonology", "bronchial", "bronchus", "respiratory tract"],
            "breast": ["breast", "mammary", "mammary gland", "mammary tissue", "breast tissue"],
            "colorectum": ["colorectum", "colorectal", "colon", "colonic", "large intestine", "large bowel", "rectal", "rectum", "rectal tissue", "anus"],
            "prostate": ["prostate", "prostatic", "prostate gland"],
            "stomach": ["stomach", "gastric", "gastrium", "gastric mucosa", "gastrointestinal tract"],
            "liver": ["liver", "hepatic", "hepatic tissue", "hepatic parenchyma"],
            "cervix": ["cervix", "cervical", "cervical canal", "uterine cervix"],
            "thyroid": ["thyroid", "thyroid gland", "thyroidal", "thyroid tissue"],
            "esophagus": ["esophagus", "oesophagus", "esophageal", "oesophageal", "esophageal mucosa", "gastroesophageal junction"],
            "ovary": ["ovary", "ovarian", "ovarian tissue", "female gonad"],
            "pancreas": ["pancreas", "pancreatic", "pancreatic tissue"],
            "bladder": ["bladder", "urinary bladder", "bladder wall", "vesical"],
            "kidney": ["kidney", "renal", "renal tissue", "renal cortex", "renal medulla"],
            "lymph nodes": ["lymph nodes", "lymphatic nodes", "nodal", "lymphoid tissue", "lymphatic system"],
            "bone marrow": ["bone marrow", "marrow", "medullary cavity", "hematopoietic tissue"],
            "skin": ["skin", "cutaneous", "dermal", "epidermis", "integumentary system"],
            "brain": ["brain", "cerebral", "cerebrum", "central nervous system", "CNS", "neural tissue", "encephalon"],
            "oral cavity": ["oral cavity", "mouth", "oral", "buccal", "oropharynx", "tongue", "floor of mouth"],
            "uterus": ["uterus", "uterine", "womb", "endometrium", "myometrium"]
        }
        tumor_region_synonyms_opts = tumor_region_synonyms_dict[tumor_region_selected_opt]
        tumor_region_synonyms = st.multiselect("Select multiple synonyms for the tumor region", 
                                       tumor_region_synonyms_opts)
        
        # INPUT 5: start date
        start_date = st.date_input("Min date cutoff for search", 
                                   min_value=datetime.date(1900, 1, 1))
        
        # INPUT 6: end date
        end_date = st.date_input("Max date cutoff for search")
        
        # INPUT 7: novel cutoff
        novel_cutoff = st.number_input("PubMed hits cutoff for novel genes", 
                                       min_value=0, 
                                       max_value=20, 
                                       value=5)
        
        # INPUT 8: novel cutoff
        recent_pubmedids_cutoff = st.number_input("How many recent pubmed ids to show", 
                                                  min_value=0, 
                                                  max_value=100, 
                                                  value=10)

    # Submit button
    submitted = st.form_submit_button("Submit")

if submitted:
    # Process all inputs after submit
    if entrez_email and entrez_api_key:
        st.write("Email:", entrez_email)
        st.write("API Key:", entrez_api_key)
    st.write("Tumor Region:", tumor_region_selected_opt)
    st.write("Selected Synonyms:", tumor_region_synonyms)
    st.write("Start Date:", start_date)
    st.write("End Date:", end_date)
    st.write("PubMed hits cutoff for novel genes:", novel_cutoff)
    st.write("How many recent pubmed ids to show:", recent_pubmedids_cutoff)
    
    st.header("Results")
    
    with st.status("Searching PubMed...", expanded=True) as status:
        results_df = validate_searched_pubmed_data(genes_list = genes_list, 
                                                    tumor_region = tumor_region_selected_opt, 
                                                    tumor_type_synonyms = tumor_region_synonyms, 
                                                    start_date = start_date, 
                                                    end_date = end_date, 
                                                    novel_cutoff = novel_cutoff, 
                                                    recent_pubmedids_cutoff = recent_pubmedids_cutoff, 
                                                    entrez_email = entrez_email, 
                                                    entrez_api_key = entrez_api_key)
        status.update(label="Done!", state="complete")
        st.write(results_df)