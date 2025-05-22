from search import validate_searched_pubmed_data
import datetime
import streamlit as st

st.title("GeneLit")

# Create two columns: left narrow, right wider
left_col, right_col = st.columns([1, 3])

with st.form("Inputs"):  # Start the form

    # Left column: Narrator text area
    with left_col:
        # INPUT 1: Gene symbols list
        genesymbols = st.text_area("Paste your gene symbols here", height=400)
        genes_list = [gene.strip() for gene in genesymbols.split("\n") if gene.strip()]
        
    # Right column: controls
    with right_col:
        # INPUT 2: tumor region
        tumor_region_opts = ["lung", "breast", "colorectum", "prostate", "stomach", "liver", "cervix", "thyroid", "esophagus", "ovary", "pancreas", "bladder", "kidney", "lymph nodes", "bone marrow", "skin", "brain", "oral cavity", "uterus"]
        tumor_region_selected_opt = st.selectbox("Select tumor region", 
                                                 tumor_region_opts)
        
        # INPUT 3: tumor type synonyms
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
        
        # INPUT 4: start date
        start_date = st.date_input("Min date cutoff for search", 
                                   min_value=datetime.date(1900, 1, 1))
        
        # INPUT 5: end date
        end_date = st.date_input("Max date cutoff for search")
        
        # INPUT 6: novel cutoff
        novel_cutoff = st.number_input("PubMed hits cutoff for novel genes", 
                                       min_value=0, 
                                       max_value=20, 
                                       value=5)
        
        # INPUT 7: novel cutoff
        recent_pubmedids_cutoff = st.number_input("How many recent pubmed ids to show", 
                                                  min_value=0, 
                                                  max_value=100, 
                                                  value=10)

    # Submit button
    submitted = st.form_submit_button("Submit")

    if submitted:
        # Process all inputs after submit
        # st.write("Gene Symbols:", genes_list)
        st.write("Tumor Region:", tumor_region_selected_opt)
        st.write("Selected Synonyms:", tumor_region_synonyms)
        st.write("Start Date:", start_date)
        st.write("End Date:", end_date)
        st.write("PubMed hits cutoff for novel genes:", novel_cutoff)
        st.write("How many recent pubmed ids to show:", recent_pubmedids_cutoff)
        
        results_df = validate_searched_pubmed_data(genes_list = genes_list, 
                                                   tumor_region = tumor_region_selected_opt, 
                                                   tumor_type_synonyms = tumor_region_synonyms, 
                                                   start_date = start_date, 
                                                   end_date = end_date, 
                                                   novel_cutoff = novel_cutoff, 
                                                   recent_pubmedids_cutoff = recent_pubmedids_cutoff)
        st.write(results_df)