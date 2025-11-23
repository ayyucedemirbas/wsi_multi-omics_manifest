import requests
import pandas as pd
from tqdm import tqdm

def gdc_request(endpoint, params):
    url = f"https://api.gdc.cancer.gov/{endpoint}"
    response = requests.get(url, params=params)
    response.raise_for_status()
    return response.json()

def get_files(experimental_strategy=None, data_category=None, project="TCGA-BRCA"):
    filters = {"op": "and", "content": []}
    if experimental_strategy:
        filters["content"].append({
            "op": "in",
            "content": {"field": "files.experimental_strategy", "value": [experimental_strategy]}
        })
    if data_category:
        filters["content"].append({
            "op": "in",
            "content": {"field": "files.data_category", "value": [data_category]}
        })
    filters["content"].append({
        "op": "in",
        "content": {"field": "cases.project.project_id", "value": [project]}
    })

    params = {
        "filters": str(filters).replace("'", '"'),
        "fields": "file_id,cases.case_id,cases.submitter_id",
        "format": "JSON",
        "size": "5000"
    }

    data = gdc_request("files", params)
    # Normalize JSON, extract 'cases' list
    df = pd.json_normalize(
        data["data"]["hits"],
        record_path=["cases"],
        meta=["file_id"],
        errors="ignore"
    )
    return df


wsi_df = get_files(experimental_strategy="Tissue Slide")
rna_df = get_files(experimental_strategy="RNA-Seq")
meth_df = get_files(data_category="DNA Methylation")
mut_df = get_files(data_category="Simple Nucleotide Variation")

def get_clinical(project="TCGA-BRCA"):
    params = {
        "filters": '{"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-BRCA"]}}',
        "fields": "case_id,submitter_id,disease_type,primary_diagnosis,ajcc_pathologic_stage,gender,age_at_diagnosis",
        "format": "JSON",
        "size": "5000"
    }
    data = gdc_request("cases", params)
    df = pd.json_normalize(data["data"]["hits"])
    return df

clinical_df = get_clinical()

wsi_ids = set(wsi_df["submitter_id"])
rna_ids = set(rna_df["submitter_id"])
meth_ids = set(meth_df["submitter_id"])
mut_ids = set(mut_df["submitter_id"])

full_patients = list(wsi_ids & rna_ids & meth_ids & mut_ids)
print(f"Total patients with all modalities: {len(full_patients)}")

manifest_rows = []

for patient in tqdm(full_patients, desc="Building manifest"):
    wsi_files = wsi_df[wsi_df["submitter_id"] == patient]["file_id"].tolist()
    rna_files = rna_df[rna_df["submitter_id"] == patient]["file_id"].tolist()
    meth_files = meth_df[meth_df["submitter_id"] == patient]["file_id"].tolist()
    mut_files = mut_df[mut_df["submitter_id"] == patient]["file_id"].tolist()
    
    clinical_info = clinical_df[clinical_df["submitter_id"] == patient]
    
    stage = clinical_info.get("ajcc_pathologic_stage", pd.Series([None])).values[0]
    diagnosis = clinical_info.get("primary_diagnosis", pd.Series([None])).values[0]
    age = clinical_info.get("age_at_diagnosis", pd.Series([None])).values[0]
    gender = clinical_info.get("gender", pd.Series([None])).values[0]

    manifest_rows.append({
        "patient_barcode": patient,
        "wsi_file_ids": ";".join(wsi_files),
        "rna_file_ids": ";".join(rna_files),
        "meth_file_ids": ";".join(meth_files),
        "mut_file_ids": ";".join(mut_files),
        "primary_diagnosis": diagnosis,
        "ajcc_pathologic_stage": stage,
        "age_at_diagnosis": age,
        "gender": gender
    })

manifest_df = pd.DataFrame(manifest_rows)

manifest_df.to_csv("TCGA_BRCA_Full_MultiOmics_WSI_Manifest.csv", index=False)
manifest_df.head()
