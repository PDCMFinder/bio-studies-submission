import os
import requests
import urllib.parse
import json
import argparse
import re

# Define the endpoint and query parameters
COLUMNS_TO_READ = [
    "external_model_id",
    "provider_name",
    "model_type",
    "histology",
    "cancer_system",
    "cancer_grade",
    "cancer_stage",
    "primary_site",
    "collection_site",
    "tumour_type",
    "patient_age",
    "patient_sex",
    "patient_ethnicity",
    "xenograft_model_specimens",
    "quality_assurance",
    "email_list",
    "data_source",
]

BASE_URL = "https://dev.cancermodels.org/api/"
SEARCH_INDEX_ENDPOINT = BASE_URL + "search_index"
MODEL_MOLECULAR_METADATA_ENDPOINT = BASE_URL + "model_molecular_metadata"
RELEASE_ENDPOINT = BASE_URL + "release_info"
PARAMS = {
    "select": ",".join(COLUMNS_TO_READ),
    "limit": 100,  # Set a limit for pagination
    "offset": 0,
}
STUDY_SECTION_ID = "s1"

# Relationship between molecular data types and the name of the tables where data is
MOLECULAR_DATA_TABLES = {
    "copy number alteration": "cna_data_table",
    "bio markers": "biomarker_data_table",
    "expression": "expression_data_table",
    "mutation": "mutation_data_table",
}

MOLECULAR_DATA_FIELDS = {
    "copy number alteration": [
        "hgnc_symbol",
        "cromosome",
        "strand",
        "log10r_cna",
        "log2r_cna",
        "seq_end_position",
        "copy_number_status",
        "gistic_value",
        "picnic_value",
    ],
    "bio markers": ["biomarker", "result"],
    "expression": [
        "hgnc_symbol",
        "rnaseq_coverage",
        "rnaseq_fpkm",
        "rnaseq_tpm",
        "rnaseq_count",
        "affy_hgea_probe_id",
        "affy_hgea_expression_value",
        "illumina_hgea_probe_id",
        "illumina_hgea_expression_value",
        "z_score",
    ],
    "mutation": [
        "hgnc_symbol",
        "amino_acid_change",
        "cromosome",
        "strand",
        "consequence",
        "read_depth",
        "allele_frequency",
        "seq_start_position",
        "ref_allele",
        "alt_allele",
        "biotype",
    ],
}


def create_folder_if_not_exists(output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


# Fetch and process data with pagination
def fetch_and_process_data(output):
    processed_data = []
    total_count = 0
    while True:
        quoted_value = urllib.parse.quote('in.("HBCx-118")', safe="(),")
        params_str = urllib.parse.urlencode(PARAMS)
        filter_str = f"external_model_id={quoted_value}"
        final_url = f"{SEARCH_INDEX_ENDPOINT}?{params_str}&{filter_str}"

        # print("final_url", final_url)

        response = requests.get(final_url)
        response.raise_for_status()
        data = response.json()

        # Process each model in the current batch
        for model in data:
            model_folder_path = (
                output + "/" + model["external_model_id"] + "_" + model["data_source"]
            )
            create_folder_if_not_exists(model_folder_path)
            study = format_model(model, model_folder_path)

            print(json.dumps(study))

        # Update metadata
        total_count += len(data)

        # Check if there are more results to fetch
        if len(data) < PARAMS["limit"]:
            break
        else:
            PARAMS["offset"] += PARAMS["limit"]

    return processed_data


def format_model(model, model_folder_path):
    study = {}
    study["accno"] = model["external_model_id"] + "_" + model["data_source"]
    # Should it br submission or another type?
    study["type"] = "submission"

    section = create_study_section(model, model_folder_path)
    study["section"] = section

    return study


def create_study_section(model, model_folder_path):
    study_section = {}
    study_section["accno"] = STUDY_SECTION_ID
    study_section["type"] = "Study"
    attributes = []
    attributes.append({"name": "Title", "value": model["histology"]})
    attributes.append({"name": "Study type", "value": model["model_type"]})
    study_section["attributes"] = attributes
    study_section = add_subsections(study_section, model, model_folder_path)
    return study_section


def add_subsections(study_section, model, model_folder_path):
    subsections = []
    organization_section = create_organization_section(model)
    subsections.append(organization_section)
    patient_tumor_subsection = create_patient_tumor_subsection(model)
    subsections.append(patient_tumor_subsection)
    if model["model_type"] == "PDX":
        pdx_model_engraftment_subsection = create_pdx_model_engraftment_subsection(
            model
        )
        subsections.append(pdx_model_engraftment_subsection)
    quality_control_subsection = create_model_quality_control_subsection(model)
    molecular_data_subsection = create_molecular_data_subsection(
        model, model_folder_path
    )
    subsections.append(molecular_data_subsection)

    subsections.append(quality_control_subsection)
    study_section["subsections"] = subsections
    return study_section


def create_organization_section(model):
    organization_section = {"accno": "o1", "type": "Organization"}
    organization_section["attributes"] = [
        {"name": "Name", "value": model["provider_name"]}
    ]
    return organization_section


def create_patient_tumor_subsection(model):
    subsection = {}
    attributes = []
    subsection["type"] = "Patient / Tumour Metadata"
    attributes.append({"name": "Patient Sex", "value": model["patient_sex"]})
    attributes.append({"name": "Patient Age", "value": model["patient_age"]})
    attributes.append(
        {"name": "Patient Ethnicity", "value": model["patient_ethnicity"]}
    )
    attributes.append({"name": "Tumour Type", "value": model["tumour_type"]})
    attributes.append({"name": "Cancer System", "value": model["cancer_system"]})
    attributes.append({"name": "Cancer Grade", "value": model["cancer_grade"]})
    attributes.append({"name": "Cancer Stage", "value": model["cancer_stage"]})
    attributes.append({"name": "Primary Site", "value": model["primary_site"]})
    attributes.append({"name": "Collection Site", "value": model["collection_site"]})

    subsection["attributes"] = attributes
    return subsection


def create_pdx_model_engraftment_subsection(model):
    rows = []
    print("model", model)
    xenograft_model_specimens = model["xenograft_model_specimens"]
    for row in xenograft_model_specimens:
        print("row", row)
        subsection_row = {}
        attributes = []

        subsection_row["type"] = "PDX model engraftment"
        attributes.append(
            {"name": "Host Strain Name", "value": row["host_strain_name"]}
        )
        attributes.append(
            {
                "name": "Host Strain Nomenclature",
                "value": row["host_strain_nomenclature"],
            }
        )

        attributes.append({"name": "Site", "value": row["engraftment_site"]})
        attributes.append({"name": "Type", "value": row["engraftment_type"]})
        attributes.append({"name": "Material", "value": row["engraftment_sample_type"]})
        attributes.append(
            {"name": "Material Status", "value": row["engraftment_sample_state"]}
        )
        attributes.append({"name": "Passage", "value": row["passage_number"]})

        subsection_row["attributes"] = attributes
        rows.append(subsection_row)

    return rows


def create_model_quality_control_subsection(model):
    rows = []

    quality_assurance = model["quality_assurance"]
    for row in quality_assurance:
        subsection_row = {}
        attributes = []

        subsection_row["type"] = "Model quality control"
        attributes.append({"name": "Technique", "value": row["validation_technique"]})
        attributes.append(
            {
                "name": "Description",
                "value": row["description"],
            }
        )

        attributes.append({"name": "Passage", "value": row["passages_tested"]})

        subsection_row["attributes"] = attributes
        rows.append(subsection_row)

    return rows


def create_molecular_data_subsection(model, model_folder_path):
    # Data will be represented as a section with files
    molecular_data_subsection = {"type": "Molecular data"}
    files = []

    url = (
        f"{MODEL_MOLECULAR_METADATA_ENDPOINT}?model_id=eq.{model['external_model_id']}"
    )

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    # Process each model in the current batch
    for model_molecular_metadata_row in data:
        sample_id = model_molecular_metadata_row["sample_id"]
        data_type = model_molecular_metadata_row["data_type"]
        platform_name = model_molecular_metadata_row["platform_name"]
        molecular_characterization_id = model_molecular_metadata_row[
            "molecular_characterization_id"
        ]

        file_path, size = fetch_molecular_data(
            model["external_model_id"],
            sample_id,
            data_type,
            platform_name,
            molecular_characterization_id,
            model_folder_path,
        )

        if size == 0:
            continue

        file = {"path": file_path.split("/", 1)[1]}
        file["type"] = "file"
        file["size"] = size

        source = model_molecular_metadata_row["source"]
        sample_type = ""
        if source == "xenograft":
            sample_type = "Engrafted Tumour"
        elif source == "patient":
            sample_type = "Patient Tumour"
        else:
            sample_type = "Tumour Cells"
        attributes = []
        if model_molecular_metadata_row["data_restricted"] == "TRUE":
            attributes.append({"name": "Request e-mail", "value": model["email_list"]})
        else:
            attributes.append({"name": "Request e-mail"})
        attributes.append({"name": "Sample ID", "value": sample_id})
        attributes.append({"name": "Sample Type", "value": sample_type})
        attributes.append(
            {
                "name": "Engrafted Tumour Passage",
                "value": model_molecular_metadata_row["xenograft_passage"],
            }
        )
        attributes.append(
            {"name": "Data Type", "value": model_molecular_metadata_row["data_type"]}
        )
        attributes.append(
            {
                "name": "Platform Used",
                "value": model_molecular_metadata_row["platform_name"],
            }
        )
        attributes.append(
            {
                "name": "Raw Data",
                "value": model_molecular_metadata_row["external_db_links"],
            }
        )
        file["attributes"] = attributes

        files.append(file)

    molecular_data_subsection["files"] = files

    return molecular_data_subsection


def clean_file_name(original_name):
    clean = re.sub(r"[\s/\\?%*:|\"<>\x7F\x00-\x1F]", "-", original_name)
    return clean


def fetch_molecular_data(
    model_name,
    sample_id,
    data_type,
    platform_name,
    molecular_characterization_id,
    model_folder_path,
):
    table_name = MOLECULAR_DATA_TABLES[data_type]

    url = f"{BASE_URL}/{table_name}?molecular_characterization_id=eq.{molecular_characterization_id}"

    file_name = clean_file_name(
        sample_id + "_" + data_type + "_" + platform_name + ".tsv"
    )

    file_path = model_folder_path + "/" + file_name

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    if data == []:
        return file_path, 0

    fields = MOLECULAR_DATA_FIELDS[data_type]
    write_molecular_data_file(data, file_path, fields)

    file_size = os.path.getsize(file_path)
    print("file size :", file_size, "bytes")
    return file_path, file_size


def write_molecular_data_file(data, file_path, fields):
    lines = []
    num_lines = 0
    # First line corresponds to the headers of the file
    lines.append("\t".join(fields))
    for row in data:
        values = []
        for field in fields:
            value = row[field]
            values.append(value)
            num_lines += 1
        lines.append("\t".join(values))

    # Create `file_path` with the content of `json` which is expected to be cna data
    with open(file_path, "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")


# Main function
def main(output):
    create_folder_if_not_exists(output)
    processed_data = fetch_and_process_data(output)

    print("done", processed_data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch, process, and formats CancerModels.org data into a structure that fits the BioStudies structurepython."
    )
    parser.add_argument(
        "--output", required=True, help="Folder where the data will be downloaded"
    )

    args = parser.parse_args()

    main(args.output)
