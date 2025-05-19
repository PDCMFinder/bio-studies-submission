from datetime import date
import os
import requests
import urllib.parse
import json
import argparse
import re

# Define the endpoint and query parameters
COLUMNS_TO_READ = [
    "pdcm_model_id",
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
    "pdx_model_publications",
    "model_relationships",
    "data_source",
]

BASE_URL = "https://dev.cancermodels.org/api/"
SEARCH_INDEX_ENDPOINT = BASE_URL + "search_index"
MODEL_MOLECULAR_METADATA_ENDPOINT = BASE_URL + "model_molecular_metadata"
DOSING_STUDIES_ENDPOINT = BASE_URL + "dosing_studies"
PATIENT_TREATMENT_ENDPOINT = BASE_URL + "patient_treatment"
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
        "chromosome",
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
        "chromosome",
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
        # TM00199
        quoted_value = urllib.parse.quote('in.("TM00199")', safe="(),")
        params_str = urllib.parse.urlencode(PARAMS)
        filter_str = f"external_model_id={quoted_value}"
        final_url = f"{SEARCH_INDEX_ENDPOINT}?{params_str}&{filter_str}"

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
    study["type"] = "submission"

    attributes = []
    attributes.append({"name": "Title", "value": model["histology"]})
    attributes.append({"name": "ReleaseDate", "value": str(date.today())})

    study["attributes"] = attributes

    section = create_study_section(model, model_folder_path)
    study["section"] = section

    return study


def create_study_section(model, model_folder_path):
    study_section = {}
    study_section["accno"] = STUDY_SECTION_ID
    study_section["type"] = "Study"

    attributes = []
    attributes.append({"name": "Model ID", "value": model["external_model_id"]})
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
    molecular_data_subsection, molecular_data_files_section = (
        create_molecular_data_subsection(model, model_folder_path)
    )
    model_treatment_subsection = create_model_treatment_subsection(model)
    patient_treatment_subsection = create_patient_treatment_subsection(model)
    publication_subsection = create_publication_subsection(model)
    related_models_subsection = create_related_models_subsection(model)
    subsections.append(quality_control_subsection)
    subsections.append(molecular_data_subsection)
    subsections.append(molecular_data_files_section)
    subsections.append(model_treatment_subsection)
    subsections.append(patient_treatment_subsection)
    subsections.append(publication_subsection)
    subsections.append(related_models_subsection)

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
    xenograft_model_specimens = model["xenograft_model_specimens"]
    for row in xenograft_model_specimens:
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


def fetch_molecular_metadata_per_model(external_model_id):
    url = f"{MODEL_MOLECULAR_METADATA_ENDPOINT}?model_id=eq.{external_model_id}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    return data


def create_molecular_data_subsection(model, model_folder_path):
    # The molecular data subsection has 2 main components:
    # 1) A table showing the molecular metadata attributes
    # 2) A section with the files containing the actual molecular data

    # Note: The files must be located as the last section. Otherwise they don't render

    molecular_metadata_sections = []
    data = fetch_molecular_metadata_per_model(model["external_model_id"])

    # molecular_data_subsection = {"type": "Molecular data"}
    files = []

    # Process each row of molecular metadata
    for model_molecular_metadata_row in data:
        # Create row for the table
        row = create_molecular_metadata_row(model_molecular_metadata_row)
        molecular_metadata_sections.append(row)

        molecular_data_file = create_molecular_data_file(
            model_molecular_metadata_row, model_folder_path
        )

        if molecular_data_file:
            files.append(molecular_data_file)

    # Add a section at the end with the files
    files_section = {}
    files_section["files"] = files

    return molecular_metadata_sections, files_section


def get_sample_type_by_source(sample_source):
    sample_type = ""

    if sample_source == "xenograft":
        sample_type = "Engrafted Tumour"
    elif sample_source == "patient":
        sample_type = "Patient Tumour"
    else:
        sample_type = "Tumour Cells"
    return sample_type


def create_molecular_metadata_row(model_molecular_metadata_row):
    molecular_metadata_row = {"type": "Molecular data"}
    attributes = []

    attributes.append(
        {"name": "Sample ID", "value": model_molecular_metadata_row["sample_id"]}
    )
    sample_type = get_sample_type_by_source(model_molecular_metadata_row["source"])
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
    external_links = model_molecular_metadata_row["external_db_links"]

    eva_raw_data_column = create_url_attribute(external_links, "ENA")
    ega_raw_data_column = create_url_attribute(external_links, "EGA")
    geo_raw_data_column = create_url_attribute(external_links, "GEO")
    dbGap_raw_data_column = create_url_attribute(external_links, "dbGAP")

    attributes.append(eva_raw_data_column)
    attributes.append(ega_raw_data_column)
    attributes.append(geo_raw_data_column)
    attributes.append(dbGap_raw_data_column)

    molecular_metadata_row["attributes"] = attributes
    return molecular_metadata_row


def create_url_attribute(external_links, name):
    link_obj = get_link_by_resource(external_links, name)
    url = None
    if link_obj:
        url = link_obj["link"]
    attribute = {"name": name}
    if url:
        attribute["value"] = name
        attribute["valqual"] = [
            {
                "name": "url",
                "value": url,
            }
        ]
    return attribute


def create_molecular_data_file(
    model_molecular_metadata_row,
    model_folder_path,
):
    sample_id = model_molecular_metadata_row["sample_id"]
    data_type = model_molecular_metadata_row["data_type"]
    platform_name = model_molecular_metadata_row["platform_name"]
    molecular_characterization_id = model_molecular_metadata_row[
        "molecular_characterization_id"
    ]
    print("data_type==>", data_type)

    model_molecular_metadata_row["sample_id"]
    # Creates the file with the molecular data
    file_path, size = fetch_molecular_data(
        sample_id,
        data_type,
        platform_name,
        molecular_characterization_id,
        model_folder_path,
    )

    if size == 0:
        return None

    file = {"path": file_path.split("/", 1)[1]}
    file["type"] = "file"
    file["size"] = size

    attributes = []

    attributes.append({"name": "Sample ID", "value": sample_id})
    sample_type = get_sample_type_by_source(model_molecular_metadata_row["source"])
    attributes.append({"name": "Sample Type", "value": sample_type})
    attributes.append(
        {
            "name": "Engrafted Tumour Passage",
            "value": model_molecular_metadata_row["xenograft_passage"],
        }
    )
    attributes.append({"name": "Data Type", "value": data_type})
    attributes.append(
        {
            "name": "Platform Used",
            "value": platform_name,
        }
    )
    file["attributes"] = attributes

    return file


def create_immune_markers_subsection(model):
    immune_markers_subsection = {"type": "Immune markers"}
    attributes = []
    url = f"{BASE_URL}/immunemarker_data_extended?model_id=eq.{model['external_model_id']}"

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    for row in data:
        additional_data = row["essential_or_additional_details"]
        marker_value = row["marker_value"]
        if additional_data:
            marker_value = f"{marker_value} ({additional_data})"
        attributes.append({"name": row["marker_name"], "value": marker_value})
    immune_markers_subsection["attributes"] = attributes


def create_model_treatment_subsection(model):
    rows = []

    url = f"{DOSING_STUDIES_ENDPOINT}?model_id=eq.{model['pdcm_model_id']}"

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    # Process each model in the current batch
    for model_treatment_metadata_row in data:
        subsection_row = {"type": "Model treatment"}
        attributes = []
        treatment_response = model_treatment_metadata_row["response"]
        processed_treatment_data_entries = process_treatment_data_entries(
            model_treatment_metadata_row["entries"]
        )
        ...
        treatment_names = processed_treatment_data_entries["names"]
        doses = processed_treatment_data_entries["doses"]
        links = processed_treatment_data_entries["links"]
        passage_range = model_treatment_metadata_row["passage_range"]

        attributes.append({"name": "Drug", "value": treatment_names})
        attributes.append({"name": "Dose", "value": doses})
        attributes.append({"name": "Passage", "value": passage_range})
        attributes.append({"name": "Response", "value": treatment_response})
        chembl_attribute = create_url_attribute(links, "ChEMBL")
        pubchem_attribute = create_url_attribute(links, "PubChem")

        attributes.append(chembl_attribute)
        attributes.append(pubchem_attribute)

        subsection_row["attributes"] = attributes

        rows.append(subsection_row)

    return rows


def create_patient_treatment_subsection(model):
    rows = []

    url = f"{PATIENT_TREATMENT_ENDPOINT}?model_id=eq.{model['pdcm_model_id']}"

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    # Process each model in the current batch
    for patient_treatment_metadata_row in data:
        subsection_row = {"type": "Patient treatment"}
        attributes = []
        treatment_response = patient_treatment_metadata_row["response"]
        processed_treatment_data_entries = process_treatment_data_entries(
            patient_treatment_metadata_row["entries"]
        )
        treatment_names = processed_treatment_data_entries["names"]
        doses = processed_treatment_data_entries["doses"]
        links = processed_treatment_data_entries["links"]

        attributes.append({"name": "Treatment", "value": treatment_names})

        attributes.append({"name": "Dose", "value": doses})

        attributes.append({"name": "Response", "value": treatment_response})
        chembl_attribute = create_url_attribute(links, "ChEMBL")
        pubchem_attribute = create_url_attribute(links, "PubChem")

        attributes.append(chembl_attribute)
        attributes.append(pubchem_attribute)

        subsection_row["attributes"] = attributes

        rows.append(subsection_row)

    return rows


def create_publication_subsection(model):
    publications_ids = model["pdx_model_publications"]
    publications_ids = publications_ids.replace(" ", "").split(",")
    publications = [get_publication_data(val) for val in publications_ids]

    rows = []

    for publication in publications:
        row = {"type": "Publications"}
        atributes = []
        pmid = publication["pmid"]

        atributes.append({"name": "Title", "value": publication["title"]})
        atributes.append({"name": "Authors", "value": publication["authorString"]})
        atributes.append({"name": "Year", "value": publication["pubYear"]})
        atributes.append({"name": "Volume", "value": publication["journalVolume"]})
        atributes.append({"name": "Issue", "value": publication["issue"]})
        atributes.append({"name": "Issn", "value": publication["journalIssn"]})

        atributes.append(
            {
                "name": "EuropePMC",
                "value": "EuropePMC",
                "valqual": [
                    {
                        "name": "url",
                        "value": f"https://europepmc.org/article/MED/{pmid}",
                    }
                ],
            }
        )

        doi = publication["doi"]
        atributes.append(
            {
                "name": "DOI",
                "value": f"DOI:{doi}",
                "valqual": [
                    {
                        "name": "url",
                        "value": f"https://dx.doi.org/{doi}",
                    }
                ],
            }
        )

        atributes.append(
            {
                "name": "PubMed",
                "value": "PubMed",
                "valqual": [
                    {
                        "name": "url",
                        "value": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}",
                    }
                ],
            }
        )

        row["attributes"] = atributes
        rows.append(row)
    return rows


def create_related_models_subsection(model):
    related_models_subsection = {"type": "Related models"}
    model_relationships = model["model_relationships"]
    links = []
    parent = model_relationships["parents"]
    if parent and parent != "null":
        parent = parent[0]["external_model_id"]
        link = {"url": parent}
        attributes = [
            {"name": "Type", "value": "BioStudies"},
            {"name": "Role", "value": "child of"},
        ]
        link["attributes"] = attributes
        links.append(link)
    children = model_relationships["children"]
    if children and children != []:
        for child in children:
            child = child["external_model_id"]
            link = {"url": child}
            attributes = [
                {"name": "Type", "value": "BioStudies"},
                {"name": "Role", "value": "parent of"},
            ]
            link["attributes"] = attributes
            links.append(link)

    related_models_subsection["links"] = links
    return related_models_subsection


def format_link(link):
    biostudies_link = {}
    attributes = []
    biostudies_link["url"] = link["url"]
    attributes.append({"name": "Description", "value": link["label"]})
    attributes.append({"name": "Type", "value": link["resource"].lower()})
    biostudies_link["attributes"] = attributes
    return biostudies_link


def get_link_by_resource(links, resource):
    if links:
        for link in links:
            if link["resource"] == resource:
                return link
    return None


def clean_file_name(original_name):
    clean = re.sub(r"[\s/\\?%*:|\"<>\x7F\x00-\x1F]", "-", original_name)
    return clean


def fetch_molecular_data(
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
    return file_path, file_size


def process_treatment_data_entries(entries: list):
    names = []
    doses = []
    links = []

    for entry in entries:
        names.append(entry["name"])
        doses.append(entry["dose"])
        if entry["external_db_links"]:
            for external_db_link in entry["external_db_links"]:
                link = {
                    "label": entry["name"],
                    "link": external_db_link["link"],
                    "resource": external_db_link["resource_label"],
                }
                links.append(link)
    return {"names": " + ".join(names), "doses": " + ".join(set(doses)), "links": links}


def write_molecular_data_file(data, file_path, fields):
    print(f"write_molecular_data_file:[BIG][{file_path}][{fields}]")
    lines = []
    num_lines = 0
    # First line corresponds to the headers of the file
    lines.append("\t".join(fields))
    for row in data:
        values = []
        for field in fields:
            value = row[field]
            if value is None:
                values.append("")
            else:
                values.append(str(value))
            num_lines += 1
        lines.append("\t".join(values))

    # Create `file_path` with the content of `json` which is expected to be cna data
    with open(file_path, "w") as f:
        for line in lines:
            f.write(line)
            f.write("\n")


def get_publication_data(pub_id):
    if pub_id == "":
        return None
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/article/MED/{pub_id.replace('PMID:', '')}?resultType=lite&format=json"
    print("URL", url)

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    result = data["result"]

    # atributes.append({"name": "PMID", "value": publication["pmid"]})
    #     atributes.append({"name": "Authors", "value": publication["authorString"]})
    #     atributes.append({"name": "Title", "value": publication["title"]})
    #     atributes.append({"name": "Year", "value": publication["pubYear"]})
    #     atributes.append({"name": "Volume", "value": publication["journalVolume"]})
    #     atributes.append({"name": "Issue", "value": publication["issue"]})
    #     atributes.append({"name": "Type", "value": publication["pubType"]})
    #     atributes.append({"name": "Issn", "value": publication["journalIssn"]})
    #     atributes.append({"name": "DOI", "value": publication["doi"]})

    return {
        "title": result["title"],
        "pubYear": result["pubYear"],
        "authorString": result["authorString"],
        "journalTitle": result["journalTitle"],
        "journalVolume": result["journalVolume"],
        "journalIssn": result["journalIssn"],
        "issue": result["issue"],
        "pubType": result["pubType"],
        "pmid": result["pmid"],
        "doi": result["doi"],
    }


# export async function getPublicationData(pubmedId: string) {
#   if (pubmedId !== "") {
#     let response = await fetch(
#       `https://www.ebi.ac.uk/europepmc/webservices/rest/article/MED/${pubmedId.replace(
#         "PMID:",
#         ""
#       )}?resultType=lite&format=json`
#     );
#     if (!response.ok) {
#       throw new Error("Network response was not ok");
#     }
#     return response
#       .json()
#       .then((d) =>
#         Object.fromEntries(
#           ["title", "pubYear", "authorString", "journalTitle", "pmid", "doi"]
#             .filter((key) => key in d.result)
#             .map((key) => [key, d.result[key]])
#         )
#       );
#   }
# }


# Main function
def main(output):
    create_folder_if_not_exists(output)
    processed_data = fetch_and_process_data(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch, process, and formats CancerModels.org data into a structure that fits the BioStudies structurepython."
    )
    parser.add_argument(
        "--output", required=True, help="Folder where the data will be downloaded"
    )

    args = parser.parse_args()

    main(args.output)
