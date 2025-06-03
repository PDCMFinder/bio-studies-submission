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
    "model_images",
    "license_name",
    "license_url",
    "data_source",
]

BASE_URL = "https://dev.cancermodels.org/api/"
SEARCH_INDEX_ENDPOINT = BASE_URL + "search_index"
CELL_MODEL_ENDPOINT = BASE_URL + "cell_model"
QUALITY_ASSURANCE = BASE_URL + "quality_assurance"
MODEL_MOLECULAR_METADATA_ENDPOINT = BASE_URL + "model_molecular_metadata"
MODEL_INFORMATION_ENDPOINT = BASE_URL + "model_information"
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
        quoted_value = urllib.parse.quote('in.("REV047")', safe="(),")
        params_str = urllib.parse.urlencode(PARAMS)
        filter_str = f"external_model_id={quoted_value}"
        final_url = f"{SEARCH_INDEX_ENDPOINT}?{params_str}&{filter_str}"

        response = requests.get(final_url)
        response.raise_for_status()
        data = response.json()

        # Process each model in the current batch
        for model in data:
            data_source = model["data_source"]
            model_folder_path = f"{output}/{data_source}/{model['external_model_id']}"
            create_folder_if_not_exists(model_folder_path)
            study = format_model(model, model_folder_path)

            json_file_name = f"{model_folder_path}/{model['external_model_id']}.json"

            with open(json_file_name, "w") as f:
                f.write(json.dumps(study, indent=2))

            # print(json.dumps(study))

        # Update metadata
        total_count += len(data)

        # Check if there are more results to fetch
        if len(data) < PARAMS["limit"]:
            break
        else:
            PARAMS["offset"] += PARAMS["limit"]

    return processed_data


def create_attribute(name, value):
    return {"name": name, "value": value if value else "N/A"}


def format_model(model, model_folder_path):
    study = {}
    # study["accno"] = model["external_model_id"] + "_" + model["data_source"]
    study["type"] = "submission"

    attributes = []

    title = f"[{model['data_source']}] [{model['model_type']}] [{model['external_model_id']}] {model['histology']}"
    attributes.append(create_attribute("Title", title))
    attributes.append(create_attribute("ReleaseDate", str(date.today())))

    study["attributes"] = attributes

    section = create_study_section(model, model_folder_path, title)
    study["section"] = section

    return study


def create_study_section(model, model_folder_path, title):
    study_section = {}
    study_section["accno"] = STUDY_SECTION_ID
    study_section["type"] = "Study"

    attributes = []
    attributes.append(create_attribute("Model ID", model["external_model_id"]))
    attributes.append(create_attribute("Title", title))
    attributes.append(create_attribute("Study type", model["model_type"]))
    attributes.append(create_attribute("Histology", model["histology"]))
    license = create_attribute("License", model["license_name"])
    license["valqual"] = [{"name": "url", "value": model["license_url"]}]

    extra_information = fetch_extra_information(model)
    if extra_information:
        model["other_model_links"] = extra_information.get("other_model_links")

    supplierLink = get_supplier(model["other_model_links"])
    if supplierLink:
        label = f"{supplierLink['link_label']} ({supplierLink['resource_label']})"
        supplier = create_attribute("Supplier", label)
        supplier["valqual"] = [{"name": "url", "value": supplierLink["link"]}]

        attributes.append(supplier)

    attributes.append(license)
    study_section["attributes"] = attributes
    study_section = add_subsections(study_section, model, model_folder_path)
    return study_section


def get_supplier(other_model_links):
    supplier = None
    if other_model_links:
        supplier = [x for x in other_model_links if x["type"] == "supplier"]
        if supplier != []:
            supplier = supplier[0]
    return supplier


def get_external_ids(other_model_links):
    external_ids = []
    if other_model_links:
        external_ids = [x for x in other_model_links if x["type"] == "external_id"]
    return external_ids


def add_section_if_not_null(current_subsections, new_subsection):
    if new_subsection is not None and new_subsection != []:
        current_subsections.append(new_subsection)


def add_subsections(study_section, model, model_folder_path):
    subsections = []
    add_section_if_not_null(subsections, create_organization_section(model))
    if model["model_type"] != "PDX":
        add_section_if_not_null(subsections, create_identifiers_subsection(model))
    add_section_if_not_null(subsections, create_patient_tumor_subsection(model))
    if model["model_type"] == "PDX":
        add_section_if_not_null(
            subsections, create_pdx_model_engraftment_subsection(model)
        )
    if model["model_type"] != "PDX":
        add_section_if_not_null(subsections, create_model_derivation_subsection(model))
    add_section_if_not_null(subsections, create_model_quality_control_subsection(model))
    molecular_metadata, molecular_data_files = create_molecular_data_subsection(
        model, model_folder_path
    )

    add_section_if_not_null(subsections, molecular_metadata)
    add_section_if_not_null(subsections, molecular_data_files)
    add_section_if_not_null(subsections, create_immune_markers_subsection(model))
    add_section_if_not_null(subsections, create_model_treatment_subsection(model))
    add_section_if_not_null(subsections, create_patient_treatment_subsection(model))
    add_section_if_not_null(subsections, create_histology_images_section(model))
    add_section_if_not_null(subsections, create_publication_subsection(model))
    add_section_if_not_null(subsections, create_related_models_subsection(model))

    study_section["subsections"] = subsections
    return study_section


def create_organization_section(model):
    organization_section = {"accno": "o1", "type": "Organization"}
    organization_section["attributes"] = [
        create_attribute("Name", f"{model['provider_name']} ({model['data_source']})")
    ]
    return organization_section


def create_identifiers_subsection(model):
    identifiers_section = []
    external_ids = get_external_ids(model["other_model_links"])
    if external_ids != []:
        for external_id in external_ids:
            row = {"type": "Identifiers"}
            attributes = []
            attributes.append(
                create_attribute("Resource", external_id["resource_label"])
            )
            link = create_attribute("Link", external_id["link_label"])
            link["valqual"] = [{"name": "url", "value": external_id["link"]}]
            attributes.append(link)
            row["attributes"] = attributes
            identifiers_section.append(row)
    return identifiers_section


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
        attributes.append(create_attribute("Host Strain Name", row["host_strain_name"]))
        attributes.append(
            create_attribute(
                "Host Strain Nomenclature", row["host_strain_nomenclature"]
            )
        )

        attributes.append(create_attribute("Site", row["engraftment_site"]))
        attributes.append(create_attribute("Type", row["engraftment_type"]))
        attributes.append(create_attribute("Material", row["engraftment_sample_type"]))
        attributes.append(
            create_attribute("Material Status", row["engraftment_sample_state"])
        )
        attributes.append(create_attribute("Passage", row["passage_number"]))

        subsection_row["attributes"] = attributes
        rows.append(subsection_row)

    return rows


def create_model_derivation_subsection(model):
    model_derivation_subsection = {"type": "Model derivation"}
    attributes = []

    url = f"{CELL_MODEL_ENDPOINT}?model_id=eq.{model['pdcm_model_id']}"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    for row in data:
        attributes.append(
            create_attribute("Growth properties", row["growth_properties"])
        )
        attributes.append(create_attribute("Growth media", row["growth_media"]))
        attributes.append(create_attribute("Plate coating", row["plate_coating"]))
        attributes.append(create_attribute("Passage", row["passage_number"]))
        attributes.append(create_attribute("Supplements", row["supplements"]))
        attributes.append(create_attribute("Contaminated", row["contaminated"]))
        attributes.append(
            create_attribute("Contamination details", row["contamination_details"])
        )
    model_derivation_subsection["attributes"] = attributes
    return [model_derivation_subsection]


def create_model_quality_control_subsection(model):
    url = f"{QUALITY_ASSURANCE}?model_id=eq.{model['pdcm_model_id']}"

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    rows = []

    for row in data:
        subsection_row = {}
        attributes = []

        subsection_row["type"] = "Model quality control"

        if model["model_type"] == "PDX":
            attributes.append(
                create_attribute("Technique", row["validation_technique"])
            )
            attributes.append(
                create_attribute(
                    "Description",
                    row["description"],
                )
            )
            attributes.append(create_attribute("Passage", row["passages_tested"]))
        else:
            attributes.append(
                create_attribute("Technique", row["validation_technique"])
            )
            attributes.append(create_attribute("Passage", row["passages_tested"]))
            attributes.append(
                create_attribute(
                    "Morphological features", row["morphological_features"]
                )
            )
            attributes.append(create_attribute("STR analysis", row["str_analysis"]))
            attributes.append(create_attribute("Model purity", row["model_purity"]))

        subsection_row["attributes"] = attributes
        rows.append(subsection_row)

    return rows


def fetch_extra_information(model):
    url = f"{MODEL_INFORMATION_ENDPOINT}?id=eq.{model['pdcm_model_id']}"
    url = url + "&select=other_model_links"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    if data and data != []:
        return data[0]
    else:
        return None


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
        create_attribute("Sample ID", model_molecular_metadata_row["sample_id"])
    )
    sample_type = get_sample_type_by_source(model_molecular_metadata_row["source"])
    attributes.append(create_attribute("Sample Type", sample_type))
    attributes.append(
        create_attribute(
            "Engrafted Tumour Passage",
            model_molecular_metadata_row["xenograft_passage"],
        )
    )
    attributes.append(
        create_attribute("Data Type", model_molecular_metadata_row["data_type"])
    )
    attributes.append(
        create_attribute(
            "Platform Used",
            model_molecular_metadata_row["platform_name"],
        )
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

    download_path = f"{model_folder_path}/molecular_data"
    create_folder_if_not_exists(download_path)

    model_molecular_metadata_row["sample_id"]
    # Creates the file with the molecular data
    file_path, size = fetch_molecular_data(
        sample_id,
        data_type,
        platform_name,
        molecular_characterization_id,
        download_path,
    )

    if size == 0:
        return None

    file = {"path": file_path.split("/", 1)[1]}
    file["type"] = "file"
    file["size"] = size

    attributes = []

    attributes.append({"name": "Sample ID", "value": sample_id})
    sample_type = get_sample_type_by_source(model_molecular_metadata_row["source"])
    attributes.append(create_attribute("Sample Type", sample_type))
    attributes.append(
        create_attribute(
            "Engrafted Tumour Passage",
            model_molecular_metadata_row["xenograft_passage"],
        )
    )
    attributes.append(create_attribute("Data Type", data_type))
    attributes.append(
        create_attribute(
            "Platform Used",
            platform_name,
        )
    )
    file["attributes"] = attributes

    return file


def create_immune_markers_subsection(model):
    immune_markers_subsection = {"type": "Immune markers"}
    subsections = []
    url = f"{BASE_URL}/immunemarker_data_extended?model_id=eq.{model['external_model_id']}"

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    if data and data != []:
        model_genomics = [row for row in data if row["marker_type"] == "Model Genomics"]

        formated_model_genomics_data = format_immune_markers_data(model_genomics)
        model_genomics_section = []

        for sample_id in formated_model_genomics_data:
            model_genomics_section_row = {}
            attributes = []
            attributes.append(create_attribute("Sample ID", sample_id))
            data_by_sample = formated_model_genomics_data[sample_id]
            for marker_name in data_by_sample:
                marker_value = data_by_sample[marker_name]
                attributes.append(create_attribute(marker_name, marker_value))
            model_genomics_section_row["attributes"] = attributes
            model_genomics_section.append(model_genomics_section_row)
        subsections.append(model_genomics_section)

        hla = [row for row in data if row["marker_type"] == "HLA type"]
        formated_hla_data = format_immune_markers_data(hla)
        hla_section = []

        for sample_id in formated_hla_data:
            hla_section_row = {"type": "HLA"}
            attributes = []
            attributes.append(create_attribute("Sample ID", sample_id))
            data_by_sample = formated_hla_data[sample_id]
            for marker_name in data_by_sample:
                marker_value = data_by_sample[marker_name]
                attributes.append(create_attribute(marker_name, marker_value))
            hla_section_row["attributes"] = attributes
            hla_section.append(hla_section_row)
        if hla_section != []:
            subsections.append(hla_section)

    else:
        return None

    immune_markers_subsection["subsections"] = subsections
    return immune_markers_subsection


def format_immune_markers_data(data):
    per_sample_data = {}
    for row in data:
        sample_id = row["sample_id"]
        if sample_id not in per_sample_data:
            per_sample_data[sample_id] = {}
        marker_name = row["marker_name"]

        if marker_name not in per_sample_data[sample_id]:
            per_sample_data[sample_id][marker_name] = ""

        marker_value = row["marker_value"]
        extra = row["essential_or_additional_details"]
        value = marker_value
        if extra:
            value = value + " (" + extra + ")"

        if per_sample_data[sample_id][marker_name] == "":
            per_sample_data[sample_id][marker_name] = value
        else:
            per_sample_data[sample_id][marker_name] = (
                per_sample_data[sample_id][marker_name] + os.linesep + value
            )
    # Sort the data so the data is sorted my marker name (on each sample_id)
    for sample_id in per_sample_data:
        per_sample_data[sample_id] = dict(sorted(per_sample_data[sample_id].items()))
    return per_sample_data


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

        attributes.append(create_attribute("Drug", treatment_names))
        attributes.append(create_attribute("Dose", doses))
        attributes.append(create_attribute("Passage", passage_range))
        attributes.append(create_attribute("Response", treatment_response))
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

        attributes.append(create_attribute("Treatment", treatment_names))

        attributes.append(create_attribute("Dose", doses))

        attributes.append(create_attribute("Response", treatment_response))
        chembl_attribute = create_url_attribute(links, "ChEMBL")
        pubchem_attribute = create_url_attribute(links, "PubChem")

        attributes.append(chembl_attribute)
        attributes.append(pubchem_attribute)

        subsection_row["attributes"] = attributes

        rows.append(subsection_row)

    return rows


def create_histology_images_section(model):
    rows = []
    data = model["model_images"]
    if not data:
        return None
    for image in data:
        row = {"type": "Histology images"}
        attributes = []
        url = image["url"]
        file_name = url.rsplit("/", 1)[1]

        attributes.append(
            {
                "name": "Url",
                "value": file_name,
                "valqual": [
                    {
                        "name": "URL",
                        "value": url,
                    }
                ],
            }
        )
        attributes.append(create_attribute("Description", image["description"]))
        attributes.append(create_attribute("Sample type", image["sample_type"]))
        attributes.append(create_attribute("Passage", image["passage"]))
        attributes.append(create_attribute("Magnification", image["magnification"]))
        attributes.append(create_attribute("Staining", image["staining"]))
        row["attributes"] = attributes
        rows.append(row)
    return rows


def create_publication_subsection(model):
    publications_ids = model["pdx_model_publications"]
    if publications_ids is None:
        return None
    publications_ids = publications_ids.replace(" ", "").split(",")

    publications = [get_publication_data(val) for val in publications_ids]

    rows = []

    for publication in publications:
        row = {"type": "Publications"}
        atributes = []
        pmid = publication["pmid"]

        atributes.append(create_attribute("Title", publication["title"]))
        atributes.append(create_attribute("Authors", publication["authorString"]))
        atributes.append(create_attribute("Year", publication["pubYear"]))
        atributes.append(create_attribute("Volume", publication["journalVolume"]))
        atributes.append(create_attribute("Issue", publication["issue"]))
        atributes.append(create_attribute("Issn", publication["journalIssn"]))

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
            create_attribute("Type", "BioStudies"),
            create_attribute("Role", "child of"),
        ]
        link["attributes"] = attributes
        links.append(link)
    children = model_relationships["children"]
    if children and children != []:
        for child in children:
            child = child["external_model_id"]
            link = {"url": child}
            attributes = [
                create_attribute("Type", "BioStudies"),
                create_attribute("Role", "parent of"),
            ]
            link["attributes"] = attributes
            links.append(link)
    if links == []:
        return None
    related_models_subsection["links"] = links
    return related_models_subsection


def format_link(link):
    biostudies_link = {}
    attributes = []
    biostudies_link["url"] = link["url"]
    attributes.append(create_attribute("Description", link["label"]))
    attributes.append(create_attribute("Type", link["resource"].lower()))
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

    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    result = data["result"]

    return {
        "title": result.get("title", "N/A"),
        "pubYear": result.get("pubYear", "N/A"),
        "authorString": result.get("authorString", "N/A"),
        "journalTitle": result.get("journalTitle", "N/A"),
        "journalVolume": result.get("journalVolume", "N/A"),
        "journalIssn": result.get("journalIssn", "N/A"),
        "issue": result.get("issue", "N/A"),
        "pubType": result.get("pubType", "N/A"),
        "pmid": result.get("pmid", "N/A"),
        "doi": result.get("doi", "N/A"),
    }


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
