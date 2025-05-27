import argparse
import json
import os
import requests

from biostudiesclient.api import Api
from biostudiesclient.auth import Auth

auth = Auth()
auth.login()
api = Api(auth)

BASE_URL = "https://wwwdev.ebi.ac.uk/"
SEARCH_URL = BASE_URL + "biostudies/api/v1/search?query={model_id}"

CHECK_IF_EXISTS = True


def process_model(data_folder, provider, model):
    print(f"process_model - {data_folder}, {provider}, {model}")
    model_path = os.path.join(data_folder, provider, model)

    model_submission_file = os.path.join(model_path, f"{model}.json")

    if not os.path.exists(model_submission_file):
        raise ValueError(f"No submission file {model_submission_file} found")

    submission = None
    with open(model_submission_file, "r") as file:
        submission = json.load(file)

    file_size = os.path.getsize(model_submission_file)
    if file_size == 0:
        raise ValueError(f"Submission file {model_submission_file} is empty")

    files_folder_name = "molecular_data"
    files_directory = os.path.join(model_path, files_folder_name)

    files = [
        os.path.join(provider, model, files_folder_name, f)
        for f in os.listdir(files_directory)
        if os.path.isfile(os.path.join(files_directory, f)) and f.__contains__(".tsv")
    ]
    for file in files:
        print("\t", file)
        print("Uploading to biostudies...")
        local_path = os.path.join(data_folder, file)
        # Uploads file to BioStudies
        api.upload_file(local_path, file)

    submit(model, submission)


def submit(model_id, submission):
    # Does it already exist?
    print(SEARCH_URL.format(model_id=model_id))
    if CHECK_IF_EXISTS:
        response = requests.get(SEARCH_URL.format(model_id=model_id))
        response.raise_for_status()
        data = response.json()
        # print(data)

        if data["totalHits"] != 0:
            hits = data["hits"]
            print("hits", hits)
            if hits != []:
                accession = hits[0]["accession"]
                print(f"Model {model_id} existed under accession {accession}")
                submission["accno"] = accession
                print("Updated submission", submission)

    response = api.create_submission(submission)

    assert response.json
    assert response.json["accno"]
    print(response.json["accno"])


def main(data_folder):
    if not os.path.exists(data_folder):
        print(f"Directory `{data_folder}` does not exist")
        return

    providers = [
        d
        for d in os.listdir(data_folder)
        if os.path.isdir(os.path.join(data_folder, d))
    ]
    for provider in providers:
        provider_path = os.path.join(data_folder, provider)
        models = [
            m
            for m in os.listdir(provider_path)
            if os.path.isdir(os.path.join(provider_path, m))
        ]
        for model in models:
            # print("Processing model", model, "(", provider, ")")
            process_model(data_folder, provider, model)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch, process, and formats CancerModels.org data into a structure that fits the BioStudies structurepython."
    )
    parser.add_argument(
        "--data_folder",
        required=True,
        help="Folder where the data to submit is located",
    )

    args = parser.parse_args()

    main(args.data_folder)
