import argparse
import json
import os
import shutil
import requests

from biostudiesclient.api import Api
from biostudiesclient.auth import Auth
import urllib

auth = Auth()
auth.login()
api = Api(auth)

BASE_URL = "https://wwwdev.ebi.ac.uk/"
SEARCH_URL = BASE_URL + "biostudies/api/v1/search?query={model_id}"

CHECK_IF_EXISTS = True

SKIP_IF_EXISTS = True

skipped = 0


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

    if os.path.exists(files_directory):
        files = [
            os.path.join(provider, model, files_folder_name, f)
            for f in os.listdir(files_directory)
            if os.path.isfile(os.path.join(files_directory, f))
            and f.__contains__(".tsv")
        ]
        for file in files:
            print("\tUploading", file)
            if " " in file:
                print("There are spaces")
                current_path = os.path.join(data_folder, provider, model)
                new_path = os.path.join("to_fix", provider, model)
                shutil.move(current_path, new_path)
                print("Moved to to_fix location")
                return False

            local_path = os.path.join(data_folder, file)
            # Uploads file to BioStudies
            print("local_path", local_path)
            print("file", file)
            api.upload_file(local_path, file)
        print("Files updated")

    submit(model, submission)
    return True


def get_biostudies_base_url_from_env():
    biostudies_api_url_dev = "https://wwwdev.ebi.ac.uk/biostudies/submissions/api"
    return os.environ.get("BIOSTUDIES_API_URL", biostudies_api_url_dev)


def get_password_from_env():
    biostudies_password_dev = "CHANGE_ME"
    return os.environ.get("BIOSTUDIES_TOKEN", biostudies_password_dev)


def exists_already(model_id):
    header = {"X-SESSION-TOKEN": get_password_from_env()}
    model_id_encoded = model_id.replace(" ", "+")
    url = (
        f"{get_biostudies_base_url_from_env()}/submissions?keywords={model_id_encoded}"
    )

    response = requests.get(url, headers=header)
    response.raise_for_status()
    data = response.json()

    accession = None

    if data and data != []:
        # Check they actually match by looking if the title contains the model id
        if f"[{model_id}]" in data[0]["title"]:
            accession = data[0]["accno"]
            print(f"Model {model_id} existed under accession {accession}")

    return accession


def submit(model_id, submission):
    # Does it already exist?
    if CHECK_IF_EXISTS:
        accession = exists_already(model_id)
        if accession:
            submission["accno"] = accession

    response = api.create_submission(submission)

    assert response.json
    assert response.json["accno"]
    print(response.json["accno"])


def main(data_folder):
    processed = 0
    skipped = 0
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
            if SKIP_IF_EXISTS:
                accession = exists_already(model)
                if accession:
                    skipped += 1
                    print("skipped", skipped)
                    current_path = os.path.join(provider_path, model)
                    new_path = os.path.join("done", provider, model)
                    print("move (dup) from", current_path, "to", new_path)
                    shutil.move(current_path, new_path)
                    continue
            ok = process_model(data_folder, provider, model)
            if ok:
                processed += 1
                if processed % 100 == 0:
                    print("Processed:", processed)
                # Move to done
                current_path = os.path.join(provider_path, model)
                new_path = os.path.join("done", provider, model)
                print("move (sub) from", current_path, "to", new_path)
                shutil.move(current_path, new_path)
    print("Ended")
    print("skipped", skipped)
    print("Processed:", processed)


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
