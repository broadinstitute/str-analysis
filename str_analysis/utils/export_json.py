import gzip
import simplejson as json
import os


def export_json(json_data, local_output_path, google_storage_dir=None):
    """Utility function for writing a json data structure to a file.

    Args:
        json_data (dict or list): The .json structure to write out.
        local_output_path (str): Local output path (where to write the json data).
        google_storage_dir (str): Optional "gs://" directory path to copy the local json file to after it's written.
    """

    if google_storage_dir and not google_storage_dir.startswith("gs://"):
        raise ValueError(f"Google Storage dir ({google_storage_dir}) doesn't start with 'gs://'")

    print(f"Writing {local_output_path}")
    open_file = gzip.open if local_output_path.endswith(".gz") else open
    with open_file(local_output_path, "wt") as f:
        json.dump(json_data, f, indent=4, ensure_ascii=True, allow_nan=False)

    if google_storage_dir:
        google_storage_path = os.path.join(google_storage_dir, os.path.basename(local_output_path))
        print(f"Copying {local_output_path} to {google_storage_path}")
        os.system(f"gsutil -m cp {local_output_path} {google_storage_path}")
