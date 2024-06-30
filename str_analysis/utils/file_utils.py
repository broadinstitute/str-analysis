import logging
logging.getLogger('asyncio').setLevel(logging.CRITICAL)

import gzip
import hailtop.fs as hfs
import io
import os
import re
import requests
import tempfile

from google.cloud import storage

gcloud_requester_pays_project = None

def set_requester_pays_project(gcloud_project):
    """Sets the requester pays project for all hailtop.fs calls"""
    global gcloud_requester_pays_project
    gcloud_requester_pays_project = gcloud_project


def open_file(path, *, download_local_copy_before_opening=False, gunzip=False, is_text_file=False):
    if path.startswith("gs://") and download_local_copy_before_opening:
        path = download_local_copy(path)

    path = os.path.expanduser(path)
    mode = "r"
    if path.startswith("gs://"):
        file = hfs.open(path, f"{mode}b", requester_pays_config=gcloud_requester_pays_project)
        if gunzip or path.endswith("gz"):
            file = gzip.GzipFile(fileobj=file, mode=mode)
    else:
        if gunzip or path.endswith("gz"):
            file = gzip.open(path, mode=mode)
        else:
            if is_text_file:
                file = open(path, f"{mode}t", encoding="utf-8")
            else:
                file = open(path, mode="rb")
            return file

    return io.TextIOWrapper(file, encoding="utf-8")


def file_exists(path):
    if path.startswith("gs://"):
        return hfs.exists(path, requester_pays_config=gcloud_requester_pays_project)

    path = os.path.expanduser(path)
    return os.path.isfile(path)


def get_file_size(path):
    if path.startswith("gs://"):
        return hfs.stat(path, requester_pays_config=gcloud_requester_pays_project).size
    else:
        return os.path.getsize(os.path.expanduser(path))


def download_local_copy(url_or_google_storage_path, verbose=False):
    """Downloads the given URL or gs:// path to a local temp file and returns the path to the local file."""

    temp_dir = tempfile.gettempdir()
    if url_or_google_storage_path.startswith("gs://"):
        path = os.path.join(temp_dir, os.path.basename(url_or_google_storage_path))
        if not os.path.isfile(path):
            if verbose:
                print(f"Downloading {url_or_google_storage_path} to {path}")
            hfs.copy(url_or_google_storage_path, path, requester_pays_config=gcloud_requester_pays_project)
    else:
        path = os.path.join(temp_dir, os.path.basename(url_or_google_storage_path))
        if not os.path.isfile(path):
            if verbose:
                print(f"Downloading {url_or_google_storage_path} to {path}")
            r = requests.get(url_or_google_storage_path)
            with open(path, "wb") as f:
                f.write(r.content)

    return path


def get_byte_range_from_google_storage(google_storage_path, start_bytes, end_bytes):
    """Downloads a byte range from a google storage path. To set a requester-pays project, call set_requester_pays_project(..)"""
    if not google_storage_path.startswith("gs://"):
        raise ValueError(f"Path {google_storage_path} must start with gs://")

    google_storage_path_match = re.match("^gs://([^/]+)/(.+)", google_storage_path)
    if not google_storage_path_match:
        raise ValueError(f"Path {google_storage_path} must be of the form gs://bucket/path/to/file")

    bucket_name, object_name = google_storage_path_match.groups()

    client = storage.Client(project=gcloud_requester_pays_project)
    bucket = client.bucket(bucket_name, user_project=gcloud_requester_pays_project)
    blob = storage.Blob(object_name, bucket)
    if not blob.exists():
        raise ValueError(f"{google_storage_path} not found")

    #print(f"Downloading {google_storage_path} [{start_bytes}-{end_bytes-1}]")
    return blob.download_as_bytes(start=start_bytes, end=end_bytes-1, raw_download=True)



