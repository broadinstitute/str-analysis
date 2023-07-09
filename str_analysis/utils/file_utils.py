import gzip
import hailtop.fs as hfs
import io
import os
import requests
import tempfile


def open_file(path, download_local_copy_before_opening=False):
    if path.startswith("gs://") and download_local_copy_before_opening:
        path = download_local_copy(path)

    path = os.path.expanduser(path)
    mode = "r"
    if path.startswith("gs://"):
        file = hfs.open(path, f"{mode}b")
        if path.endswith("gz"):
            file = gzip.GzipFile(fileobj=file, mode=mode)
    else:
        if path.endswith("gz"):
            file = gzip.open(path, mode=mode)
        else:
            file = open(path, f"{mode}t", encoding="utf-8")
            return file

    return io.TextIOWrapper(file, encoding="utf-8")


def file_exists(path):
    if path.startswith("gs://"):
        return hfs.exists(path)

    path = os.path.expanduser(path)
    return os.path.isfile(path)


def download_local_copy(url_or_google_storage_path):
    """Downloads the given URL or gs:// path to a local temp file and returns the path to the local file."""

    temp_dir = tempfile.gettempdir()
    if url_or_google_storage_path.startswith("gs://"):
        path = os.path.join(temp_dir, os.path.basename(url_or_google_storage_path))
        if not os.path.isfile(path):
            print(f"Downloading {url_or_google_storage_path} to {path}")
            hfs.copy(url_or_google_storage_path, path)
    else:
        path = os.path.join(temp_dir, os.path.basename(url_or_google_storage_path))
        if not os.path.isfile(path):
            print(f"Downloading {url_or_google_storage_path} to {path}")
            r = requests.get(url_or_google_storage_path)
            with open(path, "wb") as f:
                f.write(r.content)

    return path


