import logging
logging.getLogger('asyncio').setLevel(logging.CRITICAL)

import gzip
import hashlib
import io
import os
import re
import requests
import sys
import tempfile
import threading

from google.cloud import storage
from google.cloud.exceptions import NotFound

gcloud_requester_pays_project = None

def set_requester_pays_project(gcloud_project):
    """Sets the requester pays project for all hailtop.fs calls"""
    global gcloud_requester_pays_project
    gcloud_requester_pays_project = gcloud_project


def open_file(path, *, download_local_copy_before_opening=False, gunzip=False, is_text_file=False):
    if (path.startswith("gs://") and download_local_copy_before_opening) or path.startswith("http://") or path.startswith("https://"):
        path = download_local_copy(path)

    path = os.path.expanduser(path)
    mode = "r"
    if path.startswith("gs://"):
        try: import hailtop.fs as hfs
        except ImportError:
            print("ERROR: Hail is not installed. Please run: python3 -m pip install hail")
            sys.exit(1)

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

    if is_text_file:
        return io.TextIOWrapper(file, encoding="utf-8")
    else:
        return file


def file_exists(path):
    if path.startswith("gs://"):
        google_storage_path_match = re.match("^gs://([^/]+)/(.+)", path)
        if not google_storage_path_match:
            raise ValueError(f"Path {path} must be of the form gs://bucket/path/to/file")
        bucket_name, object_name = google_storage_path_match.groups()
        client = storage.Client(project=gcloud_requester_pays_project)
        bucket = client.bucket(bucket_name, user_project=gcloud_requester_pays_project)
        return storage.Blob(object_name, bucket).exists()

    path = os.path.expanduser(path)
    return os.path.isfile(path)


def get_file_size(path):
    if path.startswith("gs://"):
        google_storage_path_match = re.match("^gs://([^/]+)/(.+)", path)
        if not google_storage_path_match:
            raise ValueError(f"Path {path} must be of the form gs://bucket/path/to/file")
        bucket_name, object_name = google_storage_path_match.groups()
        client = storage.Client(project=gcloud_requester_pays_project)
        bucket = client.bucket(bucket_name, user_project=gcloud_requester_pays_project)
        blob = bucket.get_blob(object_name)
        if blob is None:
            raise ValueError(f"{path} not found")
        return blob.size
    else:
        return os.path.getsize(os.path.expanduser(path))


def download_local_copy(url_or_google_storage_path, verbose=False):
    """Downloads the given URL or gs:// path to a local temp file and returns the path to the local file.
    If the path is already a local file, returns it as-is.
    """

    # Handle local file paths - return as-is
    if not url_or_google_storage_path.startswith(("gs://", "http://", "https://")):
        return os.path.expanduser(url_or_google_storage_path)

    temp_dir = tempfile.gettempdir()
    # include a hash of the full source URI in the cache filename so that two different remote paths that
    # merely share a basename don't collide and silently reuse each other's cached content
    path = os.path.join(
        temp_dir,
        f"{hashlib.sha256(url_or_google_storage_path.encode()).hexdigest()[:16]}_{os.path.basename(url_or_google_storage_path)}")
    if url_or_google_storage_path.startswith("gs://"):
        if not os.path.isfile(path):
            if verbose:
                print(f"Downloading {url_or_google_storage_path} to {path}")

            # try using gsutil first (it's currently more reliable than hfs.copy)
            gcloud_requester_pays_arg = f"-u {gcloud_requester_pays_project}" if gcloud_requester_pays_project is not None else ""
            os.system(f"gsutil {gcloud_requester_pays_arg} -m cp {url_or_google_storage_path} {path}.temp")
            if not os.path.isfile(f"{path}.temp"):
                # fall back on hfs.copy
                try:
                    import hailtop.fs as hfs
                    hfs.copy(url_or_google_storage_path, f"{path}.temp", requester_pays_config=gcloud_requester_pays_project)
                except ImportError:
                    print("WARNING: Hail is not installed. Please run: python3 -m pip install hail")
                    sys.exit(1)

            os.rename(f"{path}.temp", path)
    else:
        if not os.path.isfile(path):
            if verbose:
                print(f"Downloading {url_or_google_storage_path} to {path}")
            r = requests.get(url_or_google_storage_path)
            with open(f"{path}.temp", "wb") as f:
                f.write(r.content)
            os.rename(f"{path}.temp", path)

    return path


def get_byte_range_from_google_storage(google_storage_path, start_bytes, end_bytes, client=None):
    """Downloads a byte range from a google storage path. To set a requester-pays project, call set_requester_pays_project(..).
    Pass a reusable storage.Client via `client` to avoid constructing a new one (and re-resolving credentials) on every call."""
    if not google_storage_path.startswith("gs://"):
        raise ValueError(f"Path {google_storage_path} must start with gs://")

    google_storage_path_match = re.match("^gs://([^/]+)/(.+)", google_storage_path)
    if not google_storage_path_match:
        raise ValueError(f"Path {google_storage_path} must be of the form gs://bucket/path/to/file")

    bucket_name, object_name = google_storage_path_match.groups()

    if client is None:
        client = storage.Client(project=gcloud_requester_pays_project)
    bucket = client.bucket(bucket_name, user_project=gcloud_requester_pays_project)
    blob = storage.Blob(object_name, bucket)

    #print(f"Downloading {google_storage_path} [{start_bytes}-{end_bytes-1}]")
    # download_as_bytes raises NotFound for a missing object (no separate exists() HEAD needed); translate it
    # to the ValueError this function has always raised for missing objects
    try:
        return blob.download_as_bytes(start=start_bytes, end=end_bytes-1, raw_download=True)
    except NotFound:
        raise ValueError(f"{google_storage_path} not found")


def tee_stdout_and_stderr_to_log_file(log_path):
    """Redirects stdout and stderr to both the console and a log file.

    This function creates a "tee" behavior for stdout and stderr, where all output
    is written to both the original destinations (console) and a log file simultaneously.

    Args:
        log_path (str): Path to the log file where output will be written

    Note:
        Uses daemon threads which will be terminated when the main program exits.
    """
    log_fd = os.open(log_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)

    # Save original fds
    stdout_fd = os.dup(1)
    stderr_fd = os.dup(2)

    def _tee(src_fd, dst_fds):
        while True:
            data = os.read(src_fd, 1024)
            if not data:
                break
            for fd in dst_fds:
                os.write(fd, data)

    # Create pipes
    r_out, w_out = os.pipe()
    r_err, w_err = os.pipe()

    # Redirect stdout/stderr → pipes
    os.dup2(w_out, 1)
    os.dup2(w_err, 2)

    # Start tee threads
    threading.Thread(target=_tee, args=(r_out, [stdout_fd, log_fd]), daemon=True).start()
    threading.Thread(target=_tee, args=(r_err, [stderr_fd, log_fd]), daemon=True).start()

    # Re-wrap Python streams
    sys.stdout = os.fdopen(1, 'w', buffering=1)
    sys.stderr = os.fdopen(2, 'w', buffering=1)
