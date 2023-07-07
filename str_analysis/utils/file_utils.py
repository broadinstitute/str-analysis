import gzip
import io
import os

import hailtop.fs as hfs


def open_file(path):
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
