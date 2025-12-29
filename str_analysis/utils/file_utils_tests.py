#!/usr/bin/env python3

"""Comprehensive tests for file_utils.py"""

import gzip
import io
import os
import sys
import tempfile
import unittest
from unittest import mock

from str_analysis.utils.file_utils import (
    set_requester_pays_project,
    open_file,
    file_exists,
    get_file_size,
    download_local_copy,
    get_byte_range_from_google_storage,
    tee_stdout_and_stderr_to_log_file,
)


class TestSetRequesterPaysProject(unittest.TestCase):
    """Test set_requester_pays_project function."""

    def test_set_requester_pays_project(self):
        """Test setting requester pays project."""
        import str_analysis.utils.file_utils as file_utils

        # Set project
        set_requester_pays_project("my-test-project")
        self.assertEqual(file_utils.gcloud_requester_pays_project, "my-test-project")

        # Change project
        set_requester_pays_project("another-project")
        self.assertEqual(file_utils.gcloud_requester_pays_project, "another-project")

        # Set to None
        set_requester_pays_project(None)
        self.assertIsNone(file_utils.gcloud_requester_pays_project)


class TestOpenFile(unittest.TestCase):
    """Test open_file function."""

    def setUp(self):
        """Create temporary test files."""
        self.temp_dir = tempfile.mkdtemp()

        # Create a regular text file
        self.text_file = os.path.join(self.temp_dir, "test.txt")
        with open(self.text_file, "w") as f:
            f.write("Hello World\n")

        # Create a gzipped file
        self.gz_file = os.path.join(self.temp_dir, "test.txt.gz")
        with gzip.open(self.gz_file, "wt") as f:
            f.write("Compressed Hello World\n")

        # Create a binary file
        self.bin_file = os.path.join(self.temp_dir, "test.bin")
        with open(self.bin_file, "wb") as f:
            f.write(b"\x00\x01\x02\x03")

    def test_open_regular_text_file(self):
        """Test opening a regular text file."""
        with open_file(self.text_file, is_text_file=True) as f:
            content = f.read()
            self.assertEqual(content, "Hello World\n")

    def test_open_regular_binary_file(self):
        """Test opening a regular file in binary mode."""
        with open_file(self.bin_file) as f:
            content = f.read()
            self.assertEqual(content, b"\x00\x01\x02\x03")

    def test_open_gzipped_file_auto_detect(self):
        """Test opening a gzipped file with auto-detection."""
        with open_file(self.gz_file, is_text_file=True) as f:
            content = f.read()
            self.assertEqual(content, "Compressed Hello World\n")

    def test_open_gzipped_file_explicit(self):
        """Test opening a gzipped file with explicit gunzip=True."""
        # Rename to remove .gz extension to test explicit gunzip
        no_gz_file = os.path.join(self.temp_dir, "test_compressed")
        os.rename(self.gz_file, no_gz_file)

        with open_file(no_gz_file, gunzip=True, is_text_file=True) as f:
            content = f.read()
            self.assertEqual(content, "Compressed Hello World\n")

    def test_open_file_with_tilde_expansion(self):
        """Test that ~ is properly expanded."""
        # Create a file in temp dir
        test_file = os.path.join(self.temp_dir, "tilde_test.txt")
        with open(test_file, "w") as f:
            f.write("test")

        # Mock os.path.expanduser to return our temp dir
        with mock.patch('os.path.expanduser', return_value=test_file):
            with open_file("~/tilde_test.txt", is_text_file=True) as f:
                content = f.read()
                self.assertEqual(content, "test")

    @mock.patch('str_analysis.utils.file_utils.download_local_copy')
    def test_open_http_url_downloads_first(self, mock_download):
        """Test that HTTP URLs are downloaded first."""
        mock_download.return_value = self.text_file

        with open_file("http://example.com/test.txt", is_text_file=True) as f:
            content = f.read()
            self.assertEqual(content, "Hello World\n")

        mock_download.assert_called_once_with("http://example.com/test.txt")

    @mock.patch('str_analysis.utils.file_utils.download_local_copy')
    def test_open_https_url_downloads_first(self, mock_download):
        """Test that HTTPS URLs are downloaded first."""
        mock_download.return_value = self.text_file

        with open_file("https://example.com/test.txt", is_text_file=True) as f:
            content = f.read()
            self.assertEqual(content, "Hello World\n")

        mock_download.assert_called_once_with("https://example.com/test.txt")

    @mock.patch('str_analysis.utils.file_utils.download_local_copy')
    def test_open_gs_path_with_download_flag(self, mock_download):
        """Test that gs:// paths with download flag are downloaded first."""
        mock_download.return_value = self.text_file

        with open_file("gs://bucket/test.txt", download_local_copy_before_opening=True, is_text_file=True) as f:
            content = f.read()
            self.assertEqual(content, "Hello World\n")

        mock_download.assert_called_once_with("gs://bucket/test.txt")

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


class TestFileExists(unittest.TestCase):
    """Test file_exists function."""

    def setUp(self):
        """Create temporary test file."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_file = os.path.join(self.temp_dir, "exists.txt")
        with open(self.test_file, "w") as f:
            f.write("test")

    def test_file_exists_true(self):
        """Test file_exists returns True for existing file."""
        self.assertTrue(file_exists(self.test_file))

    def test_file_exists_false(self):
        """Test file_exists returns False for non-existing file."""
        self.assertFalse(file_exists(os.path.join(self.temp_dir, "nonexistent.txt")))

    def test_file_exists_with_tilde(self):
        """Test file_exists with tilde expansion."""
        with mock.patch('os.path.expanduser', return_value=self.test_file):
            self.assertTrue(file_exists("~/exists.txt"))

    def test_file_exists_directory_returns_false(self):
        """Test file_exists returns False for directories."""
        self.assertFalse(file_exists(self.temp_dir))

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


class TestGetFileSize(unittest.TestCase):
    """Test get_file_size function."""

    def setUp(self):
        """Create temporary test files."""
        self.temp_dir = tempfile.mkdtemp()

        # Create a file with known size
        self.test_file = os.path.join(self.temp_dir, "size_test.txt")
        test_content = "A" * 1024  # 1KB
        with open(self.test_file, "w") as f:
            f.write(test_content)

    def test_get_file_size(self):
        """Test get_file_size returns correct size."""
        size = get_file_size(self.test_file)
        self.assertEqual(size, 1024)

    def test_get_file_size_empty_file(self):
        """Test get_file_size for empty file."""
        empty_file = os.path.join(self.temp_dir, "empty.txt")
        with open(empty_file, "w") as f:
            pass

        size = get_file_size(empty_file)
        self.assertEqual(size, 0)

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)


class TestDownloadLocalCopy(unittest.TestCase):
    """Test download_local_copy function."""

    @mock.patch('requests.get')
    def test_download_http_url(self, mock_get):
        """Test downloading from HTTP URL."""
        # Mock the response
        mock_response = mock.Mock()
        mock_response.content = b"HTTP content"
        mock_get.return_value = mock_response

        # Download
        local_path = download_local_copy("http://example.com/test.txt")

        # Verify
        self.assertTrue(os.path.isfile(local_path))
        self.assertTrue(local_path.endswith("test.txt"))

        with open(local_path, "rb") as f:
            content = f.read()
            self.assertEqual(content, b"HTTP content")

        # Cleanup
        os.unlink(local_path)

    @mock.patch('requests.get')
    def test_download_https_url(self, mock_get):
        """Test downloading from HTTPS URL."""
        mock_response = mock.Mock()
        mock_response.content = b"HTTPS content"
        mock_get.return_value = mock_response

        local_path = download_local_copy("https://example.com/secure.txt")

        self.assertTrue(os.path.isfile(local_path))
        self.assertTrue(local_path.endswith("secure.txt"))

        with open(local_path, "rb") as f:
            content = f.read()
            self.assertEqual(content, b"HTTPS content")

        os.unlink(local_path)

    @mock.patch('requests.get')
    def test_download_url_caching(self, mock_get):
        """Test that downloading same URL twice uses cached copy."""
        mock_response = mock.Mock()
        mock_response.content = b"Cached content"
        mock_get.return_value = mock_response

        # First download
        local_path1 = download_local_copy("http://example.com/cached.txt")
        self.assertEqual(mock_get.call_count, 1)

        # Second download - should use cached file
        local_path2 = download_local_copy("http://example.com/cached.txt")
        self.assertEqual(mock_get.call_count, 1)  # Still only called once

        self.assertEqual(local_path1, local_path2)

        os.unlink(local_path1)

    @mock.patch('requests.get')
    def test_download_with_verbose(self, mock_get):
        """Test download with verbose=True prints message."""
        mock_response = mock.Mock()
        mock_response.content = b"content"
        mock_get.return_value = mock_response

        # Capture stdout
        with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
            local_path = download_local_copy("http://example.com/file.txt", verbose=True)
            output = mock_stdout.getvalue()
            self.assertIn("Downloading", output)
            self.assertIn("http://example.com/file.txt", output)

        os.unlink(local_path)


class TestGetByteRangeFromGoogleStorage(unittest.TestCase):
    """Test get_byte_range_from_google_storage function."""

    def test_invalid_path_no_gs_prefix(self):
        """Test that non-gs:// paths raise ValueError."""
        with self.assertRaises(ValueError) as cm:
            get_byte_range_from_google_storage("/local/path", 0, 100)
        self.assertIn("must start with gs://", str(cm.exception))

    def test_invalid_path_format(self):
        """Test that invalid gs:// format raises ValueError."""
        with self.assertRaises(ValueError) as cm:
            get_byte_range_from_google_storage("gs://", 0, 100)
        self.assertIn("must be of the form", str(cm.exception))

    @mock.patch('str_analysis.utils.file_utils.storage')
    def test_valid_byte_range_download(self, mock_storage):
        """Test downloading a byte range from gs://."""
        # Mock the storage client and blob
        mock_client = mock.Mock()
        mock_bucket = mock.Mock()
        mock_blob = mock.Mock()

        mock_storage.Client.return_value = mock_client
        mock_client.bucket.return_value = mock_bucket
        mock_storage.Blob.return_value = mock_blob
        mock_blob.exists.return_value = True
        mock_blob.download_as_bytes.return_value = b"byte_range_content"

        # Test
        content = get_byte_range_from_google_storage("gs://my-bucket/path/to/file.txt", 0, 100)

        # Verify
        self.assertEqual(content, b"byte_range_content")
        mock_blob.download_as_bytes.assert_called_once_with(start=0, end=99, raw_download=True)

    @mock.patch('str_analysis.utils.file_utils.storage')
    def test_nonexistent_blob_raises_error(self, mock_storage):
        """Test that non-existent blob raises ValueError."""
        mock_client = mock.Mock()
        mock_bucket = mock.Mock()
        mock_blob = mock.Mock()

        mock_storage.Client.return_value = mock_client
        mock_client.bucket.return_value = mock_bucket
        mock_storage.Blob.return_value = mock_blob
        mock_blob.exists.return_value = False

        with self.assertRaises(ValueError) as cm:
            get_byte_range_from_google_storage("gs://my-bucket/missing.txt", 0, 100)
        self.assertIn("not found", str(cm.exception))


class TestTeeStdoutAndStderr(unittest.TestCase):
    """Test tee_stdout_and_stderr_to_log_file function."""

    def setUp(self):
        """Create temporary directory for log files."""
        self.temp_dir = tempfile.mkdtemp()
        self.log_file = os.path.join(self.temp_dir, "test.log")

    def test_tee_stdout_to_log(self):
        """Test that stdout is written to both console and log file."""
        # Save original stdout
        original_stdout = sys.stdout

        # Set up tee
        tee_stdout_and_stderr_to_log_file(self.log_file)

        # Write to stdout
        print("Test message", flush=True)

        # Give threads time to write
        import time
        time.sleep(0.1)

        # Check log file
        with open(self.log_file, "r") as f:
            log_content = f.read()
            self.assertIn("Test message", log_content)

        # Restore stdout (best effort)
        sys.stdout = original_stdout

    def test_tee_creates_log_file(self):
        """Test that log file is created."""
        self.assertFalse(os.path.exists(self.log_file))

        tee_stdout_and_stderr_to_log_file(self.log_file)

        # Log file should be created
        self.assertTrue(os.path.exists(self.log_file))

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        try:
            shutil.rmtree(self.temp_dir)
        except:
            pass  # May fail due to file descriptors


if __name__ == "__main__":
    unittest.main()
