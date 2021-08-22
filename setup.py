import os
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import unittest

with open("README.md", "rt") as fh:
    long_description = fh.read()

class CoverageCommand(build_py):
    """Run all unittests and generate a coverage report."""
    def run(self):
        os.system("python3 -m coverage run --omit='*/site-packages/*,*tests.py,setup.py,*__init__.py' ./setup.py test "
                  "&& python3 -m coverage html --include=*.py "
                  "&& open htmlcov/index.html")


class PublishCommand(build_py):
    """Publish package to PyPI"""
    def run(self):
        os.system("rm -rf dist")
        os.system("python3 setup.py sdist"
                  "&& python3 setup.py bdist_wheel"
                  "&& twine upload dist/*whl dist/*gz")


with open("requirements.txt", "rt") as f:
    requirements = [r.strip() for r in f.readlines()]


def test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('str_analysis', pattern='*tests.py')
    return test_suite


setup(
    name='str_analysis',
    version="0.2.8",
    description="Utilities  short tandem repeats (STRs)",
    install_requires=requirements,
    cmdclass={
        'coverage': CoverageCommand,
        'publish': PublishCommand,
    },
    entry_points = {
        'console_scripts': [
            'call_rfc1_canvas_alleles = str_analysis.call_rfc1_canvas_alleles:main',
            'combine_json_to_tsv = str_analysis.combine_json_to_tsv:main',
        ],
    },
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=["str_analysis", "str_analysis.utils"],
    data_files=[
        ('data', ['str_analysis/data/offtarget_regions.json.gz']),
    ],
    include_package_data=True,
    python_requires=">=3.7",
    license="MIT",
    keywords='',
    test_suite="setup.test_suite",
    #tests_require=["mock"],
    url='https://github.com/broadinstitute/str-analysis',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)
