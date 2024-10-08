import glob
import os
import unittest

from setuptools import find_packages, setup
from setuptools.command.build_py import build_py

with open("README.md", "rt") as fh:
    long_description = fh.read()

with open("requirements.txt", "rt") as f:
    requirements = [r.strip() for r in f.readlines()]


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
                  "&& python3 -m twine upload dist/*whl dist/*gz")

setup(
    name='str_analysis',
    version="1.2.10",
    description="Utilities for analyzing short tandem repeats (STRs)",
    install_requires=requirements,
    cmdclass={
        'coverage': CoverageCommand,
        'publish': PublishCommand,
    },
    entry_points = {
        'console_scripts': [
            'call_non_ref_motifs = str_analysis.call_non_ref_motifs:main',
            'filter_vcf_to_STR_variants = str_analysis.filter_vcf_to_STR_variants:main',
            'add_adjacent_loci_to_expansion_hunter_catalog = str_analysis.add_adjacent_loci_to_expansion_hunter_catalog:main',
            'split_adjacent_loci_in_expansion_hunter_catalog = str_analysis.split_adjacent_loci_in_expansion_hunter_catalog:main',
            'check_trios_for_mendelian_violations = str_analysis.check_trios_for_mendelian_violations:main',
            'simulate_str_expansions = str_analysis.simulate_str_expansions:main',
            'combine_str_catalogs = str_analysis.combine_str_catalogs:main',
            'annotate_and_filter_str_catalog = str_analysis.annotate_and_filter_str_catalog:main',
            'annotate_EHdn_locus_outliers = str_analysis.annotate_EHdn_locus_outliers:main',
            'convert_annotated_EHdn_locus_outliers_to_expansion_hunter_catalog = str_analysis.convert_annotated_EHdn_locus_outliers_to_expansion_hunter_catalog:main',
            'combine_str_json_to_tsv = str_analysis.combine_str_json_to_tsv:main',
            'combine_json_to_tsv = str_analysis.combine_json_to_tsv:main',
            'copy_EH_vcf_fields_to_json = str_analysis.copy_EH_vcf_fields_to_json:main',
            'run_reviewer = str_analysis.run_reviewer:main',
            'print_reads = str_analysis.print_reads:main',
            'make_bamlet = str_analysis.make_bamlet:main',
            'make_minicram_for_exansion_hunter = str_analysis.make_minicram_for_exansion_hunter:main',
            'filter_out_loci_with_Ns_in_flanks = str_analysis.filter_out_loci_with_Ns_in_flanks:main',
        ],
    },
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=["str_analysis", "str_analysis.utils"],
    data_files=[
        ('data', ['str_analysis/data/non_ref_motif.offtarget_regions.json.gz', 'str_analysis/data/non_ref_motif.locus_info.json']),
        ('data', glob.glob('str_analysis/data/tests/*.*')),
    ],
    include_package_data=True,
    python_requires=">=3.7",
    license="MIT",
    keywords='',
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
