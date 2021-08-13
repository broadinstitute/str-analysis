import os
import sys
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

with open("README.md", "rt") as fh:
    long_description = fh.read()

install_requires = [
    "mock",
]


class CoverageCommand(build_py):
    """Run all unittests and generate a coverage report."""
    def run(self):
        os.system("python3 -m coverage run ./setup.py test "
                  "&& python3 -m coverage html --include=*.py "
                  "&& open htmlcov/index.html")


class PublishCommand(build_py):
    """Publish package to PyPI"""
    def run(self):
        os.system("rm -rf dist")
        os.system("python3 setup.py sdist"
                  "&& python3 setup.py bdist_wheel"
                  "&& twine upload dist/*whl dist/*gz")

setup(
    name='step_pipeline',
    version="0.1",
    description="STR analysis scripts",
    install_requires=install_requires,
    cmdclass={
        'coverage': CoverageCommand,
        'publish': PublishCommand,
    },
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(where="."),
    package_dir={'': '.'},
    include_package_data=True,
    python_requires=">=3.7",
    license="MIT",
    keywords='',
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
