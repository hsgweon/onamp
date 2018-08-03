from setuptools import setup, find_packages
import os

__version__ = os.environ.get("VERSION", "0.1")

setup(
    name = "onamp",
    version = __version__,
    packages = ["onamp"],
    scripts = ["bin/onamp_createreadpairslist", "bin/onamp", "bin/onamp_dada2", "bin/onamp_reformatAssignedTaxonomy", "bin/onamp_filterASVtable"],
    description = "ONAMP (Oh not another metabarcoding pipeline!): Bioinformatics pipeline for processing NGS data... yet again.",
    long_description = "An open source tools for processing NGS amplicon data.",
    author = "Hyun Soon Gweon",
    author_email = "h.s.gweon@reading.ac.uk",
	url = "https://github.com/hsgweon/onamp",
	license = "GPA"
)
