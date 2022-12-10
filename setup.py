import pathlib
from setuptools import setup, find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

REQUIRES_PYTHON = '>=3.7.0'

setup(
    name="AGouTI",
    version="1.0.3",
    description="Annotation of Genomic and Transcriptomic Intervals using GTF/GFF files.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="ttps://github.com/zywicki-lab/agouti",
    author="Jan Kosinski",
    python_requires=REQUIRES_PYTHON,
    author_email="jankos@amu.edu.pl",
    license="GNU General Public License",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(),
    include_package_data=True,
    package_data={'': ["sample_data.tar.gz"]},
    setup_requires=["numpy", "pandas"],
    install_requires=["numpy", "pandas"],
    scripts=['agouti_pkg/agouti'],

)
