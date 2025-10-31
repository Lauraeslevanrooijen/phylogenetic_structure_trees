from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="phylo-struct-tree",
    version="0.1.0",
    author="Your Name",
    description="Create phylogenetic trees based on protein structure alignment",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "click>=8.0.0",
    ],
    entry_points={
        "console_scripts": [
            "phylo-struct-tree=phylo_struct_tree.cli:main",
        ],
    },
)
