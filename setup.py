import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
DESCRIPTION = "scTenifoldXct"

setup(
    name="scTenifoldXct",
    version="0.0.1.dev1",
    description=DESCRIPTION,
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/cailab-tamu/scTenifoldXct",
    author="Yongjian Yang, TAMU",
    author_email="yjyang027@tamu.edu",
    license="MIT",
    keywords=[
        "neural network",
        "embedding",
        "manifold-learning",
        "computational-biology",
        "single-cell",
        "cell-cell interaction",
        "gene regulatroy network",
        "visualization"
    ],
    classifiers=[
    "License :: OSI Approved :: MIT License",
    "Environment :: Console",
    "Framework :: Jupyter",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10"
    ],
    packages=find_packages(exclude=("tests*",)), # ["scTenifoldXct"]
    include_package_data=True, # MANIFEST
    entry_points={
        "console_scripts": [
            "run_Xct = scTenifoldXct.core:main",
            "run_pcNet = scTenifoldXct.pcNet:main"
        ]
    },
)