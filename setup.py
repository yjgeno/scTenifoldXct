import pathlib
from setuptools import setup, find_packages

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()
DESCRIPTION = "scTenifoldXct"
PACKAGES = find_packages(exclude=("tests*",))
exec(open('scTenifoldXct/version.py').read())

INSTALL_REQUIRES = [
        "anndata==0.8.0",
        "matplotlib~=3.5.1",
        "numpy~=1.21.6",
        "pandas~=1.4.2",
        "python_igraph~=0.9.10",
        "ray~=1.11.0",
        "scanpy==1.9.1",
        "scipy~=1.8.0",
        "statsmodels~=0.13.2",
        "torch>=1.10.2",
        "tqdm~=4.64.0",
    ]

setup(
    name="scTenifoldXct",
    version=__version__,
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
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10"
    ],
    python_requires='~=3.9.6',
    packages=PACKAGES, #["scTenifoldXct"], 
    # package_dir={"scTenifoldXct": 'scTenifoldXct'},
    # package_data={"scTenifoldXct": ['database/*.csv']},
    include_package_data=True, # MANIFEST
    install_requires=INSTALL_REQUIRES,
    entry_points={
        "console_scripts": [
            "run_Xct = scTenifoldXct.core:main",
            "run_pcNet = scTenifoldXct.pcNet:main"
        ]
    },
)