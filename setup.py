"""
for installation of package in editable mode
"""

from setuptools import setup, find_packages

setup(
    name="3D-genome-hubs",
    version="0.1.0",
    author="Akanksha Sachan",
    author_email="akanksha.11.05.07@gmail.com",
    description="A package for processing bulk 3D genome data and identifying genomic hubs of cis and trans regulatory elements responsible for gene regulation.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/akanksha-sachan/3D-genome-hubs",
    packages=find_packages(where="src"),
    package_dir={"": "src"},  # root is 'src'
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy==1.26.2",
        "scipy==1.11.4",
        "pandas==1.5.3",
        "hic-straw==1.3.1",
        "cooler==0.9.3",
        "asciitree==0.3.3",
        "astroid==3.2.0",
        "bioframe==0.5.1",
        "black==24.4.2",
        "certifi==2023.11.17",
        "charset-normalizer==3.3.2",
        "click==8.1.7",
        "cmake==3.27.7",
        "cooltools==0.6.1",
        "cramjam==2.7.0",
        "cython==3.0.5",
        "cytoolz==0.12.2",
        "dill==0.3.7",
        "exceptiongroup==1.2.0",
        "fastparquet==2023.10.1",
        "fonttools==4.45.1",
        "fsspec==2023.10.0",
        "h5py==3.10.0",
        "hic-straw==1.3.1",
        "hic2cool==1.0.1",
        "idna==3.4",
        "igraph==0.11.3",
        "imageio==2.33.0",
        "importlib-metadata==6.8.0",
        "iniconfig==2.0.0",
        "ipdb==0.13.13",
        "ipython==8.18.0",
        "isort==5.13.2",
        "js2py==0.74",
        "juicebox-notebook==0.2.1",
        "lazy-loader==0.3",
        "llvmlite==0.41.1",
        "mccabe==0.7.0",
        "memory-profiler==0.61.0",
        "multiprocess==0.70.15",
        "mypy-extensions==1.0.0",
        "networkit==10.1",
        "networkx==3.2.1",
        "numba==0.58.1",
        "pairtools==1.0.3",
        "pathspec==0.12.1",
        "pexpect==4.9.0",
        "pillow==10.1.0",
        "pluggy==1.5.0",
        "prompt-toolkit==3.0.41",
        "pyarrow==14.0.1",
        "pybigwig==0.3.22",
        "pybind11==2.12.0",
        "pyexecjs==1.5.1",
        "pyfaidx==0.7.2.2",
        "pyjsparser==2.7.1",
        "pylint==3.2.0",
        "pysam==0.22.0",
        "pytest==8.1.1",
        "python-louvain==0.16",
        "pytz==2023.3.post1",
        "pyyaml==6.0.1",
        "requests==2.31.0",
        "scikit-image==0.22.0",
        "scikit-learn==1.3.2",
        "scikit-network==0.31.0",
        "seaborn==0.13.2",
        "simplejson==3.19.2",
        "stack-data==0.6.3",
        "texttable==1.7.0",
        "tifffile==2023.9.26",
        "tomli==2.0.1",
        "tomlkit==0.12.5",
        "toolz==0.12.0",
        "tqdm==4.66.1",
        "traitlets==5.13.0",
        "trange==0.1.1",
        "typing-extensions==4.8.0",
        "tzdata==2023.3",
        "tzlocal==5.2",
        "urllib3==2.1.0",
        "wcwidth==0.2.12",
    ],
    entry_points={
        "console_scripts": [
            # CLI entry points:
            "process-2d-contacts=src.preprocessing._2Dcontacts_processing:main",
        ],
    },
)
