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
    long_description=open("README.md", encoding='utf-8').read(),
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
        "numpy>=1.19.5",
        "scipy>=1.5.4",
        "pandas>=1.1.5",
        "hicstraw>=0.0.1",  # update versions
        "cooler>=0.8.11",
    ],
    entry_points={
        "console_scripts": [
            # CLI entry points:
            "process-2d-contacts=src.preprocessing._2Dcontacts_processing:main",
        ],
    },
)