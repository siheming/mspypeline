from setuptools import setup, find_packages
from mspypeline import __version__ as version

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="mspypeline",
    version=version,
    description="PLACEHOLDER",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,  # include files specified in MANIFEST, i.e. config/
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy>=1.17.4",
        "pandas>=0.25.3",
        "scipy>=1.3.1",
        "rpy2>=2.9.4",
        "tzlocal>=2.0.0",
        "ruamel_yaml>=0.15.46",
        "matplotlib>=3.1.1",
        "matplotlib-venn>=0.11.5",
        "adjusttext>=0.7.3",
        "scikit-learn>=0.22.1",
        "plotly>=4.6.0",
    ]
)
