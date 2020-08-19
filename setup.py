from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("mspypeline/version.py", "r") as f:
    version = f.readline().split()[-1].strip('"')

setup(
    name="mspypeline",
    version=version,
    description="Package to analyze Mass Spec Data",
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
    ],
    project_urls={
        "Documentation": "https://mspypeline.readthedocs.io/en/stable/",
        "Source": "https://github.com/siheming/mspypeline",
        "Bug Tracker": "https://github.com/siheming/mspypeline/issues",
    }
)
