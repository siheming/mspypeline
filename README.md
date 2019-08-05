# README
This pipeline can be used to analyze the results of a MaxQuant analysis.

## Usage
This pipeline can be used via the command line and needs a python
installation with certain packages. A virtual environment is recommended
with all packages specified in the `requirements.txt` file. This can be
done for example via:
```bash
python3 -m venv -r requirements.txt maxquantpipeline
```
which can then be activated and deactivated via:
```bash
source maxquantpipeline  # activation
exit  # deactivation
```
When the environment is activated or the default python installation
satisfies the requirements the script can be used via:
```bash
python3 MaxQuantPipeline.py
```
To see help for the command line support type:
```bash
python3 MaxQuantPipeline.py --help
```
The arguments that can be specified when using the pipeline are:
- `--dir` the path to the directory that should be analyzed.
When this is not specified a window will open and ask to select a directory
- `--yml-file` the path to a yml file which should be used for analysis.
If the directory contains a config dir with a yml file it will be used
for analysis. Otherwise the user will be asked to select a yml file.
When this is skipped the default yml file will be used instead.
Using the default yml file can also be forced via `--yml-file default`

## Dependencies
The pipeline required multiple input files to perform the analysis. They
should be stored in a config dir on the same level as the pipeline script.
The requirements are:
- `ms_analysis_default.yml` a file which contains all defaults for the 
analysis pipeline.
- `important_protein_names.xlsx` a file which contains proteins which
should be analyzed. These impact descriptive plots and the score calculation.
- `important_receptor_names.xlsx` a file with receptors.

## Support
If additional support is required try googleing, asking a programmer or
contact me via `Simon.Heming@gmx.de`.