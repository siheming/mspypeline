# README
This pipeline can be used to analyze the results of a MaxQuant analysis.

## Requirements
It is recommended to use this pipeline with git and anaconda, which need to be installed if they aren't
already. Proxies need to be set for these tools if they are set up (like in the DKFZ).
The repository can be downloaded for example via
`git clone https://github.com/siheming/mspypeline.git` or other ways.

## Usage
This pipeline can be used via the command line and needs a python
installation with certain packages. A virtual environment is recommended
with all packages specified in the `environment.yml` file. This can be
done for example via:
```bash
conda env create python=3.7 -f environment.yml
```
which can then be activated and deactivated via:
```bash
conda activate mspypeline # activation
conda deactivate  # deactivation
```
When the environment is activated or the default python installation
satisfies the requirements the script can be used via:
```bash
python3 MSPypeline.py
```
or
```bash
python MSPypeline.py
```
If the script is started with no further arguments the first prompt will ask for the directory,
the second promp for the yml config file. If the second prompt is cancelled the default yml file is used
To see help for the command line support type:
```bash
python3 MSPypeline.py --help
```
The arguments that can be specified when using the pipeline are:
- `--dir` the path to the directory that should be analyzed.
When this is not specified a window will open and ask to select a directory
- `--yml-file` the path to a yml file which should be used for analysis.
If the directory contains a config dir with a yml file it will be used
for analysis. Otherwise the user will be asked to select a yml file.
When this is skipped the default yml file will be used instead.
Using the default yml file can also be forced via `--yml-file default`
- `--loglevel` Logging level used during run. Should be from options 
(lowest to highest): DEBUG < INFO < WARNING < ERROR.
- `--has-replicates` do the names of the experiments in the result files include technical replicates. Default is false.

## Dependencies
The pipeline required multiple input files to perform the analysis. They
should be stored in a config dir on the same level as the pipeline script.
The requirements are:
- `ms_analysis_default.yml` a file which contains all defaults for the 
analysis pipeline.
- `go_terms` a directory containing (GO-term) txt files for proteins with which
should be analyzed. This influences the enrichment analysis of the GO-term plot.
- `pathways` a directory containing (pathway) txt files for proteins with which
should be analyzed. This setting impacts descriptive plots and score calculations.

## Support
If additional support is required try googleing, asking a programmer or
contact me via `Simon.Heming@gmx.de`.