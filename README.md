# AnalysisConfigs
Repository containing analysis configurations for PocketCoffea used in the ttH(bb) working group


## Setup

Please have a look at the [Installation guide](https://pocketcoffea.readthedocs.io/en/latest/installation.html). 

The `configs` package has been created to separate the core of the framework from all the necessary configuration files
and customization code needed the different analyses. The configuration is structured as a python package to make easier
the import of customization code into the framework configuration and also to make the sharing of analysis code easier.

Once you have a `PocketCoffea` local installation, you can install the `configs` package with:

```python
pip install -e .
```

The best way is to use a virtual environment created and activated inside the usual singularity machine or conda env.

```bash
apptainer shell  -B /afs -B /cvmfs/cms.cern.ch -B /tmp  -B /eos/cms/ -B /etc/sysconfig/ngbauth-submit -B
${XDG_RUNTIME_DIR} \
--env KRB5CCNAME="FILE:${XDG_RUNTIME_DIR}/krb5cc" \
/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest

python -m venv myenv
source myenv/bin/activate 

pip install -e . # only for the first time

```

This will install the `configs` package in editable mode. 

## Common pieces
Some folders inside the `configs` package are meant to contain common pieces of the configuration that can be imported
by several PocketCoffea configurations. 

- common_datasets: the dataset definition can be common between different channels. The XS of the samples is also stored
  here centrally. This allows keeping the naming of samples common between channels if possible.
- common_parameters: parameters configurations can be grouped here. The ideal usecase would be to save here the dumped parameter
  configuration for specific channels/tests with a meaningful name and version: doing so these parameters set can be
  imported by different configuration and customized. Using the parameter set name and git version we can also have a
  stable registry of metadata. 
- common_functions: common functions to perform selections and operations that are not part of the central PocketCoffea
  repository. 
- common_cuts:  registry of Cut objects and functions that can be used by many channels


## Configurations

## ttHbb Run2
### ttHbb - semileptonic
### ttHbb - dileptonic
### ttHbb - fullyhad

## ttHbb Run3
### ttHbb - semileptonic
### ttHbb - dileptonic
### ttHbb - fullyhad
