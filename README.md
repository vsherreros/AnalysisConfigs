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

This will install the `configs` package in editable mode. 

## Configurations


