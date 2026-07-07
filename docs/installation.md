# Installation

This page describes the recommended installation paths for cryoROLE 2.0 from a
GitHub checkout.

## Requirements

cryoROLE requires Python 3.9 or newer.

Core Python dependencies are declared in `pyproject.toml`:

```text
matplotlib
numpy
pandas
scipy
```

Test dependencies are available through the optional `test` extra.

## Recommended Conda + GitHub Install

We recommend installing cryoROLE in a dedicated conda environment, then
installing the package from a GitHub checkout:

```bash
git clone https://github.com/yifancheng-ucsf/cryorole.git
cd cryorole

conda create -n cryorole python=3.10 -y
conda activate cryorole

python -m pip install .
cryorole --help
```

This installs the `cryorole` command-line entry point.

## Environment File Install

The repository also provides `environment.yml`, which creates an editable conda
environment from the repository root:

```bash
git clone https://github.com/yifancheng-ucsf/cryorole.git
cd cryorole

conda env create -f environment.yml
conda activate cryorole
cryorole --help
```

The environment file installs the package in editable mode so that local source
changes are immediately visible.

## Standard Pip Install from GitHub

If you already have a suitable Python environment, you can install directly from
the repository root:

```bash
git clone https://github.com/yifancheng-ucsf/cryorole.git
cd cryorole
python -m pip install .
cryorole --help
```

## Editable Developer Install

Use an editable install when developing cryoROLE or running tests from a source
checkout:

```bash
git clone https://github.com/yifancheng-ucsf/cryorole.git
cd cryorole
python -m pip install -e ".[test]"
cryorole --help
python -m pytest tests\core tests\normalize tests\match tests\io -q
```

On macOS/Linux, use forward slashes in the test paths:

```bash
python -m pytest tests/core tests/normalize tests/match tests/io -q
```

## Verification

After installation, check:

```bash
cryorole --help
cryorole run --help
```

For a lightweight test pass:

```bash
python -m pytest tests\core tests\normalize tests\match tests\io -q
```

## Notes

The current package version is `2.0.0a1`. Public documentation may describe this
as a 2.0 beta-preview workflow, but the package metadata remains alpha until a
formal beta package is cut.


