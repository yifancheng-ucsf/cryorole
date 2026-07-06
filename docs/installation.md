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

## Standard User Install from GitHub

Clone the repository and install from the repository root:

```bash
git clone https://github.com/yifancheng-ucsf/cryorole.git
cd cryorole
python -m pip install .
cryorole --help
```

This installs the `cryorole` command-line entry point.

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

## Conda Environment

If you prefer conda, create the environment from the repository root:

```bash
conda env create -f environment.yml
conda activate cryorole
cryorole --help
```

The environment file installs the package in editable mode so that local source
changes are immediately visible.

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

Installation has been checked in the staging workspace with:

```bash
python -m pip install . --dry-run
python -m pip install -e . --dry-run
python -m pip install -e ".[test]" --dry-run
```

These checks validate package metadata and dependency resolution in the current
environment. They cannot guarantee every future machine, Python distribution, or
operating-system configuration, so installation should still be tested after the
clean files are copied into the real GitHub repository.
