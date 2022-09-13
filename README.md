![logo](assets/logo.png)

a Python tool to parse and manipulate [JPlace files](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009), the format of [Evolutionary Placement](https://arxiv.org/abs/0911.2852) results.

![PyPI](https://img.shields.io/pypi/v/jplaceparser)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/Robaina/Jplaceparser)
[![GitHub license](https://img.shields.io/github/license/Robaina/Jplaceparser)](https://github.com/Robaina/Jplaceparser/blob/master/LICENSE)
![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4)
[![DOI](https://zenodo.org/badge/428328886.svg)](https://zenodo.org/badge/latestdoi/428328886)

<br>

## What is JPLACEparser?

JplaceParser allows not only reading Jplace files as python dictionaries but also manipulating field values and exporting back to jplace files. It contains a number of filters to remove placements that do not satisfy quality criteria. Currently, placements can be filtered by three criteria involving the [pendant and distal](https://github.com/lczech/gappa/wiki/Subcommand:-assign#automatic-ratio-example) placement length as well as the phylogenetic tree diameter:

* Filter by maximum pendant length
* Filter by maximum pendant to distal length ratio
* Filter by maximum pendant to tree diameter ratio
* Filter my minimum LWR

This is an ongoing project!

## Installation

1. ```pip install jplaceparser```

or

2. Git clone project to local directory.

   In terminal navigate to directory and enter: ```python setup.py install```

## Usage

You can find a jupyter notebook with usage examples [here](examples/examples.ipynb).

## Citation

If you use this software, please cite it as below:

Robaina-Estévez, S. (2022). JPLACEparser: a Python tool to parse and manipulate JPlace files (Version 0.0.1)[Computer software]. https://doi.org/10.5281/zenodo.7031582.

