# Liquid Cell Atlas

[![Stars](https://img.shields.io/github/stars/semenko/liquid-cell-atlas?logo=GitHub&color=yellow)](https://github.com/semenko/liquid-cell-atlas/stargazers)
[![PyPI](https://img.shields.io/pypi/v/liquid-cell-atlas.svg)](https://pypi.org/project/liquid-cell-atlas)
[![Documentation Status](https://readthedocs.org/projects/liquid-cell-atlas/badge/?version=latest)](https://scvi.readthedocs.io/en/stable/?badge=stable)
[![Coverage](https://codecov.io/gh/semenko/liquid-cell-atlas/branch/main/graph/badge.svg)](https://codecov.io/gh/semenko/liquid-cell-atlas)
[![Code
Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)


A new approach to deconvolute cell-free DNA based on methylation patterns, focused on ultra-low frequency cell populations and immune lineages.

## Subcomponents

To kick off this project, we've broken out key tasks into individual folders. Once our code matures, we'll modularize things, add tests, and combine this into a single package.

* Data Wrangling
    * External Data
        * Data aggregation & annotation
    * Internal data
        * Preprocessing pipeline (snakemake?)
        * Data aggregation & annotation
    * Data preprocessing
        * Convert to [.pat](https://github.com/nloyfer/wgbs_tools/blob/master/docs/pat_format.md) or [.beta](https://github.com/nloyfer/wgbs_tools/blob/master/docs/beta_format.md)?
* Modeling (Code)
    * Simple fingerprint models
    * Tensorflow-stored models
        * TFDF / Gradient Boosted
        * TabNet
* Testing
    * CyTOF & FACS Data
    * Stats & Graphs


## Archiecture

```mermaid
graph TD;
    A-->B;
    A-->C;
    B-->D;
    C-->D;
```


## License

Not yet determined. For now, all rights reserved. At time of publication, we may consider GPL/MIT, etc. depending on commercial licensing thoughts.
