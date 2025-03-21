# Changelog

## [0.4.0](https://github.com/kircherlab/MPRAlib/compare/v0.3.0...v0.4.0) (2025-03-21)


### ⚠ BREAKING CHANGES

* the whole code

### refactoring

* the whole code ([6f2f2b7](https://github.com/kircherlab/MPRAlib/commit/6f2f2b75ab3b1c3206527b88c629fa29fd472304))


### Features

* add filtering options for elements and variants in get_counts function ([971986b](https://github.com/kircherlab/MPRAlib/commit/971986b2593a1bf04e4c906a0000696e759a4dd5))
* full cli functionality ([6f2f2b7](https://github.com/kircherlab/MPRAlib/commit/6f2f2b75ab3b1c3206527b88c629fa29fd472304))


### Bug Fixes

* broken cli methods ([8386f83](https://github.com/kircherlab/MPRAlib/commit/8386f83fc1a9189ef9c6eb5c16343dac3ed4589c))

## [0.3.0](https://github.com/kircherlab/MPRAlib/compare/v0.2.2...v0.3.0) (2025-03-04)


### ⚠ BREAKING CHANGES

* splitting oligo and barcode counts resulting in different interface

### Features

* adding qc metric notebook ([ab17436](https://github.com/kircherlab/MPRAlib/commit/ab17436634d4b4ea83a788d3583531db96cd562f))


### Documentation

* qc metric update ([eda4ca4](https://github.com/kircherlab/MPRAlib/commit/eda4ca477ec8cc24dc88e9e7229ce9c19ca3e24f))


### Code Refactoring

* splitting oligo and barcode counts resulting in different interface ([ef09d7a](https://github.com/kircherlab/MPRAlib/commit/ef09d7a15a002748b09e771a0230fe56ab33dc7f))

## [0.2.2](https://github.com/kircherlab/MPRAlib/compare/v0.2.1...v0.2.2) (2025-02-27)


### Bug Fixes

* initialize barcode_threshold when it is not there ([4f09174](https://github.com/kircherlab/MPRAlib/commit/4f09174260cb349c5e4b89eb2c36f7b15c01ad13))

## [0.2.1](https://github.com/kircherlab/MPRAlib/compare/v0.2.0...v0.2.1) (2025-02-27)


### Features

* adding rna and dna counts correlations ([b2b4d60](https://github.com/kircherlab/MPRAlib/commit/b2b4d60258f179a68e77b90961abf85dc4b89b6b))

## [0.2.0](https://github.com/kircherlab/MPRAlib/compare/v0.1.0...v0.2.0) (2025-02-21)


### ⚠ BREAKING CHANGES

* implemented count sampling

### Features

* Control pseudo counts ([0534bb6](https://github.com/kircherlab/MPRAlib/commit/0534bb663e47e9ff77a5945340791144ad24eea5)), closes [#14](https://github.com/kircherlab/MPRAlib/issues/14)
* implement RNA and DNA count methods with sampling and filtering ([95e3be9](https://github.com/kircherlab/MPRAlib/commit/95e3be92b98677e571760c4c8efc807a4227610e))
* implemented count sampling ([bfa4da5](https://github.com/kircherlab/MPRAlib/commit/bfa4da55a73e19b069b1e44d9141b2f393e3114b))
* utils method for barcode output format ([#20](https://github.com/kircherlab/MPRAlib/issues/20)) ([81273de](https://github.com/kircherlab/MPRAlib/commit/81273de2335bc26c2b4c7fdaf1c1c6039970a726)), closes [#15](https://github.com/kircherlab/MPRAlib/issues/15)


### Bug Fixes

* issues and prepare to split up counts ([484d936](https://github.com/kircherlab/MPRAlib/commit/484d9367c723ad5004f9cb3dd59039461556d988))

## 0.1.0 (2025-02-17)


### Features

* add chromosome mapping utility and data files for hg19 and hg38 ([f037b92](https://github.com/kircherlab/MPRAlib/commit/f037b9257e5cd91d3197c21bb48035ec0fc1df9f))
* add MPRA analysis script and enhance MPRAdata class with variant mapping and counts ([00de9c0](https://github.com/kircherlab/MPRAlib/commit/00de9c0f4fca86abac764705806c1028f64bfd7f))
* add MPRA element analysis script and enhance MPRAdata class with element DNA and RNA counts ([df5e5fe](https://github.com/kircherlab/MPRAlib/commit/df5e5feaff1076fe2d88bce46a7f56cf3fc5cad3))
* add output option for MPRA data object and enhance element and variant reporting functions ([4a81736](https://github.com/kircherlab/MPRAlib/commit/4a81736157da9ca13af624d184f582b1bb2b810e))
* random sampling of barcodes ([62330d0](https://github.com/kircherlab/MPRAlib/commit/62330d04624384084ccd741deee3c78bba5628aa))


### Bug Fixes

* adjust YAML linter configuration for comments spacing ([ab6c18f](https://github.com/kircherlab/MPRAlib/commit/ab6c18f4c55fa876864b1d9f3d9f66f36ee65171))
