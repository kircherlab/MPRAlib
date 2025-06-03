# Changelog

## [0.6.5](https://github.com/kircherlab/MPRAlib/compare/v0.6.4...v0.6.5) (2025-06-03)


### Features

* adding complexity estimates ([2ef55ca](https://github.com/kircherlab/MPRAlib/commit/2ef55ca5ffd7bf7c60f1a395e7b4d4e42989ec11))
* complexity esimates in barcode level ([250c92f](https://github.com/kircherlab/MPRAlib/commit/250c92fa505d74b4d0ba955e0515dc5f9b60da46))
* complexity estimates on barcode level ([06a2173](https://github.com/kircherlab/MPRAlib/commit/06a21738d979a25c3d361b920ac70596332c347c))

## [0.6.4](https://github.com/kircherlab/MPRAlib/compare/v0.6.3...v0.6.4) (2025-05-28)


### Bug Fixes

* validate files. empty ref/alt allele for variant files must be 0 ([391703a](https://github.com/kircherlab/MPRAlib/commit/391703a80dcf5e91c2d4f5d32b8d2db02fffa82d))

## [0.6.3](https://github.com/kircherlab/MPRAlib/compare/v0.6.2...v0.6.3) (2025-05-23)


### Bug Fixes

* sampling was buggy ([abe540b](https://github.com/kircherlab/MPRAlib/commit/abe540b3f266aeb8d1a7b156e8d15ff73d6797e2))

## [0.6.2](https://github.com/kircherlab/MPRAlib/compare/v0.6.1...v0.6.2) (2025-05-22)


### Features

* cli to validate standardized file formats for MPRA using json schema ([a078ad9](https://github.com/kircherlab/MPRAlib/commit/a078ad9e0ae825b6b2ec53f403d9e69639f92ba5))
* coverage report and automatic tests using github actions ([d8139b6](https://github.com/kircherlab/MPRAlib/commit/d8139b6e834ea3e1cc8567b987e2a71ad66518c8))

## [0.6.1](https://github.com/kircherlab/MPRAlib/compare/v0.6.0...v0.6.1) (2025-04-09)


### Bug Fixes

* SyntaxError: EOL while scanning string literal ([11ce617](https://github.com/kircherlab/MPRAlib/commit/11ce6170b8a2042f47f9c9d7c7ab7f86e99db3ae))


### Documentation

* correct colab path ([8fbeb16](https://github.com/kircherlab/MPRAlib/commit/8fbeb16a81709caa7466bbd84ea43f87853f4e77))
* example data ([79cdae9](https://github.com/kircherlab/MPRAlib/commit/79cdae90e339d80e8ffd9ef1be1e97d02adc8717))
* updated mpralib notebook ([c1f53be](https://github.com/kircherlab/MPRAlib/commit/c1f53be195d8e4b668b3d03f32d4dfe92e1989e7))
* updated python notebooks ([83fdfc4](https://github.com/kircherlab/MPRAlib/commit/83fdfc47c54f9d2ea4c1a5e738d2eeffb2603b55))

## [0.6.0](https://github.com/kircherlab/MPRAlib/compare/v0.5.0...v0.6.0) (2025-04-03)


### ⚠ BREAKING CHANGES

* barcode_counts is now part available in MPRAdata as layer

### Code Refactoring

* barcode_counts is now part available in MPRAdata as layer ([c9b04ec](https://github.com/kircherlab/MPRAlib/commit/c9b04ec964ddc726a83afa54bc271657f3c8d6c9))

## [0.5.0](https://github.com/kircherlab/MPRAlib/compare/v0.4.0...v0.5.0) (2025-04-01)


### ⚠ BREAKING CHANGES

* gobal variant map. Now metadata is called "sequence design"
* dna_vs_rna and barcodes per oligo plots.
* correlation plot. breaks corrrelation due to new count type enum
* new utils sub library

### Features

* correlation plot. breaks corrrelation due to new count type enum ([f0e46a1](https://github.com/kircherlab/MPRAlib/commit/f0e46a1a2e0010080dc094a39bd192fbe4e12b90))
* dna_vs_rna and barcodes per oligo plots. ([3430538](https://github.com/kircherlab/MPRAlib/commit/34305388beeb3ef4df93dc11ca7c9677a8212227))
* gobal variant map. Now metadata is called "sequence design" ([2743e20](https://github.com/kircherlab/MPRAlib/commit/2743e20d5c891b480ec1425679af8d3c0300f132))
* normalized counts in output files ([8526f9f](https://github.com/kircherlab/MPRAlib/commit/8526f9f799547d9ef05ed416ceefb39dcea5306a))
* outlier plot ([17cf8d1](https://github.com/kircherlab/MPRAlib/commit/17cf8d1aa2160745d1fa79654babc9dfb13aff78))
* pairwise correlation plots ([bd1d118](https://github.com/kircherlab/MPRAlib/commit/bd1d118f26ee4203233df68243be68c7849a37d8))
* plot functionality ([b789c41](https://github.com/kircherlab/MPRAlib/commit/b789c416403a3bfd5880f6d909126c3924179fd6))


### Bug Fixes

* bcalm and mpralm element scripts with empty rows (Not observed) ([174657b](https://github.com/kircherlab/MPRAlib/commit/174657bc9df864b54da722b6b0f0d2c226794306))
* reporter genomic variant or element ([b250470](https://github.com/kircherlab/MPRAlib/commit/b250470774314ce2474f002edd678de37b38578d))


### Code Refactoring

* new utils sub library ([323e319](https://github.com/kircherlab/MPRAlib/commit/323e31914808dcc89ab849b19e7cc1e71332cc42))

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
