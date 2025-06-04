# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                  |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------------- | -------: | -------: | ------: | --------: |
| src/mpralib/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| src/mpralib/cli.py                    |      411 |      255 |     38% |41, 54, 67, 80, 93, 106, 119, 213-229, 266-287, 292, 320-360, 419-436, 487-556, 603-628, 679-730, 798-857, 901-988, 993, 1044-1058, 1101-1113, 1140-1147, 1167-1171, 1177 |
| src/mpralib/exception.py              |       12 |        0 |    100% |           |
| src/mpralib/mpradata.py               |      433 |       85 |     80% |27-30, 58, 71, 106, 185, 189, 203-220, 224, 252, 278, 299-304, 307, 314, 326-340, 536-547, 552-585, 589-615, 657, 663-669, 746, 756 |
| src/mpralib/utils/\_\_init\_\_.py     |        0 |        0 |    100% |           |
| src/mpralib/utils/file\_validation.py |      110 |        9 |     92% |32, 34, 41, 115-118, 130, 188 |
| src/mpralib/utils/io.py               |       97 |       10 |     90% |11-23, 67-69 |
| src/mpralib/utils/plot.py             |      108 |       34 |     69% |   117-184 |
| tests/\_\_init\_\_.py                 |        0 |        0 |    100% |           |
| tests/test\_cli.py                    |       59 |        0 |    100% |           |
| tests/test\_cli\_validate\_file.py    |       42 |        0 |    100% |           |
| tests/test\_mpradata.py               |      230 |        0 |    100% |           |
| tests/test\_utils\_io.py              |      111 |        0 |    100% |           |
| tests/test\_utils\_plot.py            |       90 |        0 |    100% |           |
|                             **TOTAL** | **1704** |  **393** | **77%** |           |


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/kircherlab/MPRAlib/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/kircherlab/MPRAlib/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2Fkircherlab%2FMPRAlib%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.