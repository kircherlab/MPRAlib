# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                  |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------------- | -------: | -------: | ------: | --------: |
| src/mpralib/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| src/mpralib/cli.py                    |      411 |      255 |     38% |41, 54, 67, 80, 93, 106, 119, 213-229, 277-298, 303, 331-371, 430-447, 516-585, 632-657, 708-759, 827-886, 930-1017, 1022, 1073-1087, 1130-1142, 1169-1176, 1196-1200, 1206 |
| src/mpralib/exception.py              |       12 |        0 |    100% |           |
| src/mpralib/mpradata.py               |      493 |       85 |     83% |46-49, 120, 171, 214, 293, 302, 310, 326-343, 347, 379, 402, 438-443, 450, 468-482, 699-710, 715-748, 754-780, 834, 848-854, 965, 975 |
| src/mpralib/utils/\_\_init\_\_.py     |        0 |        0 |    100% |           |
| src/mpralib/utils/file\_validation.py |      110 |        9 |     92% |32, 34, 41, 114-117, 129, 187 |
| src/mpralib/utils/io.py               |       97 |       10 |     90% |11-23, 64-66 |
| src/mpralib/utils/plot.py             |      108 |       34 |     69% |   117-184 |
| tests/\_\_init\_\_.py                 |        0 |        0 |    100% |           |
| tests/test\_cli.py                    |       59 |        0 |    100% |           |
| tests/test\_cli\_validate\_file.py    |       42 |        0 |    100% |           |
| tests/test\_mpradata.py               |      259 |        0 |    100% |           |
| tests/test\_utils\_io.py              |      111 |        0 |    100% |           |
| tests/test\_utils\_plot.py            |       90 |        0 |    100% |           |
|                             **TOTAL** | **1793** |  **393** | **78%** |           |


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