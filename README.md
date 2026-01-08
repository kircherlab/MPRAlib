# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                  |    Stmts |     Miss |   Cover |   Missing |
|-------------------------------------- | -------: | -------: | ------: | --------: |
| src/mpralib/\_\_init\_\_.py           |        1 |        0 |    100% |           |
| src/mpralib/cli.py                    |      440 |      206 |     53% |109, 166, 284-300, 447-461, 579-582, 655-724, 782-807, 869-922, 1005-1066, 1133-1223, 1313, 1372-1384, 1424, 1475 |
| src/mpralib/exception.py              |       12 |        0 |    100% |           |
| src/mpralib/mpradata.py               |      572 |       30 |     95% |344, 348, 357, 365, 384, 424, 482, 843-876, 884-885, 895-899, 963, 1109 |
| src/mpralib/utils/\_\_init\_\_.py     |        0 |        0 |    100% |           |
| src/mpralib/utils/file\_validation.py |      113 |        9 |     92% |32, 34, 41, 113-116, 130, 188 |
| src/mpralib/utils/io.py               |       96 |        0 |    100% |           |
| src/mpralib/utils/plot.py             |      110 |        2 |     98% |  168, 172 |
| tests/\_\_init\_\_.py                 |        0 |        0 |    100% |           |
| tests/test\_cli.py                    |       93 |        0 |    100% |           |
| tests/test\_cli\_combine.py           |       78 |        0 |    100% |           |
| tests/test\_cli\_functional.py        |      123 |        0 |    100% |           |
| tests/test\_cli\_plot.py              |       50 |        0 |    100% |           |
| tests/test\_cli\_validate\_file.py    |      111 |        0 |    100% |           |
| tests/test\_mpradata.py               |      755 |        0 |    100% |           |
| tests/test\_utils\_io.py              |      169 |        0 |    100% |           |
| tests/test\_utils\_plot.py            |       90 |        0 |    100% |           |
| **TOTAL**                             | **2813** |  **247** | **91%** |           |


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