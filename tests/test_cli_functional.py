import os
import tempfile
import pytest
from click.testing import CliRunner
from mpralib.cli import cli


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


def test_functional_group_help(runner):
    result = runner.invoke(cli, ["functional", "--help"])
    assert result.exit_code == 0
    print(result.output)
    assert "General functionality." in result.output


def test_functional_filter_outliers_help(runner):
    result = runner.invoke(cli, ["functional", "filter-outliers", "--help"])
    assert result.exit_code == 0
    print(result.output)
    assert "Usage: cli functional filter-outliers [OPTIONS]" in result.output


@pytest.fixture
def files():
    input_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")
    output_file = tempfile.NamedTemporaryFile(delete=False).name
    yield {"input": input_file, "output": output_file}
    os.remove(output_file)


def test_test_functional_filter_outliers(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter-outliers",
            "--input",
            files["input"],
            "--method",
            "global",
            "--output",
            files["output"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output"])

    # with open(files["output"], "r") as f:
    #     output_content = f.read()

    # expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")

    # with gzip.open(expected_output_file, "rt") as f:
    #     expected_content = f.read()

    # assert output_content == expected_content
