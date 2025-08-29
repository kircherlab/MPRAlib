import os
import tempfile
import gzip
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
    result = runner.invoke(cli, ["functional", "filter", "--help"])
    assert result.exit_code == 0
    print(result.output)
    assert "Usage: cli functional filter [OPTIONS]" in result.output


@pytest.fixture
def files():
    input_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")
    output_file_activity = tempfile.NamedTemporaryFile(delete=False).name
    output_file_barcode = tempfile.NamedTemporaryFile(delete=False).name
    yield {"input": input_file, "output_activity": output_file_activity, "output_barcode": output_file_barcode}
    os.remove(output_file_activity)
    os.remove(output_file_barcode)


def test_functional_filter_outliers_global(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "global",
            "--method-values",
            "{\"times_zscore\": 1000}",
            "--output-activity",
            files["output_activity"],
            "--output-barcode",
            files["output_barcode"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output_activity"])
    assert os.path.exists(files["output_barcode"])

    with open(files["output_barcode"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")
    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content

    with open(files["output_activity"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_activity.bc1.output.tsv.gz")
    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content


def test_functional_filter_outliers_large_expression(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "large_expression",
            "--output-activity",
            files["output_activity"],
            "--output-barcode",
            files["output_barcode"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output_activity"])
    assert os.path.exists(files["output_barcode"])

    with open(files["output_barcode"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")
    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content

    with open(files["output_activity"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_activity.bc1.output.tsv.gz")
    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content


def test_functional_filter_outliers_large_expression_2(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "large_expression",
            "--method-values",
            '{"times_activity": 2.0}',
            "--output-activity",
            files["output_activity"],
            "--output-barcode",
            files["output_barcode"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output_activity"])
    assert os.path.exists(files["output_barcode"])

    with open(files["output_barcode"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__), "data", "functional", "reporter_experiment_barcode.filter.large_expression.output.tsv.gz"
    )
    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content

    with open(files["output_activity"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__), "data", "functional", "reporter_experiment.filter.large_expression.output.tsv.gz"
    )
    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content


def test_functional_filter_dna_max_count_10(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "max_count",
            "--method-values",
            '{"dna_max_count": 10}',
            "--output-activity",
            files["output_activity"],
            "--output-barcode",
            files["output_barcode"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output_activity"])
    assert os.path.exists(files["output_barcode"])

    with open(files["output_barcode"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__), "data", "functional", "reporter_experiment_barcode.filter.max_count.output.tsv.gz"
    )
    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content

    with open(files["output_activity"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__), "data", "functional", "reporter_experiment.filter.max_count.output.tsv.gz"
    )
    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content
