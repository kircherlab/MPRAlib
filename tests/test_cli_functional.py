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


def test_filter_min_count_method(runner, files):
    # Should run without error and produce output files
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "min_count",
            "--method-values",
            '{"rna_min_count": 1}',
            "--bc-threshold",
            "1",
            "--output-activity",
            files["output_activity"],
            "--output-barcode",
            files["output_barcode"],
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(files["output_activity"])
    assert os.path.exists(files["output_barcode"])


def test_filter_max_count_method(runner, files):
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
            '{"dna_max_count": 10000}',
            "--output-activity",
            files["output_activity"],
            "--output-barcode",
            files["output_barcode"],
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(files["output_activity"])
    assert os.path.exists(files["output_barcode"])


def test_filter_random_method(runner, files):
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "random",
            "--method-values",
            '{"proportion": 0.5}',
            "--output-activity",
            files["output_activity"],
            "--output-barcode",
            files["output_barcode"],
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(files["output_barcode"])
    assert os.path.exists(files["output_activity"])


def test_filter_random_method_with_dict(runner, files):
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "random",
            "--method-values",
            "{'proportion': 0.5}",
            "--output-activity",
            files["output_activity"],
            "--output-barcode",
            files["output_barcode"],
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(files["output_barcode"])
    assert os.path.exists(files["output_activity"])


def test_filter_random_method_with_invalid_json(runner, files):
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "random",
            "--method-values",
            "not_a_dict",
        ],
    )
    assert result.exit_code != 0
    assert "Could not parse --method-values as dict or JSON." in result.output


def test_filter_without_method_values(runner, files):
    result = runner.invoke(
        cli,
        [
            "functional",
            "filter",
            "--input",
            files["input"],
            "--method",
            "min_count",
            "--output-activity",
            files["output_activity"],
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(files["output_activity"])


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
            '{"times_zscore": 1000}',
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
