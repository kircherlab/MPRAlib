import os
import tempfile
import pytest
from click.testing import CliRunner
import gzip
from mpralib.cli import cli


@pytest.fixture(scope="module")
def runner():

    return CliRunner()


def test_combine_group_exists(runner):
    result = runner.invoke(cli, ["combine", "--help"])
    assert result.exit_code == 0
    assert "Combine counts with other outputs." in result.output


def test_combine_get_counts_option_required(runner):
    result = runner.invoke(cli, ["combine", "get-counts"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_reporter_elements_option_required(runner):
    result = runner.invoke(cli, ["combine", "get-reporter-elements"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_reporter_variants_option_required(runner):
    result = runner.invoke(cli, ["combine", "get-reporter-variants"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_reporter_genomic_elements_option_required(runner):
    result = runner.invoke(cli, ["combine", "get-reporter-genomic-elements"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_reporter_genomic_variants_option_required(runner):
    result = runner.invoke(cli, ["combine", "get-reporter-genomic-variants"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_variant_counts_option_required(runner):
    result = runner.invoke(cli, ["combine", "get-variant-counts"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_variant_map_option_required(runner):
    result = runner.invoke(cli, ["combine", "get-variant-map"])
    assert result.exit_code != 0
    assert "Missing option '--sequence-design'" in result.output


@pytest.fixture
def files():
    base = os.path.dirname(__file__)
    output_file = tempfile.NamedTemporaryFile(delete=False).name
    yield {
        "input": {
            "barcode_counts": os.path.join(base, "data", "reporter_experiment_barcode_IGVFDS2165KBMD.input.head20000.tsv.gz"),
            "sequence_design": os.path.join(base, "data", "reporter_sequence_design.example.tsv.gz"),
        },
        "output": output_file,
    }
    os.remove(output_file)


def test_combine_get_counts_oligos(runner, files):
    result = runner.invoke(
        cli,
        [
            "combine",
            "get-counts",
            "--input",
            files["input"]["barcode_counts"],
            "--sequence-design",
            files["input"]["sequence_design"],
            "--output",
            files["output"],
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(files["output"])

    with open(files["output"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__), "data", "combine", "oligo_counts.IGVFDS2165KBMD.head20000.output.tsv.gz"
    )

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content


def test_combine_get_counts_barcodes(runner, files):
    result = runner.invoke(
        cli,
        [
            "combine",
            "get-counts",
            "--barcodes",
            "--input",
            files["input"]["barcode_counts"],
            "--sequence-design",
            files["input"]["sequence_design"],
            "--output",
            files["output"],
        ],
    )
    assert result.exit_code == 0
    assert os.path.exists(files["output"])

    with open(files["output"], "r") as f:
        output_content = f.read()
    expected_output_file = files["input"]["barcode_counts"]

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    expected_content = "ID" + expected_content[7:]
    assert output_content == expected_content


def test_combine_get_variant_map(runner, files):
    result = runner.invoke(
        cli,
        [
            "combine",
            "get-variant-map",
            "--sequence-design",
            files["input"]["sequence_design"],
            "--output",
            files["output"],
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(files["output"])

    with open(files["output"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "combine", "get_variant_map.output.tsv.gz")

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content
