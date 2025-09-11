import os
import tempfile
import pytest
from click.testing import CliRunner
from mpralib.cli import cli


@pytest.fixture(scope="module")
def runner():

    return CliRunner()


def test_combine_group_exists(runner):
    result = runner.invoke(cli, ["plot", "--help"])
    assert result.exit_code == 0
    assert "Plotting functions." in result.output


def test_combine_get_counts_option_required(runner):
    result = runner.invoke(cli, ["plot", "barcodes-per-oligo"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_reporter_elements_option_required(runner):
    result = runner.invoke(cli, ["plot", "correlation"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_reporter_variants_option_required(runner):
    result = runner.invoke(cli, ["plot", "dna-vs-rna"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_combine_get_reporter_genomic_elements_option_required(runner):
    result = runner.invoke(cli, ["plot", "outlier"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


@pytest.fixture
def files():
    base = os.path.dirname(__file__)
    output_file = tempfile.NamedTemporaryFile(delete=False).name
    yield {
        "input": {
            "barcode_counts": os.path.join(base, "data", "reporter_experiment_barcode_IGVFDS2165KBMD.input.head20000.tsv.gz"),
        },
        "output": output_file,
    }
    os.remove(output_file)


def test_plot_barcodes_per_oligo(runner, files):
    result = runner.invoke(
        cli,
        [
            "plot",
            "barcodes-per-oligo",
            "--input",
            files["input"]["barcode_counts"],
            "--output",
            files["output"],
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(files["output"])


def test_plot_correlation(runner, files):
    result = runner.invoke(
        cli,
        [
            "plot",
            "correlation",
            "--input",
            files["input"]["barcode_counts"],
            "--oligos",
            "--bc-threshold",
            10,
            "--modality",
            "activity",
            "--output",
            files["output"],
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(files["output"])


def test_plot_outlier(runner, files):
    result = runner.invoke(
        cli,
        [
            "plot",
            "outlier",
            "--input",
            files["input"]["barcode_counts"],
            "--bc-threshold",
            2,
            "--output",
            files["output"],
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(files["output"])


def test_plot_dna_vs_rna(runner, files):
    result = runner.invoke(
        cli,
        [
            "plot",
            "outlier",
            "--input",
            files["input"]["barcode_counts"],
            "--bc-threshold",
            2,
            "--output",
            files["output"],
        ],
    )

    assert result.exit_code == 0
    assert os.path.exists(files["output"])
