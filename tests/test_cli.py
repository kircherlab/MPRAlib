import pytest
import os
import tempfile
from click.testing import CliRunner
from mpralib.cli import cli, _get_chr
import gzip
import pandas as pd
from logging import Logger


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


def test_cli_group_help(runner):
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "Command line interface of MPRAlib" in result.output


@pytest.fixture
def files():
    input_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")
    output_file = tempfile.NamedTemporaryFile(delete=False).name
    yield {"input": input_file, "output": output_file}
    os.remove(output_file)


def test_barcode_activities_bc1(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "activities",
            "--input",
            files["input"],
            "--barcode-level",
            "--output",
            files["output"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output"])

    with open(files["output"], "r") as f:
        output_content = f.read()

    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()

    assert output_content == expected_content


def test_activities_bc1(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "activities",
            "--input",
            files["input"],
            "--bc-threshold",
            "1",
            "--output",
            files["output"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output"])

    with open(files["output"], "r") as f:
        output_content = f.read()

    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_activity.bc1.output.tsv.gz")

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()

    assert output_content == expected_content


def test_activities_bc10(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "activities",
            "--input",
            files["input"],
            "--bc-threshold",
            "10",
            "--output",
            files["output"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output"])

    with open(files["output"], "r") as f:
        output_content = f.read()

    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_activity.bc10.output.tsv.gz")

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()

    assert output_content == expected_content


def test_activities_bc100(runner, files):

    # Run the command
    result = runner.invoke(
        cli,
        [
            "functional",
            "activities",
            "--input",
            files["input"],
            "--bc-threshold",
            "100",
            "--output",
            files["output"],
        ],
    )

    # Check the result
    assert result.exit_code == 0
    assert os.path.exists(files["output"])

    with open(files["output"], "r") as f:
        output_content = f.read()

    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_activity.bc100.output.tsv.gz")

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()

    assert output_content == expected_content


class DummyLogger(Logger):

    def __init__(self):
        self.messages = []

    def warning(self, msg, *args, **kwargs):
        self.messages.append(msg)


@pytest.fixture
def logger():
    return DummyLogger()


def test_get_chr_found(logger):
    # Prepare a chromosome map DataFrame
    map_df = pd.DataFrame(
        {"refseq": ["NC_000001.11", "NC_000002.12"], "ucsc": ["chr1", "chr2"], "release": ["GRCh38", "GRCh38"]}
    )
    variant_id = "NC_000001.11:12345:A:T"
    result = _get_chr(map_df, variant_id, logger)
    assert result == "chr1"
    assert logger.messages == []


def test_get_chr_not_found(logger):
    map_df = pd.DataFrame(
        {"refseq": ["NC_000001.11", "NC_000002.12"], "ucsc": ["chr1", "chr2"], "release": ["GRCh38", "GRCh38"]}
    )
    variant_id = "NC_000003.13:54321:G:C"
    result = _get_chr(map_df, variant_id, logger)
    assert result is None
    assert any("Contig NC_000003.13 of SPDI NC_000003.13:54321:G:C not found" in msg for msg in logger.messages)


def test_get_chr_handles_empty_map(logger):
    map_df = pd.DataFrame(columns=["refseq", "ucsc", "release"])
    variant_id = "NC_000004.14:11111:T:A"
    result = _get_chr(map_df, variant_id, logger)
    assert result is None
    assert any("Contig NC_000004.14 of SPDI NC_000004.14:11111:T:A not found" in msg for msg in logger.messages)


def test_get_chr_with_multiple_matches(logger):
    # Should return the first match if multiple rows match
    map_df = pd.DataFrame(
        {"refseq": ["NC_000005.15", "NC_000005.15"], "ucsc": ["chr5a", "chr5b"], "release": ["GRCh38", "GRCh37"]}
    )
    variant_id = "NC_000005.15:22222:C:G"
    result = _get_chr(map_df, variant_id, logger)
    assert result in ["chr5a", "chr5b"]
    assert logger.messages == []
