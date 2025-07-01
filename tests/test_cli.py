import unittest
import os
import tempfile
from click.testing import CliRunner
from mpralib.cli import cli, _get_chr
import gzip
import pandas as pd
from logging import Logger


class TestMPRlibCLI(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

        # Create a temporary input file
        self.input_file = os.path.join(
            os.path.dirname(__file__),
            "data",
            "reporter_experiment_barcode.input.tsv.gz",
        )

    def test_barcode_activities_bc1(self):

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(delete=False) as temp_output:
            output_file = temp_output.name

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "functional",
                "activities",
                "--input",
                self.input_file,
                "--barcode-level",
                "--output",
                output_file,
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)
        self.assertTrue(os.path.exists(output_file))

        with open(output_file, "r") as f:
            output_content = f.read()

        expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")

        with gzip.open(expected_output_file, "rt") as f:
            expected_content = f.read()

        self.assertTrue(output_content == expected_content)

    def test_activities_bc1(self):

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(delete=False) as temp_output:
            output_file = temp_output.name

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "functional",
                "activities",
                "--input",
                self.input_file,
                "--bc-threshold",
                "1",
                "--output",
                output_file,
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)
        self.assertTrue(os.path.exists(output_file))

        with open(output_file, "r") as f:
            output_content = f.read()

        expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_activity.bc1.output.tsv.gz")

        with gzip.open(expected_output_file, "rt") as f:
            expected_content = f.read()

        self.assertTrue(output_content == expected_content)

    def test_activities_bc10(self):

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(delete=False) as temp_output:
            output_file = temp_output.name

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "functional",
                "activities",
                "--input",
                self.input_file,
                "--bc-threshold",
                "10",
                "--output",
                output_file,
            ],
        )

        # Check the result
        assert result.exit_code == 0
        assert os.path.exists(output_file)

        with open(output_file, "r") as f:
            output_content = f.read()

        expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_activity.bc10.output.tsv.gz")

        with gzip.open(expected_output_file, "rt") as f:
            expected_content = f.read()

        assert output_content == expected_content

    def test_activities_bc100(self):

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(delete=False) as temp_output:
            output_file = temp_output.name

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "functional",
                "activities",
                "--input",
                self.input_file,
                "--bc-threshold",
                "100",
                "--output",
                output_file,
            ],
        )

        # Check the result
        assert result.exit_code == 0
        assert os.path.exists(output_file)

        with open(output_file, "r") as f:
            output_content = f.read()

        expected_output_file = os.path.join(os.path.dirname(__file__), "data", "reporter_activity.bc100.output.tsv.gz")

        with gzip.open(expected_output_file, "rt") as f:
            expected_content = f.read()

        assert output_content == expected_content

        # Clean up
        os.remove(output_file)


class DummyLogger(Logger):

    def __init__(self):
        self.messages = []

    def warning(self, msg, *args, **kwargs):
        self.messages.append(msg)


def test_get_chr_found():
    # Prepare a chromosome map DataFrame
    map_df = pd.DataFrame({
        "refseq": ["NC_000001.11", "NC_000002.12"],
        "ucsc": ["chr1", "chr2"],
        "release": ["GRCh38", "GRCh38"]
    })
    variant_id = "NC_000001.11:12345:A:T"
    logger = DummyLogger()
    result = _get_chr(map_df, variant_id, logger)
    assert result == "chr1"
    assert logger.messages == []


def test_get_chr_not_found():
    map_df = pd.DataFrame({
        "refseq": ["NC_000001.11", "NC_000002.12"],
        "ucsc": ["chr1", "chr2"],
        "release": ["GRCh38", "GRCh38"]
    })
    variant_id = "NC_000003.13:54321:G:C"
    logger = DummyLogger()
    result = _get_chr(map_df, variant_id, logger)
    assert result is None
    assert any("Contig NC_000003.13 of SPDI NC_000003.13:54321:G:C not found" in msg for msg in logger.messages)


def test_get_chr_handles_empty_map():
    map_df = pd.DataFrame(columns=["refseq", "ucsc", "release"])
    variant_id = "NC_000004.14:11111:T:A"
    logger = DummyLogger()
    result = _get_chr(map_df, variant_id, logger)
    assert result is None
    assert any("Contig NC_000004.14 of SPDI NC_000004.14:11111:T:A not found" in msg for msg in logger.messages)


def test_get_chr_with_multiple_matches():
    # Should return the first match if multiple rows match
    map_df = pd.DataFrame({
        "refseq": ["NC_000005.15", "NC_000005.15"],
        "ucsc": ["chr5a", "chr5b"],
        "release": ["GRCh38", "GRCh37"]
    })
    variant_id = "NC_000005.15:22222:C:G"
    logger = DummyLogger()
    result = _get_chr(map_df, variant_id, logger)
    assert result in ["chr5a", "chr5b"]
    assert logger.messages == []
