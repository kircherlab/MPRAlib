import unittest
import os
import tempfile
from click.testing import CliRunner
from mpralib.cli import cli
import gzip


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
