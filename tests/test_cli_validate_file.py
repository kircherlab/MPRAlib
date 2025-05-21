import unittest
import os
from click.testing import CliRunner
from mpralib.cli import cli
from mpralib.utils.file_validation import ValidationSchema


class TestMPRlibCLIValidateFile(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

        # Create a temporary input file
        self.files = {
            ValidationSchema.REPORTER_SEQUENCE_DESIGN: os.path.join(
                os.path.dirname(__file__),
                "data",
                "reporter_sequence_design.example.tsv.gz",
            ),
            ValidationSchema.REPORTER_BARCODE_TO_ELEMENT_MAPPING: os.path.join(
                os.path.dirname(__file__),
                "data",
                "reporter_barcode_to_element_mapping.example.tsv.gz",
            ),
            ValidationSchema.REPORTER_EXPERIMENT_BARCODE: os.path.join(
                os.path.dirname(__file__),
                "data",
                "reporter_experiment_barcode.input.head101.tsv.gz",
            ),
            ValidationSchema.REPORTER_EXPERIMENT: os.path.join(
                os.path.dirname(__file__),
                "data",
                "reporter_activity.bc100.output.tsv.gz",
            ),
            ValidationSchema.REPORTER_ELEMENT: os.path.join(
                os.path.dirname(__file__),
                "data",
                "reporter_element.example.tsv.gz",
            ),
            ValidationSchema.REPORTER_VARIANT: os.path.join(
                os.path.dirname(__file__),
                "data",
                "reporter_variants.example.tsv.gz",
            ),
            ValidationSchema.REPORTER_GENOMIC_ELEMENT: os.path.join(
                os.path.dirname(__file__),
                "data",
                "reporter_genomic_element.example.bed.gz",
            ),
            ValidationSchema.REPORTER_GENOMIC_VARIANT: os.path.join(
                os.path.dirname(__file__),
                "data",
                "reporter_genomic_variant.example.bed.gz",
            )
        }

    def test_reporter_genomic_varianz(self):

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "validate-file",
                "reporter-genomic-variant",
                "--input",
                self.files[ValidationSchema.REPORTER_GENOMIC_VARIANT],
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)

    def test_reporter_genomic_element(self):

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "validate-file",
                "reporter-genomic-element",
                "--input",
                self.files[ValidationSchema.REPORTER_GENOMIC_ELEMENT],
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)

    def test_reporter_variant(self):

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "validate-file",
                "reporter-variant",
                "--input",
                self.files[ValidationSchema.REPORTER_VARIANT],
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)

    def test_reporter_element(self):

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "validate-file",
                "reporter-variant",
                "--input",
                self.files[ValidationSchema.REPORTER_ELEMENT],
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)

    def test_reporter_experiment(self):

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "validate-file",
                "reporter-experiment",
                "--input",
                self.files[ValidationSchema.REPORTER_EXPERIMENT],
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)

    def test_reporter_experiment_barcode(self):

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "validate-file",
                "reporter-experiment-barcode",
                "--input",
                self.files[ValidationSchema.REPORTER_EXPERIMENT_BARCODE],
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)

    def test_reporter_barcode_to_element_mapping(self):

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "validate-file",
                "reporter-barcode-to-element-mapping",
                "--input",
                self.files[ValidationSchema.REPORTER_BARCODE_TO_ELEMENT_MAPPING],
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)

    def test_reporter_sequence_design(self):

        # Run the command
        result = self.runner.invoke(
            cli,
            [
                "validate-file",
                "reporter-sequence-design",
                "--input",
                self.files[ValidationSchema.REPORTER_SEQUENCE_DESIGN],
            ],
        )

        # Check the result
        self.assertIs(result.exit_code, 0)
