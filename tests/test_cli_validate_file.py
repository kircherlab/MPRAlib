import os
import pytest
from click.testing import CliRunner
from mpralib.cli import cli
from mpralib.utils.file_validation import ValidationSchema


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


@pytest.fixture(scope="module")
def files():
    base = os.path.dirname(__file__)
    return {
        ValidationSchema.REPORTER_SEQUENCE_DESIGN: os.path.join(base, "data", "reporter_sequence_design.example.tsv.gz"),
        ValidationSchema.REPORTER_BARCODE_TO_ELEMENT_MAPPING: os.path.join(
            base, "data", "reporter_barcode_to_element_mapping.example.tsv.gz"
        ),
        ValidationSchema.REPORTER_EXPERIMENT_BARCODE: os.path.join(
            base, "data", "reporter_experiment_barcode.input.head101.tsv.gz"
        ),
        ValidationSchema.REPORTER_EXPERIMENT: os.path.join(base, "data", "reporter_activity.bc100.output.tsv.gz"),
        ValidationSchema.REPORTER_ELEMENT: os.path.join(base, "data", "reporter_element.example.tsv.gz"),
        ValidationSchema.REPORTER_VARIANT: os.path.join(base, "data", "reporter_variants.example.tsv.gz"),
        ValidationSchema.REPORTER_GENOMIC_ELEMENT: os.path.join(base, "data", "reporter_genomic_element.example.bed.gz"),
        ValidationSchema.REPORTER_GENOMIC_VARIANT: os.path.join(base, "data", "reporter_genomic_variant.example.bed.gz"),
        "REPORTER_GENOMIC_VARIANT_EMPTY_ALLELE": os.path.join(base, "data", "reporter_genomic_variant.example2.bed.gz"),
        "REPORTER_GENOMIC_VARIANT_FALSE": os.path.join(base, "data", "reporter_genomic_variant.example3.bed.gz"),
    }


def test_reporter_genomic_variant(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-genomic-variant",
            "--input",
            files[ValidationSchema.REPORTER_GENOMIC_VARIANT],
        ],
    )
    assert result.exit_code == 0


def test_reporter_genomic_variant_example2(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-genomic-variant",
            "--input",
            files["REPORTER_GENOMIC_VARIANT_EMPTY_ALLELE"],
        ],
    )
    assert result.exit_code == 0


def test_reporter_genomic_variant_example3(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-genomic-variant",
            "--input",
            files["REPORTER_GENOMIC_VARIANT_FALSE"],
        ],
    )
    assert result.exit_code == 1


def test_reporter_genomic_element(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-genomic-element",
            "--input",
            files[ValidationSchema.REPORTER_GENOMIC_ELEMENT],
        ],
    )
    assert result.exit_code == 0


def test_reporter_variant(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-variant",
            "--input",
            files[ValidationSchema.REPORTER_VARIANT],
        ],
    )
    assert result.exit_code == 0


def test_reporter_element(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-element",
            "--input",
            files[ValidationSchema.REPORTER_ELEMENT],
        ],
    )
    assert result.exit_code == 0


def test_reporter_experiment(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-experiment",
            "--input",
            files[ValidationSchema.REPORTER_EXPERIMENT],
        ],
    )
    assert result.exit_code == 0


def test_reporter_experiment_barcode(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-experiment-barcode",
            "--input",
            files[ValidationSchema.REPORTER_EXPERIMENT_BARCODE],
        ],
    )
    assert result.exit_code == 0


def test_reporter_barcode_to_element_mapping(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-barcode-to-element-mapping",
            "--input",
            files[ValidationSchema.REPORTER_BARCODE_TO_ELEMENT_MAPPING],
        ],
    )
    assert result.exit_code == 0


def test_reporter_sequence_design(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-sequence-design",
            "--input",
            files[ValidationSchema.REPORTER_SEQUENCE_DESIGN],
        ],
    )
    assert result.exit_code == 0
