import os
import tempfile
import pytest
from click.testing import CliRunner
from mpralib.cli import cli
from mpralib.utils.file_validation import ValidationSchema


@pytest.fixture(scope="module")
def runner():

    return CliRunner()


def test_validate_file_group_exists(runner):
    result = runner.invoke(cli, ["validate-file", "--help"])
    assert result.exit_code == 0
    assert "Validate standardized MPRA reporter formats." in result.output


def test_validate_file_reporter_sequence_design_option_required(runner):
    result = runner.invoke(cli, ["validate-file", "reporter-sequence-design"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_validate_file_reporter_barcode_to_element_mapping_option_required(runner):
    result = runner.invoke(cli, ["validate-file", "reporter-barcode-to-element-mapping"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_validate_file_reporter_experiment_barcode_option_required(runner):
    result = runner.invoke(cli, ["validate-file", "reporter-experiment-barcode"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_validate_file_reporter_experiment_option_required(runner):
    result = runner.invoke(cli, ["validate-file", "reporter-experiment"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_validate_file_reporter_element_option_required(runner):
    result = runner.invoke(cli, ["validate-file", "reporter-element"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_validate_file_reporter_variant_option_required(runner):
    result = runner.invoke(cli, ["validate-file", "reporter-variant"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_validate_file_reporter_genomic_element_option_required(runner):
    result = runner.invoke(cli, ["validate-file", "reporter-genomic-element"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


def test_validate_file_reporter_genomic_variant_option_required(runner):
    result = runner.invoke(cli, ["validate-file", "reporter-genomic-variant"])
    assert result.exit_code != 0
    assert "Missing option '--input'" in result.output


@pytest.fixture(scope="module")
def files():
    base = os.path.dirname(__file__)
    invalid_file = tempfile.NamedTemporaryFile(delete=False).name
    with open(invalid_file, "w") as f:
        f.write("badcontent\tverybad\t12\n")
    yield {
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
        "invalid_file": invalid_file,
    }
    os.remove(invalid_file)


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


def test_reporter_genomic_variant_invalid(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-genomic-variant",
            "--input",
            files["invalid_file"],
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, TypeError)


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


def test_reporter_genomic_element_invalid(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-genomic-element",
            "--input",
            files["invalid_file"],
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, TypeError)


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


def test_reporter_variant_invalid(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-variant",
            "--input",
            files["invalid_file"],
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, SystemExit)


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


def test_reporter_element_invalid(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-element",
            "--input",
            files["invalid_file"],
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, SystemExit)


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


def test_reporter_experiment_barcode_invalid(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-experiment-barcode",
            "--input",
            files["invalid_file"],
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, SystemExit)


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


def test_reporter_barcode_to_element_mapping_invalid(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-barcode-to-element-mapping",
            "--input",
            files["invalid_file"],
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, SystemExit)


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


def test_reporter_sequence_design_invalid(runner, files):
    result = runner.invoke(
        cli,
        [
            "validate-file",
            "reporter-sequence-design",
            "--input",
            files["invalid_file"],
        ],
    )
    assert result.exit_code != 0
    assert isinstance(result.exception, SystemExit)
