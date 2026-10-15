import gzip
import os
import tempfile

import pytest
from click.testing import CliRunner

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
            "barcode_counts": os.path.join(
                base,
                "data",
                "reporter_experiment_barcode_IGVFDS2165KBMD.input.head20000.tsv.gz",
            ),
            "sequence_design": os.path.join(base, "data", "reporter_sequence_design.example.tsv.gz"),
        },
        "output": output_file,
    }
    os.remove(output_file)


def test_combine_get_counts_elements(runner, files):
    result = runner.invoke(
        cli,
        [
            "combine",
            "get-counts",
            "--elements-only",
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

    with open(files["output"]) as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__),
        "data",
        "combine",
        "element_counts.IGVFDS2165KBMD.head20000.output.tsv.gz",
    )

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content


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

    with open(files["output"]) as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__),
        "data",
        "combine",
        "oligo_counts.IGVFDS2165KBMD.head20000.output.tsv.gz",
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

    with open(files["output"]) as f:
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

    with open(files["output"]) as f:
        output_content = f.read()
    expected_output_file = os.path.join(os.path.dirname(__file__), "data", "combine", "get_variant_map.output.tsv.gz")

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content


def _run_genomic_variants(runner, files, tmp_path, ref_strand, alt_strand, ref_pos="[100]", alt_pos="[100]"):
    import pandas as pd

    variant_map = pd.read_csv(
        os.path.join(os.path.dirname(__file__), "data", "combine", "get_variant_map.output.tsv.gz"),
        sep="\t",
        index_col=0,
    )
    variant_map = variant_map.groupby(level=0).first().dropna(subset=["REF", "ALT"])
    spdi, row = variant_map.index[0], variant_map.iloc[0]

    counts_file = tmp_path / "counts.tsv.gz"
    counts = pd.DataFrame(
        [
            {"barcode": f"{prefix}{i:03d}", "oligo_name": oligo, "dna_count_1": 10 + i, "rna_count_1": 20 + i}
            for prefix, oligo in (("AAAA", row["REF"]), ("CCCC", row["ALT"]))
            for i in range(3)
        ]
    )
    counts.to_csv(counts_file, sep="\t", index=False)

    design = pd.read_csv(files["input"]["sequence_design"], sep="\t", dtype=str, keep_default_na=False)
    names = design["name"].str.replace(r"[\s\[\]]", "_", regex=True)
    design.loc[names == row["REF"], "strand"] = ref_strand
    design.loc[names == row["ALT"], "strand"] = alt_strand
    design.loc[names == row["REF"], "variant_pos"] = ref_pos
    design.loc[names == row["ALT"], "variant_pos"] = alt_pos
    design_file = tmp_path / "design.tsv.gz"
    design.to_csv(design_file, sep="\t", index=False)

    stats = pd.DataFrame(
        {"logFC": [1.0], "CI.L": [0.5], "CI.R": [1.5], "P.Value": [0.01], "adj.P.Val": [0.02], "B": [1.0]},
        index=[spdi],
    )
    stats_file = tmp_path / "stats.tsv"
    stats.to_csv(stats_file, sep="\t")

    result = runner.invoke(
        cli,
        [
            "combine",
            "get-reporter-genomic-variants",
            "--input",
            str(counts_file),
            "--sequence-design",
            str(design_file),
            "--statistics",
            str(stats_file),
            "--output-reporter-genomic-variants",
            files["output"],
        ],
    )
    assert result.exit_code == 0, result.output
    with open(files["output"]) as f:
        lines = f.read().splitlines()
    assert len(lines) == 1
    return lines[0].split("\t")


@pytest.mark.parametrize(
    "ref_strand,alt_strand,expected",
    [("+", "+", "+"), ("-", "-", "-"), (".", ".", "."), ("+", "-", "."), ("-", "+", "."), ("+", ".", "."), (".", "-", ".")],
)
def test_combine_get_reporter_genomic_variants_strand(runner, files, tmp_path, ref_strand, alt_strand, expected):
    fields = _run_genomic_variants(runner, files, tmp_path, ref_strand, alt_strand)
    assert fields[5] == expected


@pytest.mark.parametrize(
    "ref_pos,alt_pos,expected",
    [
        ("[100]", "[100]", "100"),
        ("[100]", "[102]", "-1"),
        ("[100, 102]", "[100, 102]", "-1"),
        ("[100, 100]", "[100]", "100"),
    ],
)
def test_combine_get_reporter_genomic_variants_variant_pos(runner, files, tmp_path, ref_pos, alt_pos, expected):
    fields = _run_genomic_variants(runner, files, tmp_path, "+", "+", ref_pos, alt_pos)
    assert fields[16] == expected


@pytest.mark.parametrize(
    "ref_strand,alt_strand,expected",
    [("+", "+", "100"), ("-", "-", "100"), (".", ".", "100"), ("+", "-", "-1"), ("+", ".", "-1"), (".", "-", "-1")],
)
def test_combine_get_reporter_genomic_variants_variant_pos_mixed_strand(
    runner, files, tmp_path, ref_strand, alt_strand, expected
):
    fields = _run_genomic_variants(runner, files, tmp_path, ref_strand, alt_strand)
    assert fields[16] == expected
