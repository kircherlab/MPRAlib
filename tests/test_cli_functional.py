import pytest
from click.testing import CliRunner
from mpralib.cli import cli


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


def test_functional_group_help(runner):
    result = runner.invoke(cli, ["functional", "--help"])
    assert result.exit_code == 0
    assert "General functionality." in result.output
