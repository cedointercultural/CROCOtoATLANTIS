import importlib

import pytest

MODULES = [
    "croco2atlantis",
    "croco2atlantis.readers.croco",
    "croco2atlantis.geometry.bgm",
    "croco2atlantis.vertical",
    "croco2atlantis.horizontal",
    "croco2atlantis.fluxes",
    "croco2atlantis.exchanges",
    "croco2atlantis.writers.atlantis_nc",
    "croco2atlantis.validate",
    "croco2atlantis.cli",
]


@pytest.mark.parametrize("name", MODULES)
def test_module_imports(name):
    importlib.import_module(name)


def test_cli_version(capsys):
    from croco2atlantis import __version__
    from croco2atlantis.cli import main

    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out
