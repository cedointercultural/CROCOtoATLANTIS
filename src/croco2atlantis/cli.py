"""Command line entry point."""

import argparse
import sys

from . import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="croco2atlantis",
        description="Convert CROCO output into ATLANTIS hydrographic forcing files.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.parse_args(argv)
    print("croco2atlantis is under development: the conversion is not implemented yet.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
