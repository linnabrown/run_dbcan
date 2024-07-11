import sys

if sys.version_info >= (3, 8):
    # For Python 3.8 and later, use the standard library import
    from importlib.metadata import PackageNotFoundError, version
else:
    # For Python 3.7 and earlier, use the backport
    from importlib_metadata import PackageNotFoundError, version

try:
    __version__ = version("dbcan")

except PackageNotFoundError:
    try:
        from ._version import version as __version__
    except ModuleNotFoundError:
        raise RuntimeError("dbcan is not installed. Please install it with `pip install dbcan`. ")
