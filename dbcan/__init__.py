from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("dbcan")

except PackageNotFoundError:
    try:
        from ._version import version as __version__
    except ModuleNotFoundError:
        raise RuntimeError("dbcan is not installed. Please install it with `pip install dbcan`. ")
