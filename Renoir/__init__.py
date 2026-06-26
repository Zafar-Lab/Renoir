from Renoir.renoir import *
from Renoir.downstream import *

try:
    from importlib.metadata import version, PackageNotFoundError

    __version__ = version("Renoir")
except PackageNotFoundError:
    __version__ = "unknown"
