class EnhydraConfigError(ValueError):
    """Raised when a required configuration parameter is missing or invalid."""
    pass


class EnhydraIOError(IOError):
    """Raised when an input file or directory is missing or inaccessible."""
    pass


class EnhydraToolError(FileNotFoundError):
    """Raised when an external tool (mafft, trimal) is missing or not executable."""
    pass
