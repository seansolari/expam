import gzip


COMPRESSION_EXTNS = ['.tar.gz', '.tar', '.gz']
DEFAULT_MODE = "rb"
DEFAULT_OPENER = open
COMP_PARSE = {
    ".tar.gz": {"mode": "rb", "opener": gzip.open},
    ".gz": {"mode": "rb", "opener": gzip.open}
}

from .main import ExpamOptions, clear_logs, CommandGroup, PlotLogs
