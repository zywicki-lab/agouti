import sys


def eprint(*args, **kwargs):
    """Print to stderr instead of stdout
    """
    print(*args, file=sys.stderr, **kwargs)
