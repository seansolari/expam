import os
import sys

from expam.sequences import COMPRESSION_EXTNS

def die(msg):
    print("ERROR\n", msg, "\n")
    sys.exit(1)

def yield_csv(path):
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                line = line.split(",")
                yield line

def ls(path, ext=None):
    # Get list of files.
    if os.path.isfile(path):
        files = [path]
    else:
        files = [
            os.path.join(path, file)
            for file in os.listdir(path)
            if not file.startswith(".")
        ]
        
        # Filter for files not subdirectories.
        files = [file for file in files if os.path.isfile(file)]

    if ext is not None:
        # Reduce list to those that have the right file extension.
        if not isinstance(ext, str) or isinstance(ext, list):
            raise TypeError("Argument `ext` must be a string of list of strings!")

        fcheck = isfiletype if isinstance(ext, str) else isfiletypes
        files = [file for file in files if fcheck(file, ext)]

    return files

def isfiletype(path: str, ext: str, check_compressed=False):
    if check_compressed:
        path = removecompressedextension(path)

    if path.endswith(ext):
        return True
    else:
        return False

def isfiletypes(path, ext, check_compressed=False):
    if check_compressed:
        path = removecompressedextension(path)

    return any(isfiletype(path, e) for e in ext)

def removecompressedextension(path):
    for cmp_ext in COMPRESSION_EXTNS:
        path = path.rstrip(cmp_ext)

    return path

def make_path_absolute(dir, base_dir):
    if base_dir not in dir:
        dir = os.path.join(base_dir, dir)

    return dir

def is_hex(string):
    if string[0] != "#" or len(string) != 7:
        return False

    try:
        int(string[1:].upper(), 16)
    except TypeError:
        return False

    return True

def parse_int(param):
    INVALID_PARAM_MSG = ("Invalid parameter (%s), must be integer!" % str(param))
    if param is not None:
        try:
            # Convert to float.
            param = float(param)
        except ValueError:
            die(INVALID_PARAM_MSG)

        # Convert to int and see if the value changes.
        new_param = int(param)
        if new_param != param:
            die(INVALID_PARAM_MSG)
    else:
        new_param = None

    return new_param

def parse_float(param):
    INVALID_PARAM_MSG = ("Invalid parameter (%s), must be integer!" % str(param))
    if param is not None:
        try:
            param = float(param)
        except ValueError:
            die(INVALID_PARAM_MSG)

    return param

