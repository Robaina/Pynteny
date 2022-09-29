#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions and classes for general purposes
"""

from __future__ import annotations
import os
import sys
import logging
import shutil
import subprocess
import tarfile
import json
from pathlib import Path
from functools import partial
from multiprocessing import Pool

logger = logging.getLogger(__name__)


class CommandArgs():
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

        
class ConfigParser():
    """
    Base configuration class
    """
    def __init__(self, config_file: Path) -> None:
        self._config_file = Path(config_file)
        self._config = self.get_config()
    
    @classmethod
    def get_default_config(cls):
        """
        Initialize ConfigParser with default config file
        """
        return cls(cls.initialize_config_file())
    
    @staticmethod
    def initialize_config_file() -> Path:
        """
        initialize empty config file
        """
        config_file = Path(Path(Path(__file__).parent).parent) / "config.json"
        if not config_file.exists():
            config = {
                "database_dir": "",
                "upack_PGAP_database": False,
                "data_downloaded": False,
                "PGAP_database": "",
                "PGAP_meta_file": "",
                "streamlit_process": ""
            }
            with open(config_file, 'w') as f:
                json.dump(config, f, indent=4)
        return config_file
    
    def get_config_path(self) -> Path:
        """
        Show config file path
        """
        return self._config_file

    def get_config(self) -> dict:
        """
        load config file
        """
        with open(self._config_file, 'r') as file:
            config = json.loads(file.read())
        return config

    def write_config(self) -> None:
        """
        write config file
        """
        with open(self._config_file, 'w') as f:
            json.dump(self._config, f, indent=4)

    def update_config(self, key: str, value: str) -> None:
        """
        update config file
        """
        self._config[key] = value
        self.write_config()

    def get_field(self, key: str) -> str:
        """
        get field from config file
        """
        return self._config[key]


def setDefaultOutputPath(input_path: Path, tag: str = None,
                         extension: str = None,
                         only_filename: bool = False,
                         only_basename: bool = False,
                         only_dirname: bool = False) -> str:
    """
    Get default path to output file or directory
    """
    dirname = input_path.parent
    fname, ext = input_path.stem, input_path.suffix
    if extension is None:
        extension = ext
    if tag is None:
        tag = ''
    default_file = f'{fname}{tag}{extension}'
    if only_basename:
        return Path(fname)
    if only_filename:
        return Path(default_file)
    if only_dirname:
        return Path(dirname)
    else:
        return Path(os.path.abspath(os.path.join(dirname, default_file)))

def terminalExecute(command_str: str,
                    suppress_shell_output=False,
                    work_dir: Path = None,
                    return_output=False) -> subprocess.STDOUT:
    """
    Execute given command in terminal through Python
    """
    if suppress_shell_output:
        suppress_code = ">/dev/null 2>&1"
        command_str = f"{command_str} {suppress_code}"
    output = subprocess.run(
        command_str, shell=True,
        cwd=work_dir, capture_output=return_output
        )
    return output

def parallelizeOverInputFiles(callable,
                              input_list: list,
                              n_processes: int = None,
                              **callable_kwargs) -> None: 
    """
    Parallelize callable over a set of input objects using a pool 
    of workers. Inputs in input list are passed to the first argument
    of the callable.
    Additional callable named arguments may be passed.
    """
    if n_processes is None:
        n_processes = os.cpu_count - 1
    p = Pool(processes=n_processes)
    p.map(partial(callable, **callable_kwargs), input_list)
    p.close()
    p.join()

def isTarFile(tar_file: Path) -> bool:
    tar_file_str = tar_file.as_posix()
    return (
        (tar_file_str.endswith('tar.gz')) or 
        (tar_file_str.endswith('tgz')) or
        (tar_file_str.endswith('tar'))
        )

def extractTarFile(tar_file: Path, dest_dir: Path = None) -> None:
    """
    Extract tar or tar.gz files to dest_dir
    """ 
    if dest_dir is None:
        dest_dir = '.'
    if (tar_file.as_posix().endswith('tar.gz')) or (tar_file.as_posix().endswith('tgz')):
        tar = tarfile.open(tar_file, 'r:gz')
        tar.extractall(path=dest_dir)
        tar.close()
    elif tar_file.as_posix().endswith('tar'):
        tar = tarfile.open(tar_file, 'r:')
        tar.extractall(path=dest_dir)
        tar.close()
    else:
        logger.error("Input is not a tar file")
        sys.exit(1)

def listTarDir(tar_dir: Path) -> list:
    """
    List files within tar or tar.gz directory
    """
    with tarfile.open(tar_dir, 'r') as tar_obj:
        files = tar_obj.getnames()
    return files

def flattenDirectory(directory: Path) -> None:
    """
    Flatten directory, i.e. remove all subdirectories and
    copy all files to the top level directory
    """
    directory = directory.as_posix()
    for dirpath, _, filenames in os.walk(directory, topdown=False):
        for filename in filenames:
            i = 0
            source = os.path.join(dirpath, filename)
            target = os.path.join(directory, filename)
            while os.path.exists(target):
                i += 1
                file_parts = os.path.splitext(os.path.basename(filename))

                target = os.path.join(
                    directory,
                    file_parts[0] + "_" + str(i) + file_parts[1],
                )
            shutil.move(source, target)
        if dirpath != directory:
            os.rmdir(dirpath)

def isRightListNestedType(list_object: list, inner_type: type) -> bool:
    """
    check if all elements in list are of right type
    """
    return all(isinstance(x, inner_type) for x in list_object)