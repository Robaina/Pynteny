#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions and classes for general purposes
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
import sys
import tarfile
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import requests
from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


class CommandArgs:
    """Base class to hold command line arguments."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class ConfigParser:
    """Handle Pynteny configuration file."""

    def __init__(self, config_file: Path) -> None:
        self._config_file = Path(config_file)
        self._config = self.get_config()

    @classmethod
    def get_default_config(cls):
        """Initialize ConfigParser with default config file."""
        return cls(cls.initialize_config_file())

    @staticmethod
    def initialize_config_file() -> Path:
        """Initialize empty config file.

        Returns:
            Path: path to generated config file.
        """
        config_file = Path(__file__).parent / "config.json"
        if not config_file.exists():
            config = {
                "database_dir": "",
                "upack_PGAP_database": False,
                "data_downloaded": False,
                "PGAP_database": "",
                "PGAP_meta_file": "",
                "streamlit_process": "",
            }
            with open(config_file, "w", encoding="UTF-8") as f:
                json.dump(config, f, indent=4)
        return config_file

    def get_config_path(self) -> Path:
        """Show config file path."""
        return self._config_file

    def get_config(self) -> dict:
        """Load config file.

        Returns:
            dict: dict containing fields and values of config file.
        """
        with open(self._config_file, "r", encoding="UTF-8") as file:
            config = json.loads(file.read())
        return config

    def write_config(self) -> None:
        """Write config dict to file."""
        with open(self._config_file, "w", encoding="UTF-8") as f:
            json.dump(self._config, f, indent=4)

    def update_config(self, key: str, value: str) -> None:
        """Update config file

        Args:
            key (str): config file key name to be updated.
            value (str): new value.
        """
        self._config[key] = value
        self.write_config()

    def get_field(self, key: str) -> str:
        """Get field from config file.

        Args:
            key (str): key name to get the value from.

        Returns:
            str: key value.
        """
        return self._config[key]


def set_default_output_path(
    input_path: Path,
    tag: str = None,
    extension: str = None,
    only_filename: bool = False,
    only_basename: bool = False,
    only_dirname: bool = False,
) -> Path:
    """Utility function to generate a default path to output file
    or directory based on an input file name and path.

    Args:
        input_path (Path): path to input file.
        tag (str, optional): text tag to be added to file name. Defaults to None.
        extension (str, optional): change input file extension with this one. Defaults to None.
        only_filename (bool, optional): output only default filename. Defaults to False.
        only_basename (bool, optional): output only default basename (no extension). Defaults to False.
        only_dirname (bool, optional): output only path to default output directory. Defaults to False.

    Returns:
        Path: a path or name to a default output file.
    """
    input_path = Path(input_path)
    dirname = input_path.parent
    fname, ext = input_path.stem, input_path.suffix
    if extension is None:
        extension = ext
    if tag is None:
        tag = ""
    default_file = f"{fname}{tag}{extension}"
    if only_basename:
        return Path(fname)
    if only_filename:
        return Path(default_file)
    if only_dirname:
        return Path(dirname)
    else:
        return Path(os.path.abspath(os.path.join(dirname, default_file)))


def terminal_execute(
    command_str: str,
    suppress_shell_output=False,
    work_dir: Path = None,
    return_output=False,
) -> subprocess.STDOUT:
    """Execute given command in terminal through Python.

    Args:
        command_str (str): terminal command to be executed.
        suppress_shell_output (bool, optional): suppress shell output. Defaults to False.
        work_dir (Path, optional): change working directory. Defaults to None.
        return_output (bool, optional): whether to return execution output. Defaults to False.

    Returns:
        subprocess.STDOUT: subprocess output.
    """
    if suppress_shell_output:
        suppress_code = ">/dev/null 2>&1"
        command_str = f"{command_str} {suppress_code}"
    output = subprocess.run(
        command_str, shell=True, cwd=work_dir, capture_output=return_output
    )
    return output


def parallelize_over_input_files(
    callable,
    input_list: list,
    processes: int = None,
    max_tasks_per_child=10,
    **callable_kwargs,
) -> None:
    """Parallelize callable over a set of input objects using a pool
    of workers. Inputs in input list are passed to the first argument
    of the callable. Additional callable named arguments may be passed.

    Args:
        callable (_type_): function to be run.
        input_list (list): list of inputs to callable.
        n_processes (int, optional): maximum number of threads.
            Defaults to all minus one.
        max_tasks_per_child (int, optional): maximum number of tasks per child
            process before is reset. Defaults to 10.
    """
    if processes is None:
        processes = os.cpu_count - 1
    p = Pool(processes, maxtasksperchild=max_tasks_per_child)
    p.map(partial(callable, **callable_kwargs), input_list)
    p.close()
    p.join()


def download_file(url: str, output_file: Path) -> None:
    """Download file from url

    Args:
        url (str): url where file to be downloaded
        output_file (Path): path to downloaded file
    """
    output_file = Path(output_file)
    with requests.get(url, stream=True) as r:
        total_length = int(r.headers.get("Content-Length"))
        with tqdm.wrapattr(r.raw, "read", total=total_length, desc="") as raw:
            with open(output_file, "wb") as output:
                shutil.copyfileobj(raw, output)


def is_tar_file(tar_file: Path) -> bool:
    """Check whether file is tar-compressed.

    Args:
        tar_file (Path): path to file.

    Returns:
        bool: whether file is compressed or not.
    """
    tar_file = Path(tar_file)
    return tar_file.is_file() and tarfile.is_tarfile(tar_file.as_posix())


def extract_tar_file(tar_file: Path, dest_dir: Path = None) -> None:
    """Extract tar or tar.gz files to dest_dir.

    Args:
        tar_file (Path): path to tar file.
        dest_dir (Path, optional): path to destination directory
            to store the uncompressed file. Defaults to None.
    """
    tar_file = Path(tar_file)
    if dest_dir is None:
        dest_dir = "."
    if (tar_file.as_posix().endswith("tar.gz")) or (
        tar_file.as_posix().endswith("tgz")
    ):
        tar = tarfile.open(tar_file, "r:gz")
        tar.extractall(path=dest_dir)
        tar.close()
    elif tar_file.as_posix().endswith("tar"):
        tar = tarfile.open(tar_file, "r:")
        tar.extractall(path=dest_dir)
        tar.close()
    else:
        logger.error("Input is not a tar file")
        sys.exit(1)


def list_tar_dir(tar_dir: Path) -> list:
    """List files within tar or tar.gz directory.

    Args:
        tar_dir (Path): path to directory containing tar files.

    Returns:
        list: list of tar files.
    """
    tar_dir = Path(tar_dir)
    with tarfile.open(tar_dir, "r") as tar_obj:
        files = tar_obj.getnames()
    return files


def flatten_directory(directory: Path) -> None:
    """Flatten directory, i.e. remove all subdirectories and
    copy all files to the top level directory.

    Args:
        directory (Path): path to directory.
    """
    directory = Path(directory)
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


def is_right_list_nested_type(list_object: list, inner_type: type) -> bool:
    """Check if all elements in list are of the same type.

    Args:
        list_object (list): list containing elements.
        inner_type (type): type to be checked.

    Returns:
        bool: whether list contains elements of the same specified type.
    """
    return all(isinstance(x, inner_type) for x in list_object)
