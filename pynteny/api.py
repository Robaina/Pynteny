#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes to facilitate usage within Python scripts / Notebooks
"""

from pathlib import Path
from importlib import metadata

from pynteny.utils import CommandArgs
from pynteny.filter import SyntenyHits
from pynteny.subcommands import (
    synteny_search, build_database,
    download_hmms, parse_gene_ids, get_citation
    )

meta = metadata.metadata("pynteny")
__version__ = meta["Version"]
__author__ = meta["Author"]


class Command():
    """Parent class for Pynteny command
    """
    def __init__(self):
        """Parent class for Pynteny command
        """

    def _repr_html_(self):
        """Executed by Jupyter to print Author and version in html
        """
        return """
        <table>
            <tr>
                <td><strong>Pynteny version</strong></td><td>{version}</td>
            </tr><tr>
                <td><strong>Author</strong></td><td>{author}</td>
            </tr>
        </table>
        """.format(version=__version__,
                   author=__author__)

    def update(self, argname: str, value: str) -> None:
        """Update argument value in pynteny command.

        Args:
            argname (str): argument name to be updated.
            value (str): new argument value.
        """
        setattr(self._args, argname, value)
    
    @staticmethod
    def cite():
        """Display Pynteny citation
        """
        args = CommandArgs(
            version=__version__,
            author=__author__
            )
        citation = get_citation(args, silent=True)
        print(citation)


class Search(Command):
    """Search command object.
    """
    def __init__(
        self,
        data: Path,
        synteny_struc: str,
        gene_ids: bool = False,
        unordered: bool = False,
        reuse: bool = False,
        hmm_dir: Path = None,
        hmm_meta: Path = None,
        outdir: Path = None,
        prefix: str = "",
        hmmsearch_args: str = None,
        logfile: Path = None,
        processes: int = None):
        """Query sequence database for HMM hits arranged in provided synteny structure.

        Args:
            data (Path): path to input labelled database.
            synteny_struc (str): a str describing the desired synteny structure,
                structured as follows:

                '>hmm_a N_ab hmm_b bc <hmm_c'

                where N_ab corresponds to the maximum number of genes separating 
                gene found by hmm_a and gene found by hmm_b, and hmm_ corresponds 
                to the name of the hmm as provided in the keys of hmm_hits.
                More than two hmms can be concatenated. Strand location may be
                specificed by using '>' for sense and '<' for antisense. 
            gene_ids (bool, optional): whether gene symbols are used in synteny 
                string instead of HMM names. Defaults to False.
            unordered (bool, optional): whether the HMMs should be arranged in the
                exact same order displayed in the synteny_structure or in 
                any order If ordered, the filters would filter collinear rather
                than syntenic structures. Defaults to False.
            reuse (bool, optional): if True then HMMER3 won't be run again for HMMs already
            searched in the same output directory. Defaults to False.
            hmm_dir (Path, optional): path to directory containing input HMM files.
                Defaults to None, in which case the PGAP database is downloaded if not
                already so.
            hmm_meta (Path, optional): path to PGAP's metadata file. Defaults to None.
            outdir (Path, optional): path to output directory. Defaults to None.
            prefix (str, optional): prefix of output file names. Defaults to "".
            hmmsearch_args (str, optional): additional arguments to hmmsearch or hmmscan. Each
                element in the list is a string with additional arguments for each 
                input hmm (arranged in the same order), an element can also take a 
                value of None to avoid passing additional arguments for a specific 
                input hmm. A single string may also be passed, in which case the 
                same additional argument is passed to hmmsearch for all input hmms.
                Defaults to None. Defaults to None.
            logfile (Path, optional): path to log file. Defaults to None.
            processes (int, optional): maximum number of threads to be employed.
                Defaults to all minus one.
        """

        self._args = CommandArgs(
            data=Path(data),
            synteny_struc=synteny_struc,
            hmm_dir=Path(hmm_dir) if hmm_dir is not None else hmm_dir,
            hmm_meta=Path(hmm_meta) if hmm_meta is not None else hmm_meta,
            outdir=Path(outdir) if outdir is not None else outdir,
            prefix=prefix,
            hmmsearch_args=hmmsearch_args,
            gene_ids=gene_ids,
            logfile=Path(logfile) if logfile is not None else logfile,
            processes=processes,
            unordered=unordered,
            reuse=reuse
            )

    def parseGeneIDs(self, synteny_struc: str) -> str:
        """Parse gene IDs in synteny structure and find corresponding HMM names
        in provided HMM database.

        Args:
            synteny_struc (str): a str describing the desired synteny structure,
                structured as follows:

                '>hmm_a N_ab hmm_b bc <hmm_c'

                where N_ab corresponds to the maximum number of genes separating 
                gene found by hmm_a and gene found by hmm_b, and hmm_ corresponds 
                to the name of the hmm as provided in the keys of hmm_hits.
                More than two hmms can be concatenated. Strand location may be
                specificed by using '>' for sense and '<' for antisense. 

        Returns:
            str: parsed synteny structure in which gene symbols are replaced by 
                 corresponding HMM names.
        """
        args = CommandArgs(
            synteny_struc=synteny_struc,
            hmm_meta=self._args.hmm_meta,
            logfile=self._args.logfile
        )
        return parse_gene_ids(args)

    def run(self) -> SyntenyHits:
        """Run pynteny search
        """
        return synteny_search(self._args)


class Build(Command):
    """Build command object.
    """
    def __init__(
        self,
        data: Path,
        outfile: Path = None,
        logfile: Path = None,
        processes: int = None):
        """Translate nucleotide assembly file and assign contig and gene location 
           info to each identified ORF (using prodigal).

        Args:
            data (Path): _description_
            outfile (Path, optional): _description_. Defaults to None.
            logfile (Path, optional): _description_. Defaults to None.
            processes (int, optional): _description_. Defaults to None.
        """

        self._args = CommandArgs(
            data=Path(data),
            outfile=Path(outfile) if outfile is not None else outfile,
            logfile=Path(logfile) if logfile is not None else logfile,
            processes=processes
            )

    def run(self) -> None:
        """Run pynteny search
        """
        return build_database(self._args)


class Download(Command):
    """Download HMM database from NCBI.
    """
    def __init__(
        self,
        outdir: Path = None,
        logfile: Path = None,
        unpack: bool = False):
        """Download HMM database from NCBI.

        Args:
            outdir (Path, optional): path to ouput directory in which to store downloaded HMMs.
                Defaults to None.
            logfile (Path, optional): path to log file. Defaults to None.
            unpack (bool, optional): whether to unpack downloaded file. If False, then PGAP's database
                will be unpacked in each Pynteny session. Defaults to False.
        """

        self._args = CommandArgs(
            outdir=Path(outdir) if outdir is not None else outdir,
            logfile=Path(logfile) if logfile is not None else logfile,
            unpack=unpack
            )

    def run(self) -> None:
        """Run pynteny download
        """
        return download_hmms(self._args)