#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to preprocess sequence databases

1. Remove illegal characters from peptide sequences
2. Remove illegal symbols from file paths
3. Relabel fasta records and make dictionary with old labels
"""

from __future__ import annotations

import logging
import os
import sys
import tempfile
from pathlib import Path
from typing import TextIO

import pyfastx
from Bio import SeqFeature, SeqIO, SeqRecord

import pynteny.utils as utils
import pynteny.wrappers as wrappers

logger = logging.getLogger(__name__)


def remove_stop_sodon_signals(record_seq: str) -> str:
    """Remove stop codon signals from peptide sequence

    Args:
        record_seq (str): peptide sequence.

    Returns:
        str: a peptide sequence without stop codon symbols.
    """
    return record_seq.replace("*", "")


def is_legit_peptide_sequence(record_seq: str) -> bool:
    """Assert that peptide sequence only contains valid symbols.

    Args:
        record_seq (str): peptide sequence.

    Returns:
        bool: whether peptide sequence only contains legit symbols.
    """
    aas = {
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
        "*",
    }
    seq_symbols = {s.upper() for s in record_seq}
    return seq_symbols.issubset(aas)


def is_legit_DNA_sequence(record_seq: str) -> bool:
    """Assert that DNA sequence only contains valid symbols.

    Args:
        record_seq (str): nucleotide sequence.

    Returns:
        bool: whether nucleotide sequence only contains legit symbols.
    """
    nts = {"A", "G", "T", "C", "N"}
    seq_symbols = {s.upper() for s in record_seq}
    return seq_symbols.issubset(nts)


class FASTAmerger:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)

    def prepend_filename_to_record_names(self, output_dir: Path) -> None:
        """Prepend file name to each record label in fasta file

        Args:
            output_dir (Path): _description_
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        for file in self.input_dir.iterdir():
            fasta = pyfastx.Fasta(file.as_posix(), build_index=False, full_name=True)
            new_file_content = [
                f">{file.stem}_{record_name}\n{record_seq}\n"
                for record_name, record_seq in fasta
            ]
            with open(output_dir / file.name, "w+", encoding="UTF-8") as outfile:
                outfile.writelines(new_file_content)

    def merge(self, output_file: Path = None, prepend_file_name: bool = False) -> None:
        """Merge input fasta files into a one (multi)fasta file.

        Args:
            output_file (Path, optional): path to ouput merged fasta file. Defaults to None.
            prepend_file_name (bool, optional): whether to add file name as genome ID to
            each record in the result merged fasta file.

        """
        if output_file is None:
            output_file = self.input_dir / "merged.fasta"
        else:
            output_file = Path(output_file)
        logger.info("Merging FASTA files in input directory")
        # cmd_str = f'awk 1 * > {output_file}'
        # # cmd_str = "find . -maxdepth 1 -type f --exec cat {} + > " + f"{output_file}"
        cmd_str = f"printf '%s\\0' * | xargs -0 cat > {output_file}"
        if prepend_file_name:
            with tempfile.TemporaryDirectory() as tempdir:
                self.prepend_filename_to_record_names(output_dir=tempdir)
                utils.terminal_execute(cmd_str, work_dir=tempdir)
        else:
            utils.terminal_execute(cmd_str, work_dir=self.input_dir)


class FASTA:
    """Handle and process fasta files."""

    def __init__(self, input_file: Path) -> None:
        """Initialize FASTA object.

        Args:
            input_file (Path): path to input fasta file.
        """
        self._input_file = Path(input_file)

    @property
    def file_path(self) -> Path:
        """Set new path to fasta file

        Args:
            new_path (Path, optional): path to fasta file. Defaults to None.

        Returns:
            Path: path to fasta file.
        """
        return self._input_file

    @file_path.setter
    def file_path(self, new_path: Path) -> None:
        self._input_file = Path(new_path)

    @classmethod
    def from_FASTA_directory(
        cls, input_dir: Path, merged_fasta: Path = None, prepend_file_name: bool = False
    ) -> FASTA:
        """Initialize FASTA class from directory of fasta files.

        Args:
            input_dir (Path): path to input directory.
            merged_fasta (Path, optional): path to output merged fasta. Defaults to None.
            prepend_file_name (bool, optional): whether to add file name as genome ID to
                each record in the result merged fasta file.

        Returns:
            FASTA: an initialized instance of class FASTA.
        """
        input_dir = Path(input_dir)
        if merged_fasta is None:
            merged_fasta = input_dir / "merged_database.fasta"
        else:
            merged_fasta = Path(merged_fasta)
        FASTAmerger(input_dir).merge(merged_fasta, prepend_file_name)
        return cls(merged_fasta)

    def remove_duplicates(
        self,
        output_file: Path = None,
        export_duplicates: bool = False,
        point_to_new_file: bool = True,
    ) -> None:
        """Removes duplicate entries (either by sequence or ID) from fasta.

        Args:
            output_file (Path, optional): path to output fasta file. Defaults to None.
            export_duplicates (bool, optional): whether duplicated records are exported to a file. Defaults to False.
            point_to_new_file (bool, optional): whether FASTA object should point to the newly generated file. Defaults to True.

        Yields:
            None: None
        """
        if output_file is None:
            output_file = (
                Path(self._input_file.parent)
                / f"{self._input_file.stem}_noduplicates{self._input_file.suffix}"
            )
        else:
            output_file = Path(output_file)
        wrappers.run_seqkit_nodup(
            input_fasta=self._input_file,
            output_fasta=output_file,
            export_duplicates=export_duplicates,
        )
        if point_to_new_file:
            self.file_path = output_file

    def remove_corrupted_sequences(
        self,
        output_file: Path = None,
        is_peptide: bool = True,
        keep_stop_codon: bool = False,
        point_to_new_file: bool = True,
    ) -> None:
        """Filter out (DNA or peptide) sequences containing illegal characters.

        Args:
            output_file (Path, optional): path to output fasta file. Defaults to None.
            is_peptide (bool, optional): select if input is a peptide sequence, otherwise taken as nucleotide. Defaults to True.
            keep_stop_codon (bool, optional): whether to keep the stop codon in the peptide sequence. Defaults to False.
            point_to_new_file (bool, optional): whether FASTA object should point to the newly generated file. Defaults to True.
        """
        dirname = self._input_file.parent
        fname, ext = self._input_file.stem, self._input_file.suffix
        if output_file is None:
            output_file = Path(dirname) / f"{fname}_modified{ext}"
        else:
            output_file = Path(output_file)
        if is_peptide:
            isLegitSequence = is_legit_peptide_sequence
        else:
            isLegitSequence = is_legit_DNA_sequence

        fasta = pyfastx.Fasta(
            self.file_path.as_posix(), build_index=False, full_name=True
        )
        with open(output_file, "w+", encoding="UTF-8") as outfile:
            for record_name, record_seq in fasta:
                if is_peptide and (not keep_stop_codon):
                    record_seq = remove_stop_sodon_signals(record_seq)
                if isLegitSequence(record_seq):
                    outfile.write(f">{record_name}\n{record_seq}\n")
        if point_to_new_file:
            self.file_path = output_file

    def filter_by_IDs(
        self, record_ids: list, output_file: Path = None, point_to_new_file: bool = True
    ) -> None:
        """Filter records in fasta file matching provided IDs.

        Args:
            record_ids (list): list of record IDs to keep of original fasta file.
            output_file (Path, optional): path to output filtered fasta file. Defaults to None.
            point_to_new_file (bool, optional): whether FASTA object should point to the newly generated file. Defaults to True.
        """
        if output_file is None:
            output_file = (
                Path(self._input_file.parent)
                / f"{self._input_file.stem}_filtered{self._input_file.suffix}"
            )
        else:
            output_file = Path(output_file)
        with tempfile.NamedTemporaryFile(mode="w+t") as tmp_ids:
            tmp_ids.writelines("\n".join(record_ids))
            tmp_ids.flush()
            tmp_ids_path = tmp_ids.name
            cmd_str = (
                f"seqkit grep -i -f {tmp_ids_path} {self._input_file} -o {output_file}"
            )
            utils.terminal_execute(cmd_str, suppress_shell_output=True)
        if point_to_new_file:
            self.file_path = output_file

    def split_by_contigs(self, output_dir: Path = None) -> None:
        """Split large fasta file into several ones containing one contig each.

        Args:
            output_dir (Path, optional): _description_. Defaults to None.
        """
        if output_dir is None:
            output_dir = (
                Path(self._input_file.parent) / "split_" + self._input_file.name
            )
        else:
            output_dir = Path(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        contigs = pyfastx.Fasta(
            self.file_path.as_posix(), build_index=False, full_name=True
        )
        for contig_name, seq in contigs:
            output_file = (
                output_dir / f"{contig_name.split(' ')[0]}{self._input_file.suffix}"
            )
            with open(output_file, "w+", encoding="UTF-8") as outfile:
                outfile.write(f">{contig_name}\n")
                outfile.write(seq + "\n")

    def filter_by_minimum_length(
        self, min_length: int, output_file: Path = None, point_to_new_file: bool = True
    ) -> None:
        """Filter records in fasta file by minimum length.

        Args:
            min_length (int): minimal length of sequences to be kept in filtered fasta file.
            output_file (Path, optional): path to output filtered fasta file. Defaults to None.
            point_to_new_file (bool, optional): whether FASTA object should point to the newly generated file. Defaults to True.
        """
        if output_file is None:
            output_file = (
                Path(self._input_file.parent)
                / f"{self._input_file.stem}_minlength{self._input_file.suffix}"
            )
        else:
            output_file = Path(output_file)
        fasta = pyfastx.Fasta(
            self.file_path.as_posix(), build_index=False, full_name=True
        )
        with open(output_file, "w+", encoding="UTF-8") as outfile:
            for record_name, record_seq in fasta:
                if len(record_seq) >= min_length:
                    outfile.write(f">{record_name}\n{record_seq}\n")
        if point_to_new_file:
            self.file_path = output_file

    def add_prefix_to_records(
        self, prefix: str, output_file: Path = None, point_to_new_file: bool = True
    ) -> None:
        """Add prefix to sequence records in FASTA

        Args:
            prefix (str): prefix to be added.
            output_file (Path, optional): path to output filtered fasta file. Defaults to None.
            point_to_new_file (bool, optional): whether FASTA object should point to the newly generated file. Defaults to True.
        """
        if output_file is None:
            output_file = (
                Path(self._input_file.parent)
                / f"{self._input_file.stem}_prefixed{self._input_file.suffix}"
            )
        else:
            output_file = Path(output_file)
        fasta = pyfastx.Fasta(
            self.file_path.as_posix(), build_index=False, full_name=True
        )
        prefix = prefix.strip("_")
        with open(output_file, "w+", encoding="UTF-8") as outfile:
            for record_name, record_seq in fasta:
                outfile.write(f">{prefix}_{record_name}\n{record_seq}\n")
        if point_to_new_file:
            self.file_path = output_file


class LabelledFASTA(FASTA):
    """Tools to add and parse FASTA with positional info on record tags"""

    @classmethod
    def from_prodigal_output(
        cls, prodigal_faa: Path, output_file: Path = None
    ) -> LabelledFASTA:
        """Instantiate class from prodigal output file.
        Extract positional gene info from prodigal output and export to fasta file.

        Args:
            prodigal_faa (Path): path to prodigal output file containing peptide sequences
            output_file (Path, optional): path to output labelled fasta file. Defaults to None.

        Returns:
            LabelledFASTA: object containing the labelled peptide database.
        """
        number_prodigal_record_fields = 9
        prodigal_faa = Path(prodigal_faa)
        if output_file is None:
            output_file = prodigal_faa.parent / f"{prodigal_faa.stem}_longlabels.fasta"
        else:
            output_file = Path(output_file)
        data = pyfastx.Fasta(prodigal_faa.as_posix(), build_index=False, full_name=True)
        with open(output_file, "w+", encoding="UTF-8") as outfile:
            for record_name, record_seq in data:
                name_list = record_name.split(" ")
                if len(name_list) < number_prodigal_record_fields:
                    logger.error(
                        f"Invalid prodigal header format for record: {record_name}"
                    )
                    sys.exit(1)
                contig = "_".join(name_list[0].split("_")[:-1])
                gene_number = name_list[0].split("_")[-1]
                start, end = name_list[2], name_list[4]
                strand = (
                    "pos"
                    if name_list[6] == "1"
                    else ("neg" if name_list[6] == "-1" else "")
                )
                header = f">{contig}_{gene_number}__{contig}_{gene_number}_{start}_{end}_{strand}"
                outfile.write(header + "\n")
                outfile.write(record_seq + "\n")
        return cls(output_file)

    @classmethod
    def from_genbank(
        cls,
        gbk_data: Path,
        output_file: Path = None,
        prefix: str = None,
        nucleotide: bool = False,
        prepend_file_name: bool = False,
    ) -> LabelledFASTA:
        """Assign gene positional info, such as contig, gene number and loci
        to each record in genbank database and return LabelledFASTA object.

        Args:
            gbk_data (Path): path to file or directory contanining genbank files
            output_file (Path, optional): path to output labelled fasta file. Defaults to None.
            prefix (str, optional): prefix for output file. Defaults to None.
            nucleotide (bool, optional): whether records corresponds to nucleotide sequences instead of peptides. Defaults to False.
            prepend_file_name (bool, optional): whether to add file name as genome ID to
                each record in the result merged fasta file.
        Returns:
            LabelledFASTA: object containing the labelled peptide database.
        """
        gbk_data = Path(gbk_data)
        if gbk_data.is_dir():
            gbk_files = [gbk_data / f for f in gbk_data.iterdir()]
        else:
            gbk_files = [gbk_data]
        gbk_name = (
            lambda gbk_file: gbk_file.stem
            if prepend_file_name
            else lambda gbk_file: None
        )
        gbk_contigs = [
            (gbk_name(gbk_file), contig)
            for gbk_file in gbk_files
            for contig in SeqIO.parse(gbk_file, "genbank")
        ]

        if output_file is None:
            output_file = Path(gbk_files[0].parent) / f"{prefix}sequence_database.fasta"
        else:
            output_file = Path(output_file)

        with open(output_file, "w+", encoding="UTF-8") as outfile:
            for gbk_file_name, gbk_contig in gbk_contigs:
                gene_counter = 0
                for feature in gbk_contig.features:
                    if "cds" in feature.type.lower():
                        gene_counter = cls.write_record(
                            gbk_contig,
                            feature,
                            outfile,
                            gene_counter,
                            gbk_file_name,
                            nucleotide,
                        )
        return cls(output_file)

    @staticmethod
    def get_label_str(
        gbk_contig: SeqRecord, feature: SeqFeature, gene_counter: int
    ) -> str:
        name = feature.qualifiers["locus_tag"][0].replace("_", ".")
        start, end, strand = (
            str(feature.location.start),
            str(feature.location.end),
            feature.location.strand,
        )
        start = start.replace(">", "").replace("<", "")
        end = end.replace(">", "").replace("<", "")
        strand_sense = "neg" if strand == -1 else ("pos" if strand == 1 else "")
        return f">{name}__{gbk_contig.name.replace('_', '')}_{gene_counter}_{start}_{end}_{strand_sense}\n"

    @staticmethod
    def write_record(
        gbk_contig: SeqRecord,
        feature: SeqFeature,
        output_file: TextIO,
        gene_counter: int,
        prefix: str = None,
        nucleotide: bool = False,
    ) -> int:
        header = LabelledFASTA.get_label_str(gbk_contig, feature, gene_counter)
        if prefix is not None:
            header = f"{prefix}_{header}"
        if (not nucleotide) and ("translation" in feature.qualifiers):
            sequence = feature.qualifiers["translation"][0]
        elif nucleotide:
            sequence = str(feature.extract(gbk_contig).seq)
        else:
            return gene_counter
        output_file.write(header)
        output_file.write(sequence + "\n")
        gene_counter += 1
        return gene_counter


class GeneAnnotator:
    """Run prodigal on assembly, predict ORFs and extract location info"""

    def __init__(self, assembly: FASTA) -> None:
        """Initialize GeneAnnotator

        Args:
            assembly_file (FASTA): fasta object containing assembled nucleotide sequences
        """
        self._assembly = assembly

    def annotate(
        self,
        processes: int = None,
        metagenome: bool = True,
        output_file: Path = None,
        prodigal_args: str = None,
    ) -> LabelledFASTA:
        """Run prodigal on assembly and export single fasta file with peptide ORFs predictions

        Args:
            processes (int, optional): maximum number of threads. Defaults to all minus one.
            metagenome (bool, optional): whether assembled sequences correspond to metagenomic data. Defaults to True.
            output_file (Path, optional): path to output fasta file. Defaults to None.
            prodigal_args (str, optional): additional arguments to be passed to prodigal CLI. Defaults to None.

        Returns:
            LabelledFASTA: object containing the labelled peptide database.
        """
        if processes is None:
            processes = os.cpu_count() - 1
        if output_file is None:
            output_file = (
                self._assembly.file_path.parent
                / f"{self._assembly.file_path.stem}_annotated.faa"
            )
        else:
            output_file = Path(output_file)
        with tempfile.TemporaryDirectory() as contigs_dir, tempfile.TemporaryDirectory() as prodigal_dir, tempfile.NamedTemporaryFile() as temp_fasta:
            contigs_dir = Path(contigs_dir)
            prodigal_dir = Path(prodigal_dir)
            logger.info("Running prodigal on assembly data")
            self._assembly.split_by_contigs(contigs_dir)
            utils.parallelize_over_input_files(
                wrappers.run_prodigal,
                input_list=list(contigs_dir.iterdir()),
                n_processes=processes,
                output_dir=prodigal_dir,
                output_format="fasta",
                metagenome=metagenome,
                additional_args=prodigal_args,
            )
            FASTAmerger(prodigal_dir).merge(
                Path(temp_fasta.name), prepend_file_name=False
            )
            return LabelledFASTA.from_prodigal_output(
                Path(temp_fasta.name), output_file
            )


class Database:
    """_Sequence database constructor"""

    def __init__(self, data: Path) -> None:
        """Initialize Database object

        Args:
            data (Path): path to either assembly fasta file (or a
                         directory containing assembly fasta files) or
                         a genbank file containing ORF annotations (or a
                         directory containing genbank files)

        Raises:
            FileNotFoundError: if file or directory doesn't exist
        """
        self._data = Path(data)
        if not self._data.exists():
            raise FileNotFoundError(f"{self._data} does not exist")
        if self._data.is_dir():
            self._data_files = [f for f in self._data.iterdir()]
        else:
            self._data_files = [self._data]

    @staticmethod
    def is_fasta(filename: Path) -> bool:
        """Check if file is in fasta format

        Args:
            filename (Path): path to input file

        Returns:
            bool: whether the file is in fasta format
        """
        filename = Path(filename)
        if filename.exists():
            fasta = list(SeqIO.parse(str(filename), "fasta"))
            return any(fasta)
        else:
            return False

    @staticmethod
    def is_gbk(filename: Path) -> bool:
        """Check if file is in genbank format

        Args:
            filename (Path): path to input file

        Returns:
            _bool: whether the file is in genbank format
        """
        filename = Path(filename)
        if filename.exists():
            gbk = list(SeqIO.parse(str(filename), "genbank"))
            return any(gbk)
        else:
            return False

    def build(
        self,
        seq_prefix: str = None,
        prepend_file_name: bool = False,
        output_file: Path = None,
    ) -> LabelledFASTA:
        """Build database from data files.

        Args:
            prefix (str, optionall): prefix to be added to each sequence in database.
                Defaults to "".
            prepend_file_name (bool, optional): whether to add file name as genome ID to
                each record in the result merged fasta file.
            output_file (Path, optional): path to output file. Defaults to None.

        Returns:
            LabelledFASTA: object containing the labelled peptide database.
        """
        if output_file is None:
            output_file = self._data.parent / f"{self._data.stem}_labelled.faa"
        else:
            output_file = Path(output_file)
        if self.is_fasta(self._data_files[0]):
            if self._data.is_dir():
                assembly_fasta = FASTA.from_FASTA_directory(
                    self._data, prepend_file_name=prepend_file_name
                )
            else:
                assembly_fasta = FASTA(self._data)
            logger.info("Translating and annotating assembly data.")
            labelled_database = GeneAnnotator(assembly_fasta).annotate(
                output_file=output_file
            )
        elif self.is_gbk(self._data_files[0]):
            logger.info("Parsing GenBank data.")
            labelled_database = LabelledFASTA.from_genbank(
                self._data, output_file=output_file, prepend_file_name=prepend_file_name
            )
        else:
            logging.error(f"{self._data} is not a valid FASTA or genbank file")
            sys.exit(1)
        if seq_prefix is not None:
            labelled_database.add_prefix_to_records(seq_prefix, output_file)
            labelled_database = LabelledFASTA(labelled_database.file_path)
        return labelled_database
