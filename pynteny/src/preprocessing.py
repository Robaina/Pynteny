#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to preprocess sequence databases

1. Remove illegal characters from peptide sequences
2. Remove illegal symbols from file paths
3. Relabel fasta records and make dictionary with old labels
"""

from __future__ import annotations
import sys
import os
import logging
import tempfile
from pathlib import Path

from Bio import SeqIO
import pyfastx

import pynteny.src.utils as utils
import pynteny.src.wrappers as wrappers

logger = logging.getLogger(__name__)


class RecordSequence():
    """
    Tools to process nucleotide or peptide sequences
    """
    @staticmethod
    def removeStopCodonSignals(record_seq: str) -> str:
        return record_seq.replace('*', '')
    
    @staticmethod
    def isLegitPeptideSequence(record_seq: str) -> bool:
        """
        Assert that peptide sequence only contains valid symbols
        """
        aas = {
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'
        }
        seq_symbols = {s.upper() for s in record_seq}
        return seq_symbols.issubset(aas)
    
    @staticmethod
    def isLegitDNAsequence(record_seq: str) -> bool:
        """
        Assert that DNA sequence only contains valid symbols
        """
        nts = {'A', 'G', 'T', 'C', 'N'}
        seq_symbols = {s.upper() for s in record_seq}
        return seq_symbols.issubset(nts)


class FASTA():
    def __init__(self, input_file: Path) -> None:
        """
        Class to handle and process fasta files
        """
        self._input_file = Path(input_file)
        self._input_file_str = self._input_file.as_posix()

    def setFilePath(self, new_path: Path) -> None:
        """
        Set new path for fasta file
        """
        self._input_file = Path(new_path)
        self._input_file_str = self._input_file.as_posix()

    def getFilePath(self) -> Path:
        """
        Get path to fasta file
        """
        return self._input_file

    @classmethod
    def fromFASTAdirectory(cls, input_dir: Path,
                           merged_fasta: Path = None) -> FASTA:
        """
        Initialize FASTA object from directory of FASTA files
        """
        if merged_fasta is None:
            merged_fasta = input_dir / "merged_database.fasta"
        FASTA.mergeFASTAs(input_dir, merged_fasta)
        return cls(merged_fasta)

    @staticmethod
    def mergeFASTAs(input_dir: Path,
                    output_file: Path = None) -> None:
        """
        Merge input fasta files into a single fasta
        """
        if output_file is None:
            output_file = input_dir / "merged.fasta"
        logger.info(f"Merging FASTA files in input directory")
        cmd_str = f'awk 1 * > {output_file.as_posix()}'
        utils.terminalExecute(
            cmd_str,
            work_dir=input_dir,
            suppress_shell_output=False
            )
  
    def removeDuplicates(self,
                         output_file: Path = None,
                         export_duplicates: bool = False,
                         method: str = 'seqkit',
                         point_to_new_file: bool = True) -> None:
        """
        Removes duplicate entries (either by sequence or ID) from fasta.
        """
        if output_file is None:
            output_file = Path(self._input_file.parent) \
                / f"{self._input_file.stem}_noduplicates{self._input_file.suffix}"

        if 'bio' in method:
            seen_seqs, seen_ids = set(), set()
            def unique_records():
                for record in SeqIO.parse(self._input_file, 'fasta'):  
                    if (record.seq not in seen_seqs) and (record.id not in seen_ids):
                        seen_seqs.add(record.seq)
                        seen_ids.add(record.id)
                        yield record
            SeqIO.write(unique_records(), output_file, 'fasta')
        else:
            wrappers.runSeqKitNoDup(input_fasta=self._input_file, output_fasta=output_file,
                                    export_duplicates=export_duplicates)
        if point_to_new_file:
            self.setFilePath(output_file)

    def removeCorruptedSequences(self,
                                 output_file: Path = None,
                                 is_peptide: bool = True,
                                 keep_stop_codon: bool = False,
                                 point_to_new_file: bool = True) -> None:
        """
        Filter out (DNA or peptide) sequences containing illegal characters
        """
        dirname = self._input_file.parent
        fname, ext = self._input_file.stem, self._input_file.suffix
        if output_file is None:
            output_file = Path(dirname) / f'{fname}_modified{ext}'
        if is_peptide:
            isLegitSequence = RecordSequence.isLegitPeptideSequence
        else:
            isLegitSequence = RecordSequence.isLegitDNAsequence

        fasta = pyfastx.Fasta(self._input_file_str, build_index=False, full_name=True)
        with open(output_file, 'w') as outfile:
            for record_name, record_seq in fasta:
                if is_peptide and (not keep_stop_codon):
                    record_seq = RecordSequence.removeStopCodonSignals(record_seq)
                if isLegitSequence(record_seq):
                    outfile.write(f'>{record_name}\n{record_seq}\n')
        if point_to_new_file:
            self.setFilePath(output_file)

    def filterByIDs(self, record_ids: list,
                    output_file: Path = None,
                    point_to_new_file: bool = True) -> None:
        """
        Filter records in fasta file matching provided IDs
        """
        if output_file is None:
            output_file = Path(self._input_file.parent) \
                / f"{self._input_file.stem}_filtered{self._input_file.suffix}"
        with tempfile.NamedTemporaryFile(mode="w+t") as tmp_ids:
            tmp_ids.writelines("\n".join(record_ids))
            tmp_ids.flush()
            tmp_ids_path = tmp_ids.name
            cmd_str = f"seqkit grep -i -f {tmp_ids_path} {self._input_file} -o {output_file}"
            utils.terminalExecute(cmd_str, suppress_shell_output=True)
        if point_to_new_file:
            self.setFilePath(output_file)
        
    def splitByContigs(self, output_dir: Path = None) -> None:
        """
        Split large fasta file into several ones containing one contig each
        """
        if output_dir is None:
            output_dir = Path(self._input_file.parent) / "split_" + self._input_file.name
        else:
            output_dir = Path(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        contigs = pyfastx.Fasta(self._input_file_str, build_index=False, full_name=True)
        for contig_name, seq in contigs:
            outfile = output_dir / f"{contig_name.split(' ')[0]}{self._input_file.suffix}"
            with open(outfile, "w+") as file:
                file.write(f">{contig_name}\n")
                file.write(seq + "\n")

    def filterByMinimumLength(self, min_length: int,
                              output_file: Path = None,
                              point_to_new_file: bool  = True) -> None:
        """
        Filter records in fasta file by minimum length
        """
        if output_file is None:
            output_file = Path(self._input_file.parent) \
                / f"{self._input_file.stem}_minlength{self._input_file.suffix}"
        fasta = pyfastx.Fasta(self._input_file_str, build_index=False, full_name=True)
        with open(output_file, 'w') as outfile:
            for record_name, record_seq in fasta:
                if len(record_seq) >= min_length:
                    outfile.write(f'>{record_name}\n{record_seq}\n')
        if point_to_new_file:
            self.setFilePath(output_file)


class LabelledFASTA(FASTA):
    """
    Tools to add and parse FASTA with positional info on record tags
    """
    @classmethod
    def fromProdigalOutput(cls,
                           prodigal_faa: Path,
                           output_file: Path = None) -> LabelledFASTA:
        """
        Extract positional gene info from prodigal output and export to
        fasta file.
        """
        number_prodigal_record_fields = 9
        if output_file is None:
            output_file = Path(prodigal_faa.parent) / f"{prodigal_faa.stem}_longlabels.fasta"
        data = pyfastx.Fasta(prodigal_faa.as_posix(), build_index=False, full_name=True)
        with open(output_file, "w") as outfile:
            for record_name, record_seq in data:
                name_list = record_name.split(" ")
                if len(name_list) < number_prodigal_record_fields:
                    logger.error(f"Invalid prodigal header format for record: {record_name}")
                    sys.exit(1)
                contig = "_".join(name_list[0].split("_")[:-1])
                gene_number = name_list[0].split("_")[-1]
                start, end = name_list[2], name_list[4]
                strand = "pos" if name_list[6] == "1" else ("neg" if name_list[6] == "-1" else "")
                header = f">{contig}_{gene_number}__{contig}_{gene_number}_{start}_{end}_{strand}"
                outfile.write(header + "\n")
                outfile.write(record_seq + "\n")
        return cls(output_file)

    @classmethod
    def fromGenBankData(cls,
                        gbk_data: Path,
                        output_file: Path = None,
                        prefix: str = None,
                        nucleotide: bool = False) -> LabelledFASTA:
        """
        Assign gene positional info, such as contig, gene number and loci
        to each record in genbank database and return LabelledFASTA object.
        @paramms:
        gbk_data: path to either a gbk file or a directory containing gbk files
        nucleotide: if True then records are nucleotide sequences instead of peptides.
                    Note that this option will notably increase the computation time.
        """
        if gbk_data.is_dir():
            gbk_files = [gbk_data / f for f in gbk_data.listdir()]
        else:
            gbk_files = [gbk_data]
        gbk_contigs = [
            contig for gbk_file in gbk_files
            for contig in list(SeqIO.parse(gbk_file, 'genbank'))
            ]
        if output_file is None:
            output_file = Path(gbk_files.pop().parent) / f"{prefix}sequence_database.fasta"
        
        def get_label_str(gbk_contig, feature):
            name = feature.qualifiers["locus_tag"][0].replace('_', '.')
            start, end, strand = str(feature.location.start), str(feature.location.end), feature.location.strand
            start = start.replace(">", "").replace("<", "")
            end = end.replace(">", "").replace("<", "")
            strand_sense = "neg" if strand == -1 else ("pos" if strand == 1 else "")
            return f">{name}__{gbk_contig.name.replace('_', '')}_{gene_counter}_{start}_{end}_{strand_sense}\n"

        if nucleotide:
            def write_record(gbk_contig, feature, outfile, gene_counter):
                header = get_label_str(gbk_contig, feature)
                sequence = str(feature.extract(gbk_contig).seq)
                outfile.write(header)
                outfile.write(sequence + "\n")
                gene_counter += 1
                return gene_counter
        else:
            def write_record(gbk_contig, feature, outfile, gene_counter):
                if "translation" in feature.qualifiers:
                    header = get_label_str(gbk_contig, feature)
                    sequence = feature.qualifiers["translation"][0]
                    outfile.write(header)
                    outfile.write(sequence + "\n")
                    gene_counter += 1
                return gene_counter

        with open(output_file, "w") as outfile:
            for gbk_contig in gbk_contigs:
                gene_counter = 0
                for feature in gbk_contig.features:
                    if "cds" in feature.type.lower():
                        gene_counter = write_record(gbk_contig, feature, outfile, gene_counter)
        return cls(output_file)


class GeneAnnotator():
    def __init__(self, assembly_file: FASTA) -> None:
        """
        Run prodigal on assembly, predict ORFs
        and extract location info.
        """
        self._assembly_file = assembly_file

    def annotate(self, processes: int = None, metagenome: bool = True,
                 output_file: Path = None,
                 prodigal_args: str = None) -> LabelledFASTA:
        """
        Run prodigal on assembly and export single
        fasta file with peptide ORFs predictions
        """
        if processes is None:
            processes = os.cpu_count() - 1
        with tempfile.TemporaryDirectory() as contigs_dir,\
             tempfile.TemporaryDirectory() as prodigal_dir,\
             tempfile.NamedTemporaryFile() as temp_fasta:
            contigs_dir = Path(contigs_dir)
            prodigal_dir = Path(prodigal_dir)
            logger.info("Running prodigal on assembly data")
            self._assembly_file.splitByContigs(contigs_dir)
            utils.parallelizeOverInputFiles(
                wrappers.runProdigal, 
                input_list=list(contigs_dir.iterdir()),
                n_processes=processes,
                output_dir=prodigal_dir,
                output_format="fasta",
                metagenome=metagenome,
                additional_args=prodigal_args
            )
            LabelledFASTA.mergeFASTAs(
                prodigal_dir,
                output_file=Path(temp_fasta.name)
            )
            return LabelledFASTA.fromProdigalOutput(
                Path(temp_fasta.name),
                output_file
            )


class Database():
    def __init__(self, data: Path) -> None:
        """
        Initialize Database object.
        @paramms:
        data: path to either assembly fasta file (or a
              directory containing assembly fasta files) or
              a genbank file containing ORF annotations (or a
              directory containing genbank files)
        """
        self._data = Path(data)
        if not self._data.exists():
            raise FileNotFoundError(f"{self._data} does not exist")
        if self._data.is_dir():
            self._data_files = [self._data / f for f in self._data.listdir()]
        else:
            self._data_files = [self._data]
        
    @staticmethod
    def is_fasta(filename):
        if filename.exists():
            fasta = SeqIO.parse(str(filename), "fasta")
            return any(fasta)
        else:
            return False

    @staticmethod
    def is_gbk(filename):
        if filename.exists():
            gbk = SeqIO.parse(str(filename), "genbank")
            return any(gbk)
        else:
            return False

    def build(self, output_file: Path = None) -> LabelledFASTA:
        """
        Build database from data files.
        """
        if output_file is None:
            output_file = self._data.parent / f"{self._data.stem}_labelled.faa"
        if self.is_fasta(self._data_files[0]):
            logger.info("Translating and annotating assembly data.")
            if self._data.is_dir():
                assembly_fasta = FASTA.fromFASTAdirectory(self._data)
            else:
                assembly_fasta = FASTA(self._data)
            labelled_database = GeneAnnotator(
                assembly_fasta).annotate(output_file=output_file)
        elif self.is_gbk(self._data_files[0]):
            logger.info("Parsing GenBank data.")
            labelled_database = LabelledFASTA.fromGenBankData(
                self._data, output_file=output_file)
        else:
            logging.error(f"{self._data} is not a valid FASTA or genbank file")
            sys.exit(1)
        return labelled_database