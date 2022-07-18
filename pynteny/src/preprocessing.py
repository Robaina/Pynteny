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

import pynteny.src.wrappers as wrappers
from pynteny.src.utils import setDefaultOutputPath, terminalExecute

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

    @staticmethod
    def mergeFASTAs(input_dir: Path,
                    output_file: Path = None) -> None:
        """
        Merge input fasta files into a single fasta
        """
        if output_file is None:
            output_file = input_dir / "merged.fasta"
        cmd_str = f'awk 1 * > {output_file.as_posix()}'
        terminalExecute(
            cmd_str,
            work_dir=input_dir,
            suppress_shell_output=False
            )
  
    def removeDuplicates(self,
                         output_file: Path = None,
                         export_duplicates: bool = False,
                         method: str = 'seqkit') -> None:
        """
        Removes duplicate entries (either by sequence or ID) from fasta.
        """
        if output_file is None:
            output_file = setDefaultOutputPath(self._input_file, '_noduplicates')

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

    def removeCorruptedSequences(self,
                                 output_file: Path = None,
                                 is_peptide: bool = True,
                                 keep_stop_codon: bool = False) -> None:
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

    def filterByIDs(self, record_ids: list,
                    output_file: Path = None) -> None:
        """
        Filter records in fasta file matching provided IDs
        """
        if output_file is None:
            output_file = setDefaultOutputPath(self._input_file, tag="_fitered")
        with tempfile.NamedTemporaryFile(mode="w+t") as tmp_ids:
            tmp_ids.writelines("\n".join(record_ids))
            tmp_ids.flush()
            tmp_ids_path = tmp_ids.name
            cmd_str = f"seqkit grep -i -f {tmp_ids_path} {self._input_file} -o {output_file}"
            terminalExecute(cmd_str)
        
    def splitByContigs(self, output_dir: Path = None) -> None:
        """
        Split large fasta file into several ones containing one contig each
        """
        if output_dir is None:
            output_dir = Path(self._input_file.parent) / "split_" + self._input_file.name
        os.makedirs(output_dir, exist_ok=True)
        contigs = pyfastx.Fasta(self._input_file_str, build_index=False, full_name=True)
        for contig_name, seq in contigs:
            outfile = output_dir / f"{contig_name.split(' ')[0]}{self._input_file.suffix}"
            with open(outfile, "w+") as file:
                file.write(f">{contig_name}\n")
                file.write(seq + "\n")


class LabelledFASTA(FASTA):
    """
    Tools to add and parse FASTA with positional info on record tags
    """
    @classmethod
    def fromProdigalOutput(cls,
                           prodigal_faa: Path,
                           output_file: Path = None):
        """
        Extract positional gene info from prodigal output and export to
        fasta file.
        """
        number_prodigal_record_fields = 9
        if output_file is None:
            output_file = setDefaultOutputPath(prodigal_faa, tag="_longlabels")
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
    def fromGenBankFile(cls,
                        gbk_file: Path,
                        output_file: Path = None,
                        nucleotide: bool = False):
        """
        Assign gene positional info, such as contig, gene number and loci
        to each record in database
        @paramms:
        nucleotide: if True then records are nucleotide sequences instead of peptides.
                    Note that this option will notably increase the computation time.
        """
        if output_file is None:
            output_file = setDefaultOutputPath(gbk_file, extension=".fasta")
        gbk_contigs = list(SeqIO.parse(gbk_file, 'genbank'))
        
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