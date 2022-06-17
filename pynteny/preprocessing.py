#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to preprocess sequence databases

1. Remove illegal characters from peptide sequences
2. Remove illegal symbols from file paths
3. Relabel fasta records and make dictionary with old labels
"""

import os
from pathlib import Path
from typing import Self

from Bio import SeqIO
import pyfastx

import pynteny.wrappers as wrappers
from pynteny.utils import (setDefaultOutputPath,
                           terminalExecute, handle_exceptions)



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
        self._str_path = input_file.as_posix()
        self._input_file = Path(input_file)
        return None

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
        return None
  
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
        return None

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

        fasta = pyfastx.Fasta(self._input_file, build_index=False, full_name=True)
        with open(output_file, 'w') as outfile:
            for record_name, record_seq in fasta:
                if is_peptide and (not keep_stop_codon):
                    record_seq = RecordSequence.removeStopCodonSignals(record_seq)
                if isLegitSequence(record_seq):
                    outfile.write(f'>{record_name}\n{record_seq}\n')
        return None

    def filterByIDs(self, record_ids: list,
                    output_file: Path = None) -> None:
        """
        Filter records in fasta file matching provided IDs
        """
        if output_file is None:
            output_fasta = setDefaultOutputPath(self._input_file, '_fitered')
        record_ids = set(record_ids)
        fa = pyfastx.Fasta(self._input_file)
        with open(output_fasta, 'w') as fp:
            for record_id in record_ids:
                try:
                    record_obj = fa[record_id]
                    fp.write(record_obj.raw)
                except:
                    pass
        os.remove(self._input_file + ".fxi")
        return None

    def splitByContigs(self, output_dir: Path = None) -> None:
        """
        Split large fasta file into several ones containing one contig each
        """
        if output_dir is None:
            output_dir = Path(self._input_file.parent) / "split_" + self._input_file.name
        os.makedirs(output_dir, exist_ok=True)
        contigs = pyfastx.Fasta(self._input_file, build_index=False, full_name=True)
        for contig_name, seq in contigs:
            outfile = output_dir / f"{contig_name.split(' ')[0]}{self._input_file.suffix}"
            with open(outfile, "w+") as file:
                file.write(f">{contig_name}\n")
                file.write(seq + "\n")
        return None


class LabelledFASTA(FASTA):
    def __init__(self, input_file: Path) -> None:
        """
        Tools to add and parse FASTA with positional info on record tags
        """
        self._input_file = input_file
        return None
    
    @classmethod
    def fromProdigalOutput(cls,
                           prodigal_faa: Path,
                           output_file: Path = None) -> Self:
        """
        Extract positional gene info from prodigal output and export to
        fasta file.
        """
        number_prodigal_record_fields = 9
        if output_file is None:
            output_file = setDefaultOutputPath(prodigal_faa, tag="_longlabels")
        data = pyfastx.Fasta(prodigal_faa, build_index=False, full_name=True)
        with open(output_file, "w") as outfile:
            for record_name, record_seq in data:
                name_list = record_name.split(" ")
                if len(name_list) < number_prodigal_record_fields:
                    raise ValueError(f"Invalid prodigal header format for record: {record_name}")
                contig = "_".join(name_list[0].split("_")[:-1])
                gene_number = name_list[0].split("_")[-1]
                start, end = name_list[2], name_list[4]
                strand = "pos" if name_list[6] == "1" else "neg"
                header = f">{contig}_{gene_number}__{contig}_{gene_number}_{start}_{end}_{strand}"
                outfile.write(header + "\n")
                outfile.write(record_seq + "\n")
        return LabelledFASTA(output_file)

    @classmethod
    def fromGenBankFile(cls,
                        gbk_file: Path,
                        output_file: Path = None,
                        nucleotide: bool = False) -> Self:
        """
        Assign gene positional info, such as contig, gene number and loci
        to each record in database
        @paramms:
        nucleotide: if True then records are nucleotide sequences instead of peptides.
                    Note that this option will notably increase the computation time.
        """
        if output_file is None:
            output_file = setDefaultOutputPath(gbk_file, extension=".fasta")
            # output_file = Path(gbk_file.parent) / gbk_file.stem + ".fasta"
        gbk_contigs = list(SeqIO.parse(gbk_file, 'genbank'))
        
        def get_label_str(gbk_contig, feature):
            name = feature.qualifiers["locus_tag"][0].replace('_', '.')
            start, end, strand = str(feature.location.start), str(feature.location.end), feature.location.strand
            start = start.replace(">", "").replace("<", "")
            end = end.replace(">", "").replace("<", "")
            strand_sense = "neg" if strand == -1 else "pos"
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


# OLD functions

@handle_exceptions
def removeDuplicatesFromFasta(input_fasta: str,
                              output_fasta: str = None,
                              export_duplicates: bool = False,
                              method: str = 'seqkit') -> None:
    """
    Removes duplicate entries (either by sequence or ID) from fasta.
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(input_fasta, '_noduplicates')

    if 'bio' in method:
        seen_seqs, seen_ids = set(), set()
        def unique_records():
            for record in SeqIO.parse(input_fasta, 'fasta'):  
                if (record.seq not in seen_seqs) and (record.id not in seen_ids):
                    seen_seqs.add(record.seq)
                    seen_ids.add(record.id)
                    yield record

        SeqIO.write(unique_records(), output_fasta, 'fasta')

    else:
        wrappers.runSeqKitNoDup(input_fasta=input_fasta, output_fasta=output_fasta,
                                export_duplicates=export_duplicates)

@handle_exceptions
def removeCorruptedSequences(fasta_file: str,
                             output_file: str = None,
                             is_peptide: bool = True,
                             keep_stop_codon: bool = False) -> None:
    """
    Filter out (DNA or peptide) sequences containing illegal characters
    """
    dirname = os.path.dirname(fasta_file)
    basename = os.path.basename(fasta_file)
    fname, ext = os.path.splitext(basename)

    def removeStopCodonSignals(record_seq: str) -> str:
        return record_seq.replace('*', '')

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

    def isLegitDNAsequence(record_seq: str) -> bool:
        """
        Assert that DNA sequence only contains valid symbols
        """
        nts = {'A', 'G', 'T', 'C', 'N'}
        seq_symbols = {s.upper() for s in record_seq}
        return seq_symbols.issubset(nts)

    if output_file is None:
        output_file = os.path.join(dirname, f'{fname}_modified{ext}')
    else:
        output_file = os.path.abspath(output_file)
    if is_peptide:
        isLegitSequence = isLegitPeptideSequence
    else:
        isLegitSequence = isLegitDNAsequence

    fasta = pyfastx.Fasta(fasta_file, build_index=False, full_name=True)
    with open(output_file, 'w') as outfile:
        for record_name, record_seq in fasta:
            if is_peptide and (not keep_stop_codon):
                record_seq = removeStopCodonSignals(record_seq)
            if isLegitSequence(record_seq):
                outfile.write(f'>{record_name}\n{record_seq}\n')
                                

def splitFASTAbyContigs(input_fasta: str, output_dir: str = None) -> None:
    """
    Split large fasta file into several ones containing one contig each
    """
    if output_dir is None:
        output_dir = os.path.join(
            setDefaultOutputPath(input_fasta, only_dirname=True),
            "split_" + setDefaultOutputPath(input_fasta, only_basename=True)
        )
    os.makedirs(output_dir, exist_ok=True)
    base, ext = os.path.splitext(input_fasta)
    contigs = pyfastx.Fasta(input_fasta, build_index=False, full_name=True)
    for contig_name, seq in contigs:
        outfile = os.path.join(output_dir, f"{contig_name.split(' ')[0]}{ext}")
        with open(outfile, "w+") as file:
            file.write(f">{contig_name}\n")
            file.write(seq + "\n")


def mergeFASTAs(input_fastas_dir: str, output_fasta: str = None) -> None:
    """
    Merge input fasta files into a single fasta
    """
    if output_fasta is None:
        output_fasta = os.path.join(input_fastas_dir, 'merged.fasta')
    cmd_str = f'awk 1 * > {output_fasta}'
    terminalExecute(
        cmd_str,
        work_dir=input_fastas_dir,
        suppress_shell_output=False
        )


def parseProdigalOutput(prodigal_faa: str, output_file: str = None) -> None:
    """
    Extract positional gene info from prodigal output and export to
    fasta file.
    """
    number_prodigal_record_fields = 9
    if output_file is None:
        output_file = setDefaultOutputPath(prodigal_faa, tag="_longlabels")
    data = pyfastx.Fasta(prodigal_faa, build_index=False, full_name=True)
    with open(output_file, "w") as outfile:
        for record_name, record_seq in data:
            name_list = record_name.split(" ")
            if len(name_list) < number_prodigal_record_fields:
                raise ValueError(f"Invalid prodigal header format for record: {record_name}")
            contig = "_".join(name_list[0].split("_")[:-1])
            gene_number = name_list[0].split("_")[-1]
            start, end = name_list[2], name_list[4]
            strand = "pos" if name_list[6] == "1" else "neg"
            header = f">{contig}_{gene_number}__{contig}_{gene_number}_{start}_{end}_{strand}"
            outfile.write(header + "\n")
            outfile.write(record_seq + "\n")


def parseGBK(gbk_file: str, output_fasta: str = None,
                                nucleotide: bool = False) -> None:
    """
    Assign gene positional info, such as contig, gene number and loci
    to each record in database
    @paramms:
    nucleotide: if True then records are nucleotide sequences instead of peptides.
                Note that this option will notably increase the computation time.
    """
    if output_fasta is None:
        output_fasta = setDefaultOutputPath(gbk_file, extension=".fasta")
    gbk_contigs = list(SeqIO.parse(gbk_file, 'genbank'))
    
    def get_label_str(gbk_contig, feature):
        name = feature.qualifiers["locus_tag"][0].replace('_', '.')
        start, end, strand = str(feature.location.start), str(feature.location.end), feature.location.strand
        start = start.replace(">", "").replace("<", "")
        end = end.replace(">", "").replace("<", "")
        strand_sense = "neg" if strand == -1 else "pos"
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

    with open(output_fasta, "w") as outfile:
        for gbk_contig in gbk_contigs:
            gene_counter = 0
            for feature in gbk_contig.features:
                if "cds" in feature.type.lower():
                    gene_counter = write_record(gbk_contig, feature, outfile, gene_counter)