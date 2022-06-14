#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tools to preprocess sequence databases

1. Remove illegal characters from peptide sequences
2. Remove illegal symbols from file paths
3. Relabel fasta records and make dictionary with old labels
"""

import os
# import re

from Bio import SeqIO
import pyfastx

import pynteny.wrappers as wrappers
from pynteny.utils import (saveToPickleFile, setDefaultOutputPath,
                           terminalExecute, handle_exceptions)


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
        

def parseProdigalOutput(prodigal_faa: str, output_file: str = None) -> str:
    """
    Extract positional gene info from prodigal output and export to
    fasta file.
    """
    if output_file is None:
        output_file = setDefaultOutputPath(prodigal_faa, tag="_longlabels")
    data = pyfastx.Fasta(prodigal_faa, build_index=False, full_name=True)
    with open(output_file, "w") as outfile:
        for record_name, record_seq in data:
            name_list = record_name.split(" ")
            if len(name_list) < 9:
                raise ValueError(f"Invalid prodigal header format for record: {record_name}")
            contig = "_".join(name_list[0].split("_")[:-1])
            gene_number = name_list[0].split("_")[-1]
            start, end = name_list[2], name_list[4]
            strand = "pos" if name_list[6] == "1" else "neg"
            header = f">{contig}_{gene_number}__{contig}_{gene_number}_{start}_{end}_{strand}"
            outfile.write(header + "\n")
            outfile.write(record_seq + "\n")


def assignGeneLocationToRecords(gbk_file: str, output_fasta: str = None,
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