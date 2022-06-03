import os
import shutil
import warnings
from typing import OrderedDict

from Bio import SeqIO, SeqFeature, SeqRecord
import pandas as pd

from utils import setDefaultOutputPath, terminalExecute, fullPathListDir, unZipFile, mergeMultiRecordGBK



class GBK():

    def __init__(self, gbk_file: str) -> None:
        """
        Tools to parse GenBank files
        """
        print("Initializing parser...")
        # Deal with compressed GBKs
        if gbk_file.endswith('.gz'):
            unZipFile(
                input_file=gbk_file,
            )
            gbk_file = gbk_file.strip('.gz')
        
        file_handle = open(gbk_file)
        gbk_objs = list(SeqIO.parse(file_handle, 'genbank'))
        file_handle.close()

        # Deal with multi-record GBKs
        if len(gbk_objs) > 1:
            merged_gbk = setDefaultOutputPath(
                input_path=gbk_file, tag='_merged'
            )
            mergeMultiRecordGBK(
                input_file=gbk_file,
                output_file=merged_gbk
            )
            self._gbk = list(SeqIO.parse(merged_gbk, 'genbank'))[0]
            shutil.move(merged_gbk, gbk_file)
        else:
            self._gbk = gbk_objs[0]
        print("Done!")
    
    @property
    def cds(self) -> list:
        return GBKfeatureList([f for f in self._gbk.features if 'cds' in f.type.lower()])

    @property
    def meta(self) -> OrderedDict:
        return self._gbk.features[0].qualifiers

    @property
    def gbk_object(self) -> SeqRecord:
        return self._gbk

class GBKfeatureList(list):
    """
    Allow filtering cds feature objects from gbk by keywords
    """
    def __init__(self, seq=None):
        super(self.__class__, self).__init__(seq)

    def __getslice__(self, start, stop):
        return self.__class__(super(self.__class__, self).__getslice__(start, stop))

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.__class__(super(self.__class__, self).__getitem__(key))
        else:
            return super(self.__class__, self).__getitem__(key)
    
    @staticmethod
    def _text_contains_keywords(text: str, keywords: list,
                                case_insensitive: bool = True) -> bool:
        """
        Match keyword in text string
        """
        if 'any' in [k.lower() for k in keywords]:
            return True
        else:
            if case_insensitive:
                return all([key.lower() in text.lower() for key in keywords])
            else:
                return all([key in text for key in keywords])
    
    def _cds_matched(self, cds: SeqFeature, feature_keywords: dict,
                     case_insensitive: bool = True) -> bool:
        cds = cds.qualifiers
        return all(
            [(field in cds.keys() and
              self._text_contains_keywords(cds[field][0], keywords, case_insensitive))
            for field, keywords in feature_keywords.items()]
            )

    def get_by_keywords(self, keywords: dict, case_insensitive: bool = True) -> list:
        """
        Extract cds record matching keywords.
        @Arguments:
        keywords is a dictionary in which keys correspond to  cds fields 
        and values to keywords to find in each field. For instace,
        keywords = {
            'gene': ['ureC'],
            'product': ['urease', 'alpha'] 
        }
        case_insensitive: whether or not to care for case when matching keywords 
        NOTE:
        Use keywords = {
            'gene': ['any']
        }
        To get the complete list of cds entries
        """
        try:
            return [
                cds for cds in self if self._cds_matched(cds, keywords, case_insensitive)
                ]
        except Exception as e:
            raise ValueError(f'Feature not found for given keyword(s). Exception: {e}')

    def get_by_gene_id(self, gene_id: str, case_insensitive: bool = True) -> SeqFeature:
        return self.get_by_keywords({'gene': [gene_id]}, case_insensitive)

    def get_by_ec_number(self, ec: str, case_insensitive: bool = True) -> SeqFeature:
        return self.get_by_keywords({'ec_number': [ec]}, case_insensitive)