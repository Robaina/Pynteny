# Add MMP to marref labels (now only with accession ID)
import pyfastx 

from pathlib import Path
import pandas as pd 



class MARref:
    def __init__(self, mar_meta: Path):
        "Class to parse MAR database meta file"
        self._meta = pd.read_csv(mar_meta, sep="\t")

    def getTaxonomy(self, accession: str) -> str:
        """
        Get GTDB taxonomy from MARref database
        """
        try:
            taxopath = self._meta.loc[
                self._meta["acc:genbank"].str.contains(accession, case=False),
                "tax:gtdb_classification"
                ].values[0].replace(">", ";")
        except:
            taxopath = ""
        return taxopath

    def getMMPid(self, accession: str) -> str:
        """
        Get MMP id from given accession
        """
        try:
            mmp = self._meta.loc[
                self._meta["acc:genbank"].str.contains(accession, case=False),
                "id"
                ].values[0].replace("mmp.ref:", "")
        except:
            mmp = ""
        return mmp



mar = MARref("/home/robaina/Databases/MAR_database/MarRef_1.7.tsv")

output_file = "/home/robaina/Databases/MAR_database/marref_prodigal_longlabels_mmp.faa"
fasta = pyfastx.Fasta("/home/robaina/Databases/MAR_database/marref_prodigal_longlabels.faa", build_index=False, full_name=True)

with open(output_file, 'w') as outfile:
    lines = []
    for record_name, record_seq in fasta:
        mmp_id = mar.getMMPid(record_name.split("_")[0])
        lines.append(f">({mmp_id}){record_name}\n")
        lines.append(f"{record_seq}\n")
    outfile.writelines(lines)