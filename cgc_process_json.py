##########################
# to generate json file for all cgc_stardard.out file from run_dbcan
# use: python cgc_process_json.py -i cgc_standard.out -o cgc_standard.out.json
# written by Roland Madadjim in Cui's lab at Soc, UNL
# last updated: 12/09/2022
##########################

# from __future__ import print_function
import argparse
import json
import os

import numpy as np
import pandas as pd


class PrePro:
    """
    A class for preprocessing and organizing genetic data.

    Attributes
    ----------
        df (DataFrame): The initial data provided to the class.

    Args:
        data (DataFrame): The data to be processed.

    """

    def __init__(self, data):
        self.df = data

    def extract_gs(self, dataList):
        """
        Extracts gene strings from the provided data list.

        This method processes each entry in the dataList, extracting and concatenating specific gene-related information based on gene types like CAZyme, TC, TF, STP, and Null_Type.

        Args:
            dataList (list): A list of dictionaries, each containing gene and protein information.

        Returns
        -------
            str: A string representing concatenated gene information.
        """
        i = 0
        geneL = []
        # gene = list(map(lambda e: "{}".format(e["Gene_Type"]), dataList))
        # pfam = list(map(lambda e: "{}".format(e["Protein_Family"]), dataList))
        gene = ["{}".format(e["Gene_Type"]) for e in dataList]
        pfam = ["{}".format(e["Protein_Family"]) for e in dataList]
        while i < len(dataList):
            if gene[i] == "CAZyme":
                s = pfam[i]
                geneL.append(s)
            elif gene[i] == "TC":
                s = pfam[i]
                geneL.append(s)
            elif gene[i] == "TF":
                s = pfam[i]
                geneL.append(s)
            elif gene[i] == "STP":
                s = pfam[i]
                geneL.append(s)
            elif gene[i] == "Null_Type":
                s = "NA"
                geneL.append(s)
            i = i + 1
        gene_st = "-".join(geneL)
        return gene_st

    def pul_section(self):
        """
        Processes the dataframe and yields structured genetic data.

        This method groups the dataframe by 'cgc_id' and processes each group to yield a dictionary containing detailed genetic information for each 'cgc_id'.

        Yields
        ------
            dict: A dictionary containing genetic and protein information structured by 'cgc_id'.
        """
        for (cgc_id), df_pul_grouped in self.df.groupby("cgc_id"):
            datalist = list(self.cluster_section(df_pul_grouped))
            gene_str = self.extract_gs(datalist)
            yield {
                cgc_id: {
                    "changelog": [],
                    "Cluster_ID": cgc_id,
                    "Gene_String": gene_str,
                    "Contig_ID": self.df.loc[self.df["cgc_id"] == cgc_id, "contig_id"].iloc[0],
                    "ncbi_species_tax_id": [],
                    "organism_name": [],
                    "publication": [],
                    "Protein": list(self.cluster_section(df_pul_grouped)),
                    # "dbCan_Pul_accession": ID, # as string or integer?
                    # "publication": df.loc[df['ID'] == ID, 'PMID'].iloc[0],
                }
            }

    def cluster_section(self, df_pul_grouped):
        """
        Yields structured data for each gene cluster.

        This method iterates over grouped dataframe and yields a dictionary containing detailed information about each gene cluster.

        Args:
            df_pul_grouped (DataFrame): The grouped dataframe by specific gene identifiers.

        Yields
        ------
            dict: A dictionary containing information about a gene cluster.
        """
        for (
            contig_id,
            gene_type,
            protein_id,
            gene_start,
            gene_stop,
            direction,
            protein_family,
        ), _ in df_pul_grouped.groupby(
            [
                "contig_id",
                "gene_type",
                "protein_id",
                "gene_start",
                "gene_stop",
                "direction",
                "protein_family",
            ]
        ):
            yield {
                "contig_id": contig_id,
                "protein_id": protein_id,
                "Gene_Type": gene_type,
                "Gene_Start": gene_start,
                "Gene_Stop": gene_stop,
                "Strand": direction,
                "Protein_Family": protein_family,
            }

    def run_dbCan_section(self, df_puls):
        """
        Processes the dataframe and yields gene and protein information.

        This method iterates over each row of the dataframe and yields a dictionary containing gene and protein information.
        Currently, this method is a stub and needs to be implemented.

        Args:
            df_puls (DataFrame): The dataframe containing gene and protein data.

        Yields
        ------
            dict: A dictionary containing gene and protein information.
        """
        for _ in df_puls.itertuples():
            yield {
                # "Gene_Type": row.gene_type,
                # "Gene_Start": row.gene_start,
                # "Gene_Stop": row.gene_stop,
                # "Strand": row.direction,
                # "Protein_Family": row.protein_family
            }


def file_ext(choices, fname, parser):
    """
    Validates the extension of a given file name against allowed choices.

    This function checks if the file extension of `fname` is among the specified `choices`. If not, it raises an error through the provided `parser`.

    Parameters
    ----------
        choices (tuple): A tuple of allowed file extensions (without the dot).
        fname (str): The file name to be checked.
        parser (argparse.ArgumentParser): The argument parser to use for raising an error if the extension is not allowed.

    Returns
    -------
        str: The validated file name if its extension is in `choices`.

    Raises
    ------
        argparse.ArgumentError: If the file extension is not in the allowed `choices`.
    """
    ext = os.path.splitext(fname)[1][1:]
    if ext not in choices:
        parser.error("File needs to be a .out or .csv")
    return fname


class CustomEncoder(json.JSONEncoder):
    """
    Custom JSON encoder for handling numpy data types.

    This encoder extends `json.JSONEncoder` to provide custom serialization for certain data types not natively supported by the default JSON encoder. Specifically, it converts `numpy.int64` types to strings.

    Methods
    -------
        default(obj): Overrides the default method to provide custom serialization.
    """

    def default(self, obj):
        """
        Provides custom serialization for certain data types.

        Parameters
        ----------
            obj: The object to serialize.

        Returns
        -------
            The serialized object, converting `numpy.int64` to string. For other types, it relies on the superclass implementation.

        Raises
        ------
            TypeError: If the object is not a recognized type and cannot be serialized by the superclass.
        """
        if isinstance(obj, np.int64):
            return str(obj)
        return super().default(obj)


def main():
    """
    Main function to compile JSON from a 'cgc_standard.out' file.

    This function parses command-line arguments to get input and output file paths, reads the input file, processes and renames columns, and writes the processed data to a JSON file using a custom JSON encoder.

    Side Effects:
        - Reads a specified input file.
        - Writes processed data to a JSON output file.
        - Can terminate the script if the input file does not have the correct extension.
    """
    parser = argparse.ArgumentParser(description="Compiling Json from cgc_standard.out")
    parser.add_argument(
        "-i",
        required=True,
        help="path to output file (cgc_standard.out) file",
        type=lambda s: file_ext(("out", "csv"), s, parser),
    )
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    with open(args.i) as file:  ### input files
        data = pd.read_csv(file, sep="\t")
        data.rename(
            columns={
                "CGC#": "cgc_id",
                "Gene Type": "gene_type",
                "Contig ID": "contig_id",
                "Protein ID": "protein_id",
                "Gene Start": "gene_start",
                "Gene Stop": "gene_stop",
                "Direction": "direction",
                "Protein Family": "protein_family",
            },
            inplace=True,
        )
        data["gene_type"].fillna("Null_Type", inplace=True)
        data["protein_family"].fillna("0", inplace=True)
        p = PrePro(data)

    pul_list = list(p.pul_section())
    pul_dict = {}
    for sub_dict in pul_list:
        pul_dict.update(sub_dict)
    jsonPuls = json.dumps(pul_dict, indent=4, cls=CustomEncoder)

    with open(args.output, "w") as outfile:
        # with open("Json"+time.strftime("%Y%m%d%H%M%S")+".json","w") as outfile:
        outfile.write(jsonPuls)


if __name__ == "__main__":
    main()
