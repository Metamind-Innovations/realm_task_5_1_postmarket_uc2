import argparse
import glob
import json
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd


# Target genes for analysis
TARGET_GENES = ["CYP2B6", "CYP2C9", "CYP2C19", "CYP3A5", "SLCO1B1", "TPMT", "DPYD"]

SAMPLE_ID_COL = "Sample ID"
CSV_EXT = "*.csv"
CSV_SUFFIX = ".csv"


class VCFParser:
    def __init__(self):
        self.meta_info = {}
        self.header = None

    def parse_vcf(self, vcf_file):
        try:
            with open(vcf_file, "r") as f:
                lines = f.readlines()

            data_lines = []
            for line in lines:
                line = line.strip()
                if line.startswith("##"):
                    self._parse_meta_line(line)
                elif line.startswith("#CHROM"):
                    self.header = line[1:].split("\t")
                else:
                    data_lines.append(line.split("\t"))

            if not self.header or not data_lines:
                return None

            df = pd.DataFrame(data_lines, columns=self.header)
            df["Gene"] = df["INFO"].apply(self._extract_gene)

            sample_columns = self.header[9:]
            if sample_columns:
                for sample in sample_columns:
                    df[f"{sample}_GT"] = df.apply(
                        lambda row: self._extract_genotype(row["FORMAT"], row[sample]),
                        axis=1,
                    )

            return df

        except Exception as e:
            print(f"Error parsing VCF file: {e}")
            return None

    def _parse_meta_line(self, line):
        if line.startswith("##"):
            line = line[2:]
            if "=" in line:
                key, value = line.split("=", 1)
                self.meta_info[key] = value

    def _extract_gene(self, info_field):
        gene_match = re.search(r"PX=([^;]+)", info_field)
        return gene_match.group(1) if gene_match else ""

    def _extract_genotype(self, format_field, sample_field):
        format_parts = format_field.split(":")
        sample_parts = sample_field.split(":")

        if "GT" in format_parts:
            gt_index = format_parts.index("GT")
            if gt_index < len(sample_parts):
                return sample_parts[gt_index]

        return ""

    def vcf_to_csv(self, vcf_file, output_csv=None):
        if output_csv is None:
            output_csv = os.path.splitext(vcf_file)[0] + CSV_SUFFIX

        df = self.parse_vcf(vcf_file)
        if df is not None:
            df.to_csv(output_csv, index=False)
            return output_csv
        return None


def preprocess_input_data(input_dir, output_dir=None):
    if output_dir is None:
        output_dir = os.path.join(input_dir, "preprocessed")

    os.makedirs(output_dir, exist_ok=True)

    vcf_files = glob.glob(os.path.join(input_dir, "*.vcf"))
    sample_to_file = {}

    for vcf_file in vcf_files:
        filename = os.path.basename(vcf_file)
        sample_id = filename.split("_")[0]

        parser = VCFParser()
        output_csv = os.path.join(output_dir, f"{sample_id}_preprocessed.csv")

        if parser.vcf_to_csv(vcf_file, output_csv):
            sample_to_file[sample_id] = output_csv
            print(f"Converted {vcf_file} to {output_csv}")

    return sample_to_file, output_dir


def find_csv_files(input_dir):
    # Direct CSV files
    csv_files = glob.glob(os.path.join(input_dir, CSV_EXT))
    if csv_files:
        return csv_files

    # Check preprocessed directory
    preprocessed_dir = os.path.join(input_dir, "preprocessed")
    if os.path.exists(preprocessed_dir):
        return glob.glob(os.path.join(preprocessed_dir, CSV_EXT))

    return []


def main():
    parser = argparse.ArgumentParser(description="PharmCAT SHAP-based explainer")
    parser.add_argument(
        "--input_dir", required=True, help="Directory containing VCF or CSV files"
    )
    # parser.add_argument(
    #     "--phenotypes_file", required=True, help="Path to phenotypes.csv file"
    # )
    parser.add_argument(
        "--output_dir", default="pgx_shap_results", help="Output directory for results"
    )
    parser.add_argument(
        "--convert_vcf", action="store_true", help="Convert VCF files to CSV format"
    )
    # parser.add_argument(
    #     "--max_samples",
    #     type=int,
    #     default=100,
    #     help="Maximum number of samples for SHAP analysis",
    # )
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Processing input data from {args.input_dir}...")

    input_files = []
    if args.convert_vcf:
        _, csv_dir = preprocess_input_data(
            args.input_dir, os.path.join(args.output_dir, "preprocessed")
        )
        input_files = glob.glob(os.path.join(csv_dir, CSV_EXT))
        print(f"Converted VCF files to CSV format in {csv_dir}")
    else:
        input_files = find_csv_files(args.input_dir)
        if not input_files:
            raise ValueError(
                f"No CSV files found in {args.input_dir}. Use --convert_vcf to convert VCF files."
            )

    print(f"Found {len(input_files)} input files")


if __name__ == "__main__":
    main()
