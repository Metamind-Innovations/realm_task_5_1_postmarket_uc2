import argparse
import glob
import json
import os


def validate_vcf_positions(vcf_file_path):
    """
    Validates that all positions in a VCF file are within the length of their respective chromosomes.

    Args:
        vcf_file_path (str): Path to the VCF file to validate

    Returns:
        list: List of tuples containing (chrom, pos, max_length) for invalid positions
    """
    chrom_lengths = {}
    invalid_positions = []

    with open(vcf_file_path, "r") as vcf_file:
        for line in vcf_file:
            line = line.strip()

            if line.startswith("##contig=<ID="):
                contig_info = line[10:-1]
                parts = contig_info.split(",")

                chrom_id = None
                length = None

                for part in parts:
                    if part.startswith("ID="):
                        chrom_id = part[3:]
                    elif part.startswith("length="):
                        length = int(part[7:])

                if chrom_id and length:
                    chrom_lengths[chrom_id] = length

            elif line.startswith("#CHROM"):
                break

        for line in vcf_file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue

            chrom = fields[0]
            try:
                pos = int(fields[1])
                if chrom in chrom_lengths and pos > chrom_lengths[chrom]:
                    invalid_positions.append((chrom, pos, chrom_lengths[chrom]))
            except ValueError:
                continue

    return invalid_positions


def report_invalid_positions(vcf_file_path):
    """
    Checks a VCF file for invalid positions and reports them.

    Args:
        vcf_file_path (str): Path to the VCF file to validate

    Returns:
        dict: Results dictionary containing validation status and any invalid positions
    """
    invalid_positions = validate_vcf_positions(vcf_file_path)
    results = {
        "is_valid": len(invalid_positions) == 0,
        "invalid_positions": [
            {"chromosome": chrom, "position": pos, "max_length": max_length}
            for chrom, pos, max_length in invalid_positions
        ],
    }
    return results


def validate_nucleotides(vcf_file_path):
    """
    Validates that all REF and ALT fields in a VCF file contain only valid nucleotide characters (A, C, G, T, N).

    Args:
        vcf_file_path (str): Path to the VCF file to validate

    Returns:
        list: List of tuples containing (chrom, pos, field, value) for invalid nucleotide sequences
    """
    valid_nucleotides = set("ACGTN")

    invalid_sequences = []

    with open(vcf_file_path, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("#CHROM"):
                break
            elif line.startswith("#"):
                continue

        for line in vcf_file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt_field = fields[4]

            if not all(nucleotide in valid_nucleotides for nucleotide in ref.upper()):
                invalid_sequences.append((chrom, pos, "REF", ref))

            for alt in alt_field.split(","):
                if not all(
                    nucleotide in valid_nucleotides for nucleotide in alt.upper()
                ):
                    invalid_sequences.append((chrom, pos, "ALT", alt))

    return invalid_sequences


def report_invalid_nucleotides(vcf_file_path):
    """
    Checks a VCF file for invalid nucleotide sequences and reports them.

    Args:
        vcf_file_path (str): Path to the VCF file to validate

    Returns:
        dict: Results dictionary containing validation status and any invalid sequences
    """
    invalid_sequences = validate_nucleotides(vcf_file_path)
    results = {
        "is_valid": len(invalid_sequences) == 0,
        "invalid_sequences": [
            {"chromosome": chrom, "position": pos, "field": field, "value": value}
            for chrom, pos, field, value in invalid_sequences
        ],
    }
    return results


def expert_knowledge_evaluation(input_dir):
    """
    Process all VCF files in the input directory and generate validation results.

    Args:
        input_dir (str): Directory containing VCF files

    Returns:
        dict: Results dictionary containing validation results and summary statistics
    """
    vcf_files = glob.glob(os.path.join(input_dir, "*.vcf"))

    if not vcf_files:
        results = {"error": f"No VCF files found in {input_dir}"}
    else:
        results = {
            "criteria": {
                "nucleotides": {
                    "description": "Based on Hirao et. al. (https://pubs.acs.org/doi/abs/10.1021/ar200257x) allele bases can be A, C, G, T (+ N which is used for undefined bases).",
                    "valid_bases": ["A", "C", "G", "T", "N"],
                },
                "positions": {
                    "description": "For each chromosome, the position of the of the variant should be smaller than the chromosome's total length (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.39/).",
                    "source": "NCBI GRCh38 Assembly",
                },
            },
            "summary": {
                "total_files": len(vcf_files),
                "positions": {"valid": 0, "invalid": 0},
                "nucleotides": {"valid": 0, "invalid": 0},
            },
            "files": {},
        }

        for vcf_file in vcf_files:
            file_results = {
                "position_validation": report_invalid_positions(vcf_file),
                "nucleotide_validation": report_invalid_nucleotides(vcf_file),
            }
            results["files"][os.path.basename(vcf_file)] = file_results

            # Summary statistics
            if file_results["position_validation"]["is_valid"]:
                results["summary"]["positions"]["valid"] += 1
            else:
                results["summary"]["positions"]["invalid"] += 1

            if file_results["nucleotide_validation"]["is_valid"]:
                results["summary"]["nucleotides"]["valid"] += 1
            else:
                results["summary"]["nucleotides"]["invalid"] += 1

        total_files = results["summary"]["total_files"]
        if total_files > 0:
            results["summary"]["positions"]["valid_percentage"] = (
                results["summary"]["positions"]["valid"] / total_files
            ) * 100
            results["summary"]["positions"]["invalid_percentage"] = (
                results["summary"]["positions"]["invalid"] / total_files
            ) * 100
            results["summary"]["nucleotides"]["valid_percentage"] = (
                results["summary"]["nucleotides"]["valid"] / total_files
            ) * 100
            results["summary"]["nucleotides"]["invalid_percentage"] = (
                results["summary"]["nucleotides"]["invalid"] / total_files
            ) * 100

    return results


def main():
    parser = argparse.ArgumentParser(description="Validate VCF files")
    parser.add_argument(
        "--input_dir", required=True, help="Directory containing VCF files"
    )
    parser.add_argument(
        "--output_file",
        default="artifacts/expert_knowledge_results.json",
        help="Output JSON file path",
    )
    args = parser.parse_args()

    # Create artifacts directory if it doesn't exist
    os.makedirs("artifacts", exist_ok=True)

    results = expert_knowledge_evaluation(args.input_dir)

    with open(args.output_file, "w") as f:
        json.dump(results, f, indent=2)

    with open(args.output_file, "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
