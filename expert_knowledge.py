import argparse
import glob
import os


def validate_vcf_positions(vcf_file_path):
    """
    Validates that all positions in a VCF file are within the length of their respective chromosomes.

    Args:
        vcf_file_path (str): Path to the VCF file to validate

    Returns:
        list: List of tuples containing (chrom, pos, max_length) for invalid positions
    """
    # Dictionary to store chromosome lengths from the header
    chrom_lengths = {}

    # List to store invalid positions
    invalid_positions = []

    with open(vcf_file_path, "r") as vcf_file:
        # Parse header to extract chromosome lengths
        for line in vcf_file:
            line = line.strip()

            # Extract chromosome lengths from contig lines in the header
            if line.startswith("##contig=<ID="):
                # Parse the contig line to extract chromosome ID and length
                contig_info = line[10:-1]  # Remove '##contig=<' and '>'
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

            # Start processing variant lines after the header
            elif line.startswith("#CHROM"):
                break

        # Process variant lines
        for line in vcf_file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 5:  # Ensure we have at least CHROM and POS fields
                continue

            chrom = fields[0]
            try:
                pos = int(fields[1])
                # Check if position is valid
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
        bool: True if all positions are valid, False otherwise
    """
    invalid_positions = validate_vcf_positions(vcf_file_path)

    if not invalid_positions:
        print(f"All positions in {vcf_file_path} are valid.")
        return True
    else:
        print(f"Found {len(invalid_positions)} invalid positions in {vcf_file_path}:")
        for chrom, pos, max_length in invalid_positions:
            print(f"  {chrom}:{pos} exceeds maximum length of {max_length}")
        return False


def validate_nucleotides(vcf_file_path):
    """
    Validates that all REF and ALT fields in a VCF file contain only valid nucleotide characters (A, C, G, T, N).

    Args:
        vcf_file_path (str): Path to the VCF file to validate

    Returns:
        list: List of tuples containing (chrom, pos, field, value) for invalid nucleotide sequences
    """
    # Valid nucleotide characters
    valid_nucleotides = set("ACGTN")

    # List to store invalid nucleotide sequences
    invalid_sequences = []

    with open(vcf_file_path, "r") as vcf_file:
        # Skip header lines
        for line in vcf_file:
            if line.startswith("#CHROM"):
                break
            elif line.startswith("#"):
                continue

        # Process variant lines
        for line in vcf_file:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if (
                len(fields) < 5
            ):  # Ensure we have at least CHROM, POS, ID, REF, ALT fields
                continue

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt_field = fields[4]

            # Check REF field
            if not all(nucleotide in valid_nucleotides for nucleotide in ref.upper()):
                invalid_sequences.append((chrom, pos, "REF", ref))

            # Check each ALT allele (comma-separated)
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
        bool: True if all nucleotide sequences are valid, False otherwise
    """
    invalid_sequences = validate_nucleotides(vcf_file_path)

    if not invalid_sequences:
        print(f"All nucleotide sequences in {vcf_file_path} are valid.")
        return True
    else:
        print(
            f"Found {len(invalid_sequences)} invalid nucleotide sequences in {vcf_file_path}:"
        )
        for chrom, pos, field, value in invalid_sequences:
            print(f"  {chrom}:{pos} {field}={value}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Validate VCF files")
    parser.add_argument(
        "--input_dir", required=True, help="Directory containing VCF files"
    )
    args = parser.parse_args()

    # Get all VCF files in the directory
    vcf_files = glob.glob(os.path.join(args.input_dir, "*.vcf"))

    if not vcf_files:
        print(f"No VCF files found in {args.input_dir}")
        return

    print(f"Found {len(vcf_files)} VCF files to validate")

    # Process each VCF file
    for vcf_file in vcf_files:
        print(f"\nValidating {os.path.basename(vcf_file)}...")

        print("Checking positions...")
        report_invalid_positions(vcf_file)

        print("Checking nucleotides...")
        report_invalid_nucleotides(vcf_file)


if __name__ == "__main__":
    main()


# TODO: Format output so that it includes info about the criterion used as well as a source for it.
# TODO: Make output a json file.
