import argparse
import glob
import os
import re


def vcf_header_consistency(vcf_file_path):
    """
    Validates the consistency of a VCF file header and data lines.

    Checks:
    1) All data lines are TAB delimited and the number of fields in each data line
       matches the number of fields in the header line
    2) The header lines ##fileformat and #CHROM are mandatory
    3) In lines starting with ##contig it must hold that assembly=GRCh38[...]
    4) The CHROM field must be in the format "chr##"

    Args:
        vcf_file_path (str): Path to the VCF file to validate

    Returns:
        dict: Dictionary with validation results and any errors found
    """
    results = {
        "is_valid": True,
        "errors": {
            "missing_headers": [],
            "field_count_mismatch": [],
            "invalid_contig_assembly": [],
            "invalid_chrom_format": [],
        },
    }

    has_fileformat = False
    has_chrom_header = False
    expected_field_count = 0

    with open(vcf_file_path, "r") as vcf_file:
        line_number = 0

        for line in vcf_file:
            line_number += 1
            line = line.strip()

            if not line:
                continue

            # Check for fileformat header
            if line.startswith("##fileformat="):
                has_fileformat = True

            # Check contig assembly format
            elif line.startswith("##contig="):
                if "assembly=GRCh38" not in line:
                    results["errors"]["invalid_contig_assembly"].append(line_number)
                    results["is_valid"] = False

            # Check for CHROM header line and count fields
            elif line.startswith("#CHROM"):
                has_chrom_header = True
                fields = line.split("\t")
                expected_field_count = len(fields)

            elif not line.startswith("#"):
                fields = line.split("\t")

                # Check field count
                if len(fields) != expected_field_count:
                    results["errors"]["field_count_mismatch"].append(
                        f"Line {line_number}: expected {expected_field_count} fields, got {len(fields)}"
                    )
                    results["is_valid"] = False

                # Check CHROM format
                if not fields[0].startswith("chr"):
                    results["errors"]["invalid_chrom_format"].append(
                        f"Line {line_number}: '{fields[0]}' does not match format 'chr##'"
                    )
                    results["is_valid"] = False

    if not has_fileformat:
        results["errors"]["missing_headers"].append("##fileformat")
        results["is_valid"] = False

    if not has_chrom_header:
        results["errors"]["missing_headers"].append("#CHROM")
        results["is_valid"] = False

    return results


def validate_type(value, field_type):
    """
    Validates if a value matches the expected VCF field type.
    
    Args:
        value (str): The value to validate
        field_type (str): The expected type (Integer, Float, Flag, Character, String)
        
    Returns:
        bool: True if the value matches the expected type, False otherwise
    """
    if field_type == 'Integer':
        try:
            int(value)
            return True
        except ValueError:
            return False
    elif field_type == 'Float':
        try:
            float(value)
            return True
        except ValueError:
            return False
    elif field_type == 'Flag':
        return value == '' or value is None
    elif field_type == 'Character':
        return len(value) == 1
    elif field_type == 'String':
        return True
    else:
        return False


def data_type_consistency(vcf_file_path):
    """
    Validates data type consistency in INFO and FORMAT fields against their definitions in the header.
    
    Checks:
    1) For each key in the INFO column, verifies that values match the Type specified in ##INFO header
    2) For each key in the FORMAT column, verifies that values match the Type specified in ##FORMAT header
    
    Args:
        vcf_file_path (str): Path to the VCF file to validate
        
    Returns:
        dict: Dictionary with validation results and any errors found
    """
    results = {
        'is_valid': True,
        'errors': {
            'info_type_mismatch': [],
            'format_type_mismatch': [],
            'undefined_info_field': [],
            'undefined_format_field': []
        }
    }
    
    info_fields = {}
    format_fields = {}
    
    with open(vcf_file_path, 'r') as vcf_file:
        header_line = None
        column_indices = {}
        
        for line_number, line in enumerate(vcf_file, 1):
            line = line.strip()
            
            # Process header lines
            if line.startswith('##'):
                # Extract INFO field definitions
                if line.startswith('##INFO=<'):
                    # Parse the INFO field definition
                    match = re.search(r'ID=([^,]+).*Type=([^,]+)', line)
                    if match:
                        field_id, field_type = match.groups()
                        info_fields[field_id] = field_type
                
                # Extract FORMAT field definitions
                elif line.startswith('##FORMAT=<'):
                    # Parse the FORMAT field definition
                    match = re.search(r'ID=([^,]+).*Type=([^,]+)', line)
                    if match:
                        field_id, field_type = match.groups()
                        format_fields[field_id] = field_type
                
                continue
            
            # Process column header line
            if line.startswith('#CHROM'):
                header_line = line
                fields = line.split('\t')
                
                # Store column indices for easier access
                for i, field in enumerate(fields):
                    column_indices[field] = i
                
                continue
            
            # Skip empty lines or remaining header lines
            if not line or line.startswith('#'):
                continue
            
            # Process data lines
            fields = line.split('\t')
            
            # Validate INFO field
            if 'INFO' in column_indices and column_indices['INFO'] < len(fields):
                info_idx = column_indices['INFO']
                info_field = fields[info_idx]
                
                # Skip if INFO field is empty or just a dot
                if info_field and info_field != '.':
                    # Process each key-value pair in INFO
                    for item in info_field.split(';'):
                        # Handle flags (no equals sign)
                        if '=' not in item:
                            if item not in info_fields:
                                results['errors']['undefined_info_field'].append(
                                    f"Line {line_number}: Undefined INFO field '{item}'"
                                )
                                results['is_valid'] = False
                            elif info_fields[item] != 'Flag':
                                results['errors']['info_type_mismatch'].append(
                                    f"Line {line_number}: INFO field '{item}' should be a Flag"
                                )
                                results['is_valid'] = False
                            continue
                            
                        # Handle key-value pairs
                        key, value = item.split('=', 1)
                        
                        # Check if key is defined in header
                        if key not in info_fields:
                            results['errors']['undefined_info_field'].append(
                                f"Line {line_number}: Undefined INFO field '{key}'"
                            )
                            results['is_valid'] = False
                            continue
                        
                        # Get expected type
                        expected_type = info_fields[key]
                        
                        # Handle comma-separated values - validate each value individually
                        for val in value.split(','):
                            # Skip missing values (represented as '.')
                            if val == '.':
                                continue
                                
                            # Validate the value against its expected type
                            if not validate_type(val, expected_type):
                                results['errors']['info_type_mismatch'].append(
                                    f"Line {line_number}: INFO field '{key}' value '{val}' does not match type '{expected_type}'"
                                )
                                results['is_valid'] = False
            
            # Validate FORMAT field and sample data
            if 'FORMAT' in column_indices and column_indices['FORMAT'] < len(fields):
                format_idx = column_indices['FORMAT']
                format_field = fields[format_idx]
                format_keys = format_field.split(':')
                
                # Check if all FORMAT keys are defined
                for key in format_keys:
                    if key not in format_fields:
                        results['errors']['undefined_format_field'].append(
                            f"Line {line_number}: Undefined FORMAT field '{key}'"
                        )
                        results['is_valid'] = False
                
                # Validate sample data if available
                if format_idx + 1 < len(fields):
                    sample_field = fields[format_idx + 1]
                    sample_values = sample_field.split(':')
                    
                    # Check if number of values matches number of keys
                    if len(sample_values) != len(format_keys):
                        continue  # Skip validation if counts don't match
                    
                    # Validate each value against its expected type
                    for i, (key, value) in enumerate(zip(format_keys, sample_values)):
                        if key in format_fields and value != '.':
                            expected_type = format_fields[key]
                            
                            # Special handling for GT field
                            if key == 'GT':
                                # GT can be in the form 0/1, 0|1, etc.
                                if not (('/' in value and all(c.isdigit() or c == '/' for c in value)) or 
                                       ('|' in value and all(c.isdigit() or c == '|' for c in value))):
                                    results['errors']['format_type_mismatch'].append(
                                        f"Line {line_number}: FORMAT field 'GT' value '{value}' is not a valid genotype"
                                    )
                                    results['is_valid'] = False
                            # Handle other fields that might have comma-separated values
                            else:
                                # Some FORMAT fields can also have comma-separated values
                                for subval in value.split(','):
                                    if subval != '.' and not validate_type(subval, expected_type):
                                        results['errors']['format_type_mismatch'].append(
                                            f"Line {line_number}: FORMAT field '{key}' value '{subval}' does not match type '{expected_type}'"
                                        )
                                        results['is_valid'] = False
    
    return results


def check_missing_values(vcf_file_path):
    """
    Checks for missing values in key VCF columns and validates FORMAT and genotype data.

    Checks:
    1) No missing values ('.') in CHROM, POS, REF, ALT columns
    2) FORMAT column starts with 'GT'
    3) Sample column has genotype data in the form X/Y where X,Y are integers

    Args:
        vcf_file_path (str): Path to the VCF file to validate

    Returns:
        dict: Dictionary with validation results and any errors found
    """
    results = {
        "is_valid": True,
        "errors": {"missing_values": [], "invalid_format": [], "invalid_genotype": []},
    }

    with open(vcf_file_path, "r") as vcf_file:
        header_line = None
        column_indices = {}

        for line_number, line in enumerate(vcf_file, 1):
            line = line.strip()

            # Skip header lines
            if line.startswith("##"):
                continue

            # Process column header line
            if line.startswith("#CHROM"):
                header_line = line
                columns = line.split("\t")

                # Get indices of required columns
                for i, col in enumerate(columns):
                    column_indices[col] = i

                # Ensure required columns exist
                required_columns = ["#CHROM", "POS", "REF", "ALT", "FORMAT"]
                for col in required_columns:
                    if col not in column_indices:
                        results["errors"]["missing_values"].append(
                            f"Required column '{col}' not found in header"
                        )
                        results["is_valid"] = False

                continue

            # Skip if header line wasn't found
            if not header_line:
                continue

            # Process data lines
            fields = line.split("\t")

            # Check for missing values in required fields
            for col in ["#CHROM", "POS", "REF", "ALT"]:
                if col in column_indices:
                    idx = column_indices[col]
                    if idx < len(fields) and fields[idx] == ".":
                        results["errors"]["missing_values"].append(
                            f"Line {line_number}: Missing value in {col} column"
                        )
                        results["is_valid"] = False

            # Check FORMAT field starts with GT
            if "FORMAT" in column_indices:
                format_idx = column_indices["FORMAT"]
                if format_idx < len(fields):
                    format_field = fields[format_idx]
                    format_values = format_field.split(":")

                    if not format_values or format_values[0] != "GT":
                        results["errors"]["invalid_format"].append(
                            f"Line {line_number}: FORMAT field does not start with 'GT'"
                        )
                        results["is_valid"] = False

                    # Check sample genotype data
                    if format_idx + 1 < len(fields):
                        genotype_field = fields[format_idx + 1]
                        genotype_values = genotype_field.split(":")

                        if not genotype_values:
                            results["errors"]["invalid_genotype"].append(
                                f"Line {line_number}: Missing genotype data"
                            )
                            results["is_valid"] = False
                        else:
                            # Check genotype format
                            gt = genotype_values[0]
                            if not (
                                ("/" in gt and all(c.isdigit() or c == "/" for c in gt))
                                or (
                                    "|" in gt
                                    and all(c.isdigit() or c == "|" for c in gt)
                                )
                            ):
                                results["errors"]["invalid_genotype"].append(
                                    f"Line {line_number}: Invalid genotype format '{gt}'"
                                )
                                results["is_valid"] = False

    return results


def main():
    parser = argparse.ArgumentParser(description='Validate VCF files')
    parser.add_argument('--input_dir', type=str, default='data/preprocessed_vcf',
                        help='Directory containing VCF files to validate')
    args = parser.parse_args()
    
    # Find all VCF files in the input directory
    vcf_files = glob.glob(os.path.join(args.input_dir, '*.vcf'))
    
    if not vcf_files:
        print(f"No VCF files found in {args.input_dir}")
        return

    print(f"Found {len(vcf_files)} VCF files to validate")

    # Process each VCF file
    for vcf_file in vcf_files:
        print(f"\nValidating {os.path.basename(vcf_file)}...")
        
        print("Checking header consistency...")
        header_results = vcf_header_consistency(vcf_file)
        
        if header_results['is_valid']:
            print("✓ VCF file passed all header consistency checks.")
        else:
            print("✗ VCF file failed some header consistency checks:")
            for error_type, errors in header_results['errors'].items():
                if errors:
                    print(f"  {error_type.replace('_', ' ').title()}:")
                    for error in errors:
                        print(f"    - {error}")
        
        print("Checking for missing values...")
        missing_results = check_missing_values(vcf_file)
        
        if missing_results['is_valid']:
            print("✓ VCF file passed all missing value checks.")
        else:
            print("✗ VCF file failed some missing value checks:")
            for error_type, errors in missing_results['errors'].items():
                if errors:
                    print(f"  {error_type.replace('_', ' ').title()}:")
                    for error in errors:
                        print(f"    - {error}")
        
        print("Checking data type consistency...")
        type_results = data_type_consistency(vcf_file)
        
        if type_results['is_valid']:
            print("✓ VCF file passed all data type consistency checks.")
        else:
            print("✗ VCF file failed some data type consistency checks:")
            for error_type, errors in type_results['errors'].items():
                if errors:
                    print(f"  {error_type.replace('_', ' ').title()}:")
                    for error in errors:
                        print(f"    - {error}")


if __name__ == "__main__":
    main()

# TODO: Check with synthetic VCF files with different types of errors.
# TODO: Make output a json file.
