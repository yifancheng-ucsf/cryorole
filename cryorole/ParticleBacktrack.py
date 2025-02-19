import pandas as pd
import argparse
import re
import sys
import os

def read_relion_star(star):

    optics_metadata_index = {}
    optics_metadata_index_found = 0
    optics_metadata_table = []
    particles_metadata_index = {}
    particles_metadata_index_found = 0
    particles_metadata_table = []

    optics_metadata_starts = False
    particles_metadata_starts = False
    header = ""

    if not os.path.isfile(star):
        raise FileNotFoundError(f"Cannot open star file {star}.")

    with open(star, 'r') as starfile:
        for num, line in enumerate(starfile):
            sline = line.strip()

            if not sline or (sline[0] == "#"):
                header += line
                continue

            if re.match(r'data_', sline):
                if re.match(r'data_optics', sline):
                    optics_metadata_starts = True
                    particles_metadata_starts = False
                    header += line
                elif re.match(r'data_particles', sline):
                    optics_metadata_starts = False
                    particles_metadata_starts = True
                    header += line
                else:
                    raise ValueError(
                        f"Found data label not either _optics or _particles ... line {num}"
                    )
                continue

            if re.match(r'loop_', sline):
                header += line
                continue

            index_found = re.match(r'(_rln\w+).*', sline)
            if index_found:
                if optics_metadata_starts:
                    optics_metadata_index[index_found.group(1)] = optics_metadata_index_found
                    optics_metadata_index_found += 1
                    header += line
                elif particles_metadata_starts:
                    particles_metadata_index[index_found.group(1)] = particles_metadata_index_found
                    particles_metadata_index_found += 1
                    header += line
                else:
                    raise ValueError(
                        f"Found index but not in either data_optics or in data_particles ... line {num}"
                    )
                continue

            if optics_metadata_starts:
                optics_metadata_table.append(sline.split())
                header += line
            elif particles_metadata_starts:
                particles_metadata_table.append(sline.split())
            else:
                raise ValueError(
                    f"Found a possible metadata record but not in either data_optics or in data_particles ... line {num}"
                )

    particles_metadata_df = pd.DataFrame(particles_metadata_table)
    optics_metadata_df = pd.DataFrame(optics_metadata_table)

    return header, particles_metadata_df, optics_metadata_df, particles_metadata_index, optics_metadata_index

def backtrack_particles(extracted_file, star_file, output_file):

    if not os.path.isfile(extracted_file):
        print(f"Error: Extracted file {extracted_file} does not exist.")
        sys.exit(1)

    # Try reading the 'ID' column from the extracted file
    try:
        extracted_df = pd.read_csv(extracted_file, sep='\t')
    except Exception as e:
        print(f"Error reading extracted file {extracted_file}: {e}")
        sys.exit(1)

    if 'ID' not in extracted_df.columns:
        print("Error: 'ID' column not found in extracted file.")
        sys.exit(1)

    # Ensure ID is numeric
    if not pd.api.types.is_numeric_dtype(extracted_df['ID']):
        try:
            extracted_df['ID'] = pd.to_numeric(extracted_df['ID'], errors='raise')
        except ValueError:
            print("Error: 'ID' column must be numeric.")
            sys.exit(1)

    extracted_indices = extracted_df['ID'].values

    try:
        header, particles_metadata_df, _, _, _ = read_relion_star(star_file)
    except Exception as e:
        print(f"Error reading star file {star_file}: {e}")
        sys.exit(1)

    # IDs are 1-based indexing, while DataFrame .iloc uses 0-based indexing
    try:
        extracted_particles_metadata = particles_metadata_df.iloc[extracted_indices - 1]
    except IndexError:
        print("Error: One or more IDs in extracted file are out of range for the original star file.")
        sys.exit(1)

    # Write output
    try:
        with open(output_file, 'w') as f:
            f.write(header)
            for row in extracted_particles_metadata.values:
                f.write('\t'.join(map(str, row)) + '\n')
        print(f"Backtracked particles saved to {output_file}")
    except Exception as e:
        print(f"Error writing output file {output_file}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Backtrack particles from a filtered CSV file to the original star file.")
    parser.add_argument('--i', dest='input', help='Input extracted data file name', required=True)
    parser.add_argument('--s', dest='star', help='Input original star data file name', required=True)
    parser.add_argument('--o', dest='output', help='Output file name', required=True)
    args = parser.parse_args()

    extracted_file = args.input
    star_file = args.star
    output_file = args.output
    backtrack_particles(extracted_file, star_file, output_file)

if __name__ == "__main__":
    main()
