import pandas as pd
import numpy as np
import argparse
import sys

def convert_dtype(dtype, value):
    try:
        if np.issubdtype(dtype, np.number):
            return np.dtype(dtype).type(value)
        elif np.issubdtype(dtype, np.bool_):
            return value.lower() in ['true', '1']
        else:  # string or object dtype
            return value
    except ValueError:
        print(f"Error: Could not convert {value} to dtype {dtype}.")
        sys.exit(1)

def select_points_in_sphere(csv_file, center, radius, labels=None):
    # Read the CSV file into a dataframe
    if not csv_file.endswith('.csv'):
        print("Warning: The input file does not have a .csv extension, but will try to read it anyway.")

    try:
        df = pd.read_csv(csv_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: File {csv_file} not found.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: File {csv_file} is empty or not a valid CSV.")
        sys.exit(1)

    # Check for required columns
    required_cols = ['Alpha', 'Beta', 'Gamma']
    for col in required_cols:
        if col not in df.columns:
            print(f"Error: Column '{col}' not found in the file {csv_file}.")
            sys.exit(1)

    # Ensure that the required columns are numeric
    for col in required_cols:
        if not np.issubdtype(df[col].dtype, np.number):
            try:
                df[col] = pd.to_numeric(df[col], errors='raise')
            except ValueError:
                print(f"Error: Column '{col}' must be numeric.")
                sys.exit(1)

    # Compute distance from the specified center
    df['dis'] = np.sqrt((df['Alpha'] - center[0])**2 + (df['Beta'] - center[1])**2 + (df['Gamma'] - center[2])**2)

    # Build the condition for filtering
    condition = (df['dis'] <= radius)

    # If labels are provided, they come in pairs (column, value)
    if labels:
        if len(labels) % 2 != 0:
            print("Error: Labels must be provided in pairs (column value).")
            sys.exit(1)
        for i in range(0, len(labels), 2):
            column = labels[i]
            value = labels[i + 1]

            # Check if the column exists in the dataframe
            if column not in df.columns:
                print(f"Error: Column '{column}' not found in the file.")
                sys.exit(1)

            converted_value = convert_dtype(df[column].dtype, value)
            condition &= (df[column] == converted_value)

    # Apply the filtering condition
    filtered_df = df[condition].copy()

    # Drop the 'dis' column since it is only for internal calculation
    filtered_df.drop(columns=['dis'], inplace=True)

    particle_number = len(filtered_df)
    extracted_condition = f"center_{center[0]}_{center[1]}_{center[2]}_radius_{radius}"
    print(f"Extracted {particle_number} particles inside the sphere defined by {extracted_condition}")

    return filtered_df

def main():
    parser = argparse.ArgumentParser(description="Select points within a specified spherical region from a CSV file.")
    parser.add_argument('--i', dest='input', required=True, help="Path to the input CSV file.")
    parser.add_argument('--c', dest='center', nargs=3, type=float, required=True, help="Center of the sphere (Alpha Beta Gamma).")
    parser.add_argument('--r', dest='radius', type=float, required=True, help="Radius of the sphere.")
    parser.add_argument('--l', dest='labels', nargs='*', help="Pairs of column labels and values to filter on. Example: --l Class 4")
    parser.add_argument('--o', dest='output', help="Name of the output CSV file.")

    args = parser.parse_args()

    result_df = select_points_in_sphere(args.input, args.center, args.radius, args.labels)

    if args.output and result_df is not None:
        result_df.to_csv(args.output, sep='\t', index=False)
        print(f"Filtered results saved to {args.output}")
    elif result_df is not None and args.output is None:
        print("No output file specified. Results are not saved.")

if __name__ == '__main__':
    main()