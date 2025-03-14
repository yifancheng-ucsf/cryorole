import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import math
import sys

def plot_projections(data, pairs, c_column, vmin, vmax, output_prefix, limits):
    num_plots = len(pairs)
    fig, axes = plt.subplots(1, num_plots, figsize=(20, 9))

    if num_plots == 1:
        axes = [axes]  # Make it iterable

    for ax, (x_column, y_column) in zip(axes, pairs):
        sc = ax.scatter(data[x_column], data[y_column], c=data[c_column], cmap='rainbow_r', vmin=vmin, vmax=vmax, s=1,
                        rasterized=True)
        ax.set_xlabel(x_column)
        ax.set_ylabel(y_column)
        ax.set_aspect('equal')
        ax.set_title(f'{x_column}-{y_column} landscape')

        # Set axis limits if provided
        x_limits = limits.get(x_column) if limits else None
        y_limits = limits.get(y_column) if limits else None
        if x_limits:
            ax.set_xlim(x_limits)
        if y_limits:
            ax.set_ylim(y_limits)
        ax.set_aspect('equal')

    plt.tight_layout()
    cbar_ax = fig.add_axes([0.35, 0.06, 0.3, 0.03])
    cbar = fig.colorbar(sc, cax=cbar_ax, ax=axes, orientation='horizontal', fraction=0.02, pad=0.02)
    cbar.set_label(c_column)
    output_filename1 = f"{output_prefix}.pdf"
    output_filename2 = f"{output_prefix}.png"
    output_filename3 = f"{output_prefix}.svg"
    plt.savefig(output_filename1, bbox_inches='tight', dpi=300)
    plt.savefig(output_filename2, bbox_inches='tight', dpi=300)
    plt.savefig(output_filename3, bbox_inches='tight', dpi=300)
    plt.close()

def plot_histogram(data, column, bins, output_prefix):
    plt.figure(figsize=(8, 6))
    weights = (np.ones_like(data[column]) / len(data[column])) * 100  # For percentage
    plt.hist(data[column], bins=bins, weights=weights, edgecolor='black')
    plt.xlabel(column)
    plt.ylabel('Percentage')
    plt.title(f'Histogram of {column} (Percentage)')
    plt.tight_layout()
    output_filename = f"{output_prefix}_{column}_histogram.pdf"
    plt.savefig(output_filename, dpi=300)
    plt.close()
    print(f"Histogram saved to {output_filename}")

def compute_azimuth_elevation(rotvecs):
    norms = np.linalg.norm(rotvecs, axis=1, keepdims=True)
    # Avoid division by zero
    norms[norms == 0] = 1
    # Normalize the rotation vectors
    unit_vectors = rotvecs / norms
    u_x = unit_vectors[:, 0]
    u_y = unit_vectors[:, 1]
    u_z = unit_vectors[:, 2]
    # Compute azimuth (phi)
    azimuth = np.arctan2(u_y, u_x)
    # Compute elevation (theta)
    elevation = np.arcsin(u_z)
    return azimuth, elevation

def plot_azimuth_elevation(data, azimuth_column, elevation_column, c_column, vmin, vmax, output_filename):
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(data[azimuth_column], data[elevation_column], c=data[c_column], cmap='rainbow_r',
                     vmin=vmin, vmax=vmax, s=1, rasterized=True)
    plt.xlabel('Azimuth (radians)')
    plt.ylabel('Elevation (radians)')
    plt.title('Azimuth-Elevation Map')
    plt.colorbar(sc, label=c_column)
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.close()
    print(f"Azimuth-Elevation map saved to {output_filename}")

def main():
    parser = argparse.ArgumentParser(description="Visualization of Euler angles and Rotation vectors from CSV files.")
    parser.add_argument('--i', dest='input', required=True, help='Input CSV file generated by Script One.')
    parser.add_argument('--o', dest='output', default='visualization', help='Output prefix for saved plots.')
    parser.add_argument('--bins', '-b', type=int, default=50, help='Number of bins for histograms.')
    parser.add_argument('--vmin', type=float, help='Minimum value for color scale.')
    parser.add_argument('--vmax', type=float, help='Maximum value for color scale.')
    parser.add_argument('--t', dest='threshold_based_on_RND', type=float, help='Threshold for RND to filter data.')
    parser.add_argument('--p', dest='threshold_based_on_percentage', type=float,
                        help='Percentage (0 to 1) to select top data based on RND values.')
    parser.add_argument('--generate_histogram', action='store_true',
                        help="Generate the histogram of the columns.")
    parser.add_argument('--angle_limits', nargs=6, type=float, metavar=(
    'Alpha_MIN', 'Alpha_MAX', 'Beta_MIN', 'Beta_MAX', 'Gamma_MIN', 'Gamma_MAX'),
                        help='Limits for Euler angles: Euler01_min Euler01_max Euler02_min Euler02_max Euler03_min Euler03_max.',
                        default=[-180.0, 180.0, -180.0, 180.0, -180.0, 180.0])
    parser.add_argument('--rotvec_limits', nargs=6, type=float, metavar=('ROTVEC_X_MIN', 'ROTVEC_X_MAX', 'ROTVEC_Y_MIN', 'ROTVEC_Y_MAX', 'ROTVEC_Z_MIN', 'ROTVEC_Z_MAX'),
                        help='Limits for Rotation vectors: RotVec_X_min RotVec_X_max RotVec_Y_min RotVec_Y_max RotVec_Z_min RotVec_Z_max.',default=[-3.14, 3.14, -3.14, 3.14, -3.14, 3.14])
    parser.add_argument('--generate_azimuth_elevation', action='store_true',
                        help='Generate the Azimuth-Elevation map.')
    args = parser.parse_args()

    # Read the CSV file
    if not os.path.isfile(args.input):
        print(f"Error: File {args.input} does not exist.")
        sys.exit(1)

    data = pd.read_csv(args.input, sep="\t")

    required_columns = ['RND', 'Alpha', 'Beta', 'Gamma', 'RotVec_X', 'RotVec_Y', 'RotVec_Z', 'Rotation_Angle',
                        'Distance']
    for col in required_columns:
        if col not in data.columns:
            print(f"Error: Required column '{col}' not found in {args.input}.")
            sys.exit(1)

    data_original = data.copy()

    # Determine vmin and vmax for color scale
    rnd_min = data['RND'].min()
    rnd_max = data['RND'].max()

    print(f"The minimum value in column RND is {rnd_min}")
    print(f"The maxmum value in column RND is {rnd_max}")

    default_vmax = math.ceil(rnd_max * 10) / 10

    default_threshold = rnd_max / 3

    default_threshold_c = math.ceil(default_threshold * 10) / 10

    # Determine RND threshold to filter data
    rnd_threshold = None
    if args.threshold_based_on_percentage is not None:
        percentage = args.threshold_based_on_percentage
        if not (0 < percentage <= 1):
            print("Error: --threshold_based_on_percentage must be between 0 and 1.")
            sys.exit(1)
        rnd_threshold = data['RND'].quantile(1 - percentage)
        rnd_threshold = math.ceil(rnd_threshold * 10) / 10
        print(f"RND threshold based on top {percentage * 100}%: {rnd_threshold}")
    elif args.threshold_based_on_RND is not None:
        rnd_threshold = args.threshold_based_on_RND
        print(f"Using fixed RND threshold: {rnd_threshold}")
    else:
        # Default threshold: RND_MAX / 3
        rnd_threshold = default_threshold_c
        print(f"Using default RND threshold (RND_MAX / 3): {rnd_threshold}")


    data = data[data['RND'] > rnd_threshold]
    print(f"Data filtered with RND threshold > {rnd_threshold}")

    remaining_points = len(data)
    print(f"Remaining data points after filtering: {remaining_points}")

    if data.empty:
        print(f"No data points remain after filtering with RND threshold > {rnd_threshold}.")
        return

    data = data.sort_values('RND', ascending=True)

    vmin = args.vmin if args.vmin is not None else rnd_threshold
    vmax = args.vmax if args.vmax is not None else default_vmax

    print(f"Color scale vmin: {vmin}, vmax: {vmax}")

    if args.generate_histogram:
        # Plot histograms of RND, Rotation_Angle, and Distance
        histogram_threshold_output_filename = f"{args.output}_threshold_filter"
        histogram_all_output_filename = f"{args.output}_all"

        plot_histogram(data, 'RND', args.bins, histogram_threshold_output_filename)
        plot_histogram(data, 'Rotation_Angle', args.bins, histogram_threshold_output_filename)
        plot_histogram(data, 'Alpha', args.bins, histogram_threshold_output_filename)
        plot_histogram(data, 'Beta', args.bins, histogram_threshold_output_filename)
        plot_histogram(data, 'Gamma', args.bins, histogram_threshold_output_filename)
        plot_histogram(data, 'Distance', args.bins, histogram_threshold_output_filename)

        plot_histogram(data_original, 'Rotation_Angle', args.bins, histogram_all_output_filename)
        plot_histogram(data_original, 'Distance', args.bins, histogram_all_output_filename)

    # Plot projections of Euler angles
    euler_columns = ['Alpha', 'Beta', 'Gamma']
    euler_pairs = [('Alpha', 'Beta'), ('Gamma', 'Beta'), ('Alpha', 'Gamma')]

    # Apply angle limits if specified
    euler_limits = {}
    if args.angle_limits:
        euler_limits = {
            'Alpha': (args.angle_limits[0], args.angle_limits[1]),
            'Beta': (args.angle_limits[2], args.angle_limits[3]),
            'Gamma': (args.angle_limits[4], args.angle_limits[5])
        }
        print("Setting Euler angle axis limits:")
        for col in euler_columns:
            min_limit, max_limit = euler_limits[col]
            print(f" - {col}: {min_limit} to {max_limit} degrees")

    euler_output_filename = f"{args.output}_euler_projections"
    plot_projections(data, euler_pairs, 'RND', vmin, vmax, euler_output_filename, euler_limits)

    # Plot combined projections of Rotation vectors
    rotvec_columns = ['RotVec_X', 'RotVec_Y', 'RotVec_Z']
    rotvec_pairs = [('RotVec_X', 'RotVec_Y'), ('RotVec_Z', 'RotVec_Y'), ('RotVec_X', 'RotVec_Z')]

    # Apply rotation vector limits if specified
    rotvec_limits = {}
    if args.rotvec_limits:
        rotvec_limits = {
            'RotVec_X': (args.rotvec_limits[0], args.rotvec_limits[1]),
            'RotVec_Y': (args.rotvec_limits[2], args.rotvec_limits[3]),
            'RotVec_Z': (args.rotvec_limits[4], args.rotvec_limits[5])
        }
        print("Setting Rotation Vector axis limits:")
        for col in rotvec_columns:
            min_limit, max_limit = rotvec_limits[col]
            print(f" - {col}: {min_limit} to {max_limit}")

    rotvec_output_filename = f"{args.output}_rotvec_projections"
    plot_projections(data, rotvec_pairs, 'RND', vmin, vmax, rotvec_output_filename, rotvec_limits)

    if args.generate_azimuth_elevation:
        rotvecs = data[['RotVec_X', 'RotVec_Y', 'RotVec_Z']].values
        azimuth, elevation = compute_azimuth_elevation(rotvecs)
        data['Azimuth'] = azimuth
        data['Elevation'] = elevation

        # Plot the Azimuth-Elevation map
        az_el_output_filename = f"{args.output}_azimuth_elevation_map.pdf"
        plot_azimuth_elevation(data, 'Azimuth', 'Elevation', 'RND', vmin, vmax, az_el_output_filename)

        # Optionally, plot histograms of Azimuth and Elevation
        plot_histogram(data, 'Azimuth', args.bins, args.output)
        plot_histogram(data, 'Elevation', args.bins, args.output)

if __name__ == '__main__':
    main()
