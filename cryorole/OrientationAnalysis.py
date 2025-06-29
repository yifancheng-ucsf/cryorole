import numpy as np
import pandas as pd
import argparse
import re
from scipy.spatial.transform import Rotation as R
from mpi4py import MPI
from scipy.spatial import cKDTree
import gc

def read_relion_star(star):
    try:
        starfile = open(star, 'r')
    except:
        parser.error("cannot open star file {}.".format(star))

    optics_metadata_index = {}
    optics_metadata_index_found = 0
    optics_metadata_table = []
    particles_metadata_index = {}
    particles_metadata_index_found = 0
    particles_metadata_table = []
    optics_metadata_starts = False
    particles_metadata_starts = False

    for num, line in enumerate(starfile):
        sline = line.strip()

        if not sline or (sline[0] == "#"):
            continue

        if re.match(r'data_', sline):
            if re.match(r'data_optics', sline):
                optics_metadata_starts = True
                particles_metadata_starts = False
            elif re.match(r'data_particles', sline):
                optics_metadata_starts = False
                particles_metadata_starts = True
            else:
                parser.error("found data label not either _optics or _particles ... line {}".format(num))
            continue

        if re.match(r'loop_', sline):
            continue

        index_found = re.match(r'(_rln\w+).*', sline)
        if index_found:
            if optics_metadata_starts:
                optics_metadata_index[index_found.group(1)] = optics_metadata_index_found
                optics_metadata_index_found += 1
            elif particles_metadata_starts:
                particles_metadata_index[index_found.group(1)] = particles_metadata_index_found
                particles_metadata_index_found += 1
            else:
                parser.error("found index but not in either data_optics or in data_particles ... line {}".format(num))
            continue

        if True:
            if optics_metadata_starts:
                optics_metadata_table.append(sline.split())
            elif particles_metadata_starts:
                particles_metadata_table.append(sline.split())
            else:
                parser.error('found a possible metadata record but not in either data_optics or in data_particles ... line {}'.format(num))

    starfile.close()

    return particles_metadata_table, optics_metadata_table, particles_metadata_index, optics_metadata_index

def EulerAngle2RotationMatrix(rot, tilt, psi):
    alpha = np.deg2rad(rot)
    beta = np.deg2rad(tilt)
    gamma = np.deg2rad(psi)

    ca = np.cos(alpha)
    cb = np.cos(beta)
    cg = np.cos(gamma)
    sa = np.sin(alpha)
    sb = np.sin(beta)
    sg = np.sin(gamma)
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa

    R_matrix = np.array([
        [cg * cc - sg * sa, cg * cs + sg * ca, -cg * sb],
        [-sg * cc - cg * sa, -sg * cs + cg * ca, sg * sb],
        [sc, ss, cb]
    ])
    return R_matrix

def Matrix2Eulerangle(A):
    if A.shape != (3, 3):
        raise ValueError("Matrix2Eulerangle: The Euler matrix is not 3x3")

    abs_sb = np.sqrt(A[0, 2] ** 2 + A[1, 2] ** 2)

    if abs_sb > 16 * np.finfo(float).eps:
        gamma = np.arctan2(A[1, 2], -A[0, 2])
        alpha = np.arctan2(A[2, 1], A[2, 0])

        if np.abs(np.sin(gamma)) < np.finfo(float).eps:
            sign_sb = np.sign(-A[0, 2] / np.cos(gamma))
        else:
            sign_sb = np.sign(A[1, 2]) if np.sin(gamma) > 0 else -np.sign(A[1, 2])

        beta = np.arctan2(sign_sb * abs_sb, A[2, 2])
    else:
        if np.sign(A[2, 2]) > 0:
            alpha = 0
            beta = 0
            gamma = np.arctan2(-A[1, 0], A[0, 0])
        else:
            alpha = 0
            beta = np.pi
            gamma = np.arctan2(A[1, 0], -A[0, 0])

    gamma = np.degrees(gamma)
    beta = np.degrees(beta)
    alpha = np.degrees(alpha)

    euler = [alpha, beta, gamma]

    return euler

def compute_rotation_angle(rotvec):
    angle_rad = np.linalg.norm(rotvec)
    angle_deg = np.degrees(angle_rad)
    return angle_deg

def compute_rnd(rotvecs, k, epsilon=1e-7):
    coords = np.vstack(rotvecs)
    tree = cKDTree(coords)
    distances, _ = tree.query(coords, k=k+1)
    local_distances = np.mean(distances[:, 1:], axis=1)

    global_distance = np.mean(distances[:, 1])

    rnd = global_distance / (local_distances + epsilon)
    return rnd

def outlier_exclusion(rnd_initial, outlier_method, multiplier=1.5):
    if outlier_method == 'IQR':
        Q1 = np.percentile(rnd_initial, 25)
        Q3 = np.percentile(rnd_initial, 75)
        IQR = Q3 - Q1
        upper_bound = Q3 + multiplier * IQR
        outlier_mask = rnd_initial > upper_bound
    elif outlier_method == 'Z-score':
        mean_rnd = np.mean(rnd_initial)
        std_rnd = np.std(rnd_initial)
        z_scores = (rnd_initial - mean_rnd) / std_rnd
        threshold = 5
        outlier_mask = z_scores > threshold
    else:
        raise ValueError("Invalid outlier_method. Choose 'IQR' or 'Z-score'.")

    num_outliers = np.sum(outlier_mask)
    print(f"Identified {num_outliers} outliers out of {len(rnd_initial)} data points.")

    non_outlier_mask = ~outlier_mask

    return non_outlier_mask, outlier_mask

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--s1', dest='star1', required=True, help='Input the first RELION particle star file.')
    parser.add_argument('--s2', dest='star2', required=True, help='Input the second RELION particle star file.')
    parser.add_argument('--o', dest='output', required=True, help='Output root name.')
    parser.add_argument('--k', dest='k_sample', type=int, default=30, help='The nearest neighbors of each point (excluding itself, default: 30).')
    parser.add_argument('--apply_rotation', nargs=3, type=float, metavar=('ROT', 'TILT', 'PSI'),
                        help='Apply a known rotation (ROT, TILT, PSI) to all particles before processing.')
    parser.add_argument('--outlier_method', default='IQR', help="Method to identify outliers ('IQR' or 'Z-score').")
    parser.add_argument('--autoalign', action='store_true', help="If set, compute and apply Inverse Mean Rotation (IMR) to align data.")

    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Read star files on root
    if rank == 0:
        particles_metadata_table1, optics_metadata_table1, particles_metadata_index1, optics_metadata_index1 = read_relion_star(args.star1)
        particles_metadata_table2, optics_metadata_table2, particles_metadata_index2, optics_metadata_index2 = read_relion_star(args.star2)
    else:
        particles_metadata_table1 = None
        particles_metadata_table2 = None
        particles_metadata_index1 = None
        particles_metadata_index2 = None

    particles_metadata_table1 = comm.bcast(particles_metadata_table1, root=0)
    particles_metadata_table2 = comm.bcast(particles_metadata_table2, root=0)
    particles_metadata_index1 = comm.bcast(particles_metadata_index1, root=0)
    particles_metadata_index2 = comm.bcast(particles_metadata_index2, root=0)

    ParticleNumber1 = len(particles_metadata_table1)
    ParticleNumber2 = len(particles_metadata_table2)

    if ParticleNumber1 != ParticleNumber2:
        if rank == 0:
            print('Error: The number of particles in the two star files are not the same.')
        comm.Barrier()
        MPI.Finalize()
        return

    ParticleNumber = ParticleNumber1

    if args.apply_rotation:
        known_rot = args.apply_rotation
        R_known = EulerAngle2RotationMatrix(*known_rot)
    else:
        R_known = np.identity(3)

    indices = np.array_split(np.arange(ParticleNumber), size)
    local_indices = indices[rank]

    local_data_before = [None] * ParticleNumber
    local_rotation_vectors = []

    # Compute local data
    for i in local_indices:
        EulerAngles1 = [
            float(particles_metadata_table1[i][particles_metadata_index1['_rlnAngleRot']]),
            float(particles_metadata_table1[i][particles_metadata_index1['_rlnAngleTilt']]),
            float(particles_metadata_table1[i][particles_metadata_index1['_rlnAnglePsi']])
        ]

        EulerAngles2 = [
            float(particles_metadata_table2[i][particles_metadata_index2['_rlnAngleRot']]),
            float(particles_metadata_table2[i][particles_metadata_index2['_rlnAngleTilt']]),
            float(particles_metadata_table2[i][particles_metadata_index2['_rlnAnglePsi']])
        ]

        Origin1 = [
            float(particles_metadata_table1[i][particles_metadata_index1['_rlnOriginXAngst']]),
            float(particles_metadata_table1[i][particles_metadata_index1['_rlnOriginYAngst']])
        ]

        Origin2 = [
            float(particles_metadata_table2[i][particles_metadata_index2['_rlnOriginXAngst']]),
            float(particles_metadata_table2[i][particles_metadata_index2['_rlnOriginYAngst']])
        ]

        OriginX = Origin1[0] - Origin2[0]
        OriginY = Origin1[1] - Origin2[1]
        OriginX = round(OriginX, 6)
        OriginY = round(OriginY, 6)
        translationdistance = np.sqrt(OriginX ** 2 + OriginY ** 2)

        R1 = EulerAngle2RotationMatrix(EulerAngles1[0], EulerAngles1[1], EulerAngles1[2])
        R2 = EulerAngle2RotationMatrix(EulerAngles2[0], EulerAngles2[1], EulerAngles2[2])

        # Apply known rotation
        R2 = np.matmul(R2, R_known)

        # Compute Orientation_Relations
        R1_inverse = np.linalg.inv(R1)
        Orientation_Relations = np.matmul(R1_inverse, R2)

        # Convert the Orientation_Relations to rotvec and euler angles in zyx
        rotation = R.from_matrix(Orientation_Relations)
        rotvec = rotation.as_rotvec()
        rotation_angle = compute_rotation_angle(rotvec)
        euler_angles = rotation.as_euler('zyx', degrees=True)

        idx = i
        data_before = {
            'Alpha': euler_angles[0],
            'Beta': euler_angles[1],
            'Gamma': euler_angles[2],
            'RotVec_X': rotvec[0],
            'RotVec_Y': rotvec[1],
            'RotVec_Z': rotvec[2],
            'Rotation_Angle': rotation_angle,
            'Tx': OriginX,
            'Ty': OriginY,
            'Distance': translationdistance,
            'ID': i + 1
        }
        local_data_before[idx] = data_before
        local_rotation_vectors.append(rotvec)

    # Compute mean rotation vector only if autoalign is requested
    if args.autoalign:
        local_rotation_vectors_array = np.vstack(local_rotation_vectors)
        local_rotvec_sum = np.sum(local_rotation_vectors_array, axis=0)

        global_rotvec_sum = np.zeros(3)
        comm.Allreduce(local_rotvec_sum, global_rotvec_sum, op=MPI.SUM)

        N_local = len(local_rotation_vectors)
        N_total = comm.allreduce(N_local, op=MPI.SUM)

        mean_rotvec = global_rotvec_sum / N_total
        mean_rotation = R.from_rotvec(mean_rotvec)
        mean_rotation_inv = mean_rotation.inv()

        # Extract inverse mean rotation info
        inverse_mean_rotation_matrix = mean_rotation_inv.as_matrix()
        inverse_mean_euler_angles = mean_rotation_inv.as_euler('zyx', degrees=True)
        alpha_inv_mean = inverse_mean_euler_angles[0]
        beta_inv_mean = inverse_mean_euler_angles[1]
        gamma_inv_mean = inverse_mean_euler_angles[2]
        inverse_mean_rotvec = mean_rotation_inv.as_rotvec()
        inverse_mean_quat = mean_rotation_inv.as_quat()

    # Gather before-alignment data
    local_data_before = np.array(local_data_before, dtype=object)
    data_before_list = comm.gather(local_data_before, root=0)

    # Only compute after-alignment data if autoalign is requested
    if args.autoalign:
        local_data_after = [None] * ParticleNumber
        for i in local_indices:
            idx = i
            rotvec = np.array([
                local_data_before[idx]['RotVec_X'],
                local_data_before[idx]['RotVec_Y'],
                local_data_before[idx]['RotVec_Z']
            ])

            rotation = R.from_rotvec(rotvec)
            relative_rotation = mean_rotation_inv * rotation
            aligned_rotvec = relative_rotation.as_rotvec()
            rotation_angle_aligned = compute_rotation_angle(aligned_rotvec)
            euler_angles_aligned = relative_rotation.as_euler('zyx', degrees=True)

            data_after = local_data_before[idx].copy()
            data_after['Alpha'] = euler_angles_aligned[0]
            data_after['Beta'] = euler_angles_aligned[1]
            data_after['Gamma'] = euler_angles_aligned[2]
            data_after['RotVec_X'] = aligned_rotvec[0]
            data_after['RotVec_Y'] = aligned_rotvec[1]
            data_after['RotVec_Z'] = aligned_rotvec[2]
            data_after['Rotation_Angle'] = rotation_angle_aligned

            local_data_after[idx] = data_after

        local_data_after = np.array(local_data_after, dtype=object)
        data_after_list = comm.gather(local_data_after, root=0)
    else:
        data_after_list = None

    if rank == 0:
        # Combine before-alignment data
        final_data_before = [None] * ParticleNumber
        for proc_data_before in data_before_list:
            for idx, data in enumerate(proc_data_before):
                if data is not None:
                    final_data_before[idx] = data

        df_before = pd.DataFrame(final_data_before)
        df_before.sort_values('ID', inplace=True)

        # Compute RND before alignment
        k = args.k_sample
        rotvecs_before = df_before[['RotVec_X', 'RotVec_Y', 'RotVec_Z']].values
        rnd_before = compute_rnd(rotvecs_before, k)
        df_before['RND'] = rnd_before

        # Outlier exclusion
        outlier_method = args.outlier_method
        non_outlier_mask, outlier_mask = outlier_exclusion(rnd_before, outlier_method, multiplier=4)

        df_before_non_outliers = df_before[non_outlier_mask].copy()
        df_before_outliers = df_before[outlier_mask].copy()

        # Re-compute RND for non-outliers (before)
        rotvecs_before_non_outliers = df_before_non_outliers[['RotVec_X', 'RotVec_Y', 'RotVec_Z']].values
        rnd_before_non_outliers = compute_rnd(rotvecs_before_non_outliers, k)
        df_before_non_outliers['RND'] = rnd_before_non_outliers

        # Save before-alignment CSV
        output = args.output
        output_csv_before = f"{output}.csv"
        df_before_non_outliers.to_csv(output_csv_before, sep='\t', index=False)
        print(f"Data before alignment saved to {output_csv_before}")

        # Save outliers
        output_csv_outliers = f"{output}_outliers.csv"
        df_before_outliers.to_csv(output_csv_outliers, sep='\t', index=False)

        if args.autoalign:
            # Combine after-alignment data
            final_data_after = [None] * ParticleNumber
            for proc_data_after in data_after_list:
                for idx, data in enumerate(proc_data_after):
                    if data is not None:
                        final_data_after[idx] = data

            df_after = pd.DataFrame(final_data_after)
            df_after.sort_values('ID', inplace=True)

            # Compute RND after alignment for non-outliers
            df_after_non_outliers = df_after[non_outlier_mask].copy()
            rotvecs_after_non_outliers = df_after_non_outliers[['RotVec_X', 'RotVec_Y', 'RotVec_Z']].values
            rnd_after_non_outliers = compute_rnd(rotvecs_after_non_outliers, k)
            df_after_non_outliers['RND'] = rnd_after_non_outliers

            output_csv_after = f"{output}_after_alignment.csv"
            df_after_non_outliers.to_csv(output_csv_after, sep='\t', index=False)
            print(f"Data after alignment saved to {output_csv_after}")

            # Save inverse mean rotation info
            output_rotation_file = f"{args.output}_inverse_mean_rotation.txt"
            with open(output_rotation_file, 'w') as f:
                f.write("Inverse Mean Rotation Matrix:\n")
                f.write(np.array2string(inverse_mean_rotation_matrix, formatter={'float_kind': lambda x: "%.6f" % x}))
                f.write("\n\n")

                f.write("Inverse Mean Rotation (Euler Angles):\n")
                f.write(f"Alpha (α): {alpha_inv_mean:.6f} degrees\n")
                f.write(f"Beta (β): {beta_inv_mean:.6f} degrees\n")
                f.write(f"Gamma (γ): {gamma_inv_mean:.6f} degrees\n\n")

                f.write("Inverse Mean Rotation Vector:\n")
                f.write(np.array2string(inverse_mean_rotvec, formatter={'float_kind': lambda x: "%.6f" % x}))
                f.write("\n\n")

                f.write("Inverse Mean Rotation Quaternion:\n")
                f.write(np.array2string(inverse_mean_quat, formatter={'float_kind': lambda x: "%.6f" % x}))
                f.write("\n")

            print(f"Inverse Mean Rotation data saved to {output_rotation_file}")

            # Clean up large dataframes from memory
            del df_after, df_after_non_outliers, final_data_after, data_after_list

        # Clean up large dataframes and arrays from memory
        del df_before, df_before_non_outliers, df_before_outliers, final_data_before, data_before_list
        del rotvecs_before, rotvecs_before_non_outliers
        gc.collect()

    comm.Barrier()
    # MPI will finalize automatically when the script ends, freeing memory.

if __name__ == '__main__':
    main()
