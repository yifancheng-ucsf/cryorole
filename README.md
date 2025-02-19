# CryoROLE: Cryo-EM Relative Orientation LandscapE

**CryoROLE** is a python package designed to analyze the particle distribution of the relative orientation between two independently refined domains of a biological complex. The landscape is represented in angular (real) space. Detailed description of the orientational landscape can be found in the (bioRxiv link).

When two domains of a complex are independently refined to high-resolution, each particle in the dataset is assigned two set of Euler angles for both domains. Thus, there are two `.star` files in Relion (or two `.cs` files in CryoSPARC) corresponding to the same dataset in the final reconstruction. CryoROLE takes these two .star files as inputs and computes, for every particle, the relative orientation of the first domain with respect to the second. By plotting this orientations of every particle in an angular space, the program generates a particle distribution “landscape” that reflects both the range of motion between the domains and the probability of different orientations. Users can then visualize and select specific regions of the landscape to generate 3D reconstructions with defined orientation between two domains.

## Table of Contents

*   [INSTALLATION](#installation)
*   [PREREQUISITES](#prerequisites)
*   [USAGE](#usage)
    * [orientation_analysis](#orientation_analysis)
    * [landscape_projection](#landscape_projection)
    * [point_select](#point_select)
    * [particle_backtrack](#particle_backtrack)
*   [EXAMPLE\_WORKFLOW](#example_workflow)
*   [3D_Visualization_in_ChimeraX](#3d_visualization_in_chimerax)
*   [FAQ](#faq)

## INSTALLATION:

1\.  Clone the repository and cd inside

    git clone https://github.com/yifancheng-ucsf/cryorole.git
    cd cryorole


2\.  Create a conda environment with the required dependencies



    conda env create -f environment.yml

3\.  Activate the environment(always do this before using OrientationLandscape)



    conda activate cryorole

4\.  Verify the installation



    orientation_analysis --help

## Prerequisites:

#### About the input: Two Star Files

The two `.star` files must:
1, Contain the same number of particles.
2, Have the same particle ordering (i.e., particle 1 in both .star files refers to the same image).
3, Include columns such as `_rlnAngleRot`, `_rlnAngleTilt`, `_rlnAnglePsi`, `_rlnOriginXAngst`, `_rlnOriginYAngst`, etc.

*Note:* If your reconstruction is done in CryoSPARC, you can use pyem to convert the `.cs` files to `.star` files.

#### Required Python Libraries

*   numpy
*   pandas
*   matplotlib
*   scipy
*   mpi4py (for parallel execution)

#### MPI Environment (optional but recommended)

For parallel execution `orientation_analysis`, you will need an MPI implementation like OpenMPI or MPICH.

## Usage:

Below is a brief usage summary for main scripts. Use the `--help` flag for detailed arguments, default values, and optional flags.

### orientation\_analysis

1, Computes the relative orientation of each particle from two independent `.star` files.
2, Generates a CSV file containing per-particle relative orientation expressed in euler angle space (Z-Y-X）and rotation vectors, translation distance, RND scores, etc.
3, The Related\_Neighbor\_Density(RND) is used to color and filter data points. It is defined as ratio representing local density around each data point and is designed to remain comparable even across different datasets with different particle number counts.
4, Optionally compute and apply Inverse Mean Rotation (IMR) to align the orientation distribution to remove global orientation biases.

*   Basic Usage (Serial):
    ```
    orientation_analysis --s1 <star_file_1> --s2 <star_file_2> --o <output_prefix>
    ```

*   Parallel Usage(Recommended):
    ```
    mpirun -n <X> orientation_analysis --s1 <star_file_1> --s2 <star_file_2> --o <output_prefix>
    ```
    
*   Key Arguments:

> `--s1` / `--s2`: The two input .star files to compute.
> `--o`: Output prefix for the CSV files and other output.
> `--k`: Number of nearest neighbors for RND calculation (default = 30).
> `--autoalign`: If set, computes and applies IMR alignment, producing “\_before\_alignment.csv” and “\_after\_alignment.csv”.
> `--apply_rotation`: Apply a known rotation (ROT, TILT, PSI) to all particles before processing.
> `--outlier_method`: Method to identify outliers ('IQR' or 'Z-score').

### landscape\_projection

1, Visualize the CSV output from **orientation_analysis** by generating 2D projections of the orientation data (Euler angles, rotation vectors, etc.).
2, Optionally plot histograms of Alpha/Beta/Gamma, rotation angles, distances, RND scores, etc.
3, Optionally generate Azimuth-Elevation maps if requested.

*   Usage:
    ```
    landscape_projection --i <csv_file> --o <output_prefix>
    ```
*   Key Arguments:

> `--i`: Input CSV file from **orientation_analysis**.
> `--o`: Output prefix for generated plots (`.pdf`,`.png`,`.svg`).
> `--bins`: Number of bins in histograms (default = 50).
> `--t`: RND threshold for filtering (fixed value).
> `--p`: Percentage-based threshold (e.g., “keep top X% of RND”).
>  * If neigher `--t` nor `--`p is set, the default threshold is set to  (RND\_max/3).
>
> `--angle_limits`: Limits for Euler angles: Alpha\_min Alpha\_max Beta\_min Beta\_max Gamma\_min Gamma\_max, default=`[-180.0, 180.0, -180.0, 180.0, -180.0, 180.0]`.
> `--rotvec_limits`: Limits for Rotation vectors: RotVec\_X\_min RotVec\_X\_max RotVec\_Y\_min RotVec\_Y\_max RotVec\_Z\_min RotVec\_Z\_max, default=`[-3.14, 3.14, -3.14, 3.14, -3.14, 3.14]`.
> `--vmin, --vmax`: Color scale limits for scatter plots. Default would be the RND\_threshold and RND\_max.
> `--generate_histogram`: If set, plot the histogram of the columns.
> `--generate_azimuth_elevation`: If set, computes Azimuth/Elevation from rotation vectors and plots them.

### point\_select

1, Filters a CSV file to select particles whose Euler angles fall within a certain spherical region .
2, Optionally applies additional filtering based on other column-value pairs.

*   Usage:
    ```
    point_select --i <csv_file> --c <Alpha> <Beta> <Gamma> --r <radius> [--o <output_csv>]
    ```
    
*   Key Arguments:

> `--i`: Input CSV file (must have Alpha, Beta, Gamma columns).
> `--c`: Center coordinates of the “sphere” in (Alpha, Beta, Gamma) space.
> `--r`: Radius of the sphere.
> `--l`: Optional pairs of column labels and values for additional filtering.
> `--o`: Output CSV file with selected points. If not provided,  the number of matching particles is displayed in the terminal but without saved.

### particle\_backtrack

1, Maps a subset of particles (specified by their `ID` column in the CSV file) back to the original `.star` file.
2, Output a new `.star` file containing only the rows corresponding to the selected particles.

*   Usage:
    ```
    particle_backtrack --i <extracted_csv>  --s <original_star_file> --o <output_star_file>
    ```
    
*   Key Arguments:

> `--i`: The CSV file containing selected particles. Must have an `ID` column.
> `--s`: The original `.star` file with the full dataset.
> `--o`: Name of the `.star` file to write the subset to.

## EXAMPLE\_WORKFLOW:

In the folder `Sample`, two example `.star` files are provided:, `modifying_domain.star` and `condensing_domain.star`, These files correspond to separate domain refinements (as used in Nature Article <https://www.nature.com/articles/s41586-025-08782-w>). 

#### 1. Compute relative orientations relationship:
    
    mpirun -n 10 orientation_analysis --s1 modifying_domain.star --s2 condensing_domain.star --o modifying_against_condensing

This produces `modifying_against_condensing.csv` and `modifying_against_condensing.csv`(the points which have extreamly large RND have been declude, and in this dataset, it was 0).

If `--autoalign` is used, additional output `modifying_against_condensing_after_alignment.csv`.

#### 2. Visulize the orientation landscape projectons:

    landcape_projection --i modifying_against_condensing.csv --o modifying_against_condensing_visualization --p 0.4 

This produces plots with root name `modifying_against_condensing_visualization` in  SVG, PDF, and PNG formats showing the projections in euler space or rotation vector space with the top 40% RND filtered points. If `--generate_histogram` is used, the histogram of `Alpha`, `Beta`,`Gamma`, `RND` and `rotation angle` and `distance` would also generated.

You can further zoom the projections into specific ranges with the command `--angle_limits`:

    landcape_projection --i modifying_against_condensing.csv --o modifying_against_condensing_visualization_zoom --p 0.4 --angle_limits -90 90 -90 90 -90 -90 --rotvec_limits -1.57 1.57 -1.57 1.57 -1.57 1.57

#### 3. Select the points in a specified coordinate with a spherical region:

    point_select --i domain1_against_domain2.csv --c 13 0 14 --r 6  --o modifying_against_condensing_center_13_0_14_radius_6.csv

Here `(13,0,14)` in `Alpha, Beta, Gamma` space is the center and 6 is the radius, the selected points get saved to `modifying_against_condensing_center_13_0_14_radius_6.csv`.You can adjust the coordinate and radius as you need.

*Tips:*
`--o` is an optinal parameter here, If you just want to test how many particles can be extracted at a certain position and radius, `--o` is not required and the terminal will directly output information like: *Extracted 2580 particles inside the sphere defined by center\_13.0\_0.0\_14.0\_radius\_6.0*.

#### 4. Backtrack the selected points to original star file:

    python ParticleBacktrack.py --i modifying_against_condensing_center_13_0_14_radius_6.csv --s condensing_domain.star  --o modifying_against_condensing_center_13_0_14_radius_6_back_to_condensing.star

This produces a new `modifying_against_condensing_center_13_0_14_radius_6_back_to_condensing.star` file containing only the subset of particles within the specified spherical region. This subset can then be used for 3D reconstruction, local refinement, or further analysis. 

## 3D\_Visualization\_in\_ChimeraX:

To visualize the landscape in 3D, a ChimeraX Python script—developed with the help of Tom Goddard, **Points_in_ChimeraX.py**—is provided. This script registers custom commands for loading, updating, and displaying point cloud data from the CSV output of **orientation_analyisi**. The script supports two different modes for point cloud data:

*   **Euler Mode**: Uses Euler angles (`Alpha`, `Beta`, `Gamma`) with full bounds of `(-180, 180, -90, 90, -180, 180)`.
*   **Rotation Vector Mode**: Uses rotation vector components (`RotVec_X`, `RotVec_Y`, `RotVec_Z`) with full bounds of `(-3.14, 3.14, -3.14, 3.14, -3.14, 3.14)`.

### Features

*   **Load**: Load a new point cloud file with the command `points load`.
*   **Update**: Update filtering and coloring parameters on an already loaded point cloud with `points update`.
*   **Map**: Create a Gaussian-based mrc map from the point cloud with `points map`.

### Requirements

*   **ChimeraX**

### Usage

1.  **Load the script in ChimeraX:**
    Open ChimeraX and load the script by running at Command Line:

        open /path/to/Points_in_ChimeraX.py

    After loading the script in ChimeraX, the following commands will be available:

2.  **Loading a Point Cloud File**

*   **Euler Mode**
    Loads a file to display 3D point cloud using Euler angles (expects headers: `Alpha`, `Beta`, `Gamma`, `RND`):

        points load euler /path/to/yourfile.csv

    The full bounds used for mapping will be set to `(-180, 180, -90, 90, -180, 180)`.

*   **Rotation Vector Mode**
    Loads a file to display 3D point cloud using rotation vector components (expects headers: `RotVec_X`, `RotVec_Y`, `RotVec_Z`, `RND`):

        points load rotvec /path/to/yourfile.csv

    The full bounds used for mapping will be set to `(-3.14, 3.14, -3.14, 3.14, -3.14, 3.14)`.

    *Note:* You must choose one mode when loading a file. In either case, the point cloud is displayed in 3D and colored by the RND value.

3.  **Updating an Existing Point Cloud**\
    After loading a point cloud from your CSV file, you can adjust the view with the command `view` in ChimeraX, then update filtering/coloring with `points update`. For example, to hide points with `RND` values less than 0.6:


        points update #1 min_value 0.6

    Here, `#1` is the model number created when you opened the data. If you load more than one point cloud, be sure to specify the correct model. You can also update the coloring range:


        points update #1 minValue 0.6 range 0.6,2.0

    By default, the script applies a rainbow colormap over the full range of values (e.g., 0.6 to 2.0 in this case). To see all available options, type:

        usage points update

    *points update pointsModel \[palette a colormap] \[colorRange colorRange] \[minValue a number] \[maxValue a number]
    Update filtering and coloring of an existing point cloud
    colorRange: some numbers or full*

    You can specify a detailed colormap as follows:

        points update #1 minValue 0.6 palette 0.6,#0000ff10:0.8,#ff00ff10:1.0,#ffff0099:1.2,red

4.  **Creating a Map from the Point Cloud**
    Generate an mrc map from the loaded point cloud model by placing a Gaussian at each point position. The standard deviation (`sdev`) is 1 by default, and the grid spacing (`gridSpacing`) is also 1 (recommand use full bounds based on the model’s mode):

        points map #1 full_bounds true

    To see addinational options, use "usage" command:

        usage points map

    *points map pointsModel \[sdev a number] \[gridSpacing a number] \[bounds bounds] \[fullBounds true or false]
    Create a map by placing a Gaussian at each point of a point cloud
    bounds: some numbers*
    
    The map can be saved to an MRC file use:
    ```
    save test.mrc model #2
    ```
    
## FAQ

1.  **How do I ensure that the particle number and order in the two star files are identical?**
2.  **Does the sample size affect the landscape pattern?**
3.  **How can I validate that the generated landscape is correct?**
4.  **If my dataset contains multiple sub-states (e.g., from RELION 3D classification or different experimental conditions), how can I quickly obtain the landscape for each sub-state?**
5.  **Can landscapes from multiple sub-states be directly compared?**
6.  **When comparing landscapes from two similar samples refined with different initial models, what is the correct approach for comparison?**
7.  **Why does the landscape generated from my data only have a very small point in the center?**

