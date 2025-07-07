V-SYNTHES 2.0 Pipeline (Prepared by Dr. Antonina Nazarova)
This directory contains scripts necessary for executing the V-SYNTHES 2.0 workflow (located in the scripts/ folder), along with example files derived from the Enamine REAL Space library for demonstration purposes. The original V-SYNTHES algorithm was described in Sadybekov et al., 2021 (Nature), while V-SYNTHES 2.0 represents an enhanced and updated implementation detailed in Nazarova et al., 2025 (submitted). The receptor model used here is the Rho receptor, thoroughly characterized within the manuscript. The included examples represent two- and three-component reaction subsets of Enamine REAL Space. The current REAL Space library contains approximately 36 billion fully enumerated molecules across more than 164 well-established parallel synthesis reactions (for comprehensive information for your specific REAL library implemented in screening, please contact Enamine).
Execution of V-SYNTHES scripts and subsequent docking analyses require an ICM-Pro+VLS license, available from MolSoft (https://www.molsoft.com/products.html). Scripts included in this directory were validated using ICM-Pro+VLS version 3.9-4a running on Ubuntu Linux version 18.04.
Example MEL fragment files for two-component and three-component projects are provided within subdirectories under REAL_36B_files/, specifically REAL_022025_MEL_2comp.csv, REAL_022025_MEL_2comp.molt, REAL_022025_MEL_3comp.csv, and REAL_022025_MEL_3comp.molt. Complete versions of these datasets can be requested directly from Enamine LLC. Reaction libraries and synthon datasets essential for both 2- and 3-component reactions are included as REAL_36B_files/REAL_36B_022025_reactions.icb and REAL_36B_022025_syntons.molt.
To apply the V-SYNTHES 2.0 workflow using the provided examples, follow the detailed instructions outlined in Figure 2 of the manuscript, specifically described in Procedure 1 (V-SYNTHES 2.0 Automated Run). For statistical evaluations of docking poses and conformational analyses of enumerated hits, refer to Procedure 2 (RMSD Assessment between MEL fragments and Corresponding Enumerated Molecules).
Additionally, V-SYNTHES 2.0 introduces an automated algorithm, CapSelect, designed to optimize the previously manual step 2b from V-SYNTHES 1.0:
CapSelect Algorithm for V-SYNTHES 2.0
Overview
CapSelect is an innovative geometry-based algorithm that automates the identification and selection of optimal MEL fragments for enumeration within the V-SYNTHES 2.0 framework. It generates non-overlapping spherical regions around potential fragment-growth points, ensuring minimal steric interference with binding pocket residues and existing ligands, thus facilitating accurate enumeration predictions.
Key Outputs
•	Spheres: Number of viable non-overlapping spherical regions.
•	MaxMin: Maximum minimum distance from sphere surfaces to pocket residues.
•	Distance: Distances between the cap-associated sphere and centroids of other spheres.
•	CapScore: A metric evaluating available spatial volume for fragment growth.
Versions and File Generation
Two versions of CapSelect are available:
•	C++ version (CapSelect.cpp)
•	Python version (CapSelect.py)
Input files (.sdf format) are automatically generated via the script icm_generate_CapSelect_files.icm. The output file (frags_for_enum.sdf) is generated using icm_CapSelect_to_frags_for_enum.icm.
Execution Workflow
The workflow is managed by the bash script CapSelectMP.sh, which performs:
1.	Input file generation and parsing.
2.	Fragment evaluation via CapSelect (.py or .cpp).
3.	Output formatting, leveraging multi-core CPU parallelization for performance optimization.
Running CapSelect in Python
1.	Confirm Python installation (command: python --version).
2.	Place the following files into your project root directory: CapSelect.py, CapSelectMP.sh, icm_generate_CapSelect_files.icm, icm_CapSelect_to_frags_for_enum.icm.
3.	Execute the command: ./CapSelectMP.sh (specify number of cores using -c N, default is 16).
4.	The resulting file frags_for_enum.sdf will appear in root()/run/processing_file/.
Running CapSelect in C++
1.	A pre-compiled executable (CapSelect) optimized with -O2 is provided. For custom compilation, use: g++ -O2 -o CapSelect CapSelect.cpp.
2.	Place CapSelect, CapSelectMP.sh, icm_generate_CapSelect_files.icm, and icm_CapSelect_to_frags_for_enum.icm into your V-SYNTHES project directory.
3.	Run the script: ./CapSelectMP.sh (use -c N to specify the number of cores if desired).
4.	The generated output file matches that of the Python implementation.
