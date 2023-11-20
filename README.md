# V-SYNTHES-2.-_CapSelect
CapSelect Algorithm for V_SYNTHES 2.*
  Overview
CapSelect is a geometry-modular algorithm essential for automating the selection of productive MEL fragments in V_SYNTHES 2.*. It constructs a series of non-overlapping spheres, ensuring no overlap with the pocket or ligand, thereby predicting the potential location of a fully enumerated ligand upon docking.

  Key Outputs
  +Spheres: Reports the number of possible non-overlapping spheres.
  +MaxMin: Indicates the largest minimum distance from sphere meshes to the pocket surface.
  +Distance: Details the distances from the cap-associated sphere to other sphere centroids.
  +CapScore: A relative score assessing space availability for MEL fragment growth.

  Versions and File Generation
Two versions available: C++ (CapSelect.cpp) and Python (CapSelect.py).
Input files (.sdf) are auto-generated using icm_generate_CapSelect_files.icm.
Output file frags_for_enum.sdf is generated through icm_CapSelect_to_frags_for_enum.icm.

  Execution Flow
The CapSelectMP.sh bash script manages the entire workflow, including:
Generation and parsing of input files.
Processing with CapSelect.py or .cpp.
Output formatting and multi-core CPU parallel processing for efficiency.

  Running CapSelect in Python
  1.  Ensure Python is installed (python --version).
  2.  Place CapSelect.py, CapSelectMP.sh, icm_generate_CapSelect_files.icm, and icm_CapSelect_to_frags_for_enum.icm in your V_SYNTHES project root.
  3.  Run ./CapSelectMP.sh in the terminal (use -c N to specify the number of cores, default is 16).
  4.  The output file frags_for_enum.sdf will be generated in root()/run/processing_file/frags_for_enu.sdf.

  Running CapSelect in C++
  1.A pre-compiled CapSelect.exe is provided (compiled with -O2 optimization). For custom compilation: g++ -O2 -o CapSelect CapSelect.cpp.
  2.Place CapSelect, CapSelectMP.sh, icm_generate_CapSelect_files.icm, and icm_CapSelect_to_frags_for_enum.icm in your V_SYNTHES project root.
  3.Execute ./CapSelectMP.sh in the terminal (specify cores using -c N if required).
  4.Observe the same output as in the Python version.
