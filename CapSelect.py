"""
Disclaimer and Copyright Notice
"CapSelect Algorithm: Geometry-Modular Approach for MEL Fragment Selection in V_SYNTHES 2.*"

Authorship:
This code, entitled "CapSelect", is developed and maintained by Dr. Antonina L. Nazarova. 
It serves as a core algorithm for the automation and streamlining of productive MEL fragment 
selection in the V_SYNTHES 2.* version.

Intellectual Property:
The CapSelect algorithm is designed to automate and streamline the selection of productive MEL
fragments for V_SYNTHES 2.*. It employs a geometry-modular approach, centered around constructing
a series of consecutive, non-overlapping spheres that do not intersect with either the ligand or 
the pocket. This methodology effectively estimates the available space within the pocket post-MEL
binding through the arrangement of equidistant spheres. Ultimately, CapSelect provides a predictive
model of the overall geometry, forecasting the potential orientation and placement of fully 
enumerated ligands during the docking process.

Usage Rights:
The code is made available for academic, research, and non-commercial use. Any commercial use, 
duplication, modification, distribution, or replication of this software without explicit 
permission from the author is strictly prohibited.

Execution:
The CapSelect algorithm features an integrated suite of scripts, automated for ease of use. 
The CapSelectMP.sh bash script manages the workflow, starting from input file creation with 
icm_generate_CapSelect_files.icm, processing via CapSelect.py, to output generation through 
icm_CapSelect_to_frags_for_enum.icm. Enhanced by multi-core CPU parallel processing, this 
script ensures efficient and streamlined execution.

Citation Requirement:
If you use this code or its components in your research or work, proper acknowledgment is required. 
The following citation format is recommended:
Nazarova, A. L. "CapSelect Algorithm: Geometry-Modular Approach for MEL Fragment Selection in V_SYNTHES 2.*", 
[Year of Publication/Access].

"""

import math
import os
import numpy as np

def filecheck(filepath):
    """Check if a file exists."""
    try:
        with open(filepath, 'r'):
            return True
    except FileNotFoundError:
        return False


def read_sdf_frag(filename):
    print(f"Reading {filename}")

    i_res_l, x_l, y_l, z_l, elm_l, lab_l, RecConf, DockScore = [], [], [], [], [], [], [], []

    with open(filename, 'r') as myfile:
        groupId = 0
        while True:
            line = myfile.readline()
            if not line:
                break

            if line.strip() == "":
                continue

            # Skip two lines
            myfile.readline()
            myfile.readline()

            line = myfile.readline().strip()
            i_res = int(line[:3])
            i_res_l.append(i_res)

            x_group, y_group, z_group, elm_group, lab_group = [], [], [], [], []
            for _ in range(i_res):
                line = myfile.readline()
                x_group.append(float(line[:10].strip()))
                y_group.append(float(line[10:20].strip()))
                z_group.append(float(line[20:30].strip()))
                elm_group.append(line[31])
                lab_group.append(line[35])

            x_l.append(x_group)
            y_l.append(y_group)
            z_l.append(z_group)
            elm_l.append(elm_group)
            lab_l.append(lab_group)

            # Reading additional properties like RecConf and DockScore
            while True:
                line = myfile.readline()
                if not line or line.startswith("$$$$"):
                    break

                if line.startswith("> <RecConf>"):
                    RecConf_value = int(myfile.readline().strip())
                    RecConf.append(RecConf_value)
                elif line.startswith("> <Score>"):
                    DockScore_value = float(myfile.readline().strip())
                    DockScore.append(DockScore_value)

            groupId += 1
    #DockScore = np.array(DockScore)
    return x_l, y_l, z_l, elm_l, lab_l, i_res_l, RecConf, DockScore
 

def protein_atom_count(filename):
    print(f"Loading {filename}")
    atom_count = 0

    with open(filename, 'r') as myfile:
        # Skip initial lines
        for _ in range(4):
            myfile.readline()

        # Count atom lines
        for line in myfile:
            if len(line.strip()) > 30:  # Assuming lines longer than 30 characters represent atoms
                atom_count += 1
            else:
                break  # Stop counting if a line is shorter, assuming it's the end of the atom list

    print("Number of Protein Atoms:", atom_count)
    return atom_count
    

def read_sdf_protein(filename):
    print(f"Loading {filename}")

    x_p, y_p, z_p, elm_p, lab_p, i_res_p = [], [], [], [], [], []

    with open(filename, 'r') as myfile:
        # Skip header lines
        for _ in range(3):
            myfile.readline()

        # Read the atom count line
        line = myfile.readline()
        i_res_p = int(line[:4].strip())  # Extract the number of atoms

        # Read atom data
        for _ in range(i_res_p):
            line = myfile.readline()
            if len(line.strip()) < 10:  # Check if line is too short
                continue

            try:
                x = float(line[:10].strip())
                y = float(line[10:20].strip())
                z = float(line[20:30].strip())
                elm = line[31:34].strip()
                lab = line[35:].strip()

                x_p.append(x)
                y_p.append(y)
                z_p.append(z)
                elm_p.append(elm)
                lab_p.append(lab)
            except ValueError as e:
                print(f"Error reading line: {line}")
                raise e

    return x_p, y_p, z_p, elm_p, lab_p, i_res_p
    
def calculate_CapScore(i, num_lab_l, num1, min1_l, ff):
    max_label_possible = num1.shape[1]
    num_spheres = min1_l.shape[2]

    penalty_s = np.zeros(max_label_possible)
    penalty_max1_l = np.zeros(max_label_possible)
    score = np.zeros(max_label_possible)
    
    min1_l_fin = np.zeros(num_spheres)
    ff_fin = np.zeros(num_spheres)

    if num_lab_l[i] in [1, 5]:
        for j1 in range(1):
            penalty_s[j1] = (5 - num1[i][j1]) ** 2 if num1[i][j1] <= 5 else 0
            penalty_s[j1] = (penalty_s[j1] * 10) / 25

            if num1[i][j1] > 9:
                min1_value = min1_l[i][j1][9]
                if 7 <= min1_value <= 20:
                    penalty_max1_l[j1] = (7 - min1_value) ** 2
                elif min1_value > 20:
                    penalty_max1_l[j1] = 169
                else:
                    penalty_max1_l[j1] = 0
            else:
                penalty_max1_l[j1] = 0

            penalty_max1_l[j1] = (penalty_max1_l[j1] * 10) / 169
            score[j1] = 10 - penalty_s[j1] - penalty_max1_l[j1]
            CapScore_i = score[j1]
            num1_fin = num1[i][j1]
            for k in range(num1[i][j1]):
                min1_l_fin[k] = min1_l[i][j1][k]
                ff_fin[k] = ff[i][j1][k]

    if num_lab_l[i] in [10, 6, 2]:
        for j2 in range(2):
            penalty_s[j2] = (5 - num1[i][j2]) ** 2 if num1[i][j2] <= 5 else 0
            penalty_s[j2] = (penalty_s[j2] * 10) / 25

            if num1[i][j2] > 9:
                min1_value = min1_l[i][j2][9]
                if 7 <= min1_value <= 20:
                    penalty_max1_l[j2] = (7 - min1_value) ** 2
                elif min1_value > 20:
                    penalty_max1_l[j2] = 169
                else:
                    penalty_max1_l[j2] = 0
            else:
                penalty_max1_l[j2] = 0

            penalty_max1_l[j2] = (penalty_max1_l[j2] * 10) / 169
            score[j2] = 10 - penalty_s[j2] - penalty_max1_l[j2]
            
            for k1 in range(num1[i][j2]):
                min1_l_fin[k1] = min1_l[i][j2][k1]
                ff_fin[k1] = ff[i][j2][k1]

        if score[0]<score[1]:
            CapScore_i = score[1]
            num1_fin = num1[i][1]
            for k2 in range(num1[i][1]):
                min1_l_fin[k2] = min1_l[i][1][k2]
                ff_fin[k2] = ff[i][1][k2]
        else:
            CapScore_i = score[0]
            num1_fin = num1[i][0]
            for k2 in range(num1[i][0]):
                min1_l_fin[k2] = min1_l[i][0][k2]
                ff_fin[k2] = ff[i][0][k2]
    #print(f"Type of CapScore in DEF: {type(CapScore_i)}")

    return CapScore_i, min1_l_fin, ff_fin, num1_fin

def calculate_MergedScore(CapScore_i, DockScore_i, num1_fin_i):
    if CapScore_i == -100:
        MergedScore_i = -1000
    elif CapScore_i == 0 and num1_fin_i == 0:
        MergedScore_i = -1000
    elif CapScore_i == 0 and num1_fin_i == 10:
        MergedScore_i = 5.0 * np.log2(abs(DockScore_i))
    else:
        MergedScore_i = 5.0 * np.log2(abs(DockScore_i)) + 0.5 * np.log2(abs(CapScore_i))
    return MergedScore_i
    
def insert( MergedScore_i, CapScore_i, num1_fin, min1_l_fin, ff_fin):
    printstr = "> <MergedScore>\n" + str(MergedScore_i) + "\n\n" + "> <CapScore>\n" + str(CapScore_i) + "\n\n" + "> <Spheres>\n" + str(num1_fin) + "\n\n" + "> <Max(min)>\n"

    for j in range(num1_fin):
        printstr += str(min1_l_fin[j]) + (", " if j < num1_fin - 1 else "\n")

    printstr += "\n> <Distance>\n"
    for k in range(num1_fin):
        printstr += str(ff_fin[k]) + (", " if k < num1_fin - 1 else "\n")

    return printstr
    
def calculate_sphere_0(i, num_current, x_in, y_in, z_in, xll, yll, zll, num_lab_1, num_lab_2, a_s_f, a_c_f, a_s_t, a_c_t, lab_l, elm_l, x_l, y_l, z_l, i_res_l, i_res_p, x_p, y_p, z_p):
   
    r = 3.5 if ((num_lab_1 == 5 and num_current == 1) or (num_lab_2 == 5 and num_current == 2)) else 3.0
        
    min11 = float('inf')
    x_res, y_res, z_res = 0, 0, 0
    ipp = 0

    for i4 in range(72):  # Azimuth
        for i1 in range(36):  # Longitude
            z = r * a_s_f[i1] + zll[i][num_current]
            x = r * a_c_f[i1] * a_s_t[i4] + xll[i][num_current]
            y = r * a_c_f[i1] * a_c_t[i4] + yll[i][num_current]

            l_check = 1.1 if ((num_current == 1 and num_lab_1 == 5) or (num_current == 2 and num_lab_2 == 5)) else 1.3

            ip = 0
            for i3 in range(i_res_l[i]):
                if lab_l[i][i3] not in ['1', '3'] and elm_l[i][i3] != 'H':
                    r2 = math.sqrt((x - x_l[i][i3])**2 + (y - y_l[i][i3])**2 + (z - z_l[i][i3])**2)
                    r1 = math.sqrt((x_in - x_l[i][i3])**2 + (y_in - y_l[i][i3])**2 + (z_in - z_l[i][i3])**2)
                    if r2 < r or r1 < l_check:
                        ip = 1
                        break

            p_check = 2.0 if ((num_current == 1 and num_lab_1 == 5) or (num_current == 2 and num_lab_2 == 5)) else 3.0

            ip1 = 0
            for i3 in range(i_res_p):
                r2 = math.sqrt((x - x_p[i3])**2 + (y - y_p[i3])**2 + (z - z_p[i3])**2)
                if r2 < p_check:
                    ip1 = 1
                    break

            if not (ip == 1 or ip1 == 1):
                ipp = 1
                min_value = 1000.0
                for i2 in range(i_res_p):
                    r1 = math.sqrt((x - x_p[i2])**2 + (y - y_p[i2])**2 + (z - z_p[i2])**2)
                    if min_value > r1:
                        min_value = r1
                        x2 = x
                        y2 = y
                        z2 = z
                        

                if min11 > min_value:
                    min11 = min_value
                    x_res = x2
                    y_res = y2
                    z_res = z2
                    

    return ipp, min11, x_res, y_res, z_res

def calculate_sphere_non_0(i, num_current, z_res, x_res, y_res, a_s_f, a_c_f, a_s_t, a_c_t, x_l, y_l, z_l, i_res_l, lab_l, elm_l, x_p, y_p, z_p, i_res_p, num_lab_1, num_lab_2, num_lab_3, num_lab_4, x_resl, y_resl, z_resl, num):
    rrr = 5.0  # Adjusted value
    rrrr = 2.0
    xl, yl, zl = x_res, y_res, z_res

    min11 = 0.0
    ipp = 0

    for i4 in range(72):  # Azimuth
        for i1 in range(36):  # Longitude
            z = rrrr * a_s_f[i1] + zl
            x = rrrr * a_c_f[i1] * a_s_t[i4] + xl
            y = rrrr * a_c_f[i1] * a_c_t[i4] + yl

            ip = 0
            l_check = 1.1 if ((num_current == 1 and num_lab_1 == 5) or (num_current == 2 and num_lab_2 == 5)) else 1.3

            for i3 in range(i_res_l[i]):
                if lab_l[i][i3] != '1' and elm_l[i][i3] != 'H' and lab_l[i][i3] != '3':
                    r2 = math.sqrt((x - x_l[i][i3])**2 + (y - y_l[i][i3])**2 + (z - z_l[i][i3])**2)
                    r1 = math.sqrt((xl - x_l[i][i3])**2 + (yl - y_l[i][i3])**2 + (zl - z_l[i][i3])**2)
                    if r2 < rrr or r1 < l_check:
                        ip = 1
                        break

            ip1 = 0
            for i3 in range(i_res_p):
                r2 = math.sqrt((x - x_p[i3])**2 + (y - y_p[i3])**2 + (z - z_p[i3])**2)
                r1 = math.sqrt((xl - x_p[i3])**2 + (yl - y_p[i3])**2 + (zl - z_p[i3])**2)
                if r2 < 2.0 or r1 < 3.0:
                    ip1 = 1
                    break

            ip2 = 0
            if num > 1:
                for i7 in range(num):
                    r1 = math.sqrt((x - x_resl[i][num_current][i7])**2 + (y - y_resl[i][num_current][i7])**2 + (z - z_resl[i][num_current][i7])**2)
                    if r1 < 2.0:
                        ip2 = 1
                        break

            # Implementing additional conditions for num_lab_3 and num_lab_4
            ip2_1 = 0
            if num == 1:
                r2 = math.sqrt((x - xl)**2 + (y - yl)**2 + (z - zl)**2)
                
                # Conditions for num_lab_1, num_lab_2, num_lab_3, and num_lab_4
                if num_current == 1 and num_lab_1 == 5 or num_current == 2 and num_lab_2 == 5:
                    if r2 < 3.5:
                        ip2_1 = 1
                elif num_current == 1 and num_lab_3 == 1 or num_current == 2 and num_lab_4 == 1:
                    if r2 < 3.0:
                        ip2_1 = 1

            ip3 = 0
            if num > 1:
                for j2 in range(num - 1, num):
                    r2 = math.sqrt((x - x_resl[i][num_current][j2])**2 + (y - y_resl[i][num_current][j2])**2 + (z - z_resl[i][num_current][j2])**2)
                    if r2 < 3.46:
                        ip3 = 1
                        break

            # Checking all conditions together
            if not any([ip, ip1, ip2, ip3, ip2_1]):
                ipp = 1

            # Calculating minimum distance
                min_value = 1000.0
                for i2 in range(i_res_p):
                    r1 = math.sqrt((x - x_p[i2])**2 + (y - y_p[i2])**2 + (z - z_p[i2])**2)
                    if min_value > r1:
                        min_value = r1
                        x2, y2, z2 = x, y, z

                if min11 < min_value:
                    min11 = min_value
                    x_res, y_res, z_res = x2, y2, z2

    return min11, x_res, y_res, z_res, ipp



def main():
    # Constants and variables initialization
    pi = 3.14159264
    a_c_f = [math.cos(-pi / 2.0 + i * 5.0 * pi / 180.0) for i in range(36)]
    a_s_f = [math.sin(-pi / 2.0 + i * 5.0 * pi / 180.0) for i in range(36)]
    a_c_t = [math.cos(i * 2.0 * 5.0 * pi / 360.0) for i in range(72)]
    a_s_t = [math.sin(i * 2.0 * 5.0 * pi / 360.0) for i in range(72)]
    
    # Read fragment data
    x_l, y_l, z_l, elm_l, lab_l, i_res_l, RecConf, DockScore = read_sdf_frag("fragments.sdf")
    i_count_l = len(x_l)  # Number of fragments
    print("\nNumber of FRAGMENTS:", i_count_l)
    print("\n======================\n")

    xll, yll, zll, max_label= [[] for _ in range(i_count_l)], [[] for _ in range(i_count_l)], [[] for _ in range(i_count_l)], [0] * i_count_l
    CapScore = np.full((i_count_l), -100.0)
    num_lab_l = np.zeros(i_count_l)
    #max_label = np.ones(i_count_l, dtype=int)

    for i in range(i_count_l):
        if RecConf[i] == 1 and filecheck("protein1.sdf"):
            i_protein_atom_count = protein_atom_count("protein1.sdf")
            x_p, y_p, z_p, elm_p, lab_p, i_res_p = read_sdf_protein("protein1.sdf")
        elif RecConf[i] == 2 and filecheck("protein2.sdf"):
            i_protein_atom_count = protein_atom_count("protein2.sdf")
            x_p, y_p, z_p, elm_p, lab_p, i_res_p = read_sdf_protein("protein2.sdf")
        else:
            print("\nERROR: <RecConf> information invalid or file not found.")
            continue

        xl_1 = xl_2 = xl_3 = xl_4 = 0.0
        yl_1 = yl_2 = yl_3 = yl_4 = 0.0
        zl_1 = zl_2 = zl_3 = zl_4 = 0.0
        num_lab_1 = num_lab_2 = num_lab_3 = num_lab_4 = 0
        

        for i1 in range(i_res_l[i]):
            if lab_l[i][i1] == '3':  
                if num_lab_1 < 5:
                    xl_1 += x_l[i][i1]
                    yl_1 += y_l[i][i1]
                    zl_1 += z_l[i][i1]
                    num_lab_1 += 1
                else:
                    xl_2 += x_l[i][i1]
                    yl_2 += y_l[i][i1]
                    zl_2 += z_l[i][i1]
                    num_lab_2 += 1
            elif lab_l[i][i1] == '1':  
                if num_lab_3 < 1:
                    xl_3 = x_l[i][i1]
                    yl_3 = y_l[i][i1]
                    zl_3 = z_l[i][i1]
                    num_lab_3 = 1
                else:
                    xl_4 = x_l[i][i1]
                    yl_4 = y_l[i][i1]
                    zl_4 = z_l[i][i1]
                    num_lab_4 = 1

        
        if num_lab_1 == 5 and num_lab_2 == 5:
            max_label[i] = 3
            xll[i].append(xl_1 / num_lab_1)
            yll[i].append(yl_1 / num_lab_1)
            zll[i].append(zl_1 / num_lab_1)

            xll[i].append(xl_2 / num_lab_2)
            yll[i].append(yl_2 / num_lab_2)
            zll[i].append(zl_2 / num_lab_2)

            print(f"\n2_AROMATIC_CAPS: x1 = {xll[i][0]:.4f}, y1 = {yll[i][0]:.4f}, z1 = {zll[i][0]:.4f}")
            print(f"\n2_AROMATIC_CAPS: x2 = {xll[i][1]:.4f}, y2 = {yll[i][1]:.4f}, z2 = {zll[i][1]:.4f}")

        elif num_lab_3 == 1 and num_lab_4 == 1:
            max_label[i] = 3
            xll[i].append(xl_3)
            yll[i].append(yl_3)
            zll[i].append(zl_3)

            xll[i].append(xl_4)
            yll[i].append(yl_4)
            zll[i].append(zl_4)

            print(f"\n2_NON_AROMATIC_CAPS: x1 = {xll[i][0]:.4f}, y1 = {yll[i][0]:.4f}, z1 = {zll[i][0]:.4f}")
            print(f"\n2_NON_AROMATIC_CAPS: x2 = {xll[i][1]:.4f}, y2 = {yll[i][1]:.4f}, z2 = {zll[i][1]:.4f}")

        elif num_lab_1 == 5 and num_lab_3 == 1:
            max_label[i] = 3
            xll[i].append(xl_1 / num_lab_1)
            yll[i].append(yl_1 / num_lab_1)
            zll[i].append(zl_1 / num_lab_1)

            xll[i].append(xl_3)
            yll[i].append(yl_3)
            zll[i].append(zl_3)

            print(f"\nAROM_NON_AROM_CAPS: x1 = {xll[i][0]:.4f}, y1 = {yll[i][0]:.4f}, z1 = {zll[i][0]:.4f}")
            print(f"\nAROM_NON_AROM_CAPS: x2 = {xll[i][1]:.4f}, y2 = {yll[i][1]:.4f}, z2 = {zll[i][1]:.4f}")

        elif num_lab_1 == 5 and (num_lab_2 == 0 and num_lab_3 == 0 and num_lab_4 == 0):
            max_label[i] = 2
            xll[i].append(xl_1 / num_lab_1)
            yll[i].append(yl_1 / num_lab_1)
            zll[i].append(zl_1 / num_lab_1)

            print(f"\nSINGLE_CAP: x1 = {xll[i][0]:.4f}, y1 = {yll[i][0]:.4f}, z1 = {zll[i][0]:.4f}")

        elif num_lab_3 == 1 and (num_lab_1 == 0 and num_lab_2 == 0 and num_lab_4 == 0):
            max_label[i] = 2
            xll[i].append(xl_3)
            yll[i].append(yl_3)
            zll[i].append(zl_3)

            print(f"\nSINGLE_CAP: x1 = {xll[i][0]:.4f}, y1 = {yll[i][0]:.4f}, z1 = {zll[i][0]:.4f}")
        # Calculate N of CAPs

        num_lab_l[i] = num_lab_1 + num_lab_2 + num_lab_3 + num_lab_4

        # Determine max_label and num_current_1, reffering to the N of CAPs
        if num_lab_l[i] in [1, 5]:
            max_label[i] = 2
            num_current_1 = 1
            print("\nMAX_LABEL: {}".format(max_label[i]))
        elif num_lab_l[i] in [10, 6, 2]:
            max_label[i] = 3
            num_current_1 = 1
            print("\nMAX_LABEL: {}".format(max_label[i]))
            
        max_label_possible = max_label[i]-1  
        num_spheres = 10  
        
        if num_lab_1 == 0 and num_lab_2 == 0 and num_lab_3 == 0 and num_lab_4 == 0:
            CapScore[i] = -100.0
            num1_fin[i] = 0
            max_label[i] = 1
            max_label_possible = 0
            break  # Skip to the next fragment

        #Sphere-based sequence
        # Initialize the 3D NumPy arrays
        x_resl = np.zeros((i_count_l, max_label_possible, num_spheres))
        y_resl = np.zeros((i_count_l, max_label_possible, num_spheres))
        z_resl = np.zeros((i_count_l, max_label_possible, num_spheres))
        x_in = np.zeros((i_count_l, max_label_possible))
        y_in = np.zeros((i_count_l, max_label_possible))
        z_in = np.zeros((i_count_l, max_label_possible))
        
        num1 = np.zeros((i_count_l, max_label_possible), dtype=int)
        min1_l = np.zeros((i_count_l, max_label_possible, num_spheres))
        ff = np.zeros((i_count_l, max_label_possible, num_spheres))
        for num_current in range(0, max_label[i]-1):
            print(f"\nSINGLE_CAP: x1 = {xll[i][num_current]:.4f}, y1 = {yll[i][num_current]:.4f}, z1 = {zll[i][num_current]:.4f}")


        for num_current in range(0, max_label[i]-1):
            exit_for_loop = False 
            num = 0
            x_in = xll[i][num_current]
            y_in = yll[i][num_current]
            z_in = zll[i][num_current]
            x_resl[i][num_current][0] = x_in
            y_resl[i][num_current][0] = y_in
            z_resl[i][num_current][0] = z_in
            ipp, min11, x_res, y_res, z_res = calculate_sphere_0(i, num_current, x_in, y_in, z_in, xll, yll, zll, num_lab_1, num_lab_2, a_s_f, a_c_f, a_s_t, a_c_t, lab_l, elm_l, x_l, y_l, z_l, i_res_l, i_res_p, x_p, y_p, z_p)
            if ipp == 0:
                if num_lab_l[i] in [1, 5] or (num_lab_l[i] in [10, 6, 2] and num_current == 2):
                    num1[i][num_current] = num
                    min1_l[i][num_current][num] = min11
                    ff[i][num_current][num] = 0.0
                    CapScore, min1_l_fin, ff_fin, num1_fin = calculate_CapScore(i, num_lab_l, num1, min1_l, ff)
                    print(f"Fragment {i}, Cap {num_current} processed. No-spheres possible")
                    break 
                if num_lab_l[i] in [10, 6, 2] and num_current == 1:
                    continue 
            num = 1
            x_resl[i][num_current][num] = x_res
            y_resl[i][num_current][num] = y_res
            z_resl[i][num_current][num] = z_res
            min1_l[i][num_current][num] = min11
        
            while num < 10:
                min11, x_res, y_res, z_res, ipp = calculate_sphere_non_0(i, num_current, z_res, x_res, y_res, a_s_f, a_c_f, a_s_t, a_c_t, x_l, y_l, z_l, i_res_l, lab_l, elm_l, x_p, y_p, z_p, i_res_p, num_lab_1, num_lab_2, num_lab_3, num_lab_4, x_resl, y_resl, z_resl, num)
                if ipp == 0:
                    if ((num_lab_l[i] in [1, 5]) and (num == 1)) or ((num_lab_l[i] in [10, 6, 2] and num_current == 2) and (num == 1)):
                        num1[i][num_current] = 1
                        min1_l[i][num_current][num] = min1_l[i][num_current][1]
                        ff[i][num_current][num] = 0.0
                        CapScore_i, min1_l_fin, ff_fin, num1_fin = calculate_CapScore(i, num_lab_l, num1, min1_l, ff)
                        print(f"Fragment_n=1 {i}, Cap {num_current} processed. {num+1} spheres possible")
                        exit_for_loop = True 
                        break 
                    if ((num_lab_l[i] in [1, 5]) and (num > 1)) or ((num_lab_l[i] in [10, 6, 2] and num_current == 2) and (num > 1)):
                        num1[i][num_current] = num
                        CapScore_i, min1_l_fin, ff_fin, num1_fin = calculate_CapScore(i, num_lab_l, num1, min1_l, ff)
                        print(f"Fragment_n>1 {i}, Cap {num_current} processed. {num+1} spheres possible")
                        exit_for_loop = True 
                        break
                    if num_lab_l[i] in [10, 6, 2] and num_current == 1:
                        continue
                    if num == 9:
                        print(f"Fragment_n=10 {i}, Cap {num_current} processed. {num+1} spheres possible")                    
                else:
                    x_resl[i][num_current][num] = x_res
                    y_resl[i][num_current][num] = y_res
                    z_resl[i][num_current][num] = z_res
                    min1_l[i][num_current][num] = min11
                    ff[i][num_current][num] = math.sqrt((x_res - x_in)**2 + (y_res - y_in)**2 + (z_res - z_in)**2)
                    num += 1
                    #CapScore_i, min1_l_fin, ff_fin, num1_fin = calculate_CapScore(i, num_lab_l, num1, min1_l, ff)
                #print(f"Fragment_n>1 {i}, Cap {num_current} processed. {num+1} spheres possible")               
            if exit_for_loop:
                break  # Exit for loop if flag is set

    
        #CapScore_i, min1_l_fin, ff_fin, num1_fin = calculate_CapScore(i, num_lab_l, num1, min1_l, ff)
        MergedScore_i = calculate_MergedScore(CapScore_i, DockScore[i], num1_fin)
        print(f"Fragment {i}: CapScore = {CapScore_i:.4f}, MergedScore = {MergedScore_i:.4f}")             
        print(f"Type of CapScore in MAIN: {type(CapScore_i)}")

    with open("fragments.sdf", "a") as sdf_file:
        for i in range(i_count_l):
            sdf_entry = insert(MergedScore_i, CapScore_i, num1_fin, min1_l_fin, ff_fin)
            sdf_file.write(sdf_entry + "$$$$\n")

if __name__ == "__main__":
    main()