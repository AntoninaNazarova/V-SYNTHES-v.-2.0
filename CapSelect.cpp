/*
    Disclaimer and Copyright Notice
	"CapSelect Algorithm: Geometry-Modular Approach for MEL Fragment Selection in V_SYNTHES 2.*"
    Authorship:
    This code, entitled "CapSelect.cpp", is developed and maintained by Dr. Antonina L. Nazarova. 
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
	The CapSelect algorithm operates through a seamlessly integrated series of scripts, designed 
	for simplicity and ease of use. The CapSelectMP.sh bash script orchestrates the entire process:
	it initiates with input file creation using icm_generate_CapSelect_files.icm, proceeds with 
	processing by the CapSelect executable (users can compile new .exe from CapSelect.cpp using 
	the -O2 flag for optimized performance), and concludes with output file generation via 
	icm_CapSelect_to_frags_for_enum.icm. The script's utilization of multi-core CPU parallel processing.

    Citation Requirement:
    If you use this code or its components in your research or work, proper acknowledgment is required. 
    The following citation format is recommended:
    Nazarova, A. L. et al. "V-SYNTHES 2.* - The Next Generation Tool for the Screening Giga-Scale Chemical Space in Computer-Aided Drug Design"
    [Year of Publication/Access].

*/
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <functional>
#include <thread>
#include <algorithm>
#include <iomanip>


static float CapScore[50000], MergedScore[50000];
float ff_fin[50000][15];
static int num1_fin[50000];
static int i_count_l;
float min1_l_fin[50000][15];

std::string insert(int i) {
	int j;
	int k;
	std::string printstr;
	printstr = "> <MergedScore>\n"
		+ std::to_string(MergedScore[i]) + "\n"
		+ "\n"
		+ "> <CapScore>\n"
		+ std::to_string(CapScore[i]) + "\n"
		+ "\n"
		+ "> <Spheres>\n"
		+ std::to_string(num1_fin[i]) + "\n"
		+ "\n"
		+ "> <Max(min)>\n";

	for (j = 0; j < num1_fin[i]; j++)
	{
		if (j < num1_fin[i] - 1)
		{

			printstr = printstr + std::to_string(min1_l_fin[i][j]) + ", ";
		}
		else
		{
			printstr = printstr + std::to_string(min1_l_fin[i][j]) + "\n";
		}
		//
	}
	printstr = printstr + "\n" + "> <Distance>\n";
	for (k = 0; k < num1_fin[i]; k++)
	{
		if (k < num1_fin[i] - 1)
		{

			printstr = printstr + std::to_string(ff_fin[i][k]) + ", ";
		}
		else
		{
			printstr = printstr + std::to_string(ff_fin[i][k]) + "\n";
		}
		//
	}
	return printstr;
};


using namespace std;
void read_sdf_frag(string filename, float**& x_l, float**& y_l, float**& z_l, char**& elm_l, char**& lab_l, int*& i_res_l, int*& RecConf, float*& DockScore)
{
	cout << "Reading fragments.sdf" << std::endl;
	string line;

	i_res_l = new int[50000];
	x_l = new float* [50000];
	y_l = new float* [50000];
	z_l = new float* [50000];
	elm_l = new char* [50000];
	lab_l = new char* [50000];
	RecConf = new int[50000];
	DockScore = new float[50000];
	ifstream myfile(filename);
	if (myfile.is_open())
	{
		int groupId = 0;

		int i;

		while (!myfile.eof())
		{
			getline(myfile, line);
			if (line.empty()) continue;

			//skip 2 lines
			getline(myfile, line);
			getline(myfile, line);

			getline(myfile, line);
			i_res_l[groupId] = stoi(line.substr(0, 3));

			x_l[groupId] = new float[i_res_l[groupId]];
			y_l[groupId] = new float[i_res_l[groupId]];
			z_l[groupId] = new float[i_res_l[groupId]];
			elm_l[groupId] = new char[i_res_l[groupId]];
			lab_l[groupId] = new char[i_res_l[groupId]];

			for (i = 0; i < i_res_l[groupId]; i++)
			{
				getline(myfile, line);
				x_l[groupId][i] = stof(line.substr(0, 10));
				y_l[groupId][i] = stof(line.substr(10, 10));
				z_l[groupId][i] = stof(line.substr(20, 10));
				elm_l[groupId][i] = line[31];
				lab_l[groupId][i] = line[35];
			}


			//RecConf[groupId] = 1;
			groupId++;
			i_count_l = groupId;

			while (getline(myfile, line))
			{

				if (line.substr(0, 11) == "> <RecConf>")
				{
					getline(myfile, line);
					RecConf[groupId - 1] = stoi(line.substr(0, 1));
					//cout << "RecConf[groupId]" << RecConf[groupId - 1] << std::endl;
				}
				if (line.substr(0, 9) == "> <Score>")
				{
					getline(myfile, line);
					DockScore[groupId - 1] = stof(line.substr(0, 10));
					//cout << "RecConf[groupId]" << RecConf[groupId - 1] << std::endl;
				}
				if (line.substr(0, 4) == "$$$$")
					break;
			}
		}

	}
}
void protein_atom_count(string filename, int& i_protein_atom_count)
{
	i_protein_atom_count = 0;
	cout << "Loading" << filename.c_str() << std::endl;
	string line;
	ifstream myfile(filename);

	if (myfile.is_open())
	{
		int proteinId = 0;
		int atom_count;
		atom_count = 0;

		getline(myfile, line);
		//skip 2 lines
		getline(myfile, line);
		getline(myfile, line);
		getline(myfile, line);
		cout << " N of Protein Atoms_from the .sdf File" << stoi(line.substr(0, 4)) << std::endl;
		while (!myfile.eof())
		{
			getline(myfile, line);

			string str = line;

			if (str.length() > 30)
			{
				atom_count = atom_count + 1;
			}
			else
			{
				i_protein_atom_count = atom_count;
				cout << " N of Protein Atoms_calculated" << i_protein_atom_count << std::endl;
				break;
			}
		}
	}
}

void read_sdf_protein(string filename, float**& x_p, float**& y_p, float**& z_p, char**& elm_p, char**& lab_p, int*& i_res_p, int& i_protein_atom_count)
{
	cout << "Loading" << filename.c_str() << std::endl;

	string line;

	i_res_p = new int[10];

	x_p = new float* [10];
	y_p = new float* [10];
	z_p = new float* [10];
	elm_p = new char* [10];
	lab_p = new char* [10];

	ifstream myfile(filename);
	if (myfile.is_open())
	{
		int proteinId = 0;
		int i_count_p;

		int i;

		while (!myfile.eof())
		{
			getline(myfile, line);
			if (line.empty()) continue;

			//skip 2 lines
			getline(myfile, line);
			getline(myfile, line);

			getline(myfile, line);

			i_res_p[proteinId] = i_protein_atom_count;


			x_p[proteinId] = new float[i_res_p[proteinId]];
			y_p[proteinId] = new float[i_res_p[proteinId]];
			z_p[proteinId] = new float[i_res_p[proteinId]];
			elm_p[proteinId] = new char[i_res_p[proteinId]];
			lab_p[proteinId] = new char[i_res_p[proteinId]];

			for (i = 0; i < i_res_p[proteinId]; i++)
			{
				getline(myfile, line);
				x_p[proteinId][i] = stof(line.substr(0, 10));
				y_p[proteinId][i] = stof(line.substr(10, 10));
				z_p[proteinId][i] = stof(line.substr(20, 10));
				elm_p[proteinId][i] = line[31];
				lab_p[proteinId][i] = line[35];
			}
			proteinId++;
			i_count_p = proteinId + 1;
			while (getline(myfile, line))
			{
				if (line.substr(0, 4) == "$$$$")
					break;
			}
			break;
		}
	}
}



bool iterate(std::string input, std::string output, std::string delimiter)
{
	std::ifstream in(input.c_str());
	std::ofstream out(output);
	int it;
	it = 0;
	if (!in)
	{
		std::cerr << "Cannot open the File : " << input << std::endl;
		return false;
	}
	std::string line;
	while (std::getline(in, line))
	{
		if (line.substr(0, 4) == delimiter)
		{
			out << insert(it) << std::endl;
			out << delimiter << std::endl;
			it++;
		}
		else
		{
			out << line << std::endl;
		}
	}
	in.close();
	out.close();
	return true;
}

bool filecheck(const char* filename)
{
	ifstream file(filename);
	if (!file)
	{
		return false;
	}
	else
	{
		return true;
	}
}

int main()
{
	static int p_check, num_current_1, k1, k2, total, ip2_1, ip4, j2, ip3, num, i, i1, i2, i3, j1, i4, i5, ip1, ip2, i7, in, ipp;
	static int  j, ip, counter, k;
	static int  num_lab_l[50000], num_current;
	static float  m, l_check, x_in, y_in, z_in, r, r1, r3, rrr, rrrr, min11, a_c_f[500], a_c_t[500], a_s_f[500], a_s_t[500];
	static float  min, r2;
	static int cc, max_label[50000];
	static float  pi, xl, yl, zl, xl0, yl0, zl0;
	static float  x, y, z, rl, x_res, y_res, z_res, x2, y2, z2;
	static float  min1, fi, teta;
	static float  xl_1, yl_1, zl_1, xl_2, yl_2, zl_2, xl_3, yl_3, zl_3, xl_4, yl_4, zl_4;
	static int  num_lab_1, num_lab_2, num_lab_3, num_lab_4;
	float** x_l, ** y_l, ** z_l, ** x_p, ** y_p, ** z_p, * DockScore;
	char** elm_l, ** lab_l, ** elm_p, ** lab_p;
	int* i_res_l, * RecConf, * i_res_p;
	int i_protein_atom_count;



	float(*penalty_s)[15];
	penalty_s = new float[50000][15];

	float(*penalty_max1_l)[15];
	penalty_max1_l = new float[50000][15];

	float(*ff)[5][15];
	ff = new float[50000][5][15];

	float(*min1_l)[5][15];
	min1_l = new float[50000][5][15];

	float(*score)[15];
	score = new float[50000][15];

	float(*xll)[15];
	xll = new float[50000][15];

	float(*yll)[15];
	yll = new float[50000][15];

	float(*zll)[15];
	zll = new float[50000][15];

	int(*num1)[15];
	num1 = new int[50000][15];

	float(*x_resl)[5][15];
	x_resl = new float[50000][5][15];

	float(*y_resl)[5][15];
	y_resl = new float[50000][5][15];

	float(*z_resl)[5][15];
	z_resl = new float[50000][5][15];

	cout << setprecision(4) << fixed;
	cout << "    _             __                       \n" << std::endl;
	cout << "   /    _.  ._   (_    _   |   _    _  _|_ \n" << std::endl;
	cout << "   \\_  (_|  |_)  __)  (/_  |  (/_  (_   |_ \n" << std::endl;
	cout << "            |                              \n\n" << std::endl;

	read_sdf_frag("fragments.sdf", x_l, y_l, z_l, elm_l, lab_l, i_res_l, RecConf, DockScore);
	cout << "\nN of FRAGMENTS" << i_count_l << std::endl;
	//printf("\nN of FRAGMENTS \t %d", i_count_l);
	cout << "\n======================\n" << std::endl;
	//printf("\n======================\n");
	//Creating Sphere's Mesh
	pi = 3.14159264;
	for (i4 = 0; i4 < 72; i4++)	//azimut
	{
		teta = i4 * 2. * 5. * pi / 360.;
		a_c_t[i4] = cos(teta);
		a_s_t[i4] = sin(teta);
		//	printf("\n i4= %4d    teta= %4.16f   a_c_t= %4.4f    a_s_t= %4.4f ",i4,teta,a_c_t[i4],a_s_t[i4]);
	}
	for (i1 = 0; i1 < 36; i1++)	//longitude
	{
		fi = -pi / 2. + i1 * 5. * pi / 180.;
		a_c_f[i1] = cos(fi);
		a_s_f[i1] = sin(fi);
		// printf("\n fi= %4.4f   a_c_f= %4.4f    a_s_f= %4.4f ",fi,a_c_f[i1],a_s_f[i1]);
	}

	in = 0;
	k = 0;
	total = 0;
	for (i = 0; i < i_count_l; i++)	//N of fragments___________________________________
	{

		if (RecConf[i] == 1 && filecheck("protein1.sdf") == true)
		{
			protein_atom_count("protein1.sdf", i_protein_atom_count);
			read_sdf_protein("protein1.sdf", x_p, y_p, z_p, elm_p, lab_p, i_res_p, i_protein_atom_count);

		}

		else if (RecConf[i] == 2 && filecheck("protein2.sdf") == true)
		{
			protein_atom_count("protein2.sdf", i_protein_atom_count);
			read_sdf_protein("protein2.sdf", x_p, y_p, z_p, elm_p, lab_p, i_res_p, i_protein_atom_count);
		}

		else if (RecConf[i] != 2 && RecConf[i] != 1)
		{
			cout << "\nERROR: <RecConf> information invalid." << std::endl;
			break;
		}

		else if (filecheck("protein1.sdf") == false)
		{
			cout << "\nERROR: No protein file for <RecCon> 1 found!" << std::endl;
			break;
		}

		else if (filecheck("protein2.sdf") == false)
		{
			cout << "\nERROR: No protein file for <RecCon> 2 found!" << std::endl;
			break;
		}

		else {
			cout << "\nERROR: <RecConf> information invalid." << std::endl;
			break;
		}


		num = 0;
		ipp = 0;
		xl_1 = 0.;
		yl_1 = 0.;
		zl_1 = 0.;

		xl_2 = 0.;
		yl_2 = 0.;
		zl_2 = 0.;

		xl_3 = 0.;
		yl_3 = 0.;
		zl_3 = 0.;

		xl_4 = 0.;
		yl_4 = 0.;
		zl_4 = 0.;
		num_lab_1 = 0;
		num_lab_2 = 0;
		num_lab_3 = 0;
		num_lab_4 = 0;
		num_current_1 = 0;
		max_label[i] = 0;
		num_lab_l[i] = 0;

		for (i1 = 0; i1 < i_res_l[i]; i1++)	//N of fragments
		{
			if (num_lab_1 == 5)
			{
				if (lab_l[i][i1] == 51)
				{
					xl_2 = xl_2 + x_l[i][i1];
					yl_2 = yl_2 + y_l[i][i1];
					zl_2 = zl_2 + z_l[i][i1];
					num_lab_2 = num_lab_2 + 1;
				}
			}
			if (num_lab_1 < 5)
			{
				if (lab_l[i][i1] == 51)
				{
					xl_1 = xl_1 + x_l[i][i1];
					yl_1 = yl_1 + y_l[i][i1];
					zl_1 = zl_1 + z_l[i][i1];
					num_lab_1 = num_lab_1 + 1;
				}
			}
		}

		cout << "\nAROMATIC_CAP: num_lab_1=" << num_lab_1 << "    num_lab_2 = " << num_lab_2 << std::endl;

		for (i1 = 0; i1 < i_res_l[i]; i1++)	//N of fragments
		{
			if (num_lab_3 == 1)
			{
				if (lab_l[i][i1] == 49)
				{
					xl_4 = x_l[i][i1];
					yl_4 = y_l[i][i1];
					zl_4 = z_l[i][i1];
					num_lab_4 = num_lab_4 + 1;
				}
			}
			if (num_lab_3 < 1)
			{
				if (lab_l[i][i1] == 49)
				{
					xl_3 = x_l[i][i1];
					yl_3 = y_l[i][i1];
					zl_3 = z_l[i][i1];
					num_lab_3 = num_lab_3 + 1;

				}//if
			}

		}//i1
		cout << "\nNON_AROMATIC_CAP: num_lab_3=" << num_lab_3 << "    num_lab_4 = " << num_lab_4 << std::endl;

		if ((num_lab_1 == 0) && (num_lab_2 == 0) && (num_lab_3 == 0) && (num_lab_4 == 0))
		{
			CapScore[i] = -100.;//added on 09/23/2021
			num1_fin[i] = 1;
			goto pp3; //That means fragments does not have a CAP - Should go to the next fragment analysis
		}
		if ((num_lab_1 == 5) && (num_lab_2 == 5))
		{
			max_label[i] = 3;
			xll[i][1] = xl_1 / num_lab_1;
			yll[i][1] = yl_1 / num_lab_1;
			zll[i][1] = zl_1 / num_lab_1;

			xll[i][2] = xl_2 / num_lab_2;
			yll[i][2] = yl_2 / num_lab_2;
			zll[i][2] = zl_2 / num_lab_2;
			cout << "\n2_AROMATIC_CAPS: x1= " << xll[i][1] << "    y1= " << yll[i][1] << "    z1= " << zll[i][1] << std::endl;
			cout << "\n2_AROMATIC_CAPS: x2= " << xll[i][2] << "    y2= " << yll[i][2] << "    z2= " << zll[i][2] << std::endl;
		}
		if ((num_lab_3 == 1) && (num_lab_4 == 1))
		{
			max_label[i] = 3;
			xll[i][1] = xl_3 / num_lab_3;
			yll[i][1] = yl_3 / num_lab_3;
			zll[i][1] = zl_3 / num_lab_3;

			xll[i][2] = xl_4 / num_lab_4;
			yll[i][2] = yl_4 / num_lab_4;
			zll[i][2] = zl_4 / num_lab_4;
			cout << "\n2_NON_AROMATIC_CAPS: x1= " << xll[i][1] << "    y1= " << yll[i][1] << "    z1= " << zll[i][1] << std::endl;
			cout << "\n2_NON_AROMATIC_CAPS: x2= " << xll[i][2] << "    y2= " << yll[i][2] << "    z2= " << zll[i][2] << std::endl;
		}

		if ((num_lab_1 == 5) && (num_lab_3 == 1))
		{
			max_label[i] = 3;
			xll[i][1] = xl_1 / num_lab_1;
			yll[i][1] = yl_1 / num_lab_1;
			zll[i][1] = zl_1 / num_lab_1;

			xll[i][2] = xl_3 / num_lab_3;
			yll[i][2] = yl_3 / num_lab_3;
			zll[i][2] = zl_3 / num_lab_3;
			cout << "\nAROM_NON_AROM_CAPS: x1= " << xll[i][1] << "    y1= " << yll[i][1] << "    z1= " << zll[i][1] << std::endl;
			cout << "\nAROM_NON_AROM_CAPS: x2= " << xll[i][2] << "    y2= " << yll[i][2] << "    z2= " << zll[i][2] << std::endl;
		}

		if (((num_lab_1 == 5) && (num_lab_2 == 0)) && ((num_lab_3 == 0) && (num_lab_4 == 0)))
		{
			max_label[i] = 2;
			xll[i][1] = xl_1 / num_lab_1;
			yll[i][1] = yl_1 / num_lab_1;
			zll[i][1] = zl_1 / num_lab_1;
			cout << "\nSINGLE_CAP: x1= " << xll[i][1] << "    y1= " << yll[i][1] << "    z1= " << zll[i][1] << std::endl;
		}

		if (((num_lab_1 == 0) && (num_lab_2 == 0)) && ((num_lab_3 == 1) && (num_lab_4 == 0)))
		{
			max_label[i] = 2;
			xll[i][1] = xl_3 / num_lab_3;
			yll[i][1] = yl_3 / num_lab_3;
			zll[i][1] = zl_3 / num_lab_3;
			cout << "\nSINGLE_CAP: x1= " << xll[i][1] << "    y1= " << yll[i][1] << "    z1= " << zll[i][1] << std::endl;
		}

		num_lab_l[i] = num_lab_1 + num_lab_2 + num_lab_3 + num_lab_4;
		if ((num_lab_l[i] == 1) || (num_lab_l[i] == 5))
		{
			max_label[i] = 2;
			num_current_1 = 1;
		}
		if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
		{
			max_label[i] = 3;
			num_current_1 = 1;
		}
		for (cc = 1; cc < max_label[i]; cc++)
		{		
			cout << "\nCAP COORDINATES_ALL           frag=" << i << "    max_label= " << max_label[i] << "    xl= " << xll[i][cc] << "    yl= " << yll[i][cc] << "    zl= " << zll[i][cc] << std::endl;
		}

		for (num_current = 1; num_current < max_label[i]; num_current++)
		{
			num = 0;
			x_in = xll[i][num_current];
			y_in = yll[i][num_current];
			z_in = zll[i][num_current];
			x_resl[i][num_current][0] = x_in;
			y_resl[i][num_current][0] = y_in;
			z_resl[i][num_current][0] = z_in;
			
			cout << "\nCAP COORDINATES               rd=" << num_current << "    frag= " << i << "    max_label= " << max_label[i] << "    xl= " << xll[i][num_current] << "    yl= " << yll[i][num_current] << "    zl= " << zll[i][num_current] << std::endl;

			//5. ______________________ XX1,YY1,ZZ1 - based MESH generation_________
			//***************************  0-sphere generation (CAP-sphere) STEP ******************************************

			pi = 3.14159264;
			r = 3.0; 
			min11 = 0.;
			if (num_current == 1)
			{
				if (num_lab_1 == 5)
				{
					r = 3.5;
				}
			}
			if (num_current == 2)
			{
				if (num_lab_2 == 5)
				{
					r = 3.5;
				}
			}

			for (i4 = 0; i4 < 72; i4++)	//azimut
			{
				for (i1 = 0; i1 < 36; i1++)	//longitude
				{

					z = r * a_s_f[i1];
					x = r * a_c_f[i1] * a_s_t[i4];
					y = r * a_c_f[i1] * a_c_t[i4];

					//_________ MIN _______


					z = z + zll[i][num_current];
					y = y + yll[i][num_current];
					x = x + xll[i][num_current];

					l_check = 1.3;

					if (num_current == 1)
					{
						if (num_lab_1 == 5)
						{
							l_check = 1.1;
						}
					}
					if (num_current == 2)
					{
						if (num_lab_2 == 5)
						{
							l_check = 1.1;
						}
					}

					ip = 0;
					for (i3 = 0; i3 < i_res_l[i]; i3++)
					{
						if ((lab_l[i][i3] != 49) && (elm_l[i][i3] != 72) && (lab_l[i][i3] != 51))
						{
							r2 = sqrt((x - x_l[i][i3]) * (x - x_l[i][i3]) + (y - y_l[i][i3]) * (y - y_l[i][i3]) + (z - z_l[i][i3]) * (z - z_l[i][i3]));
							r1 = sqrt((x_in - x_l[i][i3]) * (x_in - x_l[i][i3]) + (y_in - y_l[i][i3]) * (y_in - y_l[i][i3]) + (z_in - z_l[i][i3]) * (z_in - z_l[i][i3]));
							if ((r2 < r) || (r1 < l_check)) { ip = 1; }
						}//if
					}//i3

					p_check = 3.0;
					if (num_current == 1)
					{
						if (num_lab_1 == 5)
						{
							p_check = 2.0;
						}
					}
					if (num_current == 2)
					{
						if (num_lab_2 == 5)
						{
							p_check = 2.0;
						}
					}
					ip1 = 0;
					for (i3 = 0; i3 < i_res_p[0]; i3++)
					{
						{
							r2 = sqrt((x - x_p[0][i3]) * (x - x_p[0][i3]) + (y - y_p[0][i3]) * (y - y_p[0][i3]) + (z - z_p[0][i3]) * (z - z_p[0][i3]));
						}
						if ((r2 < p_check)) { ip1 = 1; }
					}//i3

					if ((ip != 1) && (ip1 != 1)) { ipp = 1; }
					if ((ip == 1) || (ip1 == 1)) { goto c_lig1; }

					min = 1000.;
					for (i2 = 0; i2 < i_res_p[0]; i2++)	//
					{
						{
							r1 = sqrt((x - x_p[0][i2]) * (x - x_p[0][i2]) + (y - y_p[0][i2]) * (y - y_p[0][i2]) + (z - z_p[0][i2]) * (z - z_p[0][i2]));
						}
						if (min > r1) { min = r1; x2 = x - xl * 0.; y2 = y - yl * 0.; z2 = z - zl * 0.; }
					}//i2
					if (min11 < min) { min11 = min; x_res = x2; y_res = y2; z_res = z2; }
				c_lig1:;
				}//i1
			}//i4

					
						
						
			if (ipp == 0)
			{
				num1[i][num_current] = 0;
				min1_l[i][num_current][0] = 0.;
				ff[i][num_current][0] = 0.;
				if ((num_lab_l[i] == 2) || (num_lab_l[i] == 6) || (num_lab_l[i] == 10))
				{
					if (num_current == 2)
					{
						goto pp;
					}
					if (num_current == 1)
					{
						goto pp1;
					}
				}
				goto pp;
			}
			num = 1; 
			x_resl[i][num_current][1] = x_res;
			y_resl[i][num_current][1] = y_res;
			z_resl[i][num_current][1] = z_res;
			min1_l[i][num_current][1] = min11;

			ff[i][num_current][1] = sqrt((x_in - x_resl[i][num_current][1]) * (x_in - x_resl[i][num_current][1]) + (y_in - y_resl[i][num_current][1]) * (y_in - y_resl[i][num_current][1]) + (z_in - z_resl[i][num_current][1]) * (z_in - z_resl[i][num_current][1]));
			cout << "\nCAP                           frag=" << i+1 << "    num= " << num << "    MinMax= " << min1_l[i][num_current][1] << "    ff= " << ff[i][num_current][1] << "    x_res= " << x_resl[i][num_current][1] << "    y_res= " << y_resl[i][num_current][1] << "    z_res= " << z_resl[i][num_current][1] << std::endl;
			counter = counter + ipp;

			//***************************  NON-0 sphere generation ******************************************
		www:
			ipp = 0;
			rrr = 5.;//corrected 04072022 
			rrrr = 2.;
			zl = z_res;
			yl = y_res;
			xl = x_res;


			min11 = 0.;
			for (i4 = 0; i4 < 72; i4++)	//azimut
			{
				for (i1 = 0; i1 < 36; i1++)	//longitude
				{
					z = rrrr * a_s_f[i1];
					x = rrrr * a_c_f[i1] * a_s_t[i4];
					y = rrrr * a_c_f[i1] * a_c_t[i4];

			//_________ SEARCH of the MIN _______

					z = z + zl;
					y = y + yl;
					x = x + xl;

					ip = 0;
					l_check = 1.3;
					if (num_current == 1)
					{
						if (num_lab_1 == 5)
						{
							l_check = 1.1;
						}
					}
					if (num_current == 2)
					{
						if (num_lab_2 == 5)
						{
							l_check = 1.1;
						}
					}

					for (i3 = 0; i3 < i_res_l[i]; i3++)
					{
						if ((lab_l[i][i3] != 49) && (elm_l[i][i3] != 72) && (lab_l[i][i3] != 51))
						{
							r2 = sqrt((x - x_l[i][i3]) * (x - x_l[i][i3]) + (y - y_l[i][i3]) * (y - y_l[i][i3]) + (z - z_l[i][i3]) * (z - z_l[i][i3]));
							r1 = sqrt((xl - x_l[i][i3]) * (xl - x_l[i][i3]) + (yl - y_l[i][i3]) * (yl - y_l[i][i3]) + (zl - z_l[i][i3]) * (zl - z_l[i][i3]));
							if ((r2 < rrr) || (r1 < l_check)) { ip = 1; }

						}//if
					}//i3
					ip1 = 0;
					for (i3 = 0; i3 < i_res_p[0]; i3++)
					{
						{
							r2 = sqrt((x - x_p[0][i3]) * (x - x_p[0][i3]) + (y - y_p[0][i3]) * (y - y_p[0][i3]) + (z - z_p[0][i3]) * (z - z_p[0][i3]));
							r1 = sqrt((xl - x_p[0][i3]) * (xl - x_p[0][i3]) + (yl - y_p[0][i3]) * (yl - y_p[0][i3]) + (zl - z_p[0][i3]) * (zl - z_p[0][i3]));
						}
						if ((r2 < 2.0) || (r1 < 3.0)) { ip1 = 1; }

					}//i3

					ip2 = 0;

					if (num > 1)
					{
						for (i7 = 0; i7 < num; i7++)
						{
							r1 = sqrt((x - x_resl[i][num_current][i7]) * (x - x_resl[i][num_current][i7]) + (y - y_resl[i][num_current][i7]) * (y - y_resl[i][num_current][i7]) + (z - z_resl[i][num_current][i7]) * (z - z_resl[i][num_current][i7]));							
							if (r1 < 2.0) { ip2 = 1; }							
						}//i7
					}
					ip2_1 = 0;
					if (num == 1)
					{
						r2 = sqrt((x - x_in) * (x - x_in) + (y - y_in) * (y - y_in) + (z - z_in) * (z - z_in));//added on 09/17/21						
						if (num_current == 1)
						{							
							if (num_lab_1 == 5)
							{
								if ((r2 < 3.5)) { ip2_1 = 1; }//added on 09/17/21
							}
						}
						if (num_current == 2)
						{
							if (num_lab_2 == 5)
							{
								if ((r2 < 3.5)) { ip2_1 = 1; }//added on 09/17/21								
							}
						}

						if (num_current == 1)
						{							
							if ((num_lab_3 == 1) && (num_lab_1 != 5))
							{
								if ((r2 < 3.0)) 
								{ 
									ip2_1 = 1;								 
								}//added on 09/17/21						
							}
						}
						if (num_current == 2)
						{
							if (num_lab_4 == 1)
							{
								if ((r2 < 3.0)) 
								{ 
									ip2_1 = 1; 								
								}//added on 09/17/21
								
							}
						}
						if (num_current == 2)
						{
							if (num_lab_3 == 1)
							{
								if ((r2 < 3.0)) 
								{ 
									ip2_1 = 1;								
								}//added on 09/17/21
								
							}
						}

					}

					ip3 = 0;
					if (num > 1)
					{
						for (j2 = (num - 2); j2 < num; j2++) 
						{
							r2 = sqrt((x - x_resl[i][num_current][j2]) * (x - x_resl[i][num_current][j2]) + (y - y_resl[i][num_current][j2]) * (y - y_resl[i][num_current][j2]) + (z - z_resl[i][num_current][j2]) * (z - z_resl[i][num_current][j2]));
							if (r2 < 3.46) { ip3 = 1; } //60degree
						}
					}

					ip4 = 0;

					if ((ip != 1) && (ip1 != 1) && (ip2 != 1) && (ip3 != 1) && (ip4 != 1) && (ip2_1 != 1)) { ipp = 1; }
					if ((ip == 1) || (ip1 == 1) || (ip2 == 1) || (ip3 == 1) || (ip4 == 1) || (ip2_1 == 1)) { goto c_lig; }

					min = 1000.;
					for (i2 = 0; i2 < i_res_p[0]; i2++)	//
					{
						{
							r1 = sqrt((x - x_p[0][i2]) * (x - x_p[0][i2]) + (y - y_p[0][i2]) * (y - y_p[0][i2]) + (z - z_p[0][i2]) * (z - z_p[0][i2]));
						}
						if (min > r1) { min = r1; x2 = x - xl * 0.; y2 = y - yl * 0.; z2 = z - zl * 0.; }
					}//i2
					if (min11 < min) { min11 = min; x_res = x2; y_res = y2; z_res = z2; }

				c_lig:;
				}//i1
			}//i4
			if ((ipp == 0) && (num == 1))
			{
				num1[i][num_current] = 1;
				ff[i][num_current][0] = 0.;
				min1_l[i][num_current][num] = min1_l[i][num_current][1];//corrected 04/06/2022 (was 0)
				
				if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
				{
					if (num_current == 2)
					{

						goto pp;
					}
					if (num_current == 1)
					{
						goto pp1;
					}
				}
				goto pp;
			}



			if ((ipp == 0) && (num > 1) && (num <= 9)) 
			{
				num1[i][num_current] = num;

				if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
				{
					if (num_current == 2)
					{

						goto pp;
					}

					if (num_current == 1)
					{
						goto pp1;
					}
				}
				goto pp;

			}
			num = num + 1; 
			x_resl[i][num_current][num] = x_res;
			y_resl[i][num_current][num] = y_res;
			z_resl[i][num_current][num] = z_res;
			min1_l[i][num_current][num] = min11;

			ff[i][num_current][num] = sqrt((x_resl[i][num_current][num] - x_in) * (x_resl[i][num_current][num] - x_in) + (y_resl[i][num_current][num] - y_in) * (y_resl[i][num_current][num] - y_in) + (z_resl[i][num_current][num] - z_in) * (z_resl[i][num_current][num] - z_in));
			cout << "\n                            frag=" << i << "    num= " << num << "    MinMax= " << min1_l[i][num_current][num] << "    ff= " << ff[i][num_current][num] << "    x_res= " << x_resl[i][num_current][num] << "    y_res= " << y_resl[i][num_current][num] << "    z_res= " << z_resl[i][num_current][num] << std::endl;

			if (num > 9)
			{
				num1[i][num_current] = num;
				
				if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
				{
					if (num_current == 2)
					{

						goto pp;
					}
					if (num_current == 1)
					{

						goto pp1;
					}
				}
				goto pp;

			}//if
			goto www;

		pp1:;

		}

	pp:

	//______________________score_penalty____________________________________

		if ((num_lab_l[i] == 1) || (num_lab_l[i] == 5))
		{
			for (j1 = 1; j1 < 2; j1++)
			{
				if (num1[i][j1] <= 5)
				{
					penalty_s[i][j1] = (5 - num1[i][j1]) * (5 - num1[i][j1]);
				}
				else
				{
					penalty_s[i][j1] = 0.;
				}
			}
		}
		if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
		{
			for (j2 = 1; j2 < 3; j2++)
			{
				if (num1[i][j2] <= 5)
				{
					penalty_s[i][j2] = (5 - num1[i][j2]) * (5 - num1[i][j2]);
				}
				else
				{
					penalty_s[i][j2] = 0.;
				}
			}
		}


		if ((num_lab_l[i] == 1) || (num_lab_l[i] == 5))
		{
			for (j1 = 1; j1 < 2; j1++)
			{
				penalty_s[i][j1] = (penalty_s[i][j1] * 10.) / 25.;
			}
		}
		if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
		{
			for (j2 = 1; j2 < 3; j2++)
			{
				penalty_s[i][j2] = (penalty_s[i][j2] * 10.) / 25.;

			}
		}


		if ((num_lab_l[i] == 1) || (num_lab_l[i] == 5))
		{
			for (j1 = 1; j1 < 2; j1++)
			{
				if (num1[i][j1] > 9)
				{
					if ((min1_l[i][j1][9] >= 7.) && (min1_l[i][j1][9] <= 20.))
					{
						penalty_max1_l[i][j1] = pow((7.0 - min1_l[i][j1][9]), 2.);
					}

					if (min1_l[i][j1][9] > 20.)
					{
						penalty_max1_l[i][j1] = 169.;
					}

					if (min1_l[i][j1][9] < 7.)
					{
						penalty_max1_l[i][j1] = 0.;
					}

				}
				if (num1[i][j1] <= 9)
				{
					penalty_max1_l[i][j1] = 0.;
				}
			}
		}
		if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
		{
			for (j2 = 1; j2 < 3; j2++)
			{
				if (num1[i][j2] > 9)
				{
					if ((min1_l[i][j2][9] >= 7.) && (min1_l[i][j2][9] <= 20.))
					{
						penalty_max1_l[i][j2] = pow((7.0 - min1_l[i][j2][9]), 2.);
					}

					if (min1_l[i][j2][9] > 20.)
					{
						penalty_max1_l[i][j2] = 169.;
					}

					if (min1_l[i][j2][9] < 7.)
					{
						penalty_max1_l[i][j2] = 0.;
					}

				}
				if (num1[i][j2] <= 9)
				{
					penalty_max1_l[i][j2] = 0.;
				}
			}
		}




		if ((num_lab_l[i] == 1) || (num_lab_l[i] == 5))
		{
			for (j1 = 1; j1 < 2; j1++)
			{
				penalty_max1_l[i][j1] = (penalty_max1_l[i][j1] * 10.) / 169.;
			}
		}
		if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
		{
			for (j2 = 1; j2 < 3; j2++)
			{
				penalty_max1_l[i][j2] = (penalty_max1_l[i][j2] * 10.) / 169.;
			}
		}



		if ((num_lab_l[i] == 1) || (num_lab_l[i] == 5))
		{
			score[i][1] = 10. - penalty_s[i][1] - penalty_max1_l[i][1];
			CapScore[i] = score[i][1];
			num1_fin[i] = num1[i][1];
			for (k = 0; k < num1[i][1]; k++)
			{
				min1_l_fin[i][k] = min1_l[i][1][k];
				ff_fin[i][k] = ff[i][1][k];

			}

		}

		if ((num_lab_l[i] == 10) || (num_lab_l[i] == 6) || (num_lab_l[i] == 2))
		{
			for (j2 = 1; j2 < 3; j2++)
			{
				score[i][j2] = 10. - penalty_s[i][j2] - penalty_max1_l[i][j2];
				m = score[i][1];
				CapScore[i] = m;
				num1_fin[i] = num1[i][1];
				for (k1 = 0; k1 < num1[i][j2]; k1++)
				{
					min1_l_fin[i][k1] = min1_l[i][1][k1];
					ff_fin[i][k1] = ff[i][1][k1];
				}
				if (m < score[i][j2])
				{
					CapScore[i] = score[i][2];
					num1_fin[i] = num1[i][2];
					for (k2 = 0; k2 < num1[i][2]; k2++)
					{
						min1_l_fin[i][k2] = min1_l[i][2][k2];
						ff_fin[i][k2] = ff[i][2][k2];
					}
				}
			}
		}

		cout << "\nCOUNTER i=" << i << "    counter= " << counter << std::endl;
	pp3:;
	}

	for (i = 0; i < i_count_l; i++)	//MergedScore_function___________________________________
	{
		MergedScore[i] = 5.0 * log2(abs(DockScore[i])) + 0.5 * log2(abs(CapScore[i]));
		if (CapScore[i] == -100.)
		{
			MergedScore[i] = -1000.;
		}
		if ((CapScore[i] == 0.) && (num1_fin[i] == 0))
		{
			MergedScore[i] = -1000.;
		}
		if ((CapScore[i] == 0.) && (num1_fin[i] == 10))
		{
			MergedScore[i] = 5.0 * log2(abs(DockScore[i]));
		}
	}

	std::string input = "fragments.sdf";
	std::string output = "CapSelect.sdf";
	std::string delimiter = "$$$$";
	bool res = iterate(input, output, delimiter);
	if (res)
	{
		cout << "\nDone!CapSelect.sdf written \n" << std::endl;
	}
	else
	{
		cout << "\nError writing into the file!" << std::endl;
	};

}
