/*  *********************************************************************** 	*/
/*  Calculate dihedral angle 'phi' and 'psi' from "ATOM" record in PDB file.	*/
/*  "ATOM" record is extracted from PDB file by 'getATOM_phipsi.pl'				*/
/*																				*/
/*  'phi' and 'psi' are dihedral angles created by								*/
/*				C-N-CA-C and N-CA-C-N, respectively.							*/
/*																				*/
/*																				*/
/*  Usage   :	perl  getATOM_phipsi.pl  pdb****.ent > FILE						*/
/*		./getdihed < FILE														*/
/*																				*/
/*  Results :   res_name  res_num  'dihedral angle'  atom_name  chainID			*/
/*																				*/
/*  When you want to calculate the dihedral angle only from "CA" atoms,			*/
/*  substitute 'getATOM_phipsi.pl' for 'getATOM_ca.pl'.							*/
/*																				*/
/*  ***********************************************************************		*/


#include<stdio.h>
#include<math.h>

#define pi 3.1415926535897932384626433832795


struct coord {
	char	record[4];
	int	atom_num;
	char	atom_name[4];
	char	res_name[4];
	char	chainID[1];
	int	res_num;
	double	x;
	double	y;
	double	z;
	double	gomi1;
	double	gomi2;
	char	atom[1];
};


int main()
{
	/* Definition of variables */
	int		i, j, k, l, m,
			count = 0;

	double		r_ij_x, r_ij_y, r_ij_z,
		r_jk_x, r_jk_y, r_jk_z,
		r_kl_x, r_kl_y, r_kl_z,
		r_lm_x, r_lm_y, r_lm_z,
		op_ijk_x, op_ijk_y, op_ijk_z,
		op_jkl_x, op_jkl_y, op_jkl_z,
		op_klm_x, op_klm_y, op_klm_z,
		ip_ijkl, ip_jklm,
		ip_rij_opjkl, ip_rjk_opklm,
		abs_ijk, abs_jkl, abs_klm,
		abs_ij, abs_jk,
		value1, value2, value3, value4,
		psi, phi;

	struct coord	coord_a[20000];


	/* Read input data from pdb file */
	while (scanf("%d %lf %lf %lf %lf %lf",
		&coord_a[count].res_num,
		&coord_a[count].x,
		&coord_a[count].y,
		&coord_a[count].z,
		&coord_a[count].gomi1,
		&coord_a[count].gomi2)!=EOF)
	{
		count++;
	}

	/*
	while (scanf("%s %d %s %s %s %d %lf %lf %lf %lf %lf %s",
		&coord_a[count].record,
		&coord_a[count].atom_num,
		&coord_a[count].atom_name,
		&coord_a[count].res_name,
		&coord_a[count].chainID,
		&coord_a[count].res_num,
		&coord_a[count].x,
		&coord_a[count].y,
		&coord_a[count].z,
		&coord_a[count].gomi1,
		&coord_a[count].gomi2,
		&coord_a[count].atom)!=EOF)
	{
		count++;
	}
	*/

	/* dihedral angle */
	for (i = 1; i + 3 < count; i += 3){
		/** components of vector r_ij, r_jk, r_kl, and r_lm **/
		r_ij_x = coord_a[i].x - coord_a[i-1].x; 
		r_ij_y = coord_a[i].y - coord_a[i-1].y;
		r_ij_z = coord_a[i].z - coord_a[i-1].z;
		r_jk_x = coord_a[i+1].x - coord_a[i].x;
		r_jk_y = coord_a[i+1].y - coord_a[i].y;
		r_jk_z = coord_a[i+1].z - coord_a[i].z;
		r_kl_x = coord_a[i+2].x - coord_a[i+1].x;
		r_kl_y = coord_a[i+2].y - coord_a[i+1].y;
		r_kl_z = coord_a[i+2].z - coord_a[i+1].z;
		r_lm_x = coord_a[i+3].x - coord_a[i+2].x;
		r_lm_y = coord_a[i+3].y - coord_a[i+2].y;
		r_lm_z = coord_a[i+3].z - coord_a[i+2].z;

		/** outer products
	    (op_ijk = rij x rjk, op_jkl = rjk x rkl, and op_klm = rkl x rlm) **/

		op_ijk_x = r_ij_y*r_jk_z - r_jk_y*r_ij_z; 
		op_ijk_y = r_ij_z*r_jk_x - r_jk_z*r_ij_x; 
		op_ijk_z = r_ij_x*r_jk_y - r_jk_x*r_ij_y;

		op_jkl_x = r_jk_y*r_kl_z - r_kl_y*r_jk_z; 
		op_jkl_y = r_jk_z*r_kl_x - r_kl_z*r_jk_x; 
		op_jkl_z = r_jk_x*r_kl_y - r_kl_x*r_jk_y;

		op_klm_x = r_kl_y*r_lm_z - r_lm_y*r_kl_z; 
		op_klm_y = r_kl_z*r_lm_x - r_lm_z*r_kl_x; 
		op_klm_z = r_kl_x*r_lm_y - r_lm_x*r_kl_y;

		/** innner product (op_ijk - op_jkl, op_jkl - op_klm) and absolute values of outer products,op_ijk,op_jkl, op_klm **/
		/** components for calculation of cos(dihedral) **/

		ip_ijkl = op_ijk_x * op_jkl_x + op_ijk_y * op_jkl_y + op_ijk_z * op_jkl_z;
		ip_jklm = op_jkl_x * op_klm_x + op_jkl_y * op_klm_y + op_jkl_z * op_klm_z;
		abs_ijk = sqrt(op_ijk_x * op_ijk_x + op_ijk_y * op_ijk_y + op_ijk_z * op_ijk_z);
		abs_jkl = sqrt(op_jkl_x * op_jkl_x + op_jkl_y * op_jkl_y + op_jkl_z * op_jkl_z);
		abs_klm = sqrt(op_klm_x * op_klm_x + op_klm_y * op_klm_y + op_klm_z * op_klm_z);
     
		/** innner product (r_ij-op_jkl,rjk-op_klm) and absolute values r_ij,r_jk **/
		/** components for calculation of sin(dihedral) (denomitator is same as that of cos(dihedral)) **/

		ip_rij_opjkl = r_ij_x * op_jkl_x + r_ij_y * op_jkl_y + r_ij_z * op_jkl_z;
		ip_rjk_opklm = r_jk_x * op_klm_x + r_jk_y * op_klm_y + r_jk_z * op_klm_z;
		abs_ij = sqrt(r_ij_x * r_ij_x + r_ij_y * r_ij_y + r_ij_z * r_ij_z);
		abs_jk = sqrt(r_jk_x * r_jk_x + r_jk_y * r_jk_y + r_jk_z * r_jk_z);

		/** cos(dihedral) **/
		value1 = ip_ijkl / (abs_ijk * abs_jkl);
		phi = acos(value1);
		value2 = ip_jklm / (abs_jkl * abs_klm);
		psi = acos(value2);

		/** sin(dihedral) **/
		value3 = (ip_rij_opjkl * abs_ij) / (abs_ijk * abs_jkl);
		value4 = (ip_rjk_opklm * abs_jk) / (abs_jkl * abs_klm);
        
		/** dicision of the range of dihedral (0 to 180 or -180 to 0) **/ 
		if(value3<0) phi = -1.*phi;
		if(value4<0) psi = -1.*psi;

		/** Output data **/
		//printf("%s\t", coord_a[i].res_name);
		//printf("%d\t", coord_a[i].res_num);
		printf("%.2lf\n", 180 * phi / pi);
		//printf("  %s  - %s\n", coord_a[i].atom_name, coord_a[i+1].atom_name);
//		printf("%s\n", coord_a[i].chainID);
		//printf("%s\t", coord_a[i+1].res_name);
		//printf("%d\t", coord_a[i+1].res_num);
		//printf("%.2lf\t", 180 * psi / pi);
		//printf("  %s - %s\n", coord_a[i+1].atom_name, coord_a[i+2].atom_name);
//		printf("%s\n", coord_a[i+1].chainID);
	}
}

