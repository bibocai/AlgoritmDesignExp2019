/* sample.c */

#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>




int SetCoverLP(int *ia,int *ja, int numMatrixValueNotZero,double* ar,double* x,int numX)//ia,ja,ar描述参数矩阵 从1开始 numMatrixValueNotZero描述参数矩阵中的有效参数个数 
																							//x 是未知数 numX是x的长度 这里因为我们题目中默认是个方针
{
	glp_prob* lp;
	double z;
	//double x1, x2, x3;
s1:   lp = glp_create_prob();
	 
s2:   glp_set_prob_name(lp, "lp");
s3:   glp_set_obj_dir(lp, GLP_MIN);
s4:   glp_add_rows(lp, numX);

	for (int i = 1; i < numX+1; i++) {
		char rowName[7];
		sprintf_s(rowName, "%s%d", "p", i);
		glp_set_row_name(lp, i , rowName); //从1开始
		glp_set_row_bnds(lp, i, GLP_LO, 1.0, 0);
	}
	
	glp_add_cols(lp, numX);
	for (int i = 1; i < numX + 1; i++) {
		char colName[7];
		sprintf_s(colName, "%s%d", "x", i);
		glp_set_col_name(lp, i, colName); //从1开始
		glp_set_col_bnds(lp, i, GLP_LO, 0.0, 0);
		glp_set_obj_coef(lp, i, 1);
	}



//s5:   glp_set_row_name(lp, 1, "p");
//s6:   glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 100.0);
//s7:   glp_set_row_name(lp, 2, "q");
//s8:   glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 600.0);
//s9:   glp_set_row_name(lp, 3, "r");
//s10:  glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 300.0);
//
//s11:  glp_add_cols(lp, 3);
//s12:  glp_set_col_name(lp, 1, "x1");
//s13:  glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0);
//s14:  glp_set_obj_coef(lp, 1, 10.0);
//s15:  glp_set_col_name(lp, 2, "x2");
//s16:  glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0);
//s17:  glp_set_obj_coef(lp, 2, 6.0);
//s18:  glp_set_col_name(lp, 3, "x3");
//s19:  glp_set_col_bnds(lp, 3, GLP_LO, 0.0, 0.0);
//s20:  glp_set_obj_coef(lp, 3, 4.0);
//
//s21:  ia[1] = 1, ja[1] = 1, ar[1] = 1.0; /* a[1,1] =  1 */
//s22:  ia[2] = 1, ja[2] = 2, ar[2] = 1.0; /* a[1,2] =  1 */
//s23:  ia[3] = 1, ja[3] = 3, ar[3] = 1.0; /* a[1,3] =  1 */
//s24:  ia[4] = 2, ja[4] = 1, ar[4] = 10.0; /* a[2,1] = 10 */
//s25:  ia[5] = 3, ja[5] = 1, ar[5] = 2.0; /* a[3,1] =  2 */
//s26:  ia[6] = 2, ja[6] = 2, ar[6] = 4.0; /* a[2,2] =  4 */
//s27:  ia[7] = 3, ja[7] = 2, ar[7] = 2.0; /* a[3,2] =  2 */
//s28:  ia[8] = 2, ja[8] = 3, ar[8] = 5.0; /* a[2,3] =  5 */
//s29:  ia[9] = 3, ja[9] = 3, ar[9] = 6.0; /* a[3,3] =  6 */
s30:  glp_load_matrix(lp, numMatrixValueNotZero, ia, ja, ar);
s31:  glp_simplex(lp, NULL);
s32:  z = glp_get_obj_val(lp);
	for (int i = 0; i < numX; i++) {
		x[i] = glp_get_col_prim(lp, 1+i);
	}
//s33:  x1 = glp_get_col_prim(lp, 1);
//s34:  x2 = glp_get_col_prim(lp, 2);
//s35:  x3 = glp_get_col_prim(lp, 3);
//s36:  printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g\n",
//	z, x1, x2, x3);
  s37:  glp_delete_prob(lp);
	  return 0;


}

/* eof */
