// Copyright(c) 2018 Vissarion Fisikopoulos, Apostolos Chalkis

//Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

#ifndef SOLVE_LP_H
#define SOLVE_LP_H

#include <stdio.h>
#include "lp_lib.h"


template <typename K>
std::pair<Point,double> solveLP(std::vector<std::vector<K> > A, int d){//std::vector<std::vector<K> > A, int d){

	//typedef typename T1::FT                    K;
	//int d=P.dimension();
	//std::vector<std::vector<K> > A = P.get_matrix();
	lprec *lp;
	int Ncol=d+1, *colno = NULL, j, ret = 0, m=A.size();
	REAL *row = NULL;
	
	lp = make_lp(m, Ncol);
	
	if(lp == NULL)
		ret = 1; /* couldn't construct a new model... */
    
    REAL infinite = get_infinite(lp); /* will return 1.0e30 */
    
    if(ret == 0) {
    /* let us name our variables. Not required, but can be useful for debugging */
		
		//for (int i=0; i<d; i++){
			//set_col_name(lp, i+1, "xi");
			//set_col_name(lp, 2, "y");
		//}

		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
		row = (REAL *) malloc(Ncol * sizeof(*row));
		if((colno == NULL) || (row == NULL))
			ret = 2;
	}
	
	
	if(ret == 0) {
		set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
	}
	int i=0;
	K sum;
	while(ret==0 & i<m){
		/* construct all rows (120 x + 210 y <= 15000) */
		sum=K(0);
		for(j=0; j<d; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = A[i][j+1];
			sum += A[i][j+1]*A[i][j+1];
		}
		colno[d] = d+1; /* last column */
		row[d] = std::sqrt(sum);
		//set_bounds(lp, d, 0.0, infinite);
		

    /* add the row to lpsolve */
		if(!add_constraintex(lp, d+1, row, colno, LE, A[i][0])){
			ret = 3;
		}
		i++;
	}
	
	if(ret == 0) {
		set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
		
		// set the objective function
		for(j=0; j<d; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = 0;
			set_bounds(lp, j+1, -infinite, infinite);
		}
		colno[d] = d+1; /* last column */
		row[d] = 1.0;
		set_bounds(lp, d+1, 0.0, infinite);
		if(!set_obj_fnex(lp, d+1, row, colno)){
			ret = 4;
		}
		
	}
	
	
	if(ret == 0) {
		/* set the object direction to maximize */
		set_maxim(lp);

		/* just out of curioucity, now show the model in lp format on screen */
		/* this only works if this is a console application. If not, use write_lp and a filename */
		//write_LP(lp, stdout);
		/* write_lp(lp, "model.lp"); */

		/* I only want to see important messages on screen while solving */
		set_verbose(lp, NEUTRAL);

		/* Now let lpsolve calculate a solution */
		ret = solve(lp);
		if(ret == OPTIMAL)
			ret = 0;
		else
			ret = 5;
	}
	
	//if(ret == 0) {
	std::vector<double> temp_p(d,0);
	get_variables(lp, row);
	for(j = 0; j < d; j++){
	//printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
		temp_p[j]=double(row[j]);
	}
	
	Point xc( d , temp_p.begin() , temp_p.end() );
	double r=double(get_objective(lp));
	delete_lp(lp);
			
	return std::pair<Point,double> (xc,r);
}




template <typename K>
bool memLP_Vpoly(std::vector<std::vector<K> > &A, Point q){

	//typedef typename T1::FT                    K;
	int d=q.dimension();
	//std::vector<std::vector<K> > A = P.get_matrix();
	lprec *lp;
	int Ncol=d+1, *colno = NULL, j, ret = 0, m=A.size();
	m++;
	REAL *row = NULL;

	lp = make_lp(m, Ncol);

	if(lp == NULL)
		ret = 1; /* couldn't construct a new model... */

	REAL infinite = get_infinite(lp); /* will return 1.0e30 */

	if(ret == 0) {
		/* let us name our variables. Not required, but can be useful for debugging */

		//for (int i=0; i<d; i++){
		//set_col_name(lp, i+1, "xi");
		//set_col_name(lp, 2, "y");
		//}

		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
		row = (REAL *) malloc(Ncol * sizeof(*row));
		if((colno == NULL) || (row == NULL))
			ret = 2;
	}

	if(ret == 0) {
		set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
	}
	int i=0;
	//K sum;
	while(ret==0 & i<m-1){
		/* construct all rows (120 x + 210 y <= 15000) */
		//sum=K(0);
		for(j=0; j<d; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = A[i][j+1];
			//sum += A[i][j+1]*A[i][j+1];
		}
		colno[d] = d+1; /* last column */
		row[d] = -1.0;
		//set_bounds(lp, d, 0.0, infinite);


		/* add the row to lpsolve */
		if(!add_constraintex(lp, d+1, row, colno, LE, 0.0)){
			ret = 3;
		}
		i++;
	}
	if (ret==0){
		for(j=0; j<d; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = q[j];
			//sum += A[i][j+1]*A[i][j+1];
		}
		colno[d] = d+1; /* last column */
		row[d] = -1.0;
		//set_bounds(lp, d, 0.0, infinite);


		/* add the row to lpsolve */
		if(!add_constraintex(lp, d+1, row, colno, LE, 1.0)){
			ret = 3;
		}
	}

	//set the bounds
	if(ret == 0) {
		set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

		// set the objective function
		for(j=0; j<d; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = q[j];
			set_bounds(lp, j+1, -infinite, infinite);
		}
		colno[d] = d+1; /* last column */
		row[d] = -1.0;
		set_bounds(lp, d+1, -infinite, infinite);
		if(!set_obj_fnex(lp, d+1, row, colno)){
			ret = 4;
		}

	}

	if(ret == 0) {
		/* set the object direction to maximize */
		set_maxim(lp);

		/* just out of curioucity, now show the model in lp format on screen */
		/* this only works if this is a console application. If not, use write_lp and a filename */
		//write_LP(lp, stdout);
		/* write_lp(lp, "model.lp"); */

		/* I only want to see important messages on screen while solving */
		set_verbose(lp, NEUTRAL);

		/* Now let lpsolve calculate a solution */
		ret = solve(lp);
		if(ret == OPTIMAL)
			ret = 0;
		else
			ret = 5;
	}

	double r=double(get_objective(lp));
	delete_lp(lp);
	if(r>0.0){
		return false;
	}
	return true;


}


template <typename K>
double intersect_line_Vpoly(std::vector<std::vector<K> > &A,Point &p, Point &v, bool maxi){

	int d=v.dimension();
	//std::vector<std::vector<K> > A = P.get_matrix();
	lprec *lp;
	int m=A.size();
	m++;
	int Ncol=m, *colno = NULL, j, ret = 0;
	REAL *row = NULL;
	K res;

	lp = make_lp(d+1, Ncol);

	if(lp == NULL)
		ret = 1; /* couldn't construct a new model... */

	REAL infinite = get_infinite(lp); /* will return 1.0e30 */

	//std::cout<<"ret1 = "<<ret<<std::endl;
	if(ret == 0) {
		/* let us name our variables. Not required, but can be useful for debugging */

		//for (int i=0; i<d; i++){
		//set_col_name(lp, i+1, "xi");
		//set_col_name(lp, 2, "y");
		//}

		/* create space large enough for one row */
		colno = (int *) malloc(Ncol * sizeof(*colno));
		row = (REAL *) malloc(Ncol * sizeof(*row));
		if((colno == NULL) || (row == NULL))
			ret = 2;
	}
	//std::cout<<"ret2 = "<<ret<<std::endl;

	if(ret == 0) {
		set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
	}
	int i=0;
	//std::cout<<"ret3 = "<<ret<<std::endl;

	while(ret==0 & i<d){
		/* construct all rows (120 x + 210 y <= 15000) */
		//sum=K(0);
		for(j=0; j<m-1; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = A[j][i+1];
			//sum += A[i][j+1]*A[i][j+1];
		}
		colno[m-1] = m; /* last column */
		row[m-1] = v[i];
		//set_bounds(lp, d, 0.0, infinite);


		/* add the row to lpsolve */
		if(!add_constraintex(lp, m, row, colno, EQ, p[i])){
			ret = 3;
		}
		i++;
	}
	//std::cout<<"ret4 = "<<ret<<std::endl;
	if(ret==0){
		for(j=0; j<m-1; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = 1.0;
			//sum += A[i][j+1]*A[i][j+1];
		}
		colno[m-1] = m; /* last column */
		row[m-1] = 0.0;
		//set_bounds(lp, d, 0.0, infinite);


		/* add the row to lpsolve */
		if(!add_constraintex(lp, m, row, colno, EQ, 1.0)){
			ret = 3;
		}
	}

	//set the bounds
	//std::cout<<"ret5 = "<<ret<<std::endl;
	if(ret == 0) {
		set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */

		// set the objective function
		for(j=0; j<m-1; j++){
			colno[j] = j+1; /* j_th column */
			row[j] = 0;
			//set_bounds(lp, j+1, 0.0, infinite);
		}
		colno[m - 1] =m; /* last column */
		row[m-1] = 1.0;
		set_bounds(lp, m, -infinite, infinite);
		if(!set_obj_fnex(lp, m, row, colno)){
			ret = 4;
		}

	}

	//std::cout<<"ret6 = "<<ret<<std::endl;
	if(ret == 0) {
		/* set the object direction to maximize */
		if(maxi) {
			set_maxim(lp);
		}else{
			set_minim(lp);
		}
		//set_minim(lp);

		/* just out of curioucity, now show the model in lp format on screen */
		/* this only works if this is a console application. If not, use write_lp and a filename */
		//write_LP(lp, stdout);
		/* write_lp(lp, "model.lp"); */

		/* I only want to see important messages on screen while solving */
		set_verbose(lp, NEUTRAL);

		/* Now let lpsolve calculate a solution */
		ret = solve(lp);
		if(ret == OPTIMAL)
			ret = 0;
		else
			ret = 5;
	}
	//std::cout<<"ret7 = "<<ret<<std::endl;
	if(ret==0) {
		res = K(-get_objective(lp));
		delete_lp(lp);
	    return res;

	}
	std::cout<<K(-get_objective(lp))<<std::endl;
	exit(-1);
	return 0;

}

#endif