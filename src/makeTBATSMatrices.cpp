#include "calcBATS.h"

using namespace Rcpp ;

SEXP makeTBATSWMatrix(SEXP smallPhi_s, SEXP kVector_s, SEXP arCoefs_s, SEXP maCoefs_s, SEXP tau_s) {
	BEGIN_RCPP
	double *smallPhi, *arCoefs, *maCoefs;
	int *kVector, *tau;
	int adjustPhi = 0;
	R_len_t numSeasonal = 0, numCols = 1, p = 0, q = 0;
	int lengthSeasonal = 0;
	//Rprintf("--1\n");
	if(!Rf_isNull(smallPhi_s)) {
		smallPhi = REAL(smallPhi_s);
		adjustPhi = 1;
		numCols = numCols + 1;
	}
	//Rprintf("--2\n");
	if(!Rf_isNull(kVector_s)) {
		//Rprintf("--2-1\n");
		tau = &INTEGER(tau_s)[0];
		//Rprintf("--2-2\n");
		kVector = INTEGER(kVector_s);
		//Rprintf("--2-3\n");
		numSeasonal = LENGTH(kVector_s);
		//Rprintf("--2-4\n");
		numCols = numCols + *tau;
	}
	//Rprintf("--3\n");
	if(!Rf_isNull(arCoefs_s)) {
			arCoefs = REAL(arCoefs_s);
			p = LENGTH(arCoefs_s);
			numCols = numCols + p;
	}
	//Rprintf("--4\n");
	if(!Rf_isNull(maCoefs_s)) {
			maCoefs = REAL(maCoefs_s);
			q = LENGTH(maCoefs_s);
			numCols = numCols + q;
	}
	//Rprintf("--5\n");

	NumericMatrix wTranspose_r(1, numCols);
	arma::mat wTranspose(wTranspose_r.begin(), wTranspose_r.nrow(), wTranspose_r.ncol(), false);

	//Rprintf("--6\n");
	if(!Rf_isNull(kVector_s)) {
			wTranspose.zeros();
			//Rprintf("--7\n");

			int position = adjustPhi;

			for(R_len_t s = 0; s < numSeasonal; s++) {
				//wTranspose.submat(0,(position+1), 0, (position + kVector[s])) = arma::ones<mat>(1,  kVector[s]);
				for(int j = (position+1); j <= (position + kVector[s]); j++) {
					wTranspose(0,j) = 1;
				}
				position = position + (2 * kVector[s]);

			}


	}


	//Rprintf("--8\n");
	wTranspose(0,0) = 1;

	if(adjustPhi == 1) {
		wTranspose(0,1) = *smallPhi;
	}

	//Rprintf("--9\n");

	if(!Rf_isNull(arCoefs_s)) {
		//Rprintf("--10\n");
		for(R_len_t i = 1; i <= p; i++) {
			wTranspose(0,(adjustPhi + *tau +i)) = arCoefs[(i-1)];
		}
	}

	if(!Rf_isNull(maCoefs_s)) {
		//Rprintf("--11\n");
			for(R_len_t i = 1; i <= q; i++) {
				wTranspose(0,(adjustPhi + *tau + p + i)) = maCoefs[(i-1)];
			}

	}
	///Rprintf("--12\n");
	arma::mat w = arma::trans(wTranspose);
	//Rprintf("--13\n");
	smallPhi = 0;
	arCoefs = 0;
	maCoefs = 0;
	kVector = 0;
	//Rprintf("--14\n");
	return List::create(
			Named("w") = w,
			Named("w.transpose") = wTranspose
			);


	END_RCPP
}

SEXP makeCIMatrix(SEXP k_s, SEXP m_s) {
	BEGIN_RCPP

	double lambda, *m;
	int *k;

	k = &INTEGER(k_s)[0];
	m = &REAL(m_s)[0];
	//std::cout<<arma::math::pi()<<std::endl;
	NumericMatrix C(*k, *k);
	for(int j = 1; j<=*k; j++) {
		lambda = (2 * arma::math::pi() * j) / *m;
		C((j-1),(j-1)) = std::cos(lambda);
		//std::cout<<lambda<<std::endl;
		//std::cout<<cos(lambda)<<std::endl;
		//std::cout<<std::cos(lambda)<<std::endl;
	}
	//std::cout<<C<<std::endl;
	return wrap(C);

	END_RCPP
}

SEXP makeSIMatrix(SEXP k_s, SEXP m_s) {
	BEGIN_RCPP

	double lambda, *m;
	int *k;
	k = &INTEGER(k_s)[0];
	m = &REAL(m_s)[0];

	NumericMatrix S(*k, *k);
	for(int j = 1; j<=*k; j++) {
		lambda = (2 * arma::math::pi() * j) / *m;
		S((j-1),(j-1)) = std::sin(lambda);
		//std::cout<<lambda<<std::endl;
	}
	//std::cout<<cos(lambda)<<std::endl;
	//std::cout<<std::cos(lambda)<<std::endl;
	//std::cout<<S<<std::endl;
	return wrap(S);

	END_RCPP
}


SEXP makeAIMatrix(SEXP C_s, SEXP S_s, SEXP k_s) {
	int *k;
	k = &INTEGER(k_s)[0];

	NumericMatrix C_r(C_s);
	NumericMatrix S_r(S_s);

	arma::mat C(C_r.begin(), C_r.nrow(), C_r.ncol(), false);
	arma::mat S(S_r.begin(), S_r.nrow(), S_r.ncol(), false);
	arma::mat A((*k * 2), (*k * 2));
	A.submat(0,0, (*k -1), (*k -1)) = C;
	A.submat(0,*k, (*k -1), ((*k *2) -1)) = S;
	A.submat(*k,0, ((*k *2) -1), (*k -1)) = (-1 * S);
	A.submat(*k,*k, ((*k *2) -1), ((*k *2) -1)) = C;

	return wrap(A);

}

