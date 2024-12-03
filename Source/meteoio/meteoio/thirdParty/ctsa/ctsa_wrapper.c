// SPDX-License-Identifier: BSD-3-Clause
#include "ctsa.h"

auto_arima_object auto_arima_copy(auto_arima_object original) {
	// Allocate a new auto_arima_set on the heap
	int pqdmax[] = {original->pmax, original->dmax, original->qmax};
	int PQDmax[] = {original->Pmax, original->Dmax, original->Qmax};
	auto_arima_object copy = auto_arima_init(pqdmax, PQDmax, original->s, original->r, original->N);

	// Copy the values from the original to the copy
	copy->N = original->N;
	copy->Nused = original->Nused;
	copy->method = original->method;
	copy->optmethod = original->optmethod;
	copy->pmax = original->pmax;
	copy->dmax = original->dmax;
	copy->qmax = original->qmax;
	copy->Pmax = original->Pmax;
	copy->Dmax = original->Dmax;
	copy->Qmax = original->Qmax;
	copy->p = original->p;
	copy->d = original->d;
	copy->q = original->q;
	copy->s = original->s;
	copy->P = original->P;
	copy->D = original->D;
	copy->Q = original->Q;
	copy->r = original->r;
	copy->M = original->M;
	copy->ncoeff = original->ncoeff;
	copy->lvcov = original->lvcov;
	copy->mean = original->mean;
	copy->var = original->var;
	copy->loglik = original->loglik;
	copy->ic = original->ic;
	copy->retval = original->retval;
	copy->start = original->start;
	copy->imean = original->imean;
	copy->idrift = original->idrift;
	copy->stationary = original->stationary;
	copy->seasonal = original->seasonal;
	copy->Order_max = original->Order_max;
	copy->p_start = original->p_start;
	copy->q_start = original->q_start;
	copy->P_start = original->P_start;
	copy->Q_start = original->Q_start;
	copy->stepwise = original->stepwise;
	copy->num_models = original->num_models;
	copy->approximation = original->approximation;
	copy->verbose = original->verbose;
	copy->alpha_test = original->alpha_test;
	copy->alpha_seas = original->alpha_seas;
	copy->lambda = original->lambda;
	copy->sigma2 = original->sigma2;
	copy->aic = original->aic;
	copy->bic = original->bic;
	copy->aicc = original->aicc;

	memcpy(copy->information_criteria, original->information_criteria, 10 * sizeof(int));
	memcpy(copy->test, original->test, 10 * sizeof(int));
	memcpy(copy->type, original->type, 10 * sizeof(int));
	memcpy(copy->seas, original->seas, 10 * sizeof(int));
	int res_size;
	if (copy->P == 0 && copy->D == 0 && copy->Q==0) {
		res_size = copy->N-copy->d;
	} else {
		res_size = copy->N-copy->D*copy->s-copy->d;
	}
	memcpy(copy->res, original->res, res_size * sizeof(double));

	// If there are any members that are pointers (like phi, theta, etc.), you'll need to allocate new memory for them and copy the values
	copy->phi = (double*)malloc(original->p * sizeof(double));
	memcpy(copy->phi, original->phi, original->p * sizeof(double));
	copy->theta = (double*)malloc(original->q * sizeof(double));
	memcpy(copy->theta, original->theta, original->q * sizeof(double));
	copy->PHI = (double*)malloc(original->P * sizeof(double));
	memcpy(copy->PHI, original->PHI, original->P * sizeof(double));
	copy->THETA = (double*)malloc(original->Q * sizeof(double));
	memcpy(copy->THETA, original->THETA, original->Q * sizeof(double));
	copy->exog = (double*)malloc(original->M * sizeof(double));
	memcpy(copy->exog, original->exog, original->N*original->r * sizeof(double));
	copy->vcov = (double*)malloc(copy->lvcov * sizeof(double));

	int param_size;
	if (copy->r == 0) {
		if (copy->P == 0 && copy->D == 0 && copy->Q==0) {
			param_size = copy->p + copy->q + res_size + copy->lvcov;
		} else {
			param_size = copy->p + copy->q + copy->P + copy->Q + res_size + copy->lvcov;
		}
	} else {
		param_size = copy->p + copy->q + copy->P + copy->Q + copy->M + res_size + copy->lvcov; 
	}
	memcpy(copy->params, original->params, param_size * sizeof(double));

	return copy;
}
