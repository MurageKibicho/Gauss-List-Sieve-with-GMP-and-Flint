#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#define MPFR_PRECISION 512
#define PRINT_MPFR(name, mpfr_num){char buffer[1000];mpfr_sprintf(buffer, "%s: %.50Rf", name, mpfr_num); printf("%s\n", buffer);}
typedef struct klein_sampler_struct *KleinSampler;
struct klein_sampler_struct
{
	size_t rows;
	size_t cols;
	mpfr_t t_;
	mpfr_t *mu_;
	mpfr_t *coef_;
	mpfr_t *sPrimeSquare_;
	gmp_randstate_t rngState;
};

void MPFR_FindSquareDotProduct(mpfr_t result, size_t length, mpfr_t *a, mpfr_t *b)
{
	mpfr_t temp0;
	mpfr_init2(temp0, MPFR_PRECISION);mpfr_set_ui(temp0, 0, MPFR_RNDN );
	mpfr_set_ui(result, 0, MPFR_RNDN);
	for(size_t i = 0; i < length; i++)
	{
		mpfr_mul(temp0, a[i], b[i], MPFR_RNDN);
		mpfr_add(result, result, temp0, MPFR_RNDN);
	}
	mpfr_clear(temp0);
}

void ComputeGramSchmidt(size_t rows, size_t cols, mpfr_t *B, mpfr_t *mu, mpfr_t *c)
{
	mpfr_t bound, bound2,s,t,t1,temp0;
	mpfr_init2(bound, MPFR_PRECISION);mpfr_set_ui(bound, 0, MPFR_RNDN );
	mpfr_init2(bound2, MPFR_PRECISION);mpfr_set_ui(bound2, 0, MPFR_RNDN );
	mpfr_init2(s, MPFR_PRECISION);mpfr_set_ui(s, 0, MPFR_RNDN );
	mpfr_init2(t, MPFR_PRECISION);mpfr_set_ui(t, 0, MPFR_RNDN );
	mpfr_init2(t1, MPFR_PRECISION);mpfr_set_ui(t1, 0, MPFR_RNDN );
	mpfr_init2(temp0, MPFR_PRECISION);mpfr_set_ui(temp0, 0, MPFR_RNDN );
		
	mpfr_t *B1 = malloc(rows * cols * sizeof(mpfr_t));
	mpfr_t *b  = malloc(cols * sizeof(mpfr_t));
	mpfr_t *buf= malloc(cols * sizeof(mpfr_t));	
	
	//Copy B into B1 and init B1
	for(size_t i = 0; i < cols; i++)
	{
		for(size_t j = 0; j < rows; j++)
		{
			mpfr_init2(B1[i * rows + j], MPFR_PRECISION);	
			mpfr_set(B1[i * rows + j], B[i * rows + j], MPFR_RNDN);				
		}	
	}
	
	//Init and set b to square norm of B1
	for(size_t i = 0; i < cols; i++)
	{
		mpfr_init2(b[i], MPFR_PRECISION);mpfr_set_ui(b[i], 0, MPFR_RNDN );
		mpfr_init2(buf[i], MPFR_PRECISION);mpfr_set_ui(buf[i], 0, MPFR_RNDN );
		MPFR_FindSquareDotProduct(b[i], rows, &B1[i * rows], &B1[i * rows]);		
	}
	
	//Set precision
	double precision = 150;
	mpfr_set_ui(bound, 2, MPFR_RNDN);
	mpfr_set_ui(bound2, 2, MPFR_RNDN); 
	mpfr_pow_ui(bound, bound, 2 * (long)(0.15 * precision), MPFR_RNDN);
	mpfr_pow_ui(bound2, bound2, 2 * precision, MPFR_RNDN);
	
	for(int k = 0; k < cols; k++)
	{
		if(k > 0)
		{
			mpfr_mul(buf[0], mu[k * rows], c[0], MPFR_RNDN);
		}
		for(int j = 0; j < k; j++)
		{
			MPFR_FindSquareDotProduct(s, rows, &B1[k * rows], &B1[j * rows]);		
			mpfr_mul(t1, s, s, MPFR_RNDN);
			mpfr_mul(t1, t1, bound, MPFR_RNDN);
			
			mpfr_mul(t, b[k], b[j], MPFR_RNDN);
			
			if(mpfr_cmp(t, bound2) >= 0 && mpfr_cmp(t, t1) >= 0)
			{
				MPFR_FindSquareDotProduct(s, rows, &B[k * rows], &B[j * rows]);				
			}
			mpfr_set_ui(t1, 0, MPFR_RNDN);	
			for(int i = 0; i < j; i++)
      			{
 				mpfr_mul(t, mu[j * rows + i], buf[i], MPFR_RNDN);     			
      				mpfr_add(t1,t1,t,MPFR_RNDN);
      			}
      			mpfr_sub(t,s,t1,MPFR_RNDN);
      			mpfr_set(buf[j], t, MPFR_RNDN);
      			mpfr_div(mu[k * rows + j], t, c[j],MPFR_RNDN);		
		}
		mpfr_set_ui(s, 0, MPFR_RNDN);
		for(int j = 0; j < k; j++)
		{
			mpfr_mul(temp0,mu[k * rows + j],buf[j], MPFR_RNDN);	
			mpfr_add(s,s,temp0,MPFR_RNDN);
		}
		mpfr_sub(c[k],b[k], s,MPFR_RNDN);
	}
	
	for(size_t i = 0; i < cols; i++){mpfr_clear(b[i]);mpfr_clear(buf[i]);}
	for(size_t i = 0; i < rows * cols; i++){mpfr_clear(B1[i]);}
	mpfr_clear(bound);
	mpfr_clear(bound2);
	mpfr_clear(s);
	mpfr_clear(t);
	mpfr_clear(t1);
	mpfr_clear(temp0);
	free(B1);
	free(b);
	free(buf);
}

KleinSampler CreateSampler(mpfr_t *BIn, size_t rows, size_t cols)
{
	KleinSampler sampler = malloc(sizeof(struct klein_sampler_struct));
	gmp_randinit_default(sampler->rngState);
	sampler->rows = rows;
	sampler->cols = cols;
	mpfr_init2(sampler->t_, MPFR_PRECISION);mpfr_set_ui(sampler->t_, 0, MPFR_RNDN );
	sampler->mu_ = malloc(sampler->rows * sampler->cols * sizeof(mpfr_t));
	sampler->coef_ = malloc(sampler->rows * sizeof(mpfr_t));	
	sampler->sPrimeSquare_ = malloc(sampler->rows * sizeof(mpfr_t));	
	mpfr_t *bStarSquare = malloc(sampler->rows * sizeof(mpfr_t));	
	
	
	for(size_t i = 0; i < sampler->rows; i++)
	{
		for(size_t j = 0; j < sampler->cols; j++)
		{
			mpfr_init2(sampler->mu_[i * cols + j], MPFR_PRECISION);	
			mpfr_set_ui(sampler->mu_[i * cols + j], 0, MPFR_RNDN);
		}

		mpfr_init2(sampler->coef_[i], MPFR_PRECISION);
		mpfr_set_ui(sampler->coef_[i], 0, MPFR_RNDN );
		mpfr_init2(sampler->sPrimeSquare_[i], MPFR_PRECISION);
		mpfr_set_ui(sampler->sPrimeSquare_[i], 0, MPFR_RNDN );
		mpfr_init2(bStarSquare[i], MPFR_PRECISION);
		mpfr_set_ui(bStarSquare[i], 0, MPFR_RNDN );
	}
	
	ComputeGramSchmidt(rows, cols, BIn, sampler->mu_, bStarSquare);
	
	mpfr_t maxStarSqrNorm,sSquare;
	mpfr_init2(maxStarSqrNorm, MPFR_PRECISION);mpfr_set_ui(maxStarSqrNorm, 0, MPFR_RNDN );
	mpfr_init2(sSquare, MPFR_PRECISION);mpfr_set_ui(sSquare, 0, MPFR_RNDN );
	for(int i = 0; i < sampler->rows; i++)
	{
		if(mpfr_cmp(bStarSquare[i], maxStarSqrNorm) > 0)
		{
			mpfr_set(maxStarSqrNorm, bStarSquare[i], MPFR_RNDN);
		}
		//PRINT_MPFR("mu[i]", sampler->mu_[i]);	
		//PRINT_MPFR("bStarSquare[i]", bStarSquare[i]);	
	}
	//PRINT_MPFR("maxStarSqrNorm", maxStarSqrNorm);	
	mpfr_set_d(sampler->t_, log(sampler->rows), MPFR_RNDN);
	mpfr_mul(sSquare, maxStarSqrNorm, sampler->t_, MPFR_RNDN);
	for(int i = 0; i < sampler->rows; i++)
	{
		mpfr_div(sampler->sPrimeSquare_[i], sSquare, bStarSquare[i], MPFR_RNDN); 	
		//PRINT_MPFR("i", sampler->sPrimeSquare_[i]);
	}
	
	
	mpfr_clear(sSquare);
	mpfr_clear(maxStarSqrNorm);
	for(size_t i = 0; i < sampler->rows; i++)
	{
		mpfr_clear(bStarSquare[i]);
	}
	free(bStarSquare);
	return sampler;
}

void SampleZ(KleinSampler sampler, mpfr_t result, mpfr_t c, mpfr_t sSquare)
{
	mpfr_t s, minC, maxC, rho, x, randomX, randomRho,temp1, temp2,randNum;
	mpfr_init2(s, MPFR_PRECISION);mpfr_set_ui(s, 0, MPFR_RNDN );
	mpfr_init2(minC, MPFR_PRECISION);mpfr_set_ui(minC, 0, MPFR_RNDN );
	mpfr_init2(maxC, MPFR_PRECISION);mpfr_set_ui(maxC, 0, MPFR_RNDN );
	mpfr_init2(rho, MPFR_PRECISION);mpfr_set_ui(rho, 0, MPFR_RNDN );
	mpfr_init2(x, MPFR_PRECISION);mpfr_set_ui(x, 0, MPFR_RNDN );
	mpfr_init2(randomX, MPFR_PRECISION);mpfr_set_ui(randomX, 0, MPFR_RNDN );
	mpfr_init2(randomRho, MPFR_PRECISION);mpfr_set_ui(randomRho, 0, MPFR_RNDN );
	mpfr_init2(temp1, MPFR_PRECISION);mpfr_set_ui(temp1, 0, MPFR_RNDN );
	mpfr_init2(temp2, MPFR_PRECISION);mpfr_set_ui(temp2, 0, MPFR_RNDN );
	mpfr_init2(randNum, MPFR_PRECISION);mpfr_set_ui(randNum, 0, MPFR_RNDN );
		
	mpfr_sqrt(s, sSquare, MPFR_RNDN);
	mpfr_mul(temp1, s, sampler->t_, MPFR_RNDN);
	mpfr_sub(temp2, c, temp1, MPFR_RNDN); 
	mpfr_floor(minC, temp2);  
	
	mpfr_add(temp2, c, temp1, MPFR_RNDN); 
	mpfr_ceil(maxC, temp2); 
	
	int count = 0;
	while(true)
	{
		mpfr_urandomb(randNum, sampler->rngState); 
		mpfr_sub(temp1, maxC, minC, MPFR_RNDN);
		mpfr_mul(temp1, temp1, randNum, MPFR_RNDN); 
		mpfr_add(x, minC, temp1, MPFR_RNDN);
		//Cast x to long
		mpfr_round(x, x); 
		//PRINT_MPFR("x", x);
		mpfr_sub(temp1, x, c, MPFR_RNDN); 
		mpfr_mul(temp1, temp1, temp1, MPFR_RNDN);
		mpfr_div(temp1, temp1, sSquare, MPFR_RNDN);
		mpfr_mul_d(temp1, temp1, M_PI, MPFR_RNDN); 
		mpfr_neg(temp1, temp1, MPFR_RNDN);
		mpfr_exp(rho, temp1, MPFR_RNDN);
		
		mpfr_urandomb(randomRho, sampler->rngState); 
		if(mpfr_cmp(randomRho, rho) <= 0)
		{
			mpfr_set(result, x, MPFR_RNDN);
			break;
		}
		count += 1;
		if(count > 1000)
		{
			printf("Probable infinite loop in KleinSampler SampleZ function");
			exit(1);
		}
	}
	
	mpfr_clear(s);
	mpfr_clear(minC);
	mpfr_clear(maxC);
	mpfr_clear(rho);
	mpfr_clear(x);
	mpfr_clear(randomX);
	mpfr_clear(randomRho);
	mpfr_clear(temp1);
	mpfr_clear(temp2);
	mpfr_clear(randNum);
}

// Samples a vector in the preimage of the basis (Premultiply by basis to get lattice point)
void Sample(KleinSampler sampler, fmpz_t *v)
{
	mpfr_t temp0;
	mpz_t mpz_val;
	mpz_init(mpz_val);
	mpfr_init2(temp0, MPFR_PRECISION);mpfr_set_ui(temp0, 0, MPFR_RNDN );
	for(size_t i = 0; i < sampler->rows; i++)
	{
		mpfr_set_ui(sampler->coef_[i], 0, MPFR_RNDN);
	}
	for(int i = sampler->rows - 1; i >= 0; i--)
	{
		SampleZ(sampler, sampler->coef_[i], sampler->coef_[i], sampler->sPrimeSquare_[i]);		
		for(int j = 0; j < i; j++)
		{
			mpfr_mul(temp0, sampler->coef_[i], sampler->mu_[i * sampler->rows + j], MPFR_RNDN);
			mpfr_sub(sampler->coef_[j],sampler->coef_[j],temp0,MPFR_RNDN);
			//PRINT_MPFR("pre_mu[i]", sampler->mu_[i * sampler->rows + j]);
		}
	}
	for(size_t i = 0; i < sampler->rows; i++)
	{
		mpfr_get_z(mpz_val , sampler->coef_[i], MPFR_RNDN);
		fmpz_set_mpz(v[i], mpz_val);
		//PRINT_MPFR("pre[i]", v[i]);
	}
	mpfr_clear(temp0);
	mpz_clear(mpz_val);
}

void DestroySampler(KleinSampler sampler)
{
	if(sampler)
	{
		for(size_t i = 0; i < sampler->rows; i++)
		{
			for(size_t j = 0; j < sampler->cols; j++)
			{
				mpfr_clear(sampler->mu_[i * sampler->cols + j]);	
			}
			mpfr_clear(sampler->coef_[i]);
			mpfr_clear(sampler->sPrimeSquare_[i]);
		}
		gmp_randclear(sampler->rngState);
		mpfr_clear(sampler->t_);
		free(sampler->mu_);
		free(sampler->coef_);
		free(sampler->sPrimeSquare_);
		free(sampler);
	}
}


