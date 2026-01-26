#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"
#include "KleinSampler.h"
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fmpz.h>

#define INDEX_2D(i, j, cols) ((i) * (cols) + (j))
//clear && gcc GaussSieve.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o
typedef struct system_struct *System;
struct system_struct
{
	size_t dimension;
	size_t totalDatabaseCount;//arraySize
	size_t totalToReduceCount;//stackSize
	size_t toReduceIndex;//indexS
	size_t databaseIndex;//indexL
	size_t collisionCount;
	size_t collisionLimit;
	fmpz_t *basis;
	mpfr_t *basis_mpfr;
	fmpz_t *databasePreImage;
	fmpz_t *database;
	fmpz_t *databaseSquareRootNorm;
	fmpz_t *toReducePreImage;
	fmpz_t *toReduce;
	fmpz_t *vectorPreImage;
	fmpz_t *vector;
	fmpz_t vNorm;
};

System CreateSystem(size_t dimension, size_t totalDatabaseCount, size_t totalToReduceCount, char *matrixString[], int stringBase)
{
	System system = malloc(sizeof(struct system_struct));
	assert(totalToReduceCount >= dimension);
	
	system->dimension = dimension;
	system->totalDatabaseCount = totalDatabaseCount;
	system->totalToReduceCount = totalToReduceCount;
	system->toReduceIndex = dimension;//Start with all basis in toReduce
	system->databaseIndex = 1; //Include  0 vector
	system->collisionLimit = 0.1 * system->databaseIndex + 200;
	system->collisionCount = 0;
	
	system->basis = malloc(system->dimension * system->dimension * sizeof(fmpz_t));
	system->basis_mpfr = malloc(system->dimension * system->dimension * sizeof(mpfr_t));
	
	system->databasePreImage = malloc(system->dimension * system->totalDatabaseCount * sizeof(fmpz_t));
	system->database = malloc(system->dimension * system->totalDatabaseCount * sizeof(fmpz_t));
	system->databaseSquareRootNorm = malloc(system->dimension * system->totalDatabaseCount * sizeof(fmpz_t));
	
	system->toReducePreImage = malloc(system->dimension * system->totalToReduceCount * sizeof(fmpz_t));
	system->toReduce = malloc(system->dimension * system->totalToReduceCount * sizeof(fmpz_t));
	
	system->vectorPreImage = malloc(system->dimension * sizeof(fmpz_t));
	system->vector = malloc(system->dimension * sizeof(fmpz_t));
	
	for(size_t i = 0; i < system->dimension * system->totalDatabaseCount; i++)
	{
		fmpz_init(system->databasePreImage[i]);
		fmpz_init(system->database[i]);
		fmpz_init(system->databaseSquareRootNorm[i]);
	}
	for(size_t i = 0; i < system->dimension * system->totalToReduceCount; i++)
	{
		fmpz_init(system->toReducePreImage[i]);
		fmpz_init(system->toReduce[i]);
	}
	
	//Load basis and toReducePreimage
	for(size_t i = 0; i < system->dimension; i++)
	{
		//Init vector and preimage
		fmpz_init(system->vector[i]);
		fmpz_init(system->vectorPreImage[i]);
		for(size_t j = 0; j < system->dimension; j++)
		{
			//Load Basis
			fmpz_init(system->basis[INDEX_2D(i,j,system->dimension)]);
			mpfr_init2(system->basis_mpfr[INDEX_2D(i,j,system->dimension)], MPFR_PRECISION);
			fmpz_set_str(system->basis[INDEX_2D(i,j,system->dimension)], matrixString[INDEX_2D(i,j,system->dimension)], stringBase);
			mpfr_set_str(system->basis_mpfr[INDEX_2D(i,j,system->dimension)], matrixString[INDEX_2D(i,j,system->dimension)], stringBase, MPFR_PRECISION);			
			//Initialize toReduce with basis
			fmpz_set(system->toReduce[INDEX_2D(i,j,system->dimension)], system->basis[INDEX_2D(i,j,system->dimension)]);		
			//Set toReducePreimage
			fmpz_set_ui(system->toReducePreImage[INDEX_2D(i,j,system->dimension)], (j == i) ? 1 : 0);
		}
	}
	fmpz_init(system->vNorm);
	return system;
}

void FindSquareDotProduct(fmpz_t result, size_t length, fmpz_t *a, fmpz_t *b)
{
	fmpz_t temp0;
	fmpz_init(temp0);
	fmpz_set_ui(result, 0);
	for(size_t i = 0; i < length; i++)
	{
		fmpz_mul(temp0, a[i], b[i]);
		fmpz_add(result, result, temp0);
	}
	fmpz_clear(temp0);
}
void LoadLastVectorInToReduce(System system)
{
	//Remove vector at bottom of toReduceMatrix
	assert(system->toReduceIndex > 0);
	assert(system->toReduceIndex < system->totalToReduceCount);
	system->toReduceIndex -= 1;
	for(size_t i = 0; i < system->dimension; i++)
	{
		fmpz_set(system->vector[i], system->toReduce[INDEX_2D(system->toReduceIndex,i,system->dimension)]);
		fmpz_set(system->vectorPreImage[i], system->toReducePreImage[INDEX_2D(system->toReduceIndex,i,system->dimension)]);
		//fmpz_print(system->toReduce[INDEX_2D(system->toReduceIndex,i,system->dimension)]); printf(",");
	}

}

void DestroySystem(System system)
{
	if(system)
	{
		for(size_t i = 0; i < system->dimension * system->dimension; i++)
		{
			fmpz_clear(system->basis[i]);
			mpfr_clear(system->basis_mpfr[i]);
		}
		for(size_t i = 0; i < system->dimension * system->totalDatabaseCount; i++)
		{
			fmpz_clear(system->databasePreImage[i]);
			fmpz_clear(system->database[i]);
			fmpz_clear(system->databaseSquareRootNorm[i]);
		}
		for(size_t i = 0; i < system->dimension * system->totalToReduceCount; i++)
		{
			fmpz_clear(system->toReducePreImage[i]);
			fmpz_clear(system->toReduce[i]);
		}
		for(size_t i = 0; i < system->dimension; i++)
		{
			fmpz_clear(system->vectorPreImage[i]);
			fmpz_clear(system->vector[i]);
		}
		fmpz_clear(system->vNorm);
		free(system->basis);
		free(system->basis_mpfr);
		free(system->databasePreImage);
		free(system->database);
		free(system->databaseSquareRootNorm);
		free(system->toReducePreImage);
		free(system->toReduce);
		free(system->vectorPreImage);
		free(system->vector);
		free(system);
	}
}

void PrintSystem(System system)
{
	printf("CollisionLimit: %ld\nDatabaseIndex: %ld\n", system->collisionLimit, system->databaseIndex);
	for(size_t i = 0; i < system->dimension; i++)
	{
		for(size_t j = 0; j < system->dimension; j++)
		{
			fmpz_print(system->basis[INDEX_2D(i,j,system->dimension)]);printf(",");
		}
		printf("\n");
	}
}

void Simpledaxpy(size_t dimension, fmpz_t alpha, fmpz_t *x, fmpz_t *y)
{
	//(alpha * X[i]) + Y[i]
	for(size_t i = 0; i < dimension; i++)
	{
		fmpz_addmul(y[i], alpha, x[i]);
	}
}

bool ReduceVectors(size_t dimension, fmpz_t *vPreimage, fmpz_t *v, fmpz_t vNorm, fmpz_t *wPreimage, fmpz_t *w, fmpz_t wNorm, fmpz_t vw)		
{
	fmpz_t vwAbsDouble,quotient,negQuotient, temp0;
	fmpz_init(vwAbsDouble);fmpz_init(quotient);fmpz_init(negQuotient);fmpz_init(temp0);
	
	fmpz_mul_ui(vwAbsDouble, vw, 2);			
	fmpz_abs(vwAbsDouble, vwAbsDouble);
	
	//Reduce v with w if possible
	//printf("vwAbsDouble:");fmpz_print(vwAbsDouble);printf("\n");
	//printf("wNorm:");fmpz_print(wNorm);printf("\n");
	if(fmpz_cmp(vwAbsDouble, wNorm) > 0)
	{
		fmpz_ndiv_qr(quotient,temp0, vw, wNorm);
		fmpz_neg(negQuotient,quotient);
		//printf("Starting reduction q:");fmpz_print(quotient);printf("\n");
		//printf("vw:");fmpz_print(vw);printf("\n");
		//printf("wNorm:");fmpz_print(wNorm);printf("\n");
		
		Simpledaxpy(dimension, negQuotient, w, v);
		Simpledaxpy(dimension, negQuotient, wPreimage, vPreimage);
		//printf("daxpy v )" );for(size_t i = 0; i < dimension; i++){fmpz_print(v[i]);printf(",");}printf("\n");
		//printf("daxpy vPreimage)");for(size_t i = 0; i < dimension; i++){fmpz_print(vPreimage[i]);printf(",");}printf("\n");		
		//Update vNorm
		fmpz_mul(temp0, quotient, quotient);
		fmpz_mul(temp0, temp0, wNorm);
		fmpz_add(vNorm, vNorm, temp0);
		
		fmpz_mul(temp0, quotient, vw);
		fmpz_mul_ui(temp0, temp0, 2);
		fmpz_sub(vNorm, vNorm, temp0);
		//printf("Return true :");fmpz_print(vwAbsDouble);printf(" ");fmpz_print(wNorm);printf("\n");
		return true;
	}
	//printf("Return false :");fmpz_print(vwAbsDouble);printf(" ");fmpz_print(wNorm);printf("\n");

	fmpz_clear(vwAbsDouble);fmpz_clear(quotient);fmpz_clear(negQuotient);fmpz_clear(temp0);
	return false;
}

void AppendVectorToReduce(System system, fmpz_t *preImage, fmpz_t *vector)
{
	if(system->toReduceIndex < system->totalToReduceCount)
	{
		for(size_t i = 0; i < system->dimension; i++)
		{
			fmpz_set(system->toReduce[INDEX_2D(system->toReduceIndex,i,system->dimension)], vector[i]);
			fmpz_set(system->toReducePreImage[INDEX_2D(system->toReduceIndex,i,system->dimension)],preImage[i]);
			//fmpz_print(system->toReduce[INDEX_2D(system->toReduceIndex,i,system->dimension)]); printf(",");
		}
		system->toReduceIndex += 1;
	}
	else
	{
		printf("Can't AppendVectorToReduce\n");
	}
}

void RemoveFromDatabase(System system, int index)
{
	assert(index > -1);assert(index < system->totalDatabaseCount);
	assert(system->databaseIndex > 0);assert(system->databaseIndex < system->totalDatabaseCount);	
	//Remove element at index i from db (Replace with last in array if not last element)

	system->databaseIndex -= 1;
	for(size_t i = 0; i < system->dimension; i++)
	{
		fmpz_set(system->database[INDEX_2D(index,i,system->dimension)], system->database[INDEX_2D(system->databaseIndex,i,system->dimension)]);
		fmpz_set(system->databasePreImage[INDEX_2D(index,i,system->dimension)], system->databasePreImage[INDEX_2D(system->databaseIndex,i,system->dimension)]);
		//fmpz_print(system->toReduce[INDEX_2D(system->toReduceIndex,i,system->dimension)]); printf(",");
	}
	fmpz_set(system->databaseSquareRootNorm[index], system->databaseSquareRootNorm[system->databaseIndex]);
	
}

void AddToDatabase(System system, fmpz_t *preImage, fmpz_t *vector, fmpz_t vNorm)
{
	int destinationIndex = -1;
	if(system->databaseIndex < system->totalDatabaseCount)
	{
		destinationIndex = system->databaseIndex;
		system->databaseIndex += 1;
	}
	else
	{
		//If array has run out of space purge largest vector
		fmpz_t max;
		fmpz_init(max);
		fmpz_set_si(max, -1);
		for(int i = 0; i < system->databaseIndex; i++)
		{
			if(fmpz_cmp(system->databaseSquareRootNorm[i], max))
			{
				fmpz_set(max, system->databaseSquareRootNorm[i]);
				destinationIndex = i;
			}
		}
		fmpz_clear(max);
	}
	
	//Copy into designated index
	for(size_t i = 0; i < system->dimension; i++)
	{
		fmpz_set(system->database[INDEX_2D(destinationIndex,i,system->dimension)], vector[i]);
		fmpz_set(system->databasePreImage[INDEX_2D(destinationIndex,i,system->dimension)], preImage[i]);
	}
	fmpz_set(system->databaseSquareRootNorm[destinationIndex], vNorm);
}

void GetImageFromPreImage_Basis(size_t dimension, fmpz_t *image, fmpz_t *basis, fmpz_t *preImage)
{       //[dimension,dimension] * [dimension, 1] = [dimension,1]
	//image = Basis matrix * preImage vector 
	//Idk why but matmul is in row major
	for(size_t i = 0; i < dimension; i++)
	{
		fmpz_set_ui(image[i], 0);
		for(size_t j = 0; j < dimension; j++)
		{
			fmpz_addmul(image[i], basis[j*dimension + i], preImage[j]);
		}
	}
}

void TestSystem10()
{
	size_t dimension = 10;
	size_t totalDatabaseCount = 30000;
	size_t totalToReduceCount = 10000;
	int stringBase = 10;
	char *matrixString[] = {
        "-30", "-277", "-297", "269", "384", "20", "409", "-328", "66", "-232",
        "-86", "-313", "25", "244", "-451", "-299", "289", "-31", "-163", "-533",
        "496", "129", "-203", "348", "-248", "-354", "-453", "-415", "-114", "93",
        "-8", "-317", "-431", "-859", "-494", "-217", "-50", "231", "-105", "9",
        "196", "-146", "233", "604", "-57", "-53", "42", "115", "839", "525",
        "-185", "610", "244", "-57", "-200", "329", "653", "-37", "-625", "347",
        "65", "183", "-253", "114", "-214", "712", "-644", "490", "243", "-372",
        "-1036", "-185", "-173", "-213", "144", "-543", "-143", "-52", "396", "-179",
        "275", "89", "447", "-666", "283", "-7", "50", "-282", "632", "-559",
        "348", "430", "-871", "235", "-77", "-570", "453", "337", "-243", "55"
	};
	size_t matrixStringLength = sizeof(matrixString) / sizeof(matrixString[0]);
	assert(matrixStringLength == dimension * dimension);
	System system = CreateSystem(dimension, totalDatabaseCount, totalToReduceCount, matrixString, stringBase);
	PrintSystem(system);
	
	KleinSampler sampler = CreateSampler(system->basis_mpfr, dimension, dimension);
	int counter = 0;
	while(system->collisionCount < system->collisionLimit)
	{
		if(system->toReduceIndex > 0)
		{
			LoadLastVectorInToReduce(system);	
			//printf("\nHere (%d) Murage (%ld %ld %ld)\n",counter,system->collisionCount, system->toReduceIndex, system->databaseIndex);
			//printf("Popped vector : %ld)", system->toReduceIndex );for(size_t i = 0; i < system->dimension; i++){fmpz_print(system->vector[i]);printf(",");}printf("\n");
		}
		else
		{
			//Randomly sample a vector  (returns preimage of lattice point)
			Sample(sampler, system->vectorPreImage);
			//Convert from preimage to real coordinates (v = B * vPreimage)
			GetImageFromPreImage_Basis(system->dimension, system->vector, system->basis, system->vectorPreImage);
			//printf("\nNo, here (%d) Murage (%ld %ld %ld)\n",counter,system->collisionCount, system->toReduceIndex, system->databaseIndex);
			//printf("Sampled vectorPreImage : %ld)", system->toReduceIndex );for(size_t i = 0; i < system->dimension; i++){fmpz_print(system->vectorPreImage[i]);printf(",");}printf("\n");
			//printf("Sampled v : %ld)", system->toReduceIndex );for(size_t i = 0; i < system->dimension; i++){fmpz_print(system->vector[i]);printf(",");}printf("\n");
			//exit(1);
		}
		//printf("v: (%ld %ld)", system->collisionCount, system->toReduceIndex );for(size_t i = 0; i < system->dimension; i++){fmpz_print(system->vector[i]);printf(",");}printf("\n");
		
		//Find squared norm of vector
		FindSquareDotProduct(system->vNorm, system->dimension, system->vector, system->vector);
		//printf("Dot:");fmpz_print(system->vNorm);printf("\n");
		
		//Reduce vector
		bool reduced = false;
		for(size_t i = 0; i < system->databaseIndex; i++)
		{
			fmpz_t vwSquareNorm,wNorm;
			fmpz_init(vwSquareNorm);fmpz_init(wNorm);
			
			fmpz_t *w = &system->database[i * dimension];
			fmpz_t *wPreimage = &system->databasePreImage[i * dimension];
			fmpz_set(wNorm, system->databaseSquareRootNorm[i]);
			
			FindSquareDotProduct(vwSquareNorm, system->dimension, system->vector, w);
			//printf("w : %ld)", system->toReduceIndex );for(size_t i = 0; i < system->dimension; i++){fmpz_print(w[i]);printf(",");}printf("\n");
			//printf("{index (%ld), vw: ", i);fmpz_print(vwSquareNorm);printf("\n");
			//printf("vwSquareNorm:");fmpz_print(vwSquareNorm);printf("\n");
			//Attempt to reduce v using w, recording if at least one reduction took place		
			reduced |= ReduceVectors(system->dimension, system->vectorPreImage, system->vector, system->vNorm, wPreimage, w, wNorm, vwSquareNorm);
			//printf("Current(%d) reduced: %s\n\n", counter, reduced ? "true":"false");	
			fmpz_clear(vwSquareNorm);fmpz_clear(wNorm);
		}
		
		//Check for a collision where v was not independent of vectors in database
		if(fmpz_cmp_ui(system->vNorm, 0) == 0)
		{
			system->collisionCount += 1;
			continue;
		}
		
		//If vector was reduced then add to stack of vectors to be reduced on next round
		if(reduced)
		{
			AppendVectorToReduce(system,system->vectorPreImage, system->vector);
			if(counter == 10)
			{
				//printf("To append : %ld)", system->toReduceIndex );for(size_t i = 0; i < system->dimension; i++){fmpz_print(system->vector[i]);printf(",");}printf("\n");
				//printf("\nIndexS: (%ld)\n", system->toReduceIndex);
				//break;
			}
			continue;
		}
		if(counter == 10)
		{
			//printf("\nIndexS: (%ld)\n", system->toReduceIndex);
			//break;
		}
		// v is now fully reduced so loop through db and reduce all vectors using v
		for(int i = 0; i < system->databaseIndex; i++)
		{
			fmpz_t vwSquareNorm,wNorm;
			fmpz_init(vwSquareNorm);fmpz_init(wNorm);
			
			fmpz_t *w = &system->database[i * dimension];
			fmpz_t *wPreimage = &system->databasePreImage[i * dimension];
			fmpz_set(wNorm, system->databaseSquareRootNorm[i]);
			
			FindSquareDotProduct(vwSquareNorm, system->dimension, system->vector, w);
			//printf("vwSquareNorm:");fmpz_print(vwSquareNorm);printf("\n");
			if(ReduceVectors(system->dimension, wPreimage, w, wNorm, system->vectorPreImage, system->vector, system->vNorm, vwSquareNorm))
			{
				AppendVectorToReduce(system,wPreimage, w);
				RemoveFromDatabase(system, i);
				i -= 1;
				assert(i > -1);
				if(i < 0){break;}
			}

			fmpz_clear(vwSquareNorm);fmpz_clear(wNorm);
		}	
		AddToDatabase(system, system->vectorPreImage, system->vector, system->vNorm);	
		
		//Update collision limit based on how many vectors are in database
		system->collisionLimit = 0.2 * system->databaseIndex + 1000;
		//printf("Counter(%d) : %ld/%ld\n", counter, system->collisionCount, system->collisionLimit);
		//if(counter == 10){break;}
		counter += 1;
		
	}
	
	//Print vectors in database
	for(size_t i = 0; i < system->databaseIndex; i++)
	{
		printf("DB(%ld):",i);fmpz_print(system->databaseSquareRootNorm[i]);printf("\n");	
		//Print preimage
		printf("PreImage: ");
		for(size_t j = 0; j < system->dimension; j++)
		{
			fmpz_print(system->databasePreImage[i * system->dimension + j]);printf(",");	
		}
		printf("\n");
		printf("Lattice point: ");
		for(size_t j = 0; j < system->dimension; j++)
		{
			fmpz_print(system->database[i * system->dimension + j]);printf(",");	
		}
		printf("\n\n");
		
	}
	
	DestroySampler(sampler);
	DestroySystem(system);
}

int main()
{
	TestSystem10();
	flint_cleanup();
	return 0;
}
