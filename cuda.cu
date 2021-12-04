#include <stdio.h>
#include <string.h>
#include <cuda_runtime.h>
#include "cudaHeader.h"
#include <helper_cuda.h>

#define MAX_FIRST_TYPE_GROUPS_ROW_SIZE 9
#define MAX_FIRST_TYPE_GROUPS_COL_SIZE 5
#define MAX_SECOND_TYPE_GROUPS_ROW_SIZE 11
#define MAX_SECOND_TYPE_GROUPS_COL_SIZE 7

//COMPARE 2 LETTERS TO CHECK IF THEY ARE FROM FIRST TYPE GROUP
__device__ int search_first_type_groups(char char1, char char2)
{
    char first_type_groups[MAX_FIRST_TYPE_GROUPS_ROW_SIZE][MAX_FIRST_TYPE_GROUPS_COL_SIZE] = { "NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF" };
	int counter = 0;
	for (int i = 0; i < MAX_FIRST_TYPE_GROUPS_ROW_SIZE; i++)
	{
		for (int j = 0;j < MAX_FIRST_TYPE_GROUPS_COL_SIZE; j++)
		{
			if (first_type_groups[i][j] == char1)
				counter++;
			if (first_type_groups[i][j] == char2)
				counter++;
			if (counter == 2)
				return 1;
		}
		counter = 0;
	}
	return 0;
}

//COMPARE 2 LETTERS TO CHECK IF THEY ARE FROM SECOND TYPE GROUP
__device__ int search_second_type_groups(char char1, char char2)
{
    char second_type_groups[MAX_SECOND_TYPE_GROUPS_ROW_SIZE][MAX_SECOND_TYPE_GROUPS_COL_SIZE] = { "SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM" };
	int counter = 0;
	for (int i = 0; i < MAX_SECOND_TYPE_GROUPS_ROW_SIZE; i++)
	{
		for (int j = 0;j < MAX_SECOND_TYPE_GROUPS_COL_SIZE; j++)
		{
			if (second_type_groups[i][j] == char1)
				counter++;
			if (second_type_groups[i][j] == char2)
				counter++;
			if (counter == 2)
				return 1;
		}
		counter = 0;
	}
	return 0;
}

//RETURN THE APPROPRIATE WEIGHT
__device__ int calculate_score_between_two_characters(char c1, char c2, int* weights)
{
	if (c1 == c2)
		return weights[0];
	else if(search_first_type_groups(c1, c2) == 1)
		return (0-weights[1]);
	else if(search_second_type_groups(c1, c2) == 1)
		return (0-weights[2]);
	else 
		return (0-weights[3]);																																							
}

__global__ void calculate_score(char* cuda_seq1,char* cuda_sequance,int* cuda_res,int offset_mutant,int start_offset, int* weights)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	char seq1 = cuda_seq1[start_offset + index];
	char seq2 = cuda_sequance[index];
	
	if(index == offset_mutant)
	{	
		cuda_res[index] = 0;
	}
	else
	{
		if(index > offset_mutant)
		{
			seq2 = cuda_sequance[index - 1];
		}
		cuda_res[index] = calculate_score_between_two_characters(seq1, seq2, weights);
	}
}

int compute_on_gpu(char* seq1,char* sequences,int weights[],int offset,int offset_mutant,int length_offset_mutant, int* results)
{
	cudaError_t error = cudaSuccess;
	//ALLOCATIONS AND COPYING DATA TO CUDA
	size_t size_sequance = (strlen(sequences) + 1)*(sizeof(char));
	char* cuda_sequance;
	error = cudaMalloc((void**)&cuda_sequance, size_sequance);
	if(error != cudaSuccess)
	{
		fprintf(stderr,"failed to allocate memory %s\n",cudaGetErrorString(error));
	}
	error = cudaMemcpy(cuda_sequance,sequences,size_sequance,cudaMemcpyHostToDevice);
	if(error != cudaSuccess)
	{
		fprintf(stderr,"failed to copy memory %s\n",cudaGetErrorString(error));
	}
	
	size_t size_seq1 = (strlen(seq1) + 1)*(sizeof(char));
	char* cuda_seq1;
	error = cudaMalloc((void**)&cuda_seq1, size_seq1);
	if(error != cudaSuccess)
	{
		fprintf(stderr,"failed to allocate memory %s\n",cudaGetErrorString(error));
	}
	error = cudaMemcpy(cuda_seq1,seq1,size_seq1,cudaMemcpyHostToDevice);
	if(error != cudaSuccess)
	{
		fprintf(stderr,"failed to copy memory %s\n",cudaGetErrorString(error));
	}
	
	size_t size_res = (strlen(sequences) + 1)*(sizeof(int));
	int* cuda_res;
	error = cudaMalloc((void**)&cuda_res, size_res);
	if(error != cudaSuccess)
	{
		fprintf(stderr,"failed to allocate memory %s\n",cudaGetErrorString(error));
	}
	error = cudaMemcpy(cuda_res,results,size_res,cudaMemcpyHostToDevice);
	if(error != cudaSuccess)
	{
		fprintf(stderr,"failed to copy memory %s\n",cudaGetErrorString(error));
	}
	
	size_t size_weights = 4 * sizeof(int);
	int* cuda_w;
	error = cudaMalloc((void**)&cuda_w, size_weights);
	if(error != cudaSuccess)
	{
		fprintf(stderr,"failed to allocate memory %s\n",cudaGetErrorString(error));
	}
	error = cudaMemcpy(cuda_w,weights,size_weights,cudaMemcpyHostToDevice);
	if(error != cudaSuccess)
	{
		fprintf(stderr,"failed to copy memory %s\n",cudaGetErrorString(error));
	}
	//UPDATE THE SCORE FOR THE WORST CASE
	int max_score = strlen(sequences) *  -(weights[1] + weights[2] + weights[3]);
	int index_mutant = 0;
	
	for(int i = offset_mutant;i < (offset_mutant + length_offset_mutant);i++)
	{
		int number_element = (strlen(sequences) + 1) / 1024;
		int block_per_grid = number_element + 1;
		int thread_per_block = (strlen(sequences)+ 1) / block_per_grid;
		calculate_score<<<block_per_grid, thread_per_block>>>(cuda_seq1,cuda_sequance,cuda_res, i, offset, cuda_w);
		
		error = cudaMemcpy(results,cuda_res,size_res,cudaMemcpyDeviceToHost);
		if(error != cudaSuccess)
		{
			fprintf(stderr,"failed to copy memory %s\n",cudaGetErrorString(error));
		}

		int score = 0;
		for(int j = 0; j < strlen(sequences) + 1; j++)
		{
			score += results[j];
		}
		//UPDATE THE SCORE TO MAX SCORE
		if(score > max_score)
		{
			max_score = score;
			index_mutant = i;
		}
	}
	//SAVE THE FINAL SCORE IN A KNOWN PLACE
	results[0] = max_score;
	//FREE ALLOCATIONS IN CUDA
	if (cudaFree(cuda_sequance) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
	if (cudaFree(cuda_seq1) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
	if (cudaFree(cuda_res) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
	if (cudaFree(cuda_w) != cudaSuccess) {
        fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
	return index_mutant;
}
