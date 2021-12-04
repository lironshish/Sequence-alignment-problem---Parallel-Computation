#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <mpi.h>
#include <omp.h>
#include "cudaHeader.h"


#define WEIGHTS_NUMBER 4
#define MAX_LENGTH_SEQ1 3000
#define MAX_LENGTH_SEQ2 2000

enum ranks{ROOT};

//ANNOUNCEMENT OF FUNCTIONS
void send_initial_data_to_MPI(int weights[], char* seq,int num_proc);
void receive_initial_data_to_MPI(int rank);
void send_sequances(char* sequences, int num_proc, char* seq1, int weights[]);
void receive_sequances(char* seq1,int rank,int length_seq1, int weights[]);
int send_seq_openMP(char* seq1, char* sequences,int weights [],int offset, int rank,int *final_result_index_mutant);
int check_is_upper(char* seq, int size);

//RECEIVE THE INITIAL DATA FROM PROCESS 0 
void receive_initial_data_to_MPI(int rank)
{
	MPI_Status status;
	int length;
	int weights[WEIGHTS_NUMBER];
	MPI_Recv(weights,WEIGHTS_NUMBER,MPI_INT,ROOT,0,MPI_COMM_WORLD,&status);
	MPI_Recv(&length,1,MPI_INT,ROOT,0,MPI_COMM_WORLD,&status);
	char seq[length];
	MPI_Recv(seq,length,MPI_CHAR,ROOT,0,MPI_COMM_WORLD,&status);
    receive_sequances(seq,rank, length,weights);
}

//SEND THE INITIAL DATA FROM PROCESS 0 TO THE REST OF THE PROCESSES
void send_initial_data_to_MPI(int weights[], char* seq,int num_proc)
{
	int seq_length = strlen(seq) + 1;
	for(int  i = 1; i < num_proc; i++)
		{
			MPI_Send(weights, WEIGHTS_NUMBER, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&seq_length, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(seq, seq_length, MPI_CHAR, i, 0, MPI_COMM_WORLD);
		}	
}

//SEND THE SEQ2, LENGTH OF SEQ2 AND LENGTH OF OFFSET FROM PROCESS 0 TO THE REST OF THE PROCESSES
void send_sequances(char* sequences, int num_proc, char* seq1, int weights[])
{
	MPI_Status status;
    int seq_length = strlen(sequences) + 1;
	for(int  i = 1; i < num_proc; i++)
    {
        MPI_Send(&seq_length, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(sequences, seq_length, MPI_CHAR, i, 0, MPI_COMM_WORLD);
    }
    int length_offset_per_proc = 0;
    for(int i = 1; i <num_proc; i++)
    {
        MPI_Send(&length_offset_per_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        length_offset_per_proc += (strlen(seq1) - strlen(sequences))/(num_proc - 1);
        if((i + 1) == num_proc)
        {
            length_offset_per_proc = (strlen(seq1) - strlen(sequences));   
        }
        MPI_Send(&length_offset_per_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
	//INITIALIZE THE SCORE TO THE MINIMUM OPTION
	int final_score = (strlen(seq1) *  -(weights[1] + weights[2] + weights[3])),score,index_offset,result_index_offset,mutant,final_result_index_mutant;
	for(int i = 1; i< num_proc; i++)
	{
		//RECIEVE THE SCORE,INDEX MUTANT AND OFFSET FROM THE PROCESSES
		MPI_Recv(&score,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);    
		MPI_Recv(&index_offset,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);    
		MPI_Recv(&mutant,1,MPI_INT,i,0,MPI_COMM_WORLD,&status);   

		if(final_score < score)
		{
			final_score = score;
			result_index_offset = index_offset;
			final_result_index_mutant = mutant;
		} 	

	}
	printf("offset(n) = %d\tmutant(k) = %d\tmax score = %d\n",result_index_offset,final_result_index_mutant, final_score);
}



//RECEIVE THE SEQ2, LENGTH OF SEQ2 AND LENGTH OF OFFSET FROM PROCESS 0 
void receive_sequances(char* seq1,int rank, int length_seq1, int weights[])
{
    MPI_Status status;
    int seq_length;
	int final_result_index_mutant, result_index_mutant;
    MPI_Recv(&seq_length,1,MPI_INT,ROOT,0,MPI_COMM_WORLD,&status);
    int start_offset, end_offset,result_index_offset, score, final_score = (length_seq1 *  -(weights[1] + weights[2] + weights[3])) ;
    while(seq_length != -1)
    {
        char sequences[seq_length];
        MPI_Recv(sequences,seq_length,MPI_CHAR,ROOT,0,MPI_COMM_WORLD,&status);
        MPI_Recv(&start_offset,1,MPI_INT,ROOT,0,MPI_COMM_WORLD,&status);
        MPI_Recv(&end_offset,1,MPI_INT,ROOT,0,MPI_COMM_WORLD,&status);        
		
		for(int i = start_offset; i < end_offset; i++)
		{
        	score = send_seq_openMP(seq1,sequences,weights,i, rank, &result_index_mutant);
			if(final_score < score)
			{
				final_score = score;
				result_index_offset = i;
				final_result_index_mutant =result_index_mutant;
			}

		}
		//SEND THE SCORE,INDEX MUTANT AND OFFSET FROM THE PROCESSES TO ROOT
		MPI_Send(&final_score, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&result_index_offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&final_result_index_mutant, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		//INITIALIZE THE SCORE TO THE MINIMUM OPTION SO THAT IT DOES NOT SAVE THE RESULT OF THE LAST STRING
		final_score = (length_seq1 *  -(weights[1] + weights[2] + weights[3])) ;
		//RECEIVE A "SIGNAL" THAT WE HAVE FINISHED
        MPI_Recv(&seq_length,1,MPI_INT,ROOT,0,MPI_COMM_WORLD,&status);
    }
}



int send_seq_openMP(char* seq1, char* sequences,int weights [],int offset, int rank,int *final_result_index_mutant)
{
    int offset_mutant, length_offset_mutant,index_nutant,result_index_mutant, thread_id, score, final_score = (strlen(sequences) *  -(weights[1] + weights[2] + weights[3]));
	int results[strlen(sequences)+ 1];
    #pragma omp parallel private(thread_id, offset_mutant,length_offset_mutant, score, results, index_nutant) shared (final_score, result_index_mutant)
    {
        length_offset_mutant = (strlen(sequences)+1)/omp_get_num_threads();
        thread_id = omp_get_thread_num();
        offset_mutant = thread_id * length_offset_mutant;
        
		if(thread_id == (omp_get_num_threads() - 1))
            length_offset_mutant = strlen(sequences) - offset_mutant + 1;
        
		if(thread_id == 0)
        {
            offset_mutant = 1;
        }
        //SEND TO CUDA FOR AUXILIARY CALCULATIONS
		index_nutant = compute_on_gpu(seq1,sequences, weights,offset,offset_mutant,length_offset_mutant, results);
		score = results[0];

		//UPDATE THE SCORE TO MAX SCORE
		#pragma omp critical(final_score)
		if(final_score < score)
		{
			final_score = score;
			result_index_mutant = index_nutant;
		}
    }
	*final_result_index_mutant = result_index_mutant;
	return final_score;
}

//CHECK IF THE LETTERS ARE UPPERCASE
int check_is_upper(char* seq, int size)
{
	for(int i=0; i <size; i++)
	{
		if(!isupper(seq[i]))
			return 0;
	}
	return 1;	
}

int main(int argc, char *argv[])
{
	//VARIABLES FOR MPI INIT
 	int my_rank,num_process;
 	//MPI INIT
 	MPI_Init(&argc, &argv);
 	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
 	MPI_Comm_size(MPI_COMM_WORLD, &num_process);

	int weights[WEIGHTS_NUMBER] = { 0 };
	int num_sequences_for_comparison;

	if (my_rank==ROOT) //PROCESS 0
 	{
		//SCAN DATA
		//SCAN THE WEIGHTS VALUES
		for (int i = 0; i < WEIGHTS_NUMBER; i++)
		{
			scanf("%d", &weights[i]);
			if (weights[i] < 0) 
			{
				printf("Weight should be a positive number\n");
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
		}
    	//SCAN SEQ1
		char* seq1 = (char*)malloc(sizeof(char)*(MAX_LENGTH_SEQ1));
		scanf("%s", seq1);
		//CHECK WHETHER THE SEQUENCE IS IN CAPITAL LETTERS
		if(check_is_upper(seq1, strlen(seq1)) == 0)
		{
			printf("Sequence should be in capital letters\n");
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}

    	//SCAN ALL SEQUENCES2
		scanf("%d", &num_sequences_for_comparison);

		char** sequences = (char**)malloc(sizeof(char*)*(num_sequences_for_comparison));
		for (int i = 0; i < num_sequences_for_comparison; i++)
		{
			sequences[i] = (char*)malloc(sizeof(char)*(MAX_LENGTH_SEQ2));
			scanf("%s", sequences[i]);
			if(check_is_upper(sequences[i], strlen(sequences[i])) == 0)
			{
				printf("Sequence should be in capital letters\n");
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
		}

		double time = MPI_Wtime();
		//SEND THE WEIGHTS, SEQ1 AND LENGTH OF SEQ1 TO ALL PROCESSES WITH MPI
		send_initial_data_to_MPI(weights, seq1,num_process);
        //SEND THE SEQ2, LENGTH OF SEQ2 AND LENGTH OF OFFSET TO ALL PROCESSES WITH MPI
		for(int i = 0; i < num_sequences_for_comparison; i++)
        {
            send_sequances(sequences[i], num_process,seq1,weights);
        }
        int finish = -1;
		//SEND -1 TO KNOW THE SEQUENCES IS OVER
        for(int  i = 1; i < num_process; i++)
            MPI_Send(&finish, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

		//FREE ALLOCATIONS 
    	free(seq1);
    	for(int i =0; i <num_sequences_for_comparison; i++)
        	free(sequences[i]);
    	free(sequences);

 	} 
	else //THE OTHER PROCRSSES -> my_rank != 0
	{	
 		receive_initial_data_to_MPI(my_rank);
	}

	MPI_Finalize();
	return 0;
}