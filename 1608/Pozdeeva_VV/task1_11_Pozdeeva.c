#define _CRT_SECURE_NO_WARNINGS

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int main(int argc, char *argv[]) 
{
	int **matrix = NULL;	// матрица
	int ProcNum, ProcRank;	
	int N=0, M=0; // размерности матрицы 
	int *result = NULL; // массив с результатом	
	int *tmp =NULL;	
	MPI_Status status;
	int i, j, sum;
	double times;
	

	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	
	if(ProcRank ==0) //заполняем матрицу и передаем другим процессам нужные строки 
	{
        printf("Enter matrix sizes. Example: 5 5\n");
        scanf("%d %d", &N, &M);

		matrix =(int**) malloc(sizeof(int*)* N);			
		result = (int *) malloc(sizeof(int) *N);
		for(i =0; i< N; i++)
		{
			matrix[i] = (int*)malloc(sizeof(int)*M);
			result[i] =0;
		}		
		srand(time(0));
		rand();
		for( i=0; i<N; i++)
			for( j=0; j<M; j++)
				matrix[i][j] = rand() %100;
	}
	
    //передаем размерности всем процессам
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	
	if (ProcRank == 0)
	{
        //рассылаем строки матрицы по процессам
		for(i =0; i< N; i++)
		{
			if(i%ProcNum != 0)
			{
				MPI_Send(matrix[i], M, MPI_INT, i%ProcNum, i, MPI_COMM_WORLD);
			}
		}
	}
	times = MPI_Wtime();
    //вычисляем суммы, и возвращаем резултат в нулевой процесс
	if(ProcRank != 0)
	{
		tmp = (int*) malloc(sizeof(int)*M);
		sum = 0;
		for( i= ProcRank; i < N; i+= ProcNum)
		{
			MPI_Recv(tmp, M, MPI_INT, 0, i, MPI_COMM_WORLD, &status);
			for(j=0; j< M; j++)
			{
				sum+= tmp[j];
			}
			MPI_Send(&sum, 1, MPI_INT, 0, i, MPI_COMM_WORLD);
			
		}
		free(tmp);
		
	}
	else
	{		
		
		sum =0;
		for( i = ProcRank; i< N ; i += ProcNum)
		{
			for(j = 0; j< M; j++)
				sum += matrix[i][j];
			result[i] = sum;
		}			
	}

	printf("%d  - %20.10f\n",ProcRank, MPI_Wtime() - times);

	if(ProcRank == 0)
	{
		for(i =0; i< N; i++)
			{
				if(i%ProcNum != 0)
				{
					MPI_Recv(result + i, 1, MPI_INT, MPI_ANY_SOURCE, i , MPI_COMM_WORLD, &status);
				}
		}
        //вывод матрицы и результатов
		printf(" Matrix %dx%d \n", N, M);
		printf("----------------------------------------\n");
		for(i =0 ; i<N; i++)
		{
			for(j =0; j<M; j++) 
			{
				printf("%d ", matrix[i][j]);
			}
			printf("	sum = %d\n", result[i] );
		}		 
	}	
	
	MPI_Finalize();

    if (matrix != NULL)
    {
        for (i = 0; i < N; i++)
            free(matrix[i]);
        free(matrix);
    }
	if(result != NULL)
        free(result);

	return 0;
}
