#include <mpi.h>
#include <time.h>
#include <iostream>

using namespace std;



int main(int argc, char *argv[]) 
{
	//int **matrix = NULL;	// �������	
	int *matrix = nullptr;
	int ProcNum, ProcRank;	
	int N=0, M=0; // ����������� ������� 
	int *result = NULL; // ������ � �����������	
	int *tmp =NULL;	
	MPI_Status status;
	double times;


	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	
	if(ProcRank ==0) //��������� ������� � �������� ������ ��������� ������ ������ 
	{
        cout << "Enter matrix sizes.\nExample :\n5\n5" <<endl;
       
		cin >> N;
		cin >> M;
		
		//matrix = new int*[N];
		matrix = new int[N*M];
		result = new int[N];
 		for(int i =0; i< N; i++)
		{			
			result[i] =0;
		}		
		srand(time(0));
		for( int i=0; i<N*M; i++)
				matrix[i] = rand() %10;
		
	}
	
    //�������� ����������� ���� ���������
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	
	if (ProcRank == 0)
	{
        //��������� ������ ������� �� ���������
		for(int i =0; i< N; i++)
		{
			if(i%ProcNum != 0)
			{
				MPI_Send(matrix, M, MPI_INT, i%ProcNum, i, MPI_COMM_WORLD);
			}
		}
	}
	
    //��������� �����, � ���������� �������� � ������� �������
	if(ProcRank != 0)
	{
		tmp = (int*) malloc(sizeof(int)*M);
		int sum = 0;
		for(int i= ProcRank; i < N; i+= ProcNum)
		{
			MPI_Recv(tmp, M, MPI_INT, 0, i, MPI_COMM_WORLD, &status);
			for(int j=0; j< M; j++)
			{
				sum+= tmp[j];
			}
			MPI_Send(&sum, 1, MPI_INT, 0, i, MPI_COMM_WORLD);
			
		}
		free(tmp);
		
	}
	else
	{		
		
		int sum =0;
		for(int  i = ProcRank; i< N ; i += ProcNum)
		{
			for(int j = 0; j< M; j++)
				sum += matrix[i][j];
			result[i] = sum;
		}			
	}
		

	if(ProcRank == 0)
	{
		for(int i =0; i< N; i++)
			{
				if(i%ProcNum != 0)
				{
					MPI_Recv(result + i, 1, MPI_INT, MPI_ANY_SOURCE, i , MPI_COMM_WORLD, &status);
				}
		}
        //����� ������� � �����������
		printf(" Matrix %dx%d \n", N, M);
		printf("----------------------------------------\n");
		for(int i =0 ; i<N; i++)
		{
			for(int j =0; j<M; j++) 
			{
				printf("%d ", matrix[i][j]);
			}
			printf("	sum = %d\n", result[i] );
		}		 
	}	
	
	MPI_Finalize();
	cin >> N;

    if (matrix != NULL)
    {
        for (int i = 0; i < N; i++)
            free(matrix[i]);
        free(matrix);
    }
	if(result != NULL)
        free(result);

	return 0;
}
