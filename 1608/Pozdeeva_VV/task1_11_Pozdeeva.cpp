#include <mpi.h>
#include <time.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) 
{
	int *matrix = nullptr;
	int ProcNum, ProcRank;	
	int N=0, M=0; // размерности матрицы 
	int *result = nullptr; // массив с результатом		
	MPI_Status status;
	double times;
    int k = 0, l = 0;

	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	
	if(ProcRank ==0) //заполняем матрицу
	{
        if (argc < 3)
        {
            cout << "Enter Matrix Sizes\nN:" << endl;
            cin >> N;
            cout << "M:" << endl;
            cin >> M;
        }
        else
        {
            N = atoi(argv[1]);
            M = atoi(argv[2]);
        }
		
		matrix = new int[N*M];
        srand(time(0));
        for (int i = 0; i<N*M; i++)
            matrix[i] = rand() % 10;

		result = new int[N];
        for(int i =0; i< N; i++)
			result[i] =0;		
	}
	
    //передаем размерности всем процессам
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	
	if (ProcRank == 0)
	{
        //рассылаем строки матрицы по процессам

        k = N / ProcNum; // количество строк для каждого процесса
        l = N % ProcNum; //оставшееся количество, первому процессу		
		if (k == 0)
			k = 1;
		times = MPI_Wtime();
		
		for(int i =1; i< ProcNum; i++)
            MPI_Send(matrix + (i)*k*M, k*M, MPI_INT, i, i, MPI_COMM_WORLD);
	}
	
   

    //вычисляем суммы, и возвращаем резултат в нулевой процесс
	if(ProcRank != 0)
	{
        k = N / ProcNum;
		if (k == 0)
			k = 1;
		int *tmp = new int[M*k];
        int *res = new int[k];
        int sum;
       
        MPI_Recv(tmp, M*k, MPI_INT, 0, ProcRank, MPI_COMM_WORLD, &status);
        
        for(int i=0; i < k; i++)
		{
            sum = 0;
            for (int j = 0; j < M; j++)
                sum += tmp[j + M*i];
            res[i] = sum; 
		}	
        MPI_Send(res, k, MPI_INT, 0, ProcRank, MPI_COMM_WORLD);
        delete tmp;
        delete res;
	}
	else
	{	
        int sum;
        for (int i = 0; i < k; i++)
        {
            sum = 0;
            for (int j = 0; j<M; j++)
                sum += matrix[j+i*M];
            result[i] = sum;
        }
        
        if (l != 0)
        {            
            for (int i = 0; i < l; i++)
            {
                sum = 0;
                for (int j = 0; j < M; j++)
                    sum += matrix[k*ProcNum*M + j + i*M];
                result[i + k*ProcNum] = sum;
            }
        }
	}
	
   

	if(ProcRank == 0)
	{		
        for (int i = 1; i< ProcNum; i++)
        {
			MPI_Recv(result + i*k, k, MPI_INT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status);				
		}
		printf("Time of Proc # %d is  %.10f\n", ProcRank, MPI_Wtime() - times);
        //вывод матрицы и результатов
		printf(" Matrix %dx%d \n", N, M);
		printf("----------------------------------------\n");
		for(int i =0 ; i<N; i++)
		{
			for(int j =0; j<M; j++) 
			{
				printf("%d ", matrix[i*M+j]);
			}
			printf("	sum = %d\n", result[i] );
		}		 
		
    }	

    if (matrix != nullptr)
        delete matrix;
    if (result != nullptr)
        delete result;
    
    MPI_Finalize();
	
	return 0;
}
