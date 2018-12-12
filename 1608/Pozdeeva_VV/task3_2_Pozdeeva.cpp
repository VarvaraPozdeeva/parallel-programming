#include <iostream>
#include <time.h>
#include <mpi.h>
#include <cmath>

using namespace std;

double* CreateRandMatrix(int size)
{		
	double* tmp = new double[size*size];	
	srand(time(0));
	for(int i=0; i< size*size; i++)
		tmp[i] = (rand()%10000)/1000.0f;	 
	return tmp;
}
void MultiMatrix(double* A, double* B, double* C, int n)
{
	for(int i = 0; i < n; ++i)      
        for(int j = 0; j < n; ++j)        
            for(int k = 0; k < n; ++k)   
                C[i * n + j] += A[i * n + k] * B[k * n + j];    
}
int PrintMatrix(double* A,int N)
{
	for(int i=0; i<N*N; i+=N)
		{
			for(int j=0; j< N; j++)				
				cout<<A[i+j] << " ";
			cout << endl;
		}
	cout << endl;
	return 1;
}
int CheckEqual(double *A, double *C, int N)
{
	for(int i=0; i<N*N; i++)
		if(abs(A[i] - C[i]) > 0.000001)
			return 0;
	return 1;
}
typedef struct Grid
{
	int rank;
	int rank_grid;
	MPI_Comm Comm_grid;
	MPI_Comm Comm_Row;
	int procNum;
	int rank_row;
	int root;
	int up;
	int down;
	int q;
	int coord[2];
}GRID;

int CreateGrid(GRID * grid, int N)
{
	MPI_Comm_size(MPI_COMM_WORLD, &grid->procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &grid->rank);

	int const ndims =2;
	int dims[ndims], periods[ndims], reorder =1;
	grid->q=dims[0]=dims[1] = sqrt((double)grid->procNum);
	if((dims[0]*dims[0] != grid->procNum) || (N % dims[0] != 0))
	{
		MPI_Finalize();
		cout << " Error num of proc. It must be a perfect square!" << endl;
		return 1;
	}

	periods[0]=periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &grid->Comm_grid);

	int coords[2] = { 0, 0 };
	MPI_Cart_rank(grid->Comm_grid, coords, &grid->root);
	MPI_Cart_shift(grid->Comm_grid, 0, 1, &grid->up, &grid->down); 
	
	MPI_Comm_rank(grid->Comm_grid, &grid->rank_grid);	
	MPI_Cart_coords(grid->Comm_grid, grid->rank_grid, 2, grid->coord);

	int remain_dims[2]={0,1};
	MPI_Cart_sub(grid->Comm_grid, remain_dims, &grid->Comm_Row);
	MPI_Comm_rank(grid->Comm_Row, &grid->rank_row);

	return 0;
}
int main(int argc, char** argv)
{
    double *A = nullptr, *B = nullptr, *C = nullptr;
    double *blokA = nullptr, *blokA1 = nullptr, *blokB = nullptr, *blokC = nullptr;
	int N, sizeBlok, time;
	GRID *grid = new GRID;
	MPI_Status status;
	MPI_Datatype TypeBlok;

	if(argc>1)
		N=atoi(argv[1]);
	else
	{		
		cout << "Enter Matrix size. Example: Matrix.exe 4" << endl;		
		return 2;
	}

	MPI_Init(&argc, &argv);
	if ( CreateGrid(grid, N) == 1)
		return 1;

	sizeBlok = N/sqrt((double)grid->procNum);	
	blokA = new double[sizeBlok*sizeBlok];
	blokA1 = new double[sizeBlok*sizeBlok];
	blokB = new double[sizeBlok*sizeBlok];
	blokC = new double[sizeBlok*sizeBlok];	
	C = new double[N*N];
	for(int i=0; i<sizeBlok*sizeBlok; i++)
		blokC[i]=0;

	MPI_Datatype temptype;
	MPI_Type_vector(sizeBlok, sizeBlok, N, MPI_DOUBLE, &temptype);
	MPI_Type_create_resized(temptype, 0, sizeBlok*sizeBlok*sizeof(double), &TypeBlok);
	MPI_Type_commit(&TypeBlok);

	if(grid->rank_grid == grid->root)
	{			
		A = CreateRandMatrix(N);
		B = CreateRandMatrix(N);
		//PrintMatrix(A, N);
		//PrintMatrix(B, N);	

		for (int i = 0; i < sizeBlok; ++i)
			for (int j = 0; j < sizeBlok; j++)
			{
				blokA[i * sizeBlok + j] = A[i * N + j];
				blokB[i * sizeBlok + j] = B[i * N + j];
			}
		
		time = MPI_Wtime();
		for(int i=0; i< grid->q; i++ )
				for(int j=0; j<grid->q; j++)
					if((i!=0) || (j!=0))
					{
						int rank, coord[2] = {i,j};
						MPI_Cart_rank(grid->Comm_grid, coord, &rank);					
						MPI_Send(&A[i*sizeBlok*N + j* sizeBlok], 1, TypeBlok, rank, 0, grid->Comm_grid);
						MPI_Send(&B[i*sizeBlok*N + j* sizeBlok], 1, TypeBlok, rank, 1, grid->Comm_grid);
					}
	}
	else
	{
		MPI_Recv(blokA, sizeBlok * sizeBlok, MPI_DOUBLE, 0, 0, grid->Comm_grid, &status);
		MPI_Recv(blokB, sizeBlok * sizeBlok, MPI_DOUBLE, 0, 1, grid->Comm_grid, &status);		
	}

	for(int i=0; i<grid->q; i++)
	{
		int sourse = (grid->coord[0]+i)%(grid->q);	

		if(sourse == grid->rank_row)
			for(int j=0; j<sizeBlok*sizeBlok; j++)
				blokA1[j] = blokA[j]; 

		MPI_Bcast(blokA1, sizeBlok*sizeBlok , MPI_DOUBLE, sourse, grid->Comm_Row);
				
		MultiMatrix(blokA1, blokB, blokC, sizeBlok);		

		MPI_Sendrecv_replace(blokB, sizeBlok * sizeBlok, MPI_DOUBLE, grid->up, 1, grid->down, 1, grid->Comm_grid, &status);
	}

	if(grid->rank_grid == grid->root)
	{
		for (int i = 0; i < sizeBlok; ++i)
				for (int j = 0; j < sizeBlok; ++j)
					C[i * N + j] = blokC[i * sizeBlok + j];
		for (int i = 0; i < grid->q; ++i)
				for (int j = 0; j < grid->q; ++j)
					if ((i != 0) || (j != 0))
					{
						int source, block_coords[2] = { i, j };
						MPI_Cart_rank(grid->Comm_grid, block_coords, &source);
						MPI_Recv(&C[i * N * sizeBlok + j * sizeBlok], 1, TypeBlok, source, 4, grid->Comm_grid, &status);
					}

		printf("Time parallel = %.10f\n", MPI_Wtime() - time);
        //cout << "Matrix C (parallel)" << endl;
		//PrintMatrix(C, N);

		double *res = new double [N*N];
		for(int i=0; i<N*N; i++)
			res[i]=0;

		time = MPI_Wtime();
		MultiMatrix(A, B, res, N);
		printf("Time linear = %.10f\n", MPI_Wtime() - time);

       // cout << "Matrix C (linear) " << endl;
		//PrintMatrix(res, N);

		if(CheckEqual(C, res, N) == 1)
			cout <<"OK. Matrices are equal"<< endl;
		else 
			cout <<"ERROR. Matrices are not equal"<<endl;
        delete res;
	}

	else MPI_Send(blokC, sizeBlok * sizeBlok, MPI_DOUBLE, grid->root, 4, grid->Comm_grid);
	MPI_Type_free(&TypeBlok);
    MPI_Comm_free(&grid->Comm_grid);
    MPI_Comm_free(&grid->Comm_Row);

	MPI_Finalize();

	delete A;
	delete B;
	delete C;
	delete blokA;
	delete blokB;
	delete blokC;
	delete blokA1;
	delete grid;

	return 0;
}