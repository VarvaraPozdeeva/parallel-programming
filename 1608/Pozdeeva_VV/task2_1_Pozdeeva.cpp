#include <iostream>
#include <mpi.h>
#include <time.h>

#define PUT 101
#define GET 102
#define FINISHp 103
#define FINISHc 104

#define OK -201
#define FULL -202 // buf with res is full
#define EMPTY -203 // buf with res is empty
#define ENDED -204 // resourses ended

#define SIZE_BUF_PROD 10
#define SIZE_BUF_CONS 15
#define SIZE_BUF 5

typedef struct data
{
	int operation;
	int Rank;
	int res;

} Data;

using namespace std;

void CreateResourse(int* buf, int size, int rank)
{
	srand(time(0));
	for(int i=0; i<size; i++)
		buf[i] = rand()%100 - rank;
}
void SendResours(int* buf, int count, int rank)
{
	Data info;
	info.Rank = rank;
	info.operation = PUT;
	info.res = buf[count-1];			
	MPI_Send(&info,3,MPI_INT,0,0,MPI_COMM_WORLD);// send to 0 resourse
	
}
int PutInBuf(int *buf, int size, int res)
{
	int result;
	//cout << 0 << " PUT"<< endl;
	result = FULL;
	for(int i=0; i<size; i++)
		if(buf[i]==-100)
		{
			buf[i] = res;
			result=OK;
			break;
		}
	return result;
}
int GetFromBuf(int *buf, int size)
{
	int result;
	//cout <<  0<< " GET "<< endl;
	for(int i=0; i<5; i++)
		if(buf[i]== -100)
			result = EMPTY;
		else 
		{
			result = buf[i];
			buf[i]=-100;
			break;
		}
	return result;
}

int main(int argc, char** argv)
{
	
	int ProcNum, ProcRank;	
	MPI_Status status;

	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	

    if(ProcRank == 0) // proc with buf
	{
		
		Data info;
		int buf[SIZE_BUF];
		int result;
		int numProd = (ProcNum-1)/2 + (ProcNum-1)%2;
		int numCons = (ProcNum-1)/2;
		
		for(int i = 0; i < SIZE_BUF; i++)
			buf[i]=-100;

		double times = MPI_Wtime();
		while(true)
		{
			MPI_Recv(&info, 3, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			
			if(info.operation == PUT)
			{
				result = PutInBuf(buf, SIZE_BUF, info.res);
				MPI_Send(&result, 1, MPI_INT, info.Rank, PUT, MPI_COMM_WORLD);
			}
			else if(info.operation == FINISHp)
			{
				numProd--;
			}
			else if(info.operation == FINISHc)
			{
				numCons--;
			}
			else if(info.operation == GET)
			{
				result = GetFromBuf(buf, SIZE_BUF);
				if(result == EMPTY)
					if(numProd == 0)
						result = ENDED;	
				MPI_Send(&result, 1, MPI_INT, info.Rank, GET, MPI_COMM_WORLD);
			}

			if(numCons ==0 && numProd == 0)
				break;
			
		}
		printf("Time = %0.10f",MPI_Wtime() - times);

		cout << ProcRank<< " STOP "<< endl;
		MPI_Finalize();

    }
    else if(ProcRank%2 == 1) // producer
	{
		Data info;
		info.Rank = ProcRank;		
		int buf[SIZE_BUF_PROD];
		int result;
		int countRes = SIZE_BUF_PROD;

		CreateResourse(buf, SIZE_BUF_PROD, ProcRank);
		//cout<<"Resours proc #"<<ProcRank<<" - ";
		//for(int i=0; i<SIZE_BUF_PROD; i++)
			//cout << buf[i]<< ", ";
		//cout << endl;
		while(countRes != 0)
		{
			SendResours(buf, countRes, ProcRank);

			MPI_Recv(&result, 1 , MPI_INT , 0 , PUT , MPI_COMM_WORLD , &status);
			if(result == OK)
			{
				countRes--;
			}			
		}
		info.operation = FINISHp;
		//cout << ProcRank<< " STOP, i send  "<< countRes+sizeBuf << endl;
		MPI_Send(&info,3,MPI_INT,0,0,MPI_COMM_WORLD);		

		MPI_Finalize();
           
    }		
	else // consumer
	{ 
		int result;
		int buf[SIZE_BUF_CONS];
		int countRes =0;

		Data info;
		info.operation= GET;
		info.Rank = ProcRank;		
	
		while(countRes!= SIZE_BUF_CONS)
		{
			info.res = -1;
			MPI_Send(&info,3,MPI_INT,0,0,MPI_COMM_WORLD);
			MPI_Recv(&result, 1 , MPI_INT , 0 , GET , MPI_COMM_WORLD , &status);
			if(result != EMPTY)
				if(result != ENDED)
				{
					buf[countRes++]=result;					
				}
				else break;
		}		

		info.operation = FINISHc;
		//cout << ProcRank<< " STOP, i get  "<< countRes << endl;

		MPI_Send(&info,3,MPI_INT,0,0,MPI_COMM_WORLD);
				
		MPI_Finalize();
	}
	
	return 0;
}