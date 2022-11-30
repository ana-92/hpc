/* compute pi using Monte Carlo method */
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define CHUNKSIZE1 1000
#define CHUNKSIZE2 2000
/* message tags */
#define REQUEST 1
#define REPLY 2

int main(int argc, char *argv[])
{
	int iter, even;
	int in, out, i, iters, max, ix, iy, done, temp;
	double x, y, Pi,Pi_1, Pi_2,error, epsilon;
	int numprocs, myid, server, totalin, totalout, workerid;
	int rands_1[CHUNKSIZE1],rands_2[CHUNKSIZE2],request;
	MPI_Comm world, worker1,worker2;
	MPI_Group world_group, odd_group, even_group;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	world = MPI_COMM_WORLD;
	MPI_Comm_size(world, &numprocs);
	MPI_Comm_rank(world, &myid);
	server = numprocs - 1; /* last proc is server */

	if (myid == 0)
	{
		if (argc < 2)
		{
			fprintf(stderr, "Usage: %s epsilon\n", argv[0]);
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		sscanf(argv[1], "%lf", &epsilon);
	}
	MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Comm_group(world, &world_group);


	even = (numprocs+1)/2;
	int ranks_even[even];

	for(int i=0;i<even;i++){
		ranks_even[i]=2*i;
	}
	

	MPI_Group_incl(world_group, even, ranks_even, &even_group);
	MPI_Group_excl(world_group, even, ranks_even, &odd_group);

	MPI_Comm_create(world, world_group, &world);
	MPI_Comm_create(world, odd_group, &worker1);
	MPI_Comm_create(world, even_group, &worker2);

	MPI_Group_free(&odd_group);
	MPI_Group_free(&even_group);

	if (myid == server)
	{ /* I am the rand server */
		do
		{
			MPI_Recv(&request, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST,
					 world, &status);
			if (request)
			{
				for (i = 0; i < CHUNKSIZE2;)
				{
					rands_2[i] = random();
					if (rands_2[i] <= INT_MAX)
						i++;
				}
				MPI_Send(rands_2, CHUNKSIZE2, MPI_INT,
						 status.MPI_SOURCE, REPLY, world);
			}
		} while (request > 0);
	}
	if(myid%2!=0 && myid != server){
		request = 1;
		done = in = out = 0;
		max = INT_MAX; /* max int, for normalization */
		MPI_Send(&request, 1, MPI_INT, server, REQUEST, world);
		MPI_Comm_rank(worker1, &workerid);
		iter = 0;
		while (!done)
		{
			iter++;
			request = 1;
			MPI_Recv(rands_1, CHUNKSIZE1, MPI_INT, server, REPLY,
					 world, MPI_STATUS_IGNORE);
			for (i = 0; i < CHUNKSIZE1;)
			{
				x = (((double)rands_1[i++]) / max) * 2 - 1;
				y = (((double)rands_1[i++]) / max) * 2 - 1;
				if (x * x + y * y < 1.0)
					in++;
				else
					out++;
			}
			MPI_Allreduce(&in, &totalin, 1, MPI_INT, MPI_SUM,
						  worker1);
			MPI_Allreduce(&out, &totalout, 1, MPI_INT, MPI_SUM,
						  worker1);
			Pi_1 = (4.0 * totalin) / (totalin + totalout);
			error = fabs(Pi_1 - 3.141592653589793238462643);
			done = (error < epsilon || (totalin + totalout) > 100000000);
			request = (done) ? 0 : 1;
			if (myid == 0)
			{
				printf("\rpi worker 1= %23.20f", Pi_1);
				MPI_Send(&request, 1, MPI_INT, server, REQUEST,
						 world);
			}
			else
			{
				if (request)
					MPI_Send(&request, 1, MPI_INT, server, REQUEST,
							 world);
			}
		}
		MPI_Comm_free(&worker1);
	}
	if(myid%2==0 && myid != server)
	{ /* I am a worker process */
		request = 1;
		done = in = out = 0;
		max = INT_MAX; /* max int, for normalization */
        
		MPI_Send(&request, 1, MPI_INT, server, REQUEST, world);

		MPI_Comm_rank(worker2, &workerid);
		iter = 0;
		while (!done)
		{
			iter++;
			request = 1;
			MPI_Recv(rands_2, CHUNKSIZE2, MPI_INT, server, REPLY,
					 world, MPI_STATUS_IGNORE);
			for (i = 0; i < CHUNKSIZE2;)
			{
				x = (((double)rands_2[i++]) / max) * 2 - 1;
				y = (((double)rands_2[i++]) / max) * 2 - 1;
				if (x * x + y * y < 1.0)
					in++;
				else
					out++;
			}
			MPI_Allreduce(&in, &totalin, 1, MPI_INT, MPI_SUM,
						  worker2);
			MPI_Allreduce(&out, &totalout, 1, MPI_INT, MPI_SUM,
						  worker2);
			Pi_2 = (4.0 * totalin) / (totalin + totalout);
			error = fabs(Pi_2 - 3.141592653589793238462643);
			done = (error < epsilon || (totalin + totalout) > 100000000);
			request = (done) ? 0 : 1;
			if (myid == 0)
			{
				printf("\rpi worker 2 = %23.20f", Pi_2);
				MPI_Send(&request, 1, MPI_INT, server, REQUEST,
						 world);
			}
			else
			{
				if (request)
					MPI_Send(&request, 1, MPI_INT, server, REQUEST,
							 world);
			}
		}
		MPI_Comm_free(&worker2);
	}

	if (myid == 0)
	{
		printf("\npoints: %d\nin: %d, out: %d, <ret> to exit\n",
			   totalin + totalout, totalin, totalout);
		Pi = (4.0 * totalin) / (totalin + totalout);
		printf("\rpi global = %23.20f\n", Pi);
		getchar();
	}
	MPI_Finalize();
	return 0;
}
