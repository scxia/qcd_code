
#include<mpi.h>


MPI_Datatype MPI_position;

void create_position_MPI_datatype()
{
	
	int point_result_blocklengths[2] = { 4,2 };
	MPI_Aint point_result_displacement[2] = { 0,4 * sizeof(int) };
	MPI_Datatype point_result_types[2] = { MPI_INT,MPI_DOUBLE };
	MPI_Type_create_struct(2, point_result_blocklengths, point_result_displacement, point_result_types, &MPI_position);
	MPI_Type_commit(&MPI_position);
}




