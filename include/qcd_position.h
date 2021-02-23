


#ifndef H_QCD_POSITION
#define H_QCD_POSITION 1


typedef struct 
{
	
	int loc[4];
	double val[2];

}position_result;

extern MPI_Datatype MPI_position;

void create_position_MPI_datatype();

#endif // !H_QCD_POSITION
