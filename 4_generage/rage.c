/* ================================================================== */
/* RAGE.c - Processes a BLAST Output file, finds non-symmetric Hits   */
/* and hits to possible multidomain proteins, and Clusters into       */
/* families 							      */
/*                                                                    */ 
/* 		SERIAL VERSION 2.1a                                   */
/* 								      */
/* Anton Enright  - Computational Genomics Group (C. Ouzounis)        */
/*								      */
/* EMBL - European Bioinformatics Institute 			      */
/* ================================================================== */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>


/* DEFINITIONS - MAXIMUMS, DEFAULTS, ARCHITECTURE WORD SIZE        */

#define MAXPROTEINS 50000  /* Max Proteins in Database		   */
#define MAXHITS 1495000	  /* Max Hits 				   */
#define PROTEINSIZE 60    /* Max size of Protein Names (chars)     */
#define WORDSIZE 32       /* Word Size on This Architecture        */
#define MAXPEPTIDE 10000  /* Max Size of a peptide in Amino Acids  */ 
#define MAXDOMAINS 30	  /* Max No of Domains a Protein can Have  */

/* Function Prototypes */

void dump_matrix(int);
void dump_matrix2(int);  
void clear_matrix2(int);
void dump_similarities(int, char *);
void dump_clusters(int);
void dump_paralogues(int);
void symmetrify_matrix(int);
void multi_domain(int, int);
void reset_matrix(int);
int strptrcmp (const void *,const void *);
int strptrsearch (const void *,const void *);

int do_smith_waterman(int,int,char *,int);
void begin_cluster(int);
int cluster_matrix(int, int);

/* Smith Waterman Stuff  - Default Cutoffs / Randomizations */

int randomize = 100;
int zscore1 = 10;
int zscore2 = 3;
int temporary=0;


/* Memory Allocation Routines - Prototypes */

unsigned int **imatrix(long,long,long,long);
int *ivector(int);
void free_ivector(int *,int,int);
void free_imatrix(int **,int,int,int,int);

/* BitMatrix Stuff */

int getbit(int, int);
void setbit(int, int);
void unsetbit(int, int);
int getbit2(int, int);
void setbit2(int, int);
void unsetbit2(int, int); 

/* Multidomain Storage for Reallocation */

int mdc;
int allocated=0;
int totalhits=0;
int *multidom_index;
int *multidom;
int *domain1;
int *domain2;
int *clustered_domains;

/* Global Variable Declarations */

char database[40];
char *prssdir;
typedef struct sequencestruct SEQSTRUCT;
struct sequencestruct
        {
        char sequence_id[PROTEINSIZE];
        long sequence_offset;
	long sequence_length;
        };

SEQSTRUCT sequence[MAXPROTEINS];
unsigned int **matrix;
unsigned int **matrix2;

int *cluster_assign;
int cluster[MAXPROTEINS];
int cluster_index=0;
int curr_offset=0;
int offset_length=0;
int pos=0;
int verbose=0;
int nosym=0;
int nomd=0;
int exhaustive=0;

/* User Interface Variables - Counters-progress etc. */

int counter=0;
int progress=0;
int progress_percentage=0;


/* Main - Reads in data from BLAST input file, Builds the Matrix and */
/*        calls the symmetrify and cluster sub-routines		     */

int main(int argc, char *argv[])
{

char protein[PROTEINSIZE];
char protein1[PROTEINSIZE];
char protein2[PROTEINSIZE];
char curr_protein[PROTEINSIZE]={'\0'};
char BUFFER[250];
char *keyptr;
char **key;
SEQSTRUCT *ptr;

/* File Handles */
FILE *fp;
FILE *fp2;

int x=0;
int z=0;
int i=0;
int j=0;
int p=0;
int N=0;
int count=0;
int index1=0;
int index2=0;
int found1=0;
int found2=0;


/* ================================= */
/* Greeting And Command Line Parsing */
/* ================================= */

fprintf(stdout,"\n");
fprintf(stdout,"geneRAGE v1.02b - Serial  (c) EMBL-EBI 1999\n");
fprintf(stdout,"=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
fprintf(stdout,"Please Cite:\nAnton J. Enright & Christos A. Ouzounis:\nBioinformatics 2000 16: 451-457\n\n");

if (argc == 1)
	{
	fprintf(stdout,"Usage: rage file1 file2\n");	
	fprintf(stdout,"Where file1 is a FASTA format protein database\n  and file2 is a parsed BLAST results file\n\n");
	exit(1);
	}

        for (i=0;i < argc;i++)
                {
		                if ( (!strcmp(argv[i],"-h")) || (!strcmp(argv[i],"-help")) )
                        {
                        printf("\n\ngeneRAGE can be run as follows:\n\n");
                        printf("\'generage sequences.fasta sequences.parsed\'\n");
printf("\nOutput files will be created in the directory from which the program was called.\n");
printf("\nGeneRAGE Command Line Options\n");
printf("-----------------------------\n");
printf("-help                   Show some helpful information on commandline options\n");
printf("-verbose                Run GeneRAGE in verbose mode (Useful to see whats going on)\n");
printf("-nosym                  Turn off the symmetrification step (Not Recommended)\n");
printf("-nomd                   Turn off the multi-domain detection step (Not Recommended)\n");
printf("-exhaustive             Turn off all optimizations (Slower but more accurate)\n");
printf("-z1                     Smith-Waterman Z-Score Threshold for symmetrification (Default=10)\n");
printf("-z2                     Smith-Waterman Z-Score Threshold for multi-domain detection (Default=2)\n");
			exit(1);
                        }                                                                              
		}

if (! (prssdir=(char *) getenv("PRSSDIR")))
	{
	fprintf(stderr,"Warning: Environment variable PRSSDIR not set, looking for PRSS3 in default path\n");
	}

if ( argc < 3)
	{
	fprintf(stdout,"Error: You need to specify two files for input\n");
	exit(1);
	}

if ( argc > 3 )
	{

	for (i=3;i < argc;i++)
		{


		if ( !strcmp(argv[i],"-verbose") )
			{
			verbose=1;	
			}
	
		if ( !strcmp(argv[i],"-nomd") )
                        {
			nomd=1;
                        }
		if ( !strcmp(argv[i],"-nosym") )
                        {
                        nosym=1;
                        }     

		if ( !strcmp(argv[i],"-exhaustive") )
                        {
                        exhaustive=1;
                        }

		if ( !strcmp(argv[i],"-rand") && (i+1 < argc) )
			{
			
			if (atoi(argv[i+1]))
				{
				randomize=atoi(argv[i+1]);
				}
			else
				{
				printf("Error: One or more command line options incorrect!\n");
				exit(1);
				}
			}

		if ( !strcmp(argv[i],"-z1") && (i+1 < argc) ) 
			{

			if (atoi(argv[i+1]))
                                { 
				zscore1=(atoi(argv[i+1]));
				}
			else
                                {
                                printf("Error: One or more command line options incorrect!\n");
                                exit(1);
                                }
			}

		if ( !strcmp(argv[i],"-z2") && (i+1 < argc) )
			{
			
			if (atoi(argv[i+1]))
                                { 
				zscore2=(atoi(argv[i+1]));
				}
			else
                                {
                                printf("Error: One or more command line options incorrect!\n");
                                exit(1);
                                }
			}

		}

	}

if (  (fp=fopen(argv[1],"r")) == NULL )
	{
	fprintf(stderr,"Error: FASTA Sequences File not found\n");
        exit(1);
	}

if (  (fp2=fopen(argv[2],"r")) == NULL )
        {
        fprintf(stderr,"Error: Protein Hits File not found\n");
        exit(1);
        }

/* ============================ */
/* Setup Database Variable Name */
/* ============================ */

strcpy(database,argv[1]);

/* ================================ */
/* Read in Proteins from input file */
/* ================================ */

fprintf(stdout,"Reading in Proteins from file %s:\t",argv[1]);
while (!feof(fp))
{
curr_offset=ftell(fp);

fgets(BUFFER, 200,fp);
sscanf(BUFFER,">%s", protein);

if (  ((strcmp(protein,curr_protein)) != 0) && (feof(fp)==0) )
        {
	offset_length-= curr_offset;
	sequence[i].sequence_length=abs(offset_length);
	offset_length=curr_offset;
	/* Check for Duplicates */
	for (i=0;i<count;i++)
		{
		if (!strcmp(sequence[i].sequence_id,protein))
			{
			fprintf(stderr,"Error: Duplicated IDs found\n");
			exit(1);
			}
		}


	strcpy(sequence[count].sequence_id,protein);
	sequence[count].sequence_offset=curr_offset;
	/*printf("Sequence[%d]=%s %d\n",count,protein,sequence[count].sequence_offset);*/
        count++;
	if (count > MAXPROTEINS)
		{
		fprintf(stderr,"ERROR: MAXPROTEINS EXCEEDED PLEASE RECOMPILE (%d %d)\n",count,MAXPROTEINS);
		}
        }
strcpy(curr_protein,protein);

}
offset_length-=curr_offset;
sequence[i].sequence_length=abs(offset_length);
fclose(fp);
fprintf(stdout,"Read %d Proteins\n",count);
fflush(stdout);
fprintf(stdout,"\nSorting Input Sequences:\t");
fflush(stdout);
qsort(sequence, count, sizeof (SEQSTRUCT), strptrcmp);
fprintf(stdout,"done\n");
fflush(stdout); 

/* ========================= */
/* Set Dimensions of Matrix  */
/* ========================= */

N=count;
fprintf(stdout,"\n------------------------------------------\n");
fprintf(stdout,"Matrix Size:\t\t%d kb\nCompressed Size:\t%d kb\n",((N*N)*sizeof(int))/1024, ((N*N)*sizeof(int))/(1024*WORDSIZE));
matrix=imatrix(0,N-1,0,((N-1)/WORDSIZE)+1);
matrix2=imatrix(0,N-1,0,((N-1)/WORDSIZE)+1);
fprintf(stdout,"------------------------------------------\n");

/* ============================ */
/* Read in Hits from input file */
/* ============================ */

fprintf(stdout,"\nReading in Hits from file %s:\t",argv[2]);

x=0;

while (!feof(fp2))
{
fgets(BUFFER, 200,fp2);
sscanf(BUFFER,"%s %s", protein1, protein2);

if ( (feof(fp2)==0) )
        {
	if (x > MAXHITS)
		{
		fprintf(stderr,"ERROR: MAXHITS EXCEEDED - PLEASE RECOMPILE (%d %d)\n",x,MAXHITS);
		exit(1);
		}

	index1=0;
        index2=0; 

	keyptr=protein1;
	key=&keyptr;

	ptr=(SEQSTRUCT *)bsearch(key,sequence,count,sizeof(SEQSTRUCT),strptrsearch);
	if (ptr == NULL)
		{
		fprintf(stderr,"Error (1): Sequence %s from hits file not found in sequences file!\n",protein1);
		exit(1);
		}
	index1=((ptr-(SEQSTRUCT *)sequence));

	keyptr=protein2;
        key=&keyptr;
	ptr=(SEQSTRUCT *)bsearch(key,sequence,count,sizeof(SEQSTRUCT),strptrsearch);
	if (ptr == NULL)
                {
                fprintf(stderr,"Error (2): Sequence %s from hits file not found in sequences file!\n",protein2);
                exit(1);
                }
        index2=((ptr-(SEQSTRUCT *)sequence));

	setbit(index1,index2);
	x++;
	}
} 
fclose(fp2);

fprintf(stdout,"Read %d hits\n",x,argv[2]);
totalhits=x;
fflush(stdout);

/* ========================================== */
/* Allocate Memory for Multidomain Stuff      */
/* ========================================== */

multidom_index=(int *) calloc(N,sizeof(int));
multidom=(int *) calloc(N,sizeof(int));
domain1=(int *) calloc(N,sizeof(int));
domain2=(int *) calloc(N,sizeof(int));
clustered_domains=(int *) calloc(N,sizeof(int));
allocated=N;


/* ====================== */
/* Clear Clustering Index */
/* ====================== */

cluster_assign=(int *) calloc( (N) ,sizeof(int));
fflush(stdout);    

/*
dump_matrix(60);
*/

/* ========================================== */
/* Process Matrix and Make Symmetric	      */
/* ========================================== */

if (!nosym)
{
fprintf(stdout,"Symmetrifying Matrix:"); 
symmetrify_matrix(N);
fprintf(stdout,"done\n\n"); 
dump_similarities(N,"symmetric.out");
}

/* ========================================== */
/* Process Matrix and Locate Missing Hits     */
/* ========================================== */

if (!nomd)
{
fprintf (stdout,"Matrix Pass 2 Beginning:");
fflush(stdout);

if (exhaustive)
	{
	multi_domain(N,0);
	}
else
	{
	multi_domain(N,1);
	dump_similarities(N,"diagnostic.out"); 
	begin_cluster(N);
	dump_clusters(N);
	reset_matrix(N);
	multi_domain(N,2);
	}

dump_similarities(N,"similarities.out");
fprintf (stdout,"done\n\n");
fprintf (stdout,"\n");
}

/* ======================================================================== */
/* For Each Protein, Begin a Clustering Operation, if not already clustered */
/* Clusters recursively until no more clustering can be achieved            */
/* 									    */
/*  THIS IMPLEMENTATION IS FOR SINGLE LINKAGE CLUSTERING                    */
/* ======================================================================== */

begin_cluster(N);

/* ======================================================================== */
/* Print out the Cluster Assignments for Each Protein			    */
/* ======================================================================== */

dump_clusters(N);


/* ======================================================================== */
/* Print out the Paralogue Table                                            */
/* ======================================================================== */

dump_paralogues(N);

/* ======================================================================== */
/* Free Memory Allocated to the Matrix and Cluster List                     */
/* ======================================================================== */
/*
free_imatrix(matrix,0,N-1,0,N-1);
free_ivector(cluster_assign,0,N-1);
*/
/* ======================================================================== */
/* Our Work here is done, Exit Nicely with Return Code (0)                  */
/* ======================================================================== */

fflush(stdout);

exit(0);

}



/* ======================================================================================================== */
/*   SUB-ROUTINES ONLY BELOW THIS SECTION    								    */
/* ======================================================================================================== */


/* getbit() setbit() unsetbit()                                         */
/*                                                                      */
/* Set of Routines for fetching, setting and unsetting bits in matrices */
/*                                                                      */

int getbit(int i, int j)
{
return(matrix[i][j/WORDSIZE]&(1<<j%WORDSIZE));
}

void setbit(int i, int j)
{
matrix[i][j/WORDSIZE]|=(1<<j%WORDSIZE);
}

void unsetbit(int i, int j)
{
matrix[i][j/WORDSIZE]^=(1<<j%WORDSIZE);
}

int getbit2(int i, int j)
{
return(matrix2[i][j/WORDSIZE]&(1<<j%WORDSIZE));
}

void setbit2(int i, int j)
{
matrix2[i][j/WORDSIZE]|=(1<<j%WORDSIZE);
}

void unsetbit2(int i, int j)
{
matrix2[i][j/WORDSIZE]^=(1<<j%WORDSIZE);
}

void begin_cluster(int N)
{
int i=0;
int j=0;

  fprintf(stdout,"\nClustering Sequences:");

  if (cluster_index != 0)
	{
	printf("Resetting Clusters....\n");
	memset(cluster_assign,'\0',N*sizeof(int));
	}


  cluster_index=0;
  for (i=0;i<N;i++)
	{

        if (cluster_assign[i]==0)
                {

                for(j=0;j<MAXPROTEINS;j++)
                        {
                        cluster[j]=-1;
                        }

                pos=0;
                cluster_index++;
                cluster_matrix(i,N);
                }

	}
  fprintf(stdout,"done\n\n"); 
}

/* dump_similiarities                                                   */
/*                                                                      */
/* Dumps the Protein content of the Matrix                              */
/*                                                                      */

void dump_similarities(int N,char *simfile)
{
int i=0;
int j=0;

FILE *fp;

if (  (fp=fopen(simfile,"w")) != NULL )
        {
        printf("Outputing Symmetric Binary Relations to file: similarities.out\n");
        }

else
        {
        printf("Error: Cannot Create Output File\n");
        exit(1);
        }

for (i=0;i<N;i++)
        {
        for (j=0;j<N;j++)
                {
                if (getbit(i,j))
                        {
                        fprintf(fp,"%s\t\t%s\n",sequence[i].sequence_id,sequence[j].sequence_id);
                        }
                }
        }

fclose(fp);

}

/* dump_paralogues                                                      */
/*                                                                      */
/* Dumps GATOS Compliant Paralogue Tables (MultiDomain)                 */
/*                                                                      */ 

void dump_paralogues(N)
{

FILE *fp;
int i=0;
int j=0;
int z=0;
int p=0;

/* Array For Keeping Track of Paralogues */
int index[MAXPROTEINS];


if (  (fp=fopen("paralogues.out","w")) != NULL )
        {
        fprintf(stdout,"Outputing Paralogue Table to file: paralogues.out\n");
        }

else
        {
        fprintf(stderr,"Error: Cannot Create Output File (paralogues.out)\n");
        exit(1);
        }

clear_matrix2(N);

/* REFILL MATRIX WITH CLUSTERING DATA */

for (i=0;i<N;i++)
        {
        setbit2(i,cluster_assign[i]-1);
        }

/* Check Each Protein for an MD Hit */
for (i=0;i<N;i++)
	{

		if (multidom_index[i]==1)
                	{	

			/* FIND ITS ENTRY IN THE MULTIDOMAIN LIST */
                	for (j=0;j<mdc;j++)
                        	{
                        	if (multidom[j]==i)
                                	{
                                	clustered_domains[cluster_assign[domain1[j]]]=1;
                                	clustered_domains[cluster_assign[domain2[j]]]=1;
                                	}
                        	}

                	/* Print the Domains, then flush the index */

                	for (j=0;j<N;j++)
                        	{
                        	if (clustered_domains[j]==1)
                                	{
                                	clustered_domains[j]=0;
                               		if (!getbit2(i,j-1))
						{
						setbit2(i,j-1);
						}
                                	}

                        	}

	                }


	}

/* PRINT OUT THE PARALOGUE TABLE */	

for (i=0;i<N;i++)
        {
	fprintf(fp,"%s\t\t",sequence[i].sequence_id);

	/* FLUSH THE INDEX */

	for (j=0;j<N;j++)
		{
		index[j]=0;
		}

        for (j=0;j<N;j++)
                {
                if (getbit2(i,j))
                        {
                        for (z=0;z<N;z++)
                                {
                                if ( (getbit2(z,j)) && (z != i) )
                                        {
					if (index[z]==0)
						{
						if (p==0)
							{
							fprintf(fp,"%s",sequence[z].sequence_id);
							p=1;
							}
						else
							{
                                        		fprintf(fp,":%s",sequence[z].sequence_id);
							}
						}
					/* We Need to keep track of Which Paralogues we dump, so we don't repeat ourselves */   
					index[z]=1;

                                        }
                                }
                        }
                }
        fprintf(fp,"\n");
	p=0;
        } 


}


/* dump_clusters							*/
/*									*/
/* Dumps the content of the cluster-array			   	*/
/* 									*/

void dump_clusters(int N)
{
int i=0;
int j=0;
int z=0;
int domains=1;

/* Temporary Array for storing Domain Assignments */
unsigned int domainslist[MAXPROTEINS];

FILE *fp;

if (  (fp=fopen("clusters.out","w")) != NULL )
        {
        printf("Outputing Clusters to file: clusters.out\n");
        }

else
        {
        printf("Error: Cannot Create Output File\n");
        exit(1);
        }

for (i=0;i<MAXPROTEINS;i++)
        {
        domainslist[i]=0;
        }

for (i=0;i<N;i++)
        {
	if (multidom_index[i] != 1)
		{
		fprintf(fp,"%d\t\t%s\n",cluster_assign[i],sequence[i].sequence_id);
		}

	else
		{
		fprintf(fp,"%d\t\t%s Multi-Domain Family\n",cluster_assign[i],sequence[i].sequence_id);
		}

        if (multidom_index[i] == 1)
                {
                for (j=0;j<mdc;j++)
                        {
                        if (multidom[j] == i)
                                {
                                if (cluster_assign[domain1[j]] != cluster_assign[domain2[j]])
                                        {
                                        if ((multidom_index[domain1[j]] != 1) && (multidom_index[domain2[j]] != 1))
                                                {
						/*
                                                domainslist[cluster_assign[domain1[j]]]=1;
                                                domainslist[cluster_assign[domain2[j]]]=1;
						*/
                                                }

					if ((multidom_index[domain1[j]] != 1))
						{
						domainslist[cluster_assign[domain1[j]]]=1;
						}
					else
						{
						}

					if ((multidom_index[domain2[j]] != 1))
                                                {
						domainslist[cluster_assign[domain2[j]]]=1; 
                                                }
					else
						{
						}
                                        }
                                }
                        }

                domains=1;
                for (j=0;j<MAXPROTEINS;j++)
                        {
                        if ((domainslist[j]==1))
                        /* if ((domainslist[j]==1) && (j != cluster_assign[i]) ) */
                                {
                                fprintf(fp,"%d\t\t%s Domain %d\n",j,sequence[i].sequence_id,domains);
                                domains++;
                                }
	
                        domainslist[j]=0;
                        }

                }


        }

fclose(fp);

}

/* clear_matrix2                                                        */
/*                                                                      */
/* RESETS ALL VALUES TO ZERO					        */
/*                                                                      */
void clear_matrix2(int X)
{
int i;
int j;

for (i=0;i<X;i++)
        {
        for (j=0;j<((X/WORDSIZE)+1);j++)
                {
                matrix2[i][j]=0;
                }
        }
}

/* dump_matrix                                                          */
/*                                                                      */
/* ALLOWS THE MATRIX TO BE DUMPED TO TERM FOR DEBUGGING PURPOSES        */
/*                                                                      */

void dump_matrix(int X)
{
int i;
int j;
int checksum=0;

fprintf(stdout,"\n");
fflush(stdout);
for (i=0;i<X;i++)
        {
        for (j=0;j<X;j++)
                {
                if(getbit(i,j))
			{
			printf("1 ");
			checksum+=(i*j);
			}
		else
			{
			printf("0 ");
			}
                }
	printf("\n");
        }
printf("%d\n",checksum);
}

void dump_matrix2(int X)
{
int i;
int j;
int checksum=0;

fprintf(stdout,"\n");
fflush(stdout);

for (i=0;i<X;i++)
        {
        for (j=0;j<X;j++)
                {
                if(getbit2(i,j))
                        {
                        printf("1 ");
			checksum+=(i*j);
                        }
                else
                        {
                        printf("0 ");
                        }
                }
        printf("\n");
        }
printf("%d\n",checksum);
}

/* symmetrify_matrix 				                              	 */
/*										 */
/* Processes Matrix and Locates Unsymmetric Hits between **TWO** Proteins        */
/*                                                                     		 */

void symmetrify_matrix(int X)
{
int i;
int j;
int funkatizer=100;
int swat=0;

if (X < 100)
	{
	progress=(X/10);
	} 
else
	{
	progress=(X/100);
	}

progress_percentage=0;
counter=0;

for (i=0;i<X;i++)
        {

	if (!verbose)
		{
		counter++;
		if (counter == progress)
        		{
        		progress_percentage++;
        if (progress_percentage == 25 || progress_percentage == 50 || progress_percentage == 75 || progress_percentage == 100 )
                		{
                		fprintf(stdout,"%d",progress_percentage);
                		}
        		fprintf(stdout,".");
       			fflush (stdout);
        		counter=0;
       			}

		}

        for (j=i;j<X;j++)
                {
                if ( ((getbit(i,j)==0) && (getbit(j,i)!=0)) || ((getbit(i,j)!=0) && (getbit(j,i)==0)) )
			{

			if (verbose)
				{
				printf("\nUnsymmetric Hit Between %s %s.... ",sequence[i].sequence_id,sequence[j].sequence_id);
				}
			swat=do_smith_waterman(i,j,database,zscore1);			

		if (verbose)
			{
			if (swat==2)
				{
				printf("Correcting!\n");
				}
			else 
				{
				printf("Z-Score too Low - Removing\n");
				}
			}


			if (!getbit(i,j) && (swat==2))
				{
				setbit(i,j);
				}
			if (!getbit(j,i) && (swat==2))
                                {
				setbit(j,i);
                                }
			if (getbit(i,j) && (swat==0))
                                {
				unsetbit(i,j);
                                }
                        if (getbit(j,i) && (swat==0))
                                {
                                unsetbit(j,i);
				}

			}
                }
        }
}

/* Used by qsort Routine for Sorting */
int strptrcmp (const void *p1,const void *p2)
{
	SEQSTRUCT *a;
        SEQSTRUCT *b;

        a=(SEQSTRUCT *) p1;
        b=(SEQSTRUCT *) p2;
	return( strcmp(a->sequence_id,b->sequence_id));
}

int strptrsearch (const void *p1,const void *p2)
{
        SEQSTRUCT *a;
	char **b;

        a=(SEQSTRUCT *) p2;
	b=(char **) p1;
        return( -strcmp(a->sequence_id,*b));
} 


/* multi_domain                                                                  */
/*                                                                               */
/* Processes Matrix and Locates Inconsistencies and Multi-domain Proteins        */
/*                                                                               */

void multi_domain(int X,int passtype)
{
int i;
int j;
int z=0;
int p=0;
int hits[MAXHITS];
int *hitcount;
int *hit1;
int *hit2;
int *hitresult;
int hit_counter;
int hit_found=0;
int swat=0;                                                                                                                         
            
int flag=0;
hitcount=(int *) calloc(X,sizeof(int));
hit1=(int *) calloc(totalhits,sizeof(int));
hit2=(int *) calloc(totalhits,sizeof(int));
hitresult=(int *) calloc(totalhits,sizeof(int));



if (X < 100)
        {
        progress=(X/10);
        }
else
        {
        progress=(X/100);
        }

progress_percentage=0;
counter=0; 

for (j=0;j < X; j++)
        {
	        if (!verbose)
                {
                counter++;
                if (counter == progress)
                        {
                        progress_percentage++;
        if (progress_percentage == 25 || progress_percentage == 50 || progress_percentage == 75 || progress_percentage == 100 )
                                {
                                fprintf(stdout,"%d",progress_percentage);
                                }
                        fprintf(stdout,".");
                        fflush (stdout);
                        counter=0;
                        }

                }

	memset(hit1,'\0',totalhits*sizeof(int));
	memset(hit2,'\0',totalhits*sizeof(int));
	memset(hitresult,'\0',totalhits*sizeof(int));


        /* For Each Protein (Row) find a list of matching proteins */

        for (i=0;i< X;i++)
                {
                if (getbit(j,i))
                        {
                        hits[z]=i;
                        z++;
                        }
		/* Also Check For Multi-Domain Protein Hits */

		if (getbit2(j,i))
			{
			hits[z]=i;
			z++;
			}
                }

        /* For Each Matching Protein see if it hits all the other matches */

        for (i=0;i < z;i++)
                        {

                        for (p=i;p < z;p++)
                                {
                                if ( !getbit(hits[i],hits[p]) && !getbit2(hits[i],hits[p])) 
                                        {
					if (verbose)
					{
					if (passtype != 1)
						{
                                        	fprintf(stdout,"\nProteins (%d) %s and (%d) %s hit (%d) %s but not each other\n",hits[i],sequence[hits[i]].sequence_id,hits[p],sequence[hits[p]].sequence_id,j,sequence[j].sequence_id);
						}
					}


				if (passtype == 0)
					{
                                        /* Do A Smith-Waterman to check that no sig. similarity exists */ 
					swat=do_smith_waterman(hits[i],hits[p],database,zscore2);
					fflush(stdout);
					}

				if (passtype == 2)
					{
					printf("%s %d %d\n",sequence[j].sequence_id,cluster_assign[hits[i]],cluster_assign[hits[p]])
;
					hit_found=0;
					for (hit_counter=0;hit_counter<hitcount[i];hit_counter++)
						{
						if (hit1[hit_counter] == cluster_assign[hits[i]]) 
							{
							if (hit2[hit_counter] == cluster_assign[hits[p]])
								{
								fprintf(stdout,"Already Done This!  %d\n",hitresult[hit_counter]);
								hit_found=hit_counter+1;
								}	
							}

						if (hit1[hit_counter] == cluster_assign[hits[p]])
                                                        {
                                                        if (hit2[hit_counter] == cluster_assign[hits[i]])
                                                                {
                                                                printf("Already Done This! %d\n",hitresult[hit_counter]);
								hit_found=hit_counter+1;
                                                                }      
                                                        }
						}				

					if (hit_found)
						{
						swat=hitresult[hit_found-1];
						}
					else
						{
						swat=do_smith_waterman(hits[i],hits[p],database,zscore2);
						hitresult[hitcount[i]]=swat;
						hit1[hitcount[i]]=cluster_assign[hits[i]];
						hit2[hitcount[i]]=cluster_assign[hits[p]];
						hitcount[i]++;
						}
                                        fflush(stdout);
					}

				if (passtype == 1)
					{
					swat=0;
					}

                                        if (swat != 0)
                                                {

                                                if (!getbit(hits[i],hits[p]))
                                                        {
                                                        setbit(hits[i],hits[p]);
                                                        }
                                                if (!getbit(hits[p],hits[i]))
                                                        {
                                                        setbit(hits[p],hits[i]);
                                                        }
                                                }

                                        if (swat == 0)
                                                {
						if (verbose)
							{
							if (passtype != 1)
								{
                                              		fprintf(stdout,"\n *** Possible Multi-domain Protein %s :Domains %s %s***\n"
,sequence[j].sequence_id,sequence[hits[i]].sequence_id,sequence[hits[p]].sequence_id); 
								}
							}

                                                if (getbit(hits[i],j))
                                                        {
                                                        unsetbit(hits[i],j);
                                                        }
                                                if (getbit(j,hits[i]))
                                                        {
                                                        unsetbit(j,hits[i]);
                                                        }
                                                if (getbit(hits[p],j))
                                                        {
                                                        unsetbit(hits[p],j);
                                                        }
                                                if (getbit(j,hits[p]))
                                                        {
                                                        unsetbit(j,hits[p]);
 							}
					
                                                if (!getbit2(hits[i],j))
                                                        {
                                                        setbit2(hits[i],j);
                                                        }
                                                if (!getbit2(j,hits[i]))
                                                        {
                                                        setbit2(j,hits[i]);
                                                        }
                                                if (!getbit2(hits[p],j))
                                                        {
                                                        setbit2(hits[p],j);
                                                        }
                                                if (!getbit2(j,hits[p]))
                                                        {
                                                        setbit2(j,hits[p]);
                                                        }
		

                                                /* Add this to the MultiDomain Indexes */
						/*
                             		 	printf ("multidom_index[%d]:%s set, multidom[%d]=%d, domain1[%d]=%s, domain2[%d]=%s\
n",j,sequence[j].sequence_id,mdc,j,mdc,sequence[hits[i]].sequence_id,mdc,sequence[hits[p]].sequence_id);
						*/

						if ((passtype == 0) || (passtype == 2))
							{
							if (mdc >= allocated)
								{
								multidom=(int *) realloc(multidom,(allocated+MAXHITS)*sizeof(int));
								memset(multidom+(allocated),'\0',MAXHITS*sizeof(int));
								domain1=(int *) realloc(domain1,(allocated+MAXHITS)*sizeof(int)); 
								memset(domain1+(allocated),'\0',MAXHITS*sizeof(int));
								domain2=(int *) realloc(domain2,(allocated+MAXHITS)*sizeof(int));
								memset(domain2+(allocated),'\0',MAXHITS*sizeof(int));
								allocated=allocated+MAXHITS;
								}

                                                	multidom_index[j]=1;
                                                	multidom[mdc]=j;
                                                	domain1[mdc]=hits[i];
                                                	domain2[mdc]=hits[p];
                                                	mdc++;
							}


                                                }




                                        }
                                }
                        }                         /* Re-Zero the Hit List */

        for (i=0; i < z; i++)
                {
                hits[i]=0;
                }

        z=0;

        }   

}

/* reset_matrix									*/
/* 										*/
/* Moves bits from MD matrix back to original positions in similarity Matrix	*/
/*										*/

void reset_matrix(int X)
{
int i=0;
int j=0;

for (i=0;i<X;i++)
	{
	for (j=0;j<X;j++)
		{
		if (getbit2(i,j))
			{
			unsetbit2(i,j);
			setbit(i,j);
			}
		}
	}


}

/* do_smith_waterman 								*/
/*										*/
/* Calls PRSS3 to do a Smith-Waterman Between Two Protein Sequences		*/
/* Executes the SW_CLUMP script(perl5) wrapper for a Randomized Ssearch 	*/
/*										*/

int do_smith_waterman(int query1,int query2, char *database, int zscore)
{
char command_string[100];
char BUFFER[MAXPEPTIDE]={'\0'};

int shuffs=0;
double mu=0;
double variance=0;
double Zscore=0;
double stdev=0;
int score=0;

int swreturn=1;
int search_complete=0;
int indexer=0;
int record1=0;
int record2=0;
int temp=0;

int index1=0;
int index2=0;
int i=0;

char zscorecut[10]="";
char random[10]="";
char substring[1];

FILE *filein;
FILE *fileout1;
FILE *fileout2;
FILE *fileres;

/* First Work out the Line Numbers for each Input Sequence */

index1=sequence[query1].sequence_offset;
index2=sequence[query2].sequence_offset;

if (  (filein=fopen(database,"r")) != NULL )
        {
        }
else
        {
        fprintf(stderr,"Error: FASTA formatted sequence file containing query sequences not found\n");
        exit(1);
        }

if (  (fileout1=fopen("rage1.temp","w")) != NULL )
        {
        }
else
        {
        fprintf(stderr,"Error: Cannot create temporary output file for Smith-Waterman Analysis (1) \n");
        exit(1);
        }

if (  (fileout2=fopen("rage2.temp","w")) != NULL )
        {
        }
else
        {
        fprintf(stderr,"Error: Cannot create temporary output file for Smith-Waterman Analysis (2) \n");
        if (errno)
                { 
		perror("Error returned from system call");
                errno=0;
                }

        exit(1);
        }

        /* Use Direct File I/O and the Index Offset Array to quickly Fetch Sequences */

        fseek(filein,index1,SEEK_SET);
        fread(BUFFER, (sequence[query1].sequence_length),1,filein);
        fprintf(fileout1,"%s",BUFFER);
        fclose(fileout1);

        /* EMPTY THE BUFFER */

        for (i=0;i<MAXPEPTIDE;i++)
                {
                BUFFER[i]=0;
                }

        fseek(filein,index2,SEEK_SET);
        fread(BUFFER, sequence[query2].sequence_length,1,filein);
        fprintf(fileout2,"%s",BUFFER);
        fclose(fileout2);
        fclose(filein);

	 /* EMPTY THE BUFFER */

        for (i=0;i<MAXPEPTIDE;i++)
                {
                BUFFER[i]=0;
                } 

/* Prepare a Command String to Call PRSS3 */

sprintf(random,"%d", randomize);
strcpy(command_string,"");
if (prssdir != NULL)
	{
	strcat(command_string,prssdir);
	strcat(command_string,"/prss33 -H -d ");
	}
else
	{
	strcat(command_string,"prss33 -H -d ");	
	}
strcat(command_string,random);
strcat(command_string," -q -z 0 rage1.temp rage2.temp > prss.out");
strcat(command_string," ");

if (verbose)
        {
        fprintf(stdout,"Running Smith-Waterman\n");
        fflush(stdout);
        }

swreturn=(system(command_string));

if (swreturn != 0)
        {
        fprintf(stderr,"Error: Smith-Waterman analysis failed\n");
        exit(1);
        }

if (errno)
        {
        perror("Error returned from system call");
        errno=0;
	        }


/* OPEN THE PRSS3 RESULTS FILE AND CALCULATE Z-SCORE */

if (  (fileres=fopen("prss.out","r")) != NULL )
        {
        }
else
        {
        fprintf(stderr,"Error: Cannot open output file from Smith-Waterman Analysis\n");
        exit(1);
        }

while ( (!feof(fileres)) )
        {
        fgets(BUFFER, 200,fileres);
        //sscanf(BUFFER," unscaled statistics: mu= %lf  var=%lf",&mu,&variance);
        //sscanf(BUFFER," unshuffled s-w score: %d;",&score);
        sscanf(BUFFER, "Statistics: (shuffled [%d]) Unscaled statistics: mu= %lf var=%lf", &shuffs, &mu, &variance);
        sscanf(BUFFER, "Smith-Waterman score: %d;", &score);
        }

/* CLOSE Results File */
fclose(fileres);

/* Calculate Z-Score */
stdev=sqrt(variance);
Zscore=((score-mu) / stdev); 

if (verbose)
        {
        fprintf(stdout,"mu=%lf variance=%lf score=%d stdev=%lf Z-Score=%lf\n",mu,variance,score,stdev,Zscore);
        fflush(stdout);
        }

if (Zscore > zscore)
        {
        swreturn=2;
        }

/* Remove Temporary Files */
remove("rage1.temp");
remove("rage2.temp");
remove("prss.out");

return swreturn;

}

int cluster_matrix(int X,int N)
{
int i=0;
int j=0;
int found=0;

printf(".");


/* 
printf("pos =%d\n",pos);
printf("Adding %s to pos[%d] Cluster:%d\n",sequence[X]sequence_id,pos,cluster_index);
*/

cluster[pos]=X;
cluster_assign[X]=cluster_index;
pos++;

/*
printf("Clustering for Protein %s\n",sequence[X].sequence_id);
printf("protein %s hits:",sequence[X].sequence_id);

for (i=0;i<N;i++)
        {
	if ( getbit(X,i))
                {
		printf("%s ",sequence[i].sequence_id); 
		}
	}

printf("\n"); 

*/


for (i=0;i<N;i++)
	{
	
	if ( getbit(X,i))
		{
		/* printf("Hits to %s\n",sequence[i].sequence_id); */
		found=0;		
		for(j=0;j<pos;j++)
			{
			if (cluster[j] == i)
				{
				found=1;
				}

			}
			if (found != 1)
				{
	/*			printf("Calling cluster_matrix(%d,%d)\n",i,N); */
				cluster_matrix(i,N);
				}
		}

	}

return(0);

}

unsigned int **imatrix(long nrl, long nrh, long ncl, long nch)
{
#define NR_END 1 

        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        unsigned int **m;

        /* allocate pointers to rows */

        m=(unsigned int **) calloc((size_t)(nrow+NR_END),sizeof(unsigned int*));
        if (!m) 
		{
		printf("allocation failure 1 in matrix()");
		}

        m += NR_END;
        m -= nrl;


        /* allocate rows and set pointers to them */

        m[nrl]=(unsigned int *) calloc((size_t)(nrow*ncol+NR_END),sizeof(unsigned int));

        if (!m[nrl]) 
		{	
		printf("allocation failure 2 in matrix()");
		}

        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) 
		{
		m[i]=m[i-1]+ncol;
		}

        /* return pointer to array of pointers to rows */
        return m;
}

int *ivector(int nh)
{
int *v;

printf("ivector called for %d\n",nh); 
v=(int *) calloc((nh+1),sizeof(int));

if (!v)
	{
	printf("allocation failure in vector()");
	return NULL;
	}

return v;
}

void free_ivector(int *v, int nl, int nh)
{
free((char*) (v+nl-1));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
free((unsigned int*) (m[nrl]+ncl-1));
free((unsigned int*) (m+nrl-1));
}
