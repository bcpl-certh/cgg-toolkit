/* ================================================================== */
/* diffuse.c - Program for detection of gene fusions events between   */
/* genomes                                                            */
/*                                                                    */ 
/* 		SERIAL VERSION 1.0a                                   */
/* 								      */
/* Anton Enright  - Computational Genomics Group (C. Ouzounis)        */
/*								      */
/* EMBL - European Bioinformatics Institute 			      */
/* ================================================================== */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<errno.h>

/* DEFINITIONS - MAXIMUMS, DEFAULTS, ARCHITECTURE WORD SIZE        */

#define MAXPROTEINS 200000 /* Max Proteins in Database		   */
#define MAXHITS 1000000	  /* Max Hits 				   */
#define PROTEINSIZE 60    /* Max size of Protein Names (chars)     */
#define WORDSIZE 32       /* Word Size on This Architecture        */
#define MAXPEPTIDE 25000  /* Max Size of a peptide in Amino Acids  */

/* Function Prototypes */

void dump_matrix(int);
void dump_matrix2(int);  
void symmetrify_matrix(int);
void detect_fusions(int, int);
void dump_similarities(int);

int do_smith_waterman(int, int,char *, int);

/* Smith Waterman Stuff  - Default Cutoffs / Randomizations */

int randomize = 100;
int zscore1 = 10;
int zscore2 = 10;


/* Memory Allocation Routines - Prototypes */

unsigned int **imatrix(long,long,long,long);
int **matrix(int, int);
void free_imatrix(unsigned int **,int,int,int,int);

/* BitMatrix Stuff */

int getbit(int, int);
void setbit(int, int);
void unsetbit(int, int);
int getbit2(int, int);
void setbit2(int, int);
void unsetbit2(int, int); 

/* New Variable to store SWAT results */
int **swatscores;

/* Global Variable Declarations */

char query_database[100];
char reference_database[100];

char sequence1[MAXPROTEINS][PROTEINSIZE];
char sequence2[MAXPROTEINS][PROTEINSIZE];

long seq_index1[MAXPROTEINS];
long seq_index2[MAXPROTEINS];

unsigned int **query_query_matrix;
unsigned int **query_ref_matrix;

int verbose=0;
int symmetrify=1;
int outfile=0;
int symout=0;
int interspecies=0;
int storeswat=0;
int substrindex=0;

/* User Interface Variables - Counters-progress etc. */

int counter=0;
int total_fusion_count=0;
int progress=0;
int progress_percentage=0;
FILE *fileout;
FILE *symmetric_out;

/* Main - Reads in data from BLAST input file, Builds the Matrix and */
/*        calls the symmetrify and cluster sub-routines		     */

int main(int argc, char *argv[])
{

char protein1[PROTEINSIZE];
char protein2[PROTEINSIZE];
char curr_protein[PROTEINSIZE];
char BUFFER[250];

FILE *fp;
FILE *fp2;
FILE *fp3;
FILE *fp4;

int x=0;
int z=0;
int i=0;
int j=0;
int p=0;
int N=0;
int count=0;
int sequence1_count=0;
int sequence2_count=0;
int index1=0;
int index2=0;
int found1=0;
int found2=0;

long curr_offset=0;

/* ================================= */
/* Greeting And Command Line Parsing */
/* ================================= */

fprintf(stdout,"\n");
fprintf(stdout,"difFUSE v1.0 - Serial  (c) EMBL-EBI 1999\n");
fprintf(stdout,"=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");
fprintf(stdout,"Anton Enright: enright@ebi.ac.uk	\n\n");

if ( argc == 1)
	{
	fprintf(stdout,"Usage: diffuse file1 file2 file3 file4\n");
	fprintf(stdout,"\nFor help: diffuse -help\n");
	exit(0);
	}

if ( !strcmp(argv[1],"-help")) 
	{
	fprintf(stdout,"Usage: diffuse file1 file2 file3 file4\n\n");	
	fprintf(stdout,"Where:\nfile1 is a FASTA file containing the query sequences\n");
	fprintf(stdout,"file2 is a FASTA file containing the reference sequences\n");
	fprintf(stdout,"file3 is a parsed similarity file for the query vs query similarities\n");
	fprintf(stdout,"file4 is a parsed similarity file fot the query vs reference similarities\n");
	fprintf(stdout,"\nOptional Parameters\n");
	fprintf(stdout,"===================\n");
	fprintf(stdout,"-z1 N\t\t Set symmetrification Smith-Waterman cut-off to N (Default=10)\n");
	fprintf(stdout,"-z2 N\t\t Set fusion detection Smith-Waterman cut-off to N (Default=10)\n");
	fprintf(stdout,"-rand X\t\t Do X randomizations during Smith-Waterman (Default=100)\n");
	fprintf(stdout,"-symout file\t\t Only symmetrify query matrix to dump to file\n");
	fprintf(stdout,"-nosym\t\t Do not symmetrify the query matrix\n");
	fprintf(stdout,"-verbose\t Verbose Output\n");
	fprintf(stdout,"-outfile file\t Output results to this file (non-verbose)\n");
	fprintf(stdout,"-interspecies\t Use Interspecies fusion mode\n");
	fprintf(stdout,"-storeswat\t Never re-run the same SWAT (needs more memory)\n");
	fprintf(stdout,"\n");
	exit(0);
	}

if ( argc < 5)
	{
	fprintf(stdout,"Error: You need to specify at least four files for input\n");
	exit(1);
	}

if ( argc > 5 )
	{

	for (i=5;i < argc;i++)
		{
		if ( !strcmp(argv[i],"-verbose") && (!outfile) )
			{
			verbose=1;	
			}

		if ( !strcmp(argv[i],"-nosym") )
			{
			symmetrify=0;
			}

		if ( !strcmp(argv[i],"-interspecies") )
                        {
			interspecies=1;
			}

		if ( !strcmp(argv[i],"-storeswat") )
                        {
			storeswat=1;
			}


		if ( !strcmp(argv[i],"-rand") && (i+1 < argc) )
			{
			
			if (atoi(argv[i+1]))
				{
				randomize=atoi(argv[i+1]);
				}
			else
				{
				fprintf(stderr,"Error: One or more command line options incorrect!\n");
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
                                fprintf(stderr,"Error: One or more command line options incorrect!\n");
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
                                fprintf(stderr,"Error: One or more command line options incorrect!\n");
                                exit(1);
                                }
			}

		if ( !strcmp(argv[i],"-outfile") && (i+1 < argc) )
                        {

                        if ((fileout=fopen(argv[i+1],"w")) != NULL )
                                {
				verbose=0;
				outfile=1;	
                                }
                        else
                                {
                                fprintf(stderr,"Error: Output File cannot be created\n");
                                exit(1);
                                }
                        }

		if ( !strcmp(argv[i],"-symout") && (i+1 < argc) )
                        {

                        if ((symmetric_out=fopen(argv[i+1],"w")) != NULL )
                                {
				symout=1;
				fprintf(stdout,"Symmetric Query Hits will be stored in %s\n",argv[i+1]);
				fflush(stdout);
                                }
                        else
                                {
                                fprintf(stderr,"Error: Output File for symmetric matrix cannot be created\n");
                                exit(1);
                                }
                        }
		

		}

	}


if (  (fp=fopen(argv[1],"r")) != NULL )
        {
        /* printf("Successful\n"); */
        }

else
        {
        fprintf(stderr,"Error: FASTA formatted sequence file containing query sequences not found\n");
        exit(1);
        }

if (  (fp2=fopen(argv[2],"r")) != NULL )
        {
        /* printf("Successful\n"); */
        }

else
        {
	fprintf(stderr,"Error: FASTA formatted sequence file containing reference sequences not found\n");
        exit(1);
        }

if (  (fp3=fopen(argv[3],"r")) != NULL )
        {
        /* printf("Successful\n"); */
        }

else
        {
        fprintf(stderr,"Error: Parsed Blast File containing Query vs Query hits not found: %s\n",argv[3]);
        exit(1);
        }

if (  (fp4=fopen(argv[4],"r")) != NULL )
        {
        /* printf("Successful\n"); */
        }

else
        {
        fprintf(stderr,"Error: Parsed Blast Containing Query vs Reference Hits not found: %s\n",argv[4]);
        exit(1);
        }

/* READ IN QUERY PROTEINS */

fprintf(stdout,"\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
fprintf(stdout,"Reading in Query Proteins from file %s........",argv[1]);
fflush(stdout);

/* Set Current Protein to Nothing and Zero the line counter*/
strcpy(curr_protein,"");
i=0;

while (!feof(fp))
{

/* Store File offset Value */
curr_offset=ftell(fp);

fgets(BUFFER, 200,fp);
sscanf(BUFFER,">%s", protein1);

if (  ((strcmp(protein1,curr_protein)) != 0) && (feof(fp)==0))
        {
        strcpy(sequence1[count],protein1);
	seq_index1[count]=curr_offset;
        count++;
        }

strcpy(curr_protein,protein1);
i++;
}
fclose(fp);
sequence1_count=count;


/* READ IN REFERENCE PROTEINS FROM FILE */

fprintf(stdout,"Read %d Proteins\n",sequence1_count);
fflush(stdout);

/* ReZero the Counters */
count=0;
i=0;

fprintf(stdout,"Reading in Reference Proteins from file %s........",argv[2]);

while (!feof(fp2))
{
fgets(BUFFER, 200,fp2);
sscanf(BUFFER,">%s", protein1);

if (  ((strcmp(protein1,curr_protein)) != 0) && (feof(fp2)==0) )
        {
        strcpy(sequence2[count],protein1);
	seq_index2[count]=ftell(fp2);
        count++;
        }

strcpy(curr_protein,protein1);

i++;
}
fclose(fp2);

sequence2_count=count;

fprintf(stdout,"Read %d Proteins\n",sequence2_count);
fflush(stdout);

count=0; 


/* ============================ */
/* Setup Database Variable Name */
/* ============================ */

strcpy(query_database,argv[1]);
strcpy(reference_database,argv[2]);

fprintf(stdout,"Query Database is %s\n",query_database);
fprintf(stdout,"Reference Database is %s\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n",reference_database);

/* ========================= */
/* Set Dimensions of Matrix  */
/* ========================= */

fprintf(stdout,"\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
fprintf(stdout,"Mallocing %d kb for QxQ Matrix,\n",((sequence1_count * sequence1_count)*sizeof(int))/WORDSIZE);
fprintf(stdout,"Mallocing %d kb for QxR Matrix......",((sequence1_count * sequence2_count)*sizeof(int))/WORDSIZE);
query_query_matrix=imatrix(0,sequence1_count-1,0,((sequence1_count-1)/WORDSIZE)+1);
query_ref_matrix=imatrix(0,sequence1_count-1,0,((sequence2_count-1)/WORDSIZE)+1);
fprintf(stdout,"done\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");

/* ========================================= */
/* Set All Values in Lookup Matrices to Zero */
/* ========================================= */ 

for (i=0;i<sequence1_count;i++)
        {
        for (j=0;j<((sequence1_count/WORDSIZE)+1);j++)
                {
                query_query_matrix[i][j]=0;
                }

        }

for (i=0;i<sequence1_count;i++)
        {
        for (j=0;j<((sequence2_count/WORDSIZE)+1);j++)
                {
                query_ref_matrix[i][j]=0;
                }

        }



/* ============================================ */
/* Read in Query vs Query Hits from input files */
/* ============================================ */

fprintf(stdout,"\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n"); 
fprintf(stdout,"Reading in similarity information:\n",argv[3]);
fflush(stdout);

while (!feof(fp3))
        {
        fgets(BUFFER, 200,fp3);
        sscanf(BUFFER,"%s %s", protein1, protein2);
        found1=found2=0;
        index1=index2=0;
        for (i=0;i<sequence1_count;i++)
                {
                if (  ((strcmp(protein1,sequence1[i])) == 0) && found1==0 )
                        {
                        index1=i;
                        found1=1;
                        }

                if (  ((strcmp(protein2,sequence1[i])) == 0) && found2==0 )
                        { 
			index2=i;
                        found2=1;
                        }
                }

        if (found1==0 || found2==0)
                {
                fprintf(stderr,"Error: One or More Sequences Not Found in File %s\n",argv[3]);

		if (found1 ==0)
                        {
                        fprintf(stderr,"Sequence: %s\n",protein1);
                        }
                if (found2 ==0)
                        {
                        fprintf(stderr,"Sequence: %s\n",protein2);
                        } 


                exit(1);
                }
        count++;
	setbit(index1,index2);
        }

fprintf(stdout,"Read %d Query vs Query Hits\n",count - 1);
fflush(stdout);

/* ========================================= */
/* Read in Query vs Reference Hits from file */
/* ========================================= */ 

fflush(stdout);

count=0;

while (!feof(fp4))
        {
        fgets(BUFFER, 200,fp4);
        sscanf(BUFFER,"%s %s", protein1, protein2);

        /* FIND SEQUENCE NUMBERS */

        found1=found2=0;
        index1=index2=0;
        for (i=0;i<sequence1_count;i++)
                {
                if (  ((strcmp(protein1,sequence1[i])) == 0) && found1==0 )
                        {
                        index1=i;
                        found1=1;
                        }
                }

        for (i=0;i<sequence2_count;i++)
                {
                if (   ((strcmp(protein2,sequence2[i])) == 0) && found2==0 )
                        {
                        found2=1;
                        index2=i;
			}
                }

        if ( (found1 == 0) || (found2 ==0) )
                {
                fprintf(stderr,"Error: 1 or more sequences not found in file %s->",argv[4]);

		if (found1 ==0)
			{
			fprintf(stderr,"Sequence: %s\n",protein1);
			}
		if (found2 ==0)
			{
			fprintf(stderr,"Sequence: %s\n",protein1);
			}

                exit(0);
                }

	setbit2(index1,index2);
        count++;
        }

fprintf(stdout,"Read %d Query vs Reference Hits\n",count - 1 );
fprintf(stdout,"=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n"); 
fflush(stdout);

/* ========================================== */
/* Process Matrix and Make Symmetric	      */
/* ========================================== */

if (interspecies)
	{
	for (i=0;i<sequence1_count;i++)
		{
		for (j=0;j<sequence1_count;j++)
			{
			substrindex=0;
			for (z=0;z<4;z++)
				{
				if (sequence1[i][z] == sequence1[j][z])
					{
					substrindex++;
					}
				}
			if (substrindex == 4)
				{
				setbit(i,j);
				}
		

			}
		}
	}


/* TEST CODE FOR STORING SWAT RESULTS */
if (storeswat)
	{
	printf("Allocating Memory for storage of SWAT results\n");
	swatscores=(int **)matrix(sequence1_count,sequence1_count);
	}

if (symmetrify)
	{
	fprintf(stdout,"\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
	fprintf(stdout,"Symmetrifying Query Matrix:\n");
	symmetrify_matrix(sequence1_count);
	fprintf(stdout,"done\n");
	fprintf(stdout,"=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");

	if (symout)
		{
		dump_similarities(sequence1_count);
		exit(0);
		}
	}

fflush(stdout);

/* ========================================== */
/* Process Matrix and Locate Missing Hits     */
/* ========================================== */

fprintf(stdout,"\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
fprintf(stdout,"Detecting Fusion Events:\n\n");
detect_fusions(sequence1_count,sequence2_count);
fprintf(stdout,"=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");

/* ======================================================================== */
/* Free Memory Allocated to the Matrix and Cluster List                     */
/* ======================================================================== */

free_imatrix(query_query_matrix,0,sequence1_count-1,0,((sequence1_count-1)/WORDSIZE)+1);
free_imatrix(query_ref_matrix,0,sequence1_count-1,0,((sequence2_count-1)/WORDSIZE)+1); 


/* ======================================================================== */
/* Our Work here is done, Exit Nicely with Return Code (0)                  */
/* ======================================================================== */

fprintf(stdout,"\nRun Complete - Have a nice day\n");
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
return(query_query_matrix[i][j/WORDSIZE]&(1<<j%WORDSIZE));
}

void setbit(int i, int j)
{
query_query_matrix[i][j/WORDSIZE]|=(1<<j%WORDSIZE);
}

void unsetbit(int i, int j)
{
query_query_matrix[i][j/WORDSIZE]^=(1<<j%WORDSIZE);
}

int getbit2(int i, int j)
{
return(query_ref_matrix[i][j/WORDSIZE]&(1<<j%WORDSIZE));
}

void setbit2(int i, int j)
{
query_ref_matrix[i][j/WORDSIZE]|=(1<<j%WORDSIZE);
}

void unsetbit2(int i, int j)
{
query_ref_matrix[i][j/WORDSIZE]^=(1<<j%WORDSIZE);
}

/* dump_matrix                                                          */
/*                                                                      */
/* ALLOWS THE MATRIX TO BE DUMPED TO TERM FOR DEBUGGING PURPOSES        */
/*                                                                      */

void dump_matrix(int X)
{
int i;
int j;

fflush(stdout);
for (i=0;i<X;i++)
        {
        for (j=0;j<X;j++)
                {
                if(getbit(i,j))
			{
			printf("1 ");
			}
		else
			{
			printf("0 ");
			}
                }
	printf("\n");
        }
}

void dump_matrix2(int X)
{
int i;
int j;

for (i=0;i<X;i++)
        {
        for (j=0;j<X;j++)
                {
                if(getbit2(i,j))
                        {
                        printf("1 ");
                        }
                else
                        {
                        printf("0 ");
                        }
                }
        printf("\n");
        }
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
                		fprintf(stdout,"%d%%",progress_percentage);
				fflush(stdout);
                		}
			else
				{
					if (progress_percentage < 100)
						{
        					fprintf(stdout,".");
						}
				}

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
				printf("\nUnsymmetric Hit Between %s %s.... ",sequence1[i],sequence1[j]);
				}

			swat=0;
			swat=do_smith_waterman(i,j,query_database,zscore1);			

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

/* detect_fusions()                                                              */
/*                                                                               */
/* Processes Matrix and Locates Possible Fusion Events                           */
/*                                                                               */

void detect_fusions(int X, int Y)
{
int i;
int j;
int z=0;
int p=0;
int hits[MAXHITS];
int swat=0;                                                                                                                                     


if (Y < 100)
        {
        progress=(Y/10);
        }
else
        {
        progress=(Y/100);
        }

progress_percentage=0;
counter=0;

for (j=0;j < Y; j++)
        {

		if (outfile)
                {
                counter++;
                if (counter == progress)
                        {
                        progress_percentage++;
        if (progress_percentage == 25 || progress_percentage == 50 || progress_percentage == 75 || progress_percentage == 100 )
                                {
                                fprintf(stdout,"%d%%",progress_percentage);
                                }

			else 
				{
				if (progress_percentage < 100)
                                                {
                                                fprintf(stdout,".");
                                                }

				}

                        fflush (stdout);
                        counter=0;
                        }

                }


        /* For Each Query Protein (column) find a list of matching proteins */

        for (i=0;i< X;i++)
                {
                if (getbit2(i,j))
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
                                if (!getbit(hits[i],hits[p]))
                                        {

					fflush(stdout);
						if (verbose)
						{
						fprintf(stdout,"-----------------------------------------------------------------------------------\n");
                       		                fprintf(stdout,"proteins (%d) %s and (%d) %s hit (%d) %s but not each other\n",hits[i],sequence1[hits[i]],hits[p],sequence1[hits[p]],j,sequence2[j]);
						}

					if (storeswat)
						{
						printf("Checking: %d %d=%d\n",hits[i],hits[p],swatscores[hits[i]][hits[p]]);
						if (swatscores[hits[i]][hits[p]]==0)
							{
							swat=do_smith_waterman(hits[i],hits[p],query_database,zscore2);
							swatscores[hits[i]][hits[p]]=swat+1;
							}
						else
							{
							swat=swatscores[hits[i]][hits[p]]-1;
							printf("Already done this Smith-Waterman - taking previous result\n");
							}
						}
					else
						{
            					swat=do_smith_waterman(hits[i],hits[p],query_database,zscore2);
						}

					 if (swat == 0)
                                		{
						if (!outfile)
						{	
                                		fprintf(stdout,"Possible fusion event: of %s with %s to %s\n", sequence1[hits[i]], sequence1[hits[p]],sequence2[j]);
						}
				
						if (outfile)
							{
							fprintf(fileout,"Possible fusion event: of %s with %s to %s\n", sequence1[hits[i]],sequence1[hits[p]],sequence2[j]);
							fflush(fileout);
							}


						total_fusion_count++;

						fflush(stdout);
                                		}    

                                        }
                                }
                        } 
	
        for (i=0; i < z; i++)
                {
                hits[i]=0;
                }

        z=0;

        }   

fprintf(stdout,"\n\nTotal Fusions Detected: %d\n",total_fusion_count);

// WARNING - Disabling it temporarily until a better one is found
//fflush(fileout);
//fclose(fileout); 

}

/* do_smith_waterman 								*/
/*										*/
/* Calls PRSS3 to do a Smith-Waterman Between Two Protein Sequences		*/
/* Executes the SW_CLUMP script(perl5) wrapper for a Randomized Ssearch 	*/
/*										*/

int do_smith_waterman(int query1,int query2, char *database, int zscore)
{
char command_string[100];
char BUFFER[MAXPEPTIDE];

double shuffs=0;
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
index1=seq_index1[query1];
index2=seq_index1[query2];

if (  (filein=fopen(database,"r")) != NULL )
	{
	}
else
	{
	fprintf(stderr,"Error: FASTA formatted sequence file containing query sequences not found\n");
        exit(1);
	}

if (  (fileout1=fopen("diffuse1.temp","w")) != NULL )
        {
        }
else
        {
        fprintf(stderr,"Error: Cannot create temporary output file for Smith-Waterman Analysis (1) \n");
        exit(1);
        }

if (  (fileout2=fopen("diffuse2.temp","w")) != NULL )
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
        fread(&BUFFER, (seq_index1[query1+1]-seq_index1[query1]),1,filein);
        fprintf(fileout1,"%s",BUFFER);
	fclose(fileout1);

	/* EMPTY THE BUFFER */

	for (i=0;i<MAXPEPTIDE;i++)
		{
		BUFFER[i]=0;
		}

        fseek(filein,index2,SEEK_SET);
        fread(&BUFFER, seq_index1[query2+1]-seq_index1[query2],1,filein);
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
strcat(command_string,"./ssearch35 -E 10000 -H -B -d 1 -z 21 -Z ");
strcat(command_string,random);
strcat(command_string, " -k ");
strcat(command_string, random);
strcat(command_string," diffuse1.temp diffuse2.temp > prss.out");
strcat(command_string," ");
printf("\n%s\n", command_string);

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
remove("diffuse1.temp");
remove("diffuse2.temp");
remove("prss.out");

return swreturn;
}

unsigned int **imatrix(long nrl, long nrh, long ncl, long nch)
{
#define NR_END 1 

        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        unsigned int **m;

        /* allocate pointers to rows */

        m=(unsigned int **) malloc((size_t)((nrow+NR_END)*sizeof(unsigned int*)));
        if (!m) 
		{
		printf("allocation failure 1 in matrix()");
		}

        m += NR_END;
        m -= nrl;


        /* allocate rows and set pointers to them */

        m[nrl]=(unsigned int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(unsigned int)));

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


int **matrix(int nrh,int nch)
{
int i=0;
int **m;

m=(int **) calloc (nrh,sizeof(int *));
if (!m)
	{
	printf("Error in memory allocation!\n");
	exit(1);
	}

for (i=0;i<nrh;i++)
	{
	m[i]= (int *) calloc (nch,sizeof(int));
	if (!m[i])
		{
		printf("Error in memory allocation!\n");
		exit(1);
		}
	}
return m;
}

void dump_similarities(int N)
{
int i=0;
int j=0;

for (i=0;i<N;i++)
        {
        for (j=0;j<N;j++)
                {
                if (getbit(i,j))
                        {
                        fprintf(symmetric_out,"%s\t%s\n",sequence1[i],sequence1[j]);
                        }
                }
        }

fclose(symmetric_out);

}

void free_imatrix(unsigned int **m,int nrl,int nrh,int ncl,int nch)
{
free((unsigned int*) (m[nrl]+ncl-1));
free((unsigned int*) (m+nrl-1));
}
