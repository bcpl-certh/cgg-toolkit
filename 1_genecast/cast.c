/*
******************************************************************
CAST an algorithm that masks composition biased regions in 
protein sequences 

Authors: Vassilis Promponas, Anton Enright, David Kreil.
Please Cite:	
 	Promponas V.J., Enright A.J., Tsoka S., Kreil D., Leroy C., 
	Hamodrakas S., Sander C., Ouzounis C.A.; 
	Bioinformatics 16(10) 915-922 (2000).

Copyright (C) 2002  European Bioinformatics Institute, (EMBL-EBI)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
	
	
***************************************************************************
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXSEQLEN 30000
int init(char SequenceFile[]);
int read_matrix(char *);
int user_input=0;
int counter=0;
void detectBiased(char seqStr[]);
void checkAA( int resIndex, char seqStr[]);
void dump_matrix(void);
void print_banner(void);

/* 
 ***********************************************************************************************
 This is the Blosum62 matrix. Its rows & cols follow the same order as allAA 

 WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!WARNING!!

 If it is replaced by another matrix the order should be  the same 
   A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,  F,  P,  S,  T,  W,  Y,  V,  X,  B,  Z,  *

 The values corresponding to X against other residues is calculated,according 
 to the review article by Altschul S.F., et. al (1994), Nature Genetics, vol. 6, 
 as the mean value of the residue-pair scores in the same row.(all 20 values=>including A-A etc) 
 ***********************************************************************************************
*/

float Matrix[24][24]=
{
{  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -0.8, -2, -1, -4},  
{ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -0.9, -1,  0, -4}, 
{ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -0.9,  3,  0, -4}, 
{ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -1.3,  4,  1, -4}, 
{  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -1.55,-3, -3, -4}, 
{ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -0.6,  0,  3, -4}, 
{ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2, -0.85, 1,  4, -4}, 
{  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1.7, -1, -2, -4}, 
{ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -0.85, 0,  0, -4}, 
{ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -1.35,-3, -3, -4}, 
{ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -1.25,-4, -3, -4}, 
{ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -0.85, 0,  1, -4}, 
{ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -0.65,-3, -1, -4},  
{ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -1.25,-3, -3, -4}, 
{ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -1.6, -2, -1, -4}, 
{  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -0.55, 0,  0, -4}, 
{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -0.7, -1, -1, -4}, 
{ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -1.6, -4, -3, -4}, 
{ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -0.8, -3, -2, -4}, 
{  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -1.1, -3, -2, -4}		     
};

float scores[MAXSEQLEN];

int Cutoff = 40 ;
char allAA[24]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X','B','Z','*'};
char titleStr[1000]={'\0'}; 

FILE* fileOut;
FILE* fileStat;
char output_file[200];
int from =0;
int to = 0;

/* Command-Line Flags */

int output=0;
int nostderr=0;
int dostderr=0;
int debug=0;

int main(int argc, char *argv[] )
{

char  SequenceFile[100]= {'\0'};	  
int nbOpt = 0;
int result = 0;

#ifdef WINMAC
  if ( (argc == 1) || (( (argv[1])[0] == '-' ) &&  ( (argv[1])[1] != 'h' )))
  {
  user_input=1;
  print_banner();
  printf ("Please Cite:\n");
  printf ("\n\tPromponas V.J., Enright A.J., Tsoka S., Kreil D.P.,\n");
  printf ("\tLeroy C., Hamodrakas S., Sander C., Ouzounis C.A;\n");
  printf ("\n\tBioinformatics 16(10) 915-922 (2000)\n\n");
  printf ("No sequence file entered...\n");
  printf ("To get some type %s -help from the command line\n\n",argv[0]); 
  printf("\nYou can run CAST by dragging a FASTA file onto the cast application icon\n\n");
  printf("Alternatively run cast from the command-line: cast somfile.fasta\n\n");
  printf("Please type the name of a file to be filtered: ");
  scanf("%s",argv[1]);
  }
#endif

 if ( ((argc == 1) || (( (argv[1])[0] == '-' ) &&  ( (argv[1])[1] != 'h' ))) && (!user_input))
 {
  print_banner();
  printf ("Please Cite:\n");
  printf ("\n\tPromponas V.J., Enright A.J., Tsoka S., Kreil D.P.,\n");
  printf ("\tLeroy C., Hamodrakas S., Sander C., Ouzounis C.A;\n");
  printf ("\n\tBioinformatics 16(10) 915-922 (2000)\n\n");
  printf ("No sequence file entered...\n");
  printf ("To get some help try %s -help from the command line\n\n",argv[0]); 
  return -1;
 }
else
 {


strcpy(SequenceFile,argv[1]);
#ifdef WINMAC
print_banner();
strcpy(output_file,SequenceFile);
strcat(output_file,"-filtered");
printf("Filtering File: %s\n",SequenceFile);
printf("Output File   : %s\n",output_file);
fileOut=fopen(output_file,"w");
#endif


	for(nbOpt=1; nbOpt < argc; nbOpt++)
	{
        if (strcmp(argv[nbOpt] , "-thr")==0)
         {
         Cutoff = atoi(argv[++nbOpt]);
         } 
  	if ( strcmp(argv[nbOpt], "-help")==0 || strcmp(argv[nbOpt], "-h")==0)
  	 {
	    print_banner();
            printf (" Usage: %s SequenceFile [options]\n\n",argv[0]);
            printf ("                 -help    ... print this text\n");
            printf ("                 -thr t   ... set the threshold score for reported regions\n");
            printf ("                              default is %d\n", Cutoff);
            printf ("                              t should be an integer number\n");
            printf ("                 -stat    ... outputs statistics information to file cast.stat\n");
	    printf ("                 -matrix  ... use different mutation matrix (.mat) file\n");
	    printf ("                 -verbose ... verbose mode prints filtering information to standard output\n");
	    printf ("                 -stderr  ... verbose mode prints filtering information to standard error\n\n");

	    return 1;
	 }	 

	if (strcmp(argv[nbOpt], "-stat") == 0) 
	 {
	 output=1;
	 }

	if (strcmp(argv[nbOpt], "-verbose") == 0)
	 {
	 nostderr=1;
	 dostderr=0;
	 }

	if (strcmp(argv[nbOpt], "-stderr") == 0)
         {
         nostderr=0;
         dostderr=1;
         } 

	if (strcmp(argv[nbOpt], "-debug") == 0)
         {
	 debug=1;
         }

	if (strcmp(argv[nbOpt], "-matrix") == 0)
         {

	 if (argc <= nbOpt+1)
		{
		fprintf(stderr,"Error: Mutation matrix file not specified\n");
		exit(1);
		}

	 /* Read In User Specified .mat File */

	 read_matrix(argv[nbOpt+1]);
         }

	}	

	result = init(SequenceFile);
	if (result == -1)
	{
	 fprintf(stderr,"Error: Program not terminated correctly\n");
	 return -1;
	}
 }

#ifdef WINMAC
fclose(fileOut);
printf("\nCAST v1.0 Complete\n");
#endif
return 1;
}
/*
 *****************************************
 ***************END OF MAIN***************
 *****************************************  
*/



int init( char SequenceFile[])
{
char seqStr[MAXSEQLEN]={'\0'};
char c;

FILE *fileIn;

int i=0;

if (!(fileIn = fopen(SequenceFile, "r")))
	{
	fprintf(stderr,"Error: Sequence input file not found or not readable\n");
	exit(1);
	}

if (output == 1)
	{
	fileStat = fopen("cast.stat", "w");
	}


	/******** 24/01/99  - Removed Original File I/O Section *******************/
	/******** Replaced with shorter more efficient, bug-free code  ************/
	/********						       ************/
	/********		       Anton Enright - anton@ebi.ac.uk ************/
	/**************************************************************************/
        /******** This version works on the concept that anything          ********/
  	/******** between a '>' and a '>' or a '>' and an EOF is sequence  ********/
	/******** After removing title string, spaces and carriage returns ********/
	/**************************************************************************/

	while ((c=fgetc(fileIn)) != EOF)
		{

		if (c=='>')
			{
			ungetc(c,fileIn);
			i=0;
			while ( (c=fgetc(fileIn)) != '\n')
				{
				titleStr[i]=c;
				i++;
				}	
			titleStr[i]='\0';

			#ifdef WINMAC
			counter++;
			fprintf(stdout,".");
			if (counter > 60)
			 {
			 fprintf(stdout,"\n");
			 counter=0;
			 } 
			fflush(stdout);
			fprintf(fileOut,"%s\n",titleStr);
			#else
			printf("%s\n",titleStr);
			#endif

			i=0;
			while ( ((c=fgetc(fileIn)) != EOF ) && (c != '>') )
				{
					if ((c != '\n') && (c != ' ') && (c != '\r'))
						{

						/* 15/12/2000			*/
						/*				*/
						/* Convert Lower Case to Upper  */
						/*				*/
						/* Corrected by Anton		*/
		
						if (((c >= 'a') && (c <= 'z')))
                                                        {
							c += 'A' - 'a'; 
                                                        }    

						if (((c >= 'A') && (c <= 'Z')) || (c == '*'))
							{
							seqStr[i]=c;
							i++;	
							}
						else
							{
					fprintf(stderr,"Error: Sequence file contains non-standard or lowercase characters\n");
							exit(1);
							}
						}
				}
			seqStr[i]='\0';
			ungetc(c,fileIn);
			detectBiased(seqStr);
			}

		}
	if (output) { fclose(fileStat); }
return(1);
}

/*
 *******************************************
 ***************END OF INIT*****************
 *******************************************
*/


void checkAA( int resIndex, char seqStr[])
{
int i=0;
int pos = 0 ;
int toto;
char *cc;
int s = 0 ;
int score =0;
int foundaa =0;
from=0;
to=0;


/*
 **************************************************
 count the scores....this part matches to run_SW_VA 
 **************************************************
*/

/* forward */

cc =(char *) strchr(allAA,seqStr[0]);
toto = cc-allAA;
s = Matrix[resIndex][toto];
scores[0] = s;
if ( s  < 0 ) s  = 0;

for(pos = 1; pos < strlen(seqStr); pos++)
{
	/* 15/12/2000								 */
	/* New Code Added By Anton To Stop Crashes with Non-Standard Amino Acids */
	/* 									 */
	/*									 */
 if ((seqStr[pos] != 'J') && (seqStr[pos] != 'O') && (seqStr[pos] != 'U'))
	{
	 cc = (char *) strchr(allAA,seqStr[pos]);
 	toto = cc-allAA;
 	s += Matrix[resIndex][toto];
 	scores[pos] = s;
 	if ( s  < 0 ) s  = 0;
	}
	
}

/* backwards */
cc = (char *) strchr(allAA,seqStr[strlen(seqStr) - 1]);
toto = cc-allAA;
s  = Matrix[resIndex][toto];
if ( s  < 0 ) s  = 0;
for ( pos = strlen(seqStr) - 2 ; pos >= 0; pos--)
{

	/* 15/12/2000                                                            */
        /* New Code Added By Anton To Stop Crashes with Non-Standard Amino Acids */
        /*                                                                       */
        /*                                                                       */
	if ((seqStr[pos] != 'J') && (seqStr[pos] != 'O') && (seqStr[pos] != 'U')) 
	{
	cc = (char *) strchr(allAA,seqStr[pos]);
	toto = cc-allAA;
	scores[pos] += s ;
	s  += Matrix[resIndex][toto];
	if ( s  < 0 ) s  = 0;
	}
}

/*
 **************************************************** 
*/
score = 0;
 for(i=0;i<strlen(seqStr);i++)
 {
 if (score < scores[i])
   {
   score = scores[i];
   from = i;
   }
 }
 for(i=from+1;i<strlen(seqStr);i++)
 {
 if (scores[i] < score) break;
 }
 to = i-1;


  /* 10/04/2001                                                            */
  /* New Code Added By Anton To Stop the Non-Existent AA Loop of Death     */
  /*                                                                       */
  /*                                                                       */ 
  for (i=from;i<=to;i++)
        {
        if (seqStr[i] == allAA[resIndex])
                {
                foundaa=1;
                }
        }

 if (!foundaa)
        {
        scores[from]=0;
        }    

}

/*
 ***********************************************
 *******************END OF checkAA**************
 ***********************************************
*/


void detectBiased(char seqStr[])
{
int counter=0;
int Mscore =0;
int score =0;
int i = 0 ;
int j;
int mfrom=0;
int mto=0;
char maa = ' ';


while(1)
{
maa = ' ';
Mscore =0;
 for(counter=0;counter<20;counter++)
 {

checkAA(counter,seqStr);
score = scores[from];

  if(score >= Cutoff)
  {
   if( Mscore < score )
   {
   maa = allAA[counter];
   mfrom = from;
   mto = to;
   Mscore = score;
   }
  }
 }
if(maa != ' ')
{

if (output == 1)
	{
	fprintf(fileStat,"%s: %c-rich region from %d to %d corrected with score %d\n",titleStr,maa,mfrom+1,mto+1,Mscore);
	}

if (nostderr)
	{
	fprintf(stdout,"%c-rich region from %d to %d corrected with score %d\n",maa,mfrom+1,mto+1,Mscore);
	}

if (dostderr)
	{
	fprintf(stderr,"%c-rich region from %d to %d corrected with score %d\n",maa,mfrom+1,mto+1,Mscore);
	}
 for(i=mfrom;i<=mto;i++)
 {

 if( seqStr[i]==maa ) 
	{
	seqStr[i]='X';
	}
 }
}
 
 i=0;
 j=1;
 
while( (seqStr[i]!='\0') && (maa == ' ') )
{

#ifdef WINMAC
fprintf(fileOut,"%c",seqStr[i]);
#else
fprintf(stdout,"%c",seqStr[i]);
#endif

i++;
if (j==60)
 {
	#ifdef WINMAC
	fprintf(fileOut,"\n");	
	#else
	fprintf(stdout,"\n");
	#endif
  j=0;
 }
j++;
}
if(j!=1) 
	{	
	#ifdef WINMAC
	fprintf(fileOut,"\n");
	#else
	fprintf(stdout,"\n");
	#endif
	}

if(maa == ' ') break;

}

}

/*
 ******************************************************
 ****************END OF detectBiased******************* 
 ******************************************************
*/


/*
 ******************************************************
 **************** read_matrix() ***********************
 ******************************************************
 ******************************************************
 ******* Read in a user defined mutation matrix *******
 ******************************************************
 Anton Enright - Added 21/1/00
*/

/* read_matrix() function, to allow the user to use a
   different mutation matrix file - Anton Enright 01/00
*/

int read_matrix(char *matrix_file)
{

FILE *matrixin;
char temp;
float total_score=0;

int i=0;
int j=0;
int z=0;
int found=0;
int aa_count=0;
int total_aa=0;

int aa_index=0;

/* Initialize the AminoAcid Index */

for (i=0;i<24;i++)
	{
	allAA[i]=0;
	}

/* Clear the Mutation Matrix */

for (i=0;i<24;i++)
	{
	for (j=0;j<24;j++)
		{
		Matrix[i][j]=0;
		}
	}

if (  (matrixin=fopen(matrix_file,"r")) != NULL )
	{
	while( (temp=fgetc(matrixin)) !=  EOF)
		{

		/* Ignore Comments Lines in .mat file */

		if (temp=='#')
			{
			while ( (temp=fgetc(matrixin)) != '\n') {}
			ungetc(temp,matrixin);
			}
	

		/* Read in the First Line of Amino Acids */	

		if (temp == ' ')
			{
			aa_count=0;
			while ( (temp=fgetc(matrixin)) != '\n') 
				{
				if (temp != ' ') 
					{
					if ( (temp != '*') && ((temp < 'A') || (temp > 'Z')) )
						{
						fprintf(stderr,"Error: Matrix File Format not Correct (%c)\n",temp);
						exit(1);
						}

					allAA[aa_count]=temp;
					aa_count++;
					total_aa=aa_count;
					}

				}
			
			
			/* Check to see if the 'B' 'Z' 'X' and '*' columns are present */
			/* If not, then add them to the end */

			found=0;
                        for (j=0;j<aa_count;j++)
                                 {
                                 if (allAA[j] == 'B')
                                           {
                                           found=1;
                                           }
                                 }

                                 if (found==0)
                                        {
                                        allAA[aa_count]='B';
					aa_count++;
                                        }

			found=0;
                        for (j=0;j<aa_count;j++)
                                 {
                                 if (allAA[j] == 'Z')
                                           {
                                           found=1;
                                           }
                                 }

                                 if (found==0)
                                        {
                                        allAA[aa_count]='Z';
                                        aa_count++;
                                        }

			found=0;
                        for (j=0;j<aa_count;j++)
                                 {
                                 if (allAA[j] == 'X')
                                           {
                                           found=1;
                                           }
                                 }

                                 if (found==0)
                                        {
                                        allAA[aa_count]='X';
                                        aa_count++;
                                        }


			found=0;
                        for (j=0;j<24;j++)
                                 {
                                 if (allAA[j] == '*')
                                           {
                                           found=1;
                                           }
                                 }

                                 if (found==0)
                                        {

					/* Set the '*' Values Automatically */

                                        allAA[aa_count]='*';

					for (j=0;j<24;j++)
						{
						Matrix[j][aa_count]=-4;
						Matrix[aa_count][j]=-4;
						}
					Matrix[aa_count][aa_count]=1;
                                        }


			}

		else
			{
			/* Read in Data Line for Each Amino Acid (ROW) */

			for (i=0;i<total_aa;i++)
				{

				 if ( (temp != '*') && ((temp < 'A') || (temp > 'Z')) )
                                                {
                                                fprintf(stderr,"Error: Matrix File Format not Correct (%c)\n",temp);
                                                exit(1);
                                                }          

				if (temp == allAA[i])
					{
					aa_index=i;
					}
				}			

			for (i=0;i<total_aa;i++)
				{
				fscanf(matrixin,"%f",&Matrix[aa_index][i]);
				}

			while ( (temp=fgetc(matrixin)) != '\n') {}	
			}
	
		}
	}
else
	{
	fprintf(stderr,"Error: Matrix file not found or not readable\n");
	exit(1);
	}


/* Correct The Matrix 'X' Values                                                                   */
/* The values corresponding to X against other residues is calculated,according                    */
/* to the review article by Altschul S.F., et. al (1994), Nature Genetics, vol. 6,                 */
/* as the mean value of the residue-pair scores in the same row.(all 20 values=>including A-A etc) */

for (i=0;i<20;i++)
	{

	/* First, Total Up the Scores */

	total_score=0;
	for (j=0;j<20;j++)
		{
		total_score+=Matrix[i][j];
		}

	/* Calculate the Mean */

	total_score=(total_score/20); 

	/* Fix the Matrix */

	for (z=0;z<24;z++)
		{
		if (allAA[z] == 'X')
			{
			Matrix[i][z]=total_score;
			Matrix[z][i]=total_score;
			}
		}
	}

return(1);
}

/* Dump_matrix() - Function for debugging, error-checking the Mutation Matrix */
/* Anton Enright - 01/00                                                      */

void dump_matrix(void)
{
int i=0;
int j=0;

	for (i=0;i<24;i++)
		{
		printf("%c ",allAA[i]);
		}
	printf("\n");

	for (i=0;i<24;i++)
		{
		printf("%c ",allAA[i]);
		for (j=0;j<24;j++)
			{
			printf("%1.1f ",Matrix[i][j]);
			}
		printf("\n");
		}


}

void print_banner(void)
{
printf ("\nCAST v1.0 - Compositional Bias Filtering Algorithm\n\n");
printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
/* Added GPL License Statement */
printf("CAST v1.0, Copyright (C) 2002, European Bioinformatics Institute, EMBL-EBI.\n");
printf("CAST comes with ABSOLUTELY NO WARRANTY; for details please see the file\n");
printf("LICENSE.  This is free software, and you are welcome to redistribute\n");
printf("it under certain conditions; please see the file LICENSE for details.\n");
printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
}
