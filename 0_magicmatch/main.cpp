#include <cstdlib>

#include "OS.h"
#include "CmdLineOption.h"
#include "FastaFileDB.h"
#include "CMagicException.h"


int main( int argc, char* argv[] )
{

	CmdLineOption cmdLine(argc, argv);

	cmdLine.addUsage("Usage: MagicMatch [-help show usage] [-f file1 file2 ...] [-fe file1 file2 ...] [-t translate to ID-KEY || -c find identicals between files]");
	cmdLine.addUsage("	-f		fasta file names (at least one file is required)");
	cmdLine.addUsage("	-fe		previously encoded fasta files ");
	cmdLine.addUsage("	-t		translate the sequences of the files");
	cmdLine.addUsage("	-c		map between ids of all the files");
	cmdLine.addUsage("	-h		show usage");


	
	try{
		//check to show usage 
		int compare = cmdLine.GetArgumentCount("-h");
		if(compare != -1){
			throw(new CMagicException("", NO_ERROR));
		}

		int i;
		FastaFileDB* fastaFileDB;
		
		//getting number of fasta files to process
		int num_of_fasta_files = cmdLine.GetArgumentCount("-f");

		//getting number of fasta files to process
		int num_of_encoded_files = cmdLine.GetArgumentCount("-fe");

		int num_of_all_files=0;

		//if no files throw exception
		if(num_of_fasta_files == -1 && num_of_encoded_files == -1){
			throw(new CMagicException("No files, please specify at least one file\n", USAGE_ERROR));
		}

		num_of_all_files += num_of_fasta_files<0?0:num_of_fasta_files;
		num_of_all_files += num_of_encoded_files<0?0:num_of_encoded_files;

		//allocate FastaFileDB to process all the files
		fastaFileDB = new FastaFileDB[num_of_all_files];

		int translate = cmdLine.GetArgumentCount("-t");
		if(translate != -1){
			//process files for translation
			for(i=0; i<num_of_fasta_files; ++i){
				fastaFileDB[i].analyzeStoreFile(cmdLine.GetArgument("-f",i).c_str(), FastaFileDB::RUN_TRANSLATE);
			}

		}else{
			int compare = cmdLine.GetArgumentCount("-c");
			if(compare == -1){
				throw(new CMagicException("Please specify either -compare || -translate \n", USAGE_ERROR));
			}
			//process files for comparison
			cmdLine.GetArgument("-c",0).c_str();
			//if only one file exist throw exception (nothing to compare to)
			if(num_of_all_files == 1){
				throw(new CMagicException("No files to compare, please specify at least another file\n", USAGE_ERROR));
			}
			for(i=0; i<num_of_fasta_files; ++i){
				fastaFileDB[i].analyzeStoreFile(cmdLine.GetArgument("-f",i).c_str(), FastaFileDB::RUN_COMPARE);
			}
			for(int j=i; j<num_of_all_files; ++j){
				fastaFileDB[j].analyzeStoreFile(cmdLine.GetArgument("-fe",j-i).c_str());
			}

			for(i=0; i<num_of_all_files; ++i){
				for(int j=i+1; j<num_of_all_files; j++){
					fastaFileDB[i].compareIds(&fastaFileDB[j]);
				}
			}
		}
		
	}catch(CMagicException* e){
		if(e->getErrorID() == USAGE_ERROR){
			cout<<"USAGE ERROR!!\n\n";
		}else if(e->getErrorID() == FILE_ERROR){
			cout<<"FILE PROCESSING ERROR!!\n\n";
		}

		cout<<e->getErrorMsg();

		if(e->getErrorID() == USAGE_ERROR || e->getErrorID() == NO_ERROR){
			cout<<cmdLine.getUsage();
		}

		exit(1);
	}
	return 0; 
}
