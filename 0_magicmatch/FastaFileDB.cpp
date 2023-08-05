

#include "FastaFileDB.h"
#include "md5class.h"

#include <cstdio>
#include <cstdlib>

#include <sys/stat.h>


FastaFileDB::FastaFileDB()
: m_pValid_chars("ABCDEFGHIKLMNPQRSTUVWXYZ-"), // Valid IUPAC 1-letter codes
  m_pLowercase("abcdefghijklmnopqrstuvwxyz"),
  m_pSequence(new char[MAX_SEQ_NUM]),
  m_numOfRecords(0),
  m_lineCount(0),
  m_pMd5(new CMD5()),
  m_fileStoreDB(0),
  m_pFileBuffer(0)
{
	m_fileStoreDB = new FastaDB[MAX_RECORD_NUM];
}

void FastaFileDB::analyzeStoreFile(const char* i_fileName)
{
	if(m_fileStoreDB==0){
		throw(new CMagicException("Not enough memory...", MEMORY_ERROR));
	}
	readFile(i_fileName);
	analyze();
	delete m_pFileBuffer;
	m_pFileBuffer=0;
}


void FastaFileDB::analyzeStoreFile(const char* i_fileName, EExecute i_execute)
{
	readFile(i_fileName);
	analyze(i_execute);
	delete m_pFileBuffer;
	m_pFileBuffer=0;
}

void FastaFileDB::readFile(const char* i_fileName)
{
	FILE *fp;
	struct stat file_stat;
	const  char *file_mode = "rb";
	long no_of_items_read, item_size, no_items =1;

	// File open
	if ((fp = fopen(i_fileName, file_mode)) == NULL)
	{
		throw(new CMagicException("Cannot open file ... quitting\n", FILE_ERROR));
	}

	// File read
	stat (i_fileName, &file_stat);
	item_size = file_stat.st_size;

	m_pFileBuffer = new char[item_size+1];
	if(m_pFileBuffer==0){
		throw("Not enough memory...", MEMORY_ERROR);
	}

	no_of_items_read = fread (	m_pFileBuffer,
								item_size,
								no_items,
								fp );

	if (no_of_items_read != no_items)
	{
		throw(new CMagicException("File read error ... quitting\n", FILE_ERROR));
	}

	// File close
	if (fclose(fp))
	{
		throw(new CMagicException("Close file error ... quitting\n", FILE_ERROR));
	}

	m_pFileBuffer[item_size]='\0';
}


void FastaFileDB::analyze(EExecute i_execute)
{
	unsigned int localIdlength;
	FastaDB* pDB = m_fileStoreDB;
	char* pData =  m_pFileBuffer;
	while( (pData=strchr(pData, '>')) != NULL){
		m_numOfRecords++;

		pDB->m_pRecord = pData;
		pDB->m_pSeq    = strchr(pDB->m_pRecord, '\n') + 1;

		localIdlength = (unsigned int)(pDB->m_pSeq - pData - 1);
		localIdlength = (localIdlength>MAX_ID_LENGTH) ? MAX_ID_LENGTH:localIdlength;

		memcpy (pDB->m_id, pData+1, localIdlength);
		pDB->m_id[localIdlength - 1] = '\0';
		pData += localIdlength;
		pDB++;

		if(m_numOfRecords == MAX_RECORD_NUM-1){
			char buffer[128];
			sprintf(buffer,"Max number of record reached please split the file so that the max number of records will not exceed %d !!\n", MAX_RECORD_NUM);
			throw(new CMagicException(buffer, LIMIT_ERROR));
		}
	}

	char* pSeq_front;
	char* pSeq_back;
	char* last_char;

	pDB = m_fileStoreDB;

	for (unsigned int i=0; i<m_numOfRecords; i++){
		m_pSequence[0] = '\0';
		pSeq_front = pSeq_back = pDB->m_pSeq;
		pDB++;

		while(pSeq_front < pDB->m_pRecord && i<m_numOfRecords-1){

			// Filter out formatting chars (newline, tab & whitespace) and 
			// '*' char which are used to terminate PIR FASTA files
			pSeq_front = strpbrk (pSeq_back, "\n\t *");     

			if (*pSeq_front == '\n'){ m_lineCount++; }

			strncat (m_pSequence, pSeq_back, pSeq_front - pSeq_back);
			pSeq_front++;
			pSeq_back = pSeq_front;
		}	

		// Extract final sequence record
		if (i == m_numOfRecords-1)
		{
			last_char = pSeq_back + strlen(pSeq_back);
			pSeq_front = pSeq_back;

			while (pSeq_front < last_char)
			{
				// Filter out formatting chars (newline, tab & whitespace)
				pSeq_front = strpbrk (pSeq_back, "\n\t *");    

				if (*pSeq_front == '\n'){ m_lineCount++; }

				strncat (m_pSequence, pSeq_back, pSeq_front-pSeq_back);
				pSeq_front++;
				pSeq_back = pSeq_front;
			}
			i++;
		}

		// Convert any lowercase chars in sequence to uppercase
		pSeq_front = m_pSequence;

		while (pSeq_front != '\0')
		{
			pSeq_front = strpbrk (pSeq_front, m_pLowercase);

			if (pSeq_front)
			{
				if (pSeq_front < m_pSequence+strlen(m_pSequence))
				{
					*(pSeq_front)-= 0x20;
				}
				pSeq_front++;
			}
		}

		// Validate sequence data ... report line/posn of non-valid chars		
		unsigned int len = strspn (m_pSequence, m_pValid_chars);

		if (len != strlen(m_pSequence))
		{
			char buffer[1024];
			sprintf(buffer,"Invalid char found ...protein id:%s char:'%c' at line %d, posn %d\n", (pDB-1)->m_id, m_pSequence[len], (int)(m_lineCount+i+1), (int)len+1);
			throw(new CMagicException(buffer, FILE_ERROR));
		}

		m_pMd5->setPlainText(m_pSequence);

		strcpy((pDB-1)->m_checksum, m_pMd5->getMD5Digest());
		
	
		if(i_execute == RUN_COMPARE){
			m_fastaDBmap.insert(key_pair((pDB-1)->m_checksum,(pDB-1)->m_id));
			m_fastaDBmap.size();
		}

		if(i_execute == RUN_TRANSLATE){
			cout<<(pDB-1)->m_checksum<<"\t"<<(pDB-1)->m_id<<endl;
		}
	}
}


void FastaFileDB::analyze()
{
	int i;
	char* pData =  m_pFileBuffer;
	while(*pData != '\0'){
		if(*pData == '\n'){
			pData++;
			continue;
		}
		i=0;
		for(pData ; *pData != '\t'; pData++){
			m_fileStoreDB->m_checksum[i++] = *pData;
		}
		m_fileStoreDB->m_checksum[i]='\0';

		i=0;


		for(pData ; *pData != '\n' || *pData=='\0' ; pData++){
			m_fileStoreDB->m_id[i++] = *pData;
		}
		m_fileStoreDB->m_id[i]='\0';
		m_fastaDBmap.insert(key_pair(m_fileStoreDB->m_checksum,m_fileStoreDB->m_id));
		m_fileStoreDB++;
	}
}

void FastaFileDB::compareIds(const FastaFileDB* i_fastaFileDB) const
{
	FastaDBmap::const_iterator it = i_fastaFileDB->m_fastaDBmap.begin();
	FastaDBmap::const_iterator it_found;
	for(it; it!=i_fastaFileDB->m_fastaDBmap.end(); it++){

		for (it_found = m_fastaDBmap.lower_bound(it->first); 
				it_found != m_fastaDBmap.upper_bound(it->first) && m_fastaDBmap.end() != it_found; 
				++it_found) {

			cout<<it_found->second<<"\t"<<it->second<<endl; 
		}
	}
}

void FastaFileDB::show() const
{
	FastaDBmap::const_iterator it = m_fastaDBmap.begin();
	for(it; it!=m_fastaDBmap.end(); it++){
			cout<<it->first<<"\t"<<it->second<<endl; 
	}
}
