
#ifndef _FASTA_FILE_DB_H_
#define _FASTA_FILE_DB_H_

#include <cstring>

#include "md5class.h"
#include "CMagicException.h"
#include "OS.h"


//const unsigned int MAX_ID_LENGTH = 920;
const unsigned int MAX_ID_LENGTH = 600;
const unsigned int MAX_RECORD_NUM = 3200000;
const unsigned int MAX_SEQ_NUM = 64000;


struct FastaDB {
	char*	m_pRecord;
	char	m_id[MAX_ID_LENGTH];
	char*	m_pSeq;
	char	m_checksum[33];

	FastaDB(){};
};

struct compare {
    bool operator() (const char* i_pStr1, const char* i_pStr2) const {
        return strcmp(i_pStr1,i_pStr2)>0;
    }
};


typedef multimap<const char*, const char*, compare> FastaDBmap;
typedef pair<const char*, const char*> key_pair;

class FastaFileDB {

private:
	FastaDB*		m_fileStoreDB;
	CMD5*			m_pMd5;
	FastaDBmap		m_fastaDBmap;

	char*			m_pFileBuffer;
	char*			m_pSequence;
	const char*		m_pValid_chars;
	const char*		m_pLowercase;
	unsigned int	m_numOfRecords;
	unsigned int	m_lineCount;


public:
	enum	EExecute	{ RUN_COMPARE, RUN_TRANSLATE };

	FastaFileDB();
	virtual ~FastaFileDB() {};

	void	analyzeStoreFile(const char* i_fileName , EExecute i_execute);

	void	analyzeStoreFile(const char* i_fileName);

	void	show() const;

	void	compareIds(const FastaFileDB* i_fastaFileDB) const;

protected:
	void	readFile(const char* i_fileName);
	void	analyze(EExecute i_execute);
	void	analyze();
};

#endif
