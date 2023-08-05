#ifndef _CMAGIC_EXCEPTION_H_
#define _CMAGIC_EXCEPTION_H_

enum MagicError {USAGE_ERROR, FILE_ERROR, NO_ERROR, LIMIT_ERROR, MEMORY_ERROR};

class CMagicException
{
private:
	char*		m_errorMsg;
	MagicError	m_eid;
public:

	//constructor
	CMagicException(const char*, MagicError i_eid);

	//destructor
	virtual ~CMagicException();

	//get error id
	MagicError	getErrorID() const { return m_eid; }

	//get error msg
	const char*	getErrorMsg() const{ return m_errorMsg; }

private:
	//empy constructor not allowd
	CMagicException();
};

#endif
