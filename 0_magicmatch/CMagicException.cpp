
#include "CMagicException.h"
#include <string.h>

CMagicException::CMagicException(const char* i_errorMsg, MagicError i_eid)
: m_eid(i_eid),
  m_errorMsg(0)
{
	if(i_errorMsg != 0){
		m_errorMsg = new char[strlen(i_errorMsg)];
		strcpy(m_errorMsg, i_errorMsg);
	}

}


CMagicException::~CMagicException()
{
	if(m_errorMsg != NULL){
		delete m_errorMsg;
	}
}
