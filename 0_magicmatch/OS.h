

#ifndef _OS_H_
#define _OS_H_

#ifdef WIN32 //WINDOWS

#pragma warning(disable:4786)

#include <map>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

#else //UNIX

using namespace std;

#include <map>
#include <string>
#include <vector>
#include <iostream>

#define vector std::vector
#define string std::string
#endif


#endif
