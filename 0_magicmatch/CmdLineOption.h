#ifndef _CMD_LINE_OPTION_H_
#define _CMD_LINE_OPTION_H_


// tell the compiler to shut up
#pragma warning(disable:4786)

//#include <iostream> // you may need this
#include "OS.h"

// handy little container for our argument vector
struct CCmdParam
{
   vector<string> m_strings;
};

// this class is actually a map of strings to vectors
typedef map<string, CCmdParam> _CCmdLine;

// the command line parser class
class CmdLineOption : public _CCmdLine
{
private:
	string m_usage;
public:
   /*------------------------------------------------------
      int CCmdLine::CmdLineOption(int argc, char **argv)

      parse the command line into switches and arguments.

      returns number of switches found
   ------------------------------------------------------*/
   CmdLineOption(int argc, char **argv);

   /*------------------------------------------------------
      bool CCmdLine::HasSwitch(const char *pSwitch)

      was the switch found on the command line ?

      ex. if the command line is : app.exe -a p1 p2 p3 -b p4 -c -d p5

      call                          return
      ----                          ------
      cmdLine.HasSwitch("-a")       true
      cmdLine.HasSwitch("-z")       false
   ------------------------------------------------------*/   
   bool        HasSwitch(const char *pSwitch);

   /*------------------------------------------------------

      string CCmdLine::GetSafeArgument(const char *pSwitch, int iIdx, const char *pDefault)

      fetch an argument associated with a switch . if the parameter at
      index iIdx is not found, this will return the default that you
      provide.

      example :
  
      command line is : app.exe -a p1 p2 p3 -b p4 -c -d p5

      call                                      return
      ----                                      ------
      cmdLine.GetSafeArgument("-a", 0, "zz")    p1
      cmdLine.GetSafeArgument("-a", 1, "zz")    p2
      cmdLine.GetSafeArgument("-b", 0, "zz")    p4
      cmdLine.GetSafeArgument("-b", 1, "zz")    zz

   ------------------------------------------------------*/

   string  GetSafeArgument(const char *pSwitch, int iIdx, const char *pDefault);

   /*------------------------------------------------------

      string CCmdLine::GetArgument(const char *pSwitch, int iIdx)

      fetch a argument associated with a switch. throws an exception 
      of (int)0, if the parameter at index iIdx is not found.

      example :
  
      command line is : app.exe -a p1 p2 p3 -b p4 -c -d p5

      call                             return
      ----                             ------
      cmdLine.GetArgument("-a", 0)     p1
      cmdLine.GetArgument("-b", 1)     throws (int)0, returns an empty string

   ------------------------------------------------------*/
   string  GetArgument(const char *pSwitch, int iIdx); 

   /*------------------------------------------------------
      int CCmdLine::GetArgumentCount(const char *pSwitch)

      returns the number of arguments found for a given switch.

      returns -1 if the switch was not found

   ------------------------------------------------------*/
   int         GetArgumentCount(const char *pSwitch);


   void			addUsage(const char* i_usage);

   string		getUsage() { return m_usage; }

protected:
   /*------------------------------------------------------

   protected member function
   test a parameter to see if it's a switch :

   switches are of the form : -x
   where 'x' is one or more characters.
   the first character of a switch must be non-numeric!

   ------------------------------------------------------*/
   bool        IsSwitch(const char *pParam);

};

#endif
