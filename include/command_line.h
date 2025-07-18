// $Id: command_line.h 13 2018-12-15 00:34:28Z n7dr $

// Released under the GNU Public License, version 2
//   see: https://www.gnu.org/licenses/gpl-2.0.html

// Principal author: N7DR

// Copyright owners:
//    N7DR

/*! \file   command_line.h

    API for managing the command line
*/

#ifndef COMMANDLINEH
#define COMMANDLINEH

#include <string>

// -----------  command_line  ----------------

/*! \class  command_line
    \brief  Class that implements management of the command line
*/

class command_line
{
protected:
  unsigned int _argc;       ///< Number of arguments

// these can't be const because of the = operator
  char** _argv;             ///< Pointers to arguments

  std::string* _arg;        ///< Pointers to arguments

/// internal initialisation function
  void _init(void);

public:

/*! \brief          Constructor
    \param  argc    number of arguments
    \param  argv    pointer to array of individual arguments
*/
inline command_line(int argc, char** argv) :
    _argc(argc),
    _argv((char**)argv)
  { _init(); }

/*!	\brief	        Copy constructor
	\param	obj 	object to be copied
*/
inline command_line(const command_line& cl) :
    _argc(cl._argc),
    _argv(cl._argv)
  { _init(); }

/*!	\brief	Destructor
*/
  inline ~command_line(void)
    { delete [] _arg; }

/// command_line = command_line
  void operator=(const command_line&);

/*! \brief      Return parameter number (wrt 1)
    \param  n   parameter number
    \return     the value of the <i>n</i>th parameter

    If the value of <i>n</i> does not correspond to a parameter that was actually present, this functions throws an x_command_line_invalid_parameter()
*/
  std::string parameter(const unsigned int n) const;

/*! \brief      Return parameter number (wrt 1)
    \param  n   parameter number
    \return     the value of the <i>n</i>th parameter

    If the value of <i>n</i> does not correspond to a parameter that was actually present, this functions throws an x_command_line_invalid_parameter()
*/
  inline std::string operator[](const unsigned int n) const
    { return parameter(n); }

/*! \brief      Obtain the name of the program
    \return     the name of the program
*/
  inline std::string program_name(void) const
    { return _arg[0]; }

/*! \brief      Obtain the base name of the program
    \return     the base name of the program (i.e., with no "/" characters)
*/
  std::string base_program_name(void) const;

/*! \brief      Obtain the number of parameters passed to the program
    \return     the number of parameters
*/
  inline int n_parameters(void) const
    { return (_argc - 1); }

/*! \brief  Convert the entire command line to lower case
*/
  void tolower(void);

/*! \brief  Convert the entire command line to lower case
*/
  inline void to_lower(void)
    { tolower(); }

/*! \brief  Convert the entire command so that the case matches exactly what was originally passed to the program
*/
  void tooriginal(void);

/*! \brief  Convert the entire command so that the case matches exactly what was originally passed to the program
*/
  inline void to_original(void)
    { tooriginal(); }

/*! \brief  Convert the entire command line to upper case
*/
  void toupper(void);

/*! \brief  Convert the entire command line to upper case
*/
  inline void to_upper(void)
    { toupper(); }

/*! \brief      Convert a particular parameter to lower case
    \param  n   parameter number to convert (wrt 0)
*/
  void tolower(const unsigned int n);

/*! \brief      Convert a particular parameter to lower case
    \param  n   parameter number to convert (wrt 0)
*/
  inline void to_lower(const unsigned int n)
    { tolower(n); }

/*! \brief      Convert a particular parameter to its original case
    \param  n   parameter number to convert (wrt 0)
*/
  void tooriginal(const unsigned int n);

/*! \brief      Convert a particular parameter to its original case
    \param  n   parameter number to convert (wrt 0)
*/
  inline void to_original(const unsigned int n)
    { tooriginal(n); }

/*!     \brief  Convert a particular parameter to upper case
        \param  n  Parameter number to convert (wrt 0) 
*/
  void toupper(const unsigned int n);

/*!     \brief      Is a particular value present?
        \param  v   value for which to look
        \return     whether the value corresponding to <i>v</i> is present
        
        A "value" is something like a parameter to a -xxx option. If, for example, value_present("-xxx") is TRUE, it means that -xxx is present, and a value follows it  
*/
  bool value_present(const std::string_view v) const;
  
/*!     \brief      Return a particular value
        \param  v   value to return
        \return     the value
        
        A "value" is something like a parameter to a -xxx option. If, for example, the command line contains "-xxx burble", then value("-xxx") will return "burble"  
*/
  std::string value(const std::string_view v) const;
  
/*!     \brief      Is a particular parameter present?
        \param  p   parameter for which to look
        \return     whether the parameter <i>p</i> is present
        
        A "parameter" is an actual parameter that appears on the command line.
*/
  bool parameter_present(const std::string_view p) const;

/*! \brief      Return a particular value if it's present
    \param  v   value to return
    \return     the value, if it's present; otherwise, string()

    A "value" is something like a parameter to a -xxx option. If, for example, the command line contains "-xxx burble", then value("-xxx") will return "burble"
*/
  inline std::string value_if_present(const std::string_view v) const
    { return (value_present(v) ? value(v) : std::string { } ); }
};

// ---------------------------  exceptions  ----------------------

/*! \class  x_command_line_invalid_parameter
    \brief  Trivial class for exceptions in command line processing
*/

class x_command_line_invalid_parameter
{
};

#endif    // !COMMANDLINEH
