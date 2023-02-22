// ProcessOptions.h: interface for the ProcessFragmentFile class.
//
//////////////////////////////////////////////////////////////////////
/*
  Copyright by Mohamed I. Abouelhoda (C) 2003
  =====================================
  You may use, copy and distribute this file freely as long as you
   - do not change the file,
   - leave this copyright notice in the file,
   - do not make any profit with the distribution of this file
   - give credit where credit is due
  You are not allowed to copy or distribute this file otherwise
  The commercial usage and distribution of this file is prohibited
  Please report bugs and suggestions to <mibrahim@informatik.uni-ulm.de>
  
*/



#if !defined(AFX_PROCESSOPTIONS_H__84CD68B1_169D_4118_8470_D28E438BAC88__INCLUDED_)
#define AFX_PROCESSOPTIONS_H__84CD68B1_169D_4118_8470_D28E438BAC88__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class ProcessOptions  
{
public:
	//////////////////// Attributes /////////////////////////////////
  double* optionsarray;
  long* fileindexarray;
	//////////////////// Operations ////////////////////////////////
 public:
	ProcessOptions();
	~ProcessOptions();
	int ParseOptions(long argc, char **argv);
	int getnumerical(char*, long* value);
	int getnumfloat(char* token, double* value);
	int getnumint(char* token, long* value);
	int cmpdecimal(char* token, int type);
};

#endif // !defined(AFX_PROCESSOPTIONS_H__84CD68B1_169D_4118_8470_D28E438BAC88__INCLUDED_)






