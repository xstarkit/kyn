#include <stdio.h>

/*
   Some routines defined in xspec and used in KY local models xs*.c
   they are needed when we use routines outside xspec...
*/

// This function should return string with path to KBHtablesNN.fits
// We define the string to be "./" (the tables are in working dir)
char* FGMODF(void){
  return("./");
}

// This function should return string with path to KBHtablesNN.fits
// defined in KYDIR
// We define the string to be ""
char* FGMSTR(char* dname){
  return("");
}

// This function should define and set some variable inside xspec
void FPMSTR(const char* value1, const char* value2){
  return;
}

// Next subroutine should write a text to the terminal
int xs_write(char* text,int  idest){
  fprintf(stdout,"%s\n",text);
  return(0);
}
