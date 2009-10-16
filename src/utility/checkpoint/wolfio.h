/* (c) Wolfgang Tichy 2002 */

/* wolfio.c contains a collection a functions, which facilitate reading
   from files:
   
int fgotonext(FILE *in, const char *label);
  moves filepointer forward to the place after label 

int fgetparameter(FILE *in, const char *par, char *str);
  reads value of parameter par into string str 
  right now the filepointer points to the char after the parameter value 
  (this should not be the case, fgetpos and fsetpos cause seg faults!!!)
  (yet we should do this without moving the filepointer!)  

int extract_after_EQ(char *str);
  extracts everything in str behind the equal sign (=) and writes it back into str 

int extrstr_before_after_EQ(const char *str, char *before, char *after);
  writes LHS of str into before and RHS of str into after
  
int extrstr_before_after(const char *str, char *before, char *after, char z);
   writes all of str before char z into before and all after z into after
   returns the position of z in str

int fscanline(FILE *in,char *str);
  read a whole line into str   
  note that sscanf(str, const char *format, ...); 
  can then read strings or floats or ... from str
  
int find_before_after(const char *str, char *before, char *after, const char *z);
  writes all of str before string z into before and all after z into after
  returns the position of z in str or EOF if z in not found

int pfind_before_after(const char *str,int p,char *before,char *after,const char *z);
  writes all of str (starting at position p) before string z into before 
  and all after z into after
  returns the position of z in str or EOF if z in not found
  
int sscan_word_at_p(const char *str, int p, char *word);
  read everything (starting at position p) in str before the next
  space, \t or \n, write it into word,
  returns the position of the space, \t or \n,
  returns EOF if end of string reached and no word found	

int fscan_str_using_getc(FILE *in, char *str);
  fscan_str_using_getc is a replacement for reading one string
  with fscanf(in, "%s", str); but using only fgetc 

int fscanf1(FILE *in, char *fmt, char *str);
  read one string as with fscanf but using fscan_str_using_getc 
*/

/*
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
*/

/* Functions */
int fgotonext(FILE *in, const char *label);
int fgetparameter(FILE *in, const char *par, char *str);
int extract_after_EQ(char *str);
int extrstr_before_after_EQ(const char *str, char *before, char *after);
int fscanline(FILE *in,char *str);
int extrstr_before_after(const char *str, char *before, char *after, char z);
int find_before_after(const char *str, char *before, char *after, const char *z);
int pfind_before_after(const char *str,int p,char *before,char *after,const char *z);
int sscan_word_at_p(const char *str, int p, char *word);
int fscan_str_using_getc(FILE *in, char *str);
int fscanf1(FILE *in, char *fmt, char *str);
