/* (c) Wolfgang Tichy 2002 */

/* wolfio.c contains a collection a functions, which facilitate reading
   from files:
   
int fgotonext(FILE *in, const char *label);
  moves filepointer forward to the place after label 

int fgetparameter(FILE *in, const char *par, char *str);
  reads value of next parameter par into string str,
  the filepointer points to the char after the parameter value 

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
*/


#define MAXSTRINGLEN 65535   

/* use this to use wolfio fscanf1 */
#define FSCANF1 fscanf1
/* otherwise use: #define FSCANF1 fscanf */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


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


/* move filepointer forward to the place after label */
int fgotonext(FILE *in, const char *label)
{
 int z;
 char str[MAXSTRINGLEN]; 
 
 while((z=FSCANF1(in,"%s",str)) !=EOF)
 {
   if(strcmp(str,label) ==0) break;
   else continue;
 }
 return z;
}



/* find parameter par in file in */
/* fgetparameter writes the parameter value into str
   it returns EOF if the parameter is not found       */
int fgetparameter(FILE *in, const char *par, char *str)
{
 int z;
 size_t parlength=strlen(par);
 

 while((z=FSCANF1(in,"%s",str)) !=EOF)
 {
   if(strncmp(par, str , parlength) ==0)
   { 
     /* check if str is of form par */
     if(parlength==strlen(str)) 
     { 
       /* read str after par */
       FSCANF1(in,"%s",str);
       if( (strlen(str)==1) && (str[0] == '=') ) 
       { 
         FSCANF1(in,"%s",str);  break;
       }
       if( (strlen(str)>1) && (str[0] == '=') )
       {
         extract_after_EQ(str);  break;
       }
       else   
         break;
     }
     else
     {
       /* check if str is of form par=value */
       if( strlen(str)>parlength+1 )
       {
         extract_after_EQ(str);   break;
       }
       else 
       {
         FSCANF1(in,"%s",str);  break;
       }
     }
   }
   else continue;
 }

 if(z==EOF) str[0]=0;
 return z;
}


/* read everything in str behind '=' and write it again into str */
int extract_after_EQ(char *str)
{
 char *before;
 char *after;
 int i;
 
 before = (char *) calloc( strlen(str)+2, sizeof(char) );
 after =  (char *) calloc( strlen(str)+2, sizeof(char) );
 
 i=extrstr_before_after_EQ(str,before,after);
 strcpy(str,after);
 
 free(before);
 free(after);
 return i;
}


/* read everything in str before and after '=' write it into before and after */
int extrstr_before_after_EQ(const char *str, char *before, char *after)
{
  int i,EQpos;
  
  before[0]=0;
  after[0]=0;
  
  for(i=0; (i<=strlen(str))&&(str[i]!= '=') ;i++)
  {
    before[i]=str[i];
  }
  EQpos=i;
  before[EQpos]=0;

  for(i=EQpos+1; i<=strlen(str) ;i++)
  {
    after[i-EQpos-1]=str[i];
  }
  if(EQpos>=strlen(str)) return EOF;
  return EQpos;
}


/* fscanline : ganze Zeile (bis \n) als string einelesen */
/* N=fscanline(fp,str) liefert : 
      str = String aus Datei (mit Zeiger fp) bis zum Zeilenende
      N   = EOF wenn Dateiende erreicht, ansonsten N=Laenge von str */
/* Aus Datei "in" ganze Zeile (bis \n) als String str einelesen  */
int fscanline(FILE *in,char *str)
{
 int z,strlength=0;
 char ch[3];

 str[0]=0;
 for( z=fgetc(in); (z!='\n')&&(z!=EOF) ;z=fgetc(in) , strlength++ )  
 {
  ch[0]=z; ch[1]=0;
  strcat(str,ch); 
 }
 if(z==EOF && strlength==0 ) return EOF;
 return strlength;
}


/* read everything in str before and after char z,
    write it into before and after     */
int extrstr_before_after(const char *str, char *before, char *after, char z)
{
  int i,zpos;
  
  before[0]=0;
  after[0]=0;
  
  for(i=0; (i<=strlen(str))&&(str[i]!= z) ;i++)
  {
    before[i]=str[i];
  }
  zpos=i;
  before[zpos]=0;

  for(i=zpos+1; i<=strlen(str) ;i++)
  {
    after[i-zpos-1]=str[i];
  }
  if(zpos>=strlen(str)) return EOF;
  return zpos;
}


/* read everything in str before and after string z,
    write it into before and after     */
int find_before_after(const char *str, char *before, char *after, const char *z)
{
  int i,zpos;
  int strl=strlen(str);
  int zl=strlen(z);
  
  before[0]=0;
  after[0]=0;
  
  for(i=0; (i<=strl)&&(strncmp(str+i,z,zl)!= 0) ;i++)
  {    
    before[i]=str[i];
  }
  zpos=i;
  before[zpos]=0;

  for(i=zpos+zl; i<=strl ;i++)
  {
    after[i-zpos-zl]=str[i];
  }
  if(zpos+zl>strl) return EOF;
  return zpos;
}


/* read everything (starting at position p) in str before and after string z,
    write it into before and after     */
int pfind_before_after(const char *str,int p,char *before,char *after,const char *z)
{
  int i,zpos;
  int strl=strlen(str);
  int zl=strlen(z);
  
  before[0]=0;
  after[0]=0;
  if(p>=strl) return EOF;
  
  for(i=p; (i<=strl)&&(strncmp(str+i,z,zl)!= 0) ;i++)
  {    
    before[i-p]=str[i];
  }
  zpos=i;
  before[zpos-p]=0;

  for(i=zpos+zl; i<=strl ;i++)
  {
    after[i-zpos-zl]=str[i];
  }
  if(zpos+zl>strl) return EOF;
  return zpos;
}



/* read everything (starting at position p) in str before the next
   space, \t or \n, write it into word,
   return the position of the space, \t or \n, 
   return EOF if end of string reached and no word found		*/
int sscan_word_at_p(const char *str, int p, char *word)
{
 int i;
 int strl=strlen(str);
 const char *pstr;
 
 if(p>=strl) return EOF;
 if(p<0)     return EOF;
 
 for(i=0;i+p<strl;i++)
   if(str[p+i]!=' ' && str[p+i]!='\t' && str[p+i]!='\n') break;

 if(p+i>=strl) return EOF;
 pstr=str+p+i;
 
 word[0]=0;
 sscanf(pstr,"%s",word);
 return strlen(word)+p+i;
}


/* fscan_str_using_getc is a replacement for reading one string
   with fscanf(in, "%s", str); but using only fgetc */
int fscan_str_using_getc(FILE *in, char *str)
{
  fpos_t pos[2];
  int ch;
  int i=0;

  do
  {
    fgetpos(in, pos);
    ch = fgetc(in);
    if( isspace(ch) && i==0 ) continue;
    if( ch==EOF     && i==0 ) break;
    if( isspace(ch) || ch==EOF)
    {
      str[i]=0;
      break;
    }
    str[i] = ch;
    i++;
  } while(ch!=EOF); 
  fsetpos(in, pos); /* go back to before last read char */

  if(i==0) return EOF; /* return EOF if no string found */
  return i;            /* return string length */
}

/* read one string as with fscanf but using fscan_str_using_getc */
int fscanf1(FILE *in, char *fmt, char *str)
{
  int ret;

  ret = fscan_str_using_getc(in, str);
  if( strcmp(fmt, "%s")!=0 ) 
    printf("fscanf1: wolfio error, format string is not %%s\n");
  if(ret==EOF) return EOF;
  return 1;
}
