/* Parse output of diff and check for differences in floating point numbers */
/* (c) Wolfgang Tichy 2002 */

/*
floatdiff can find the relative and absolute differences in diff output
which looks e.g. like this:
--------------------------------------------------------------------
diff L8_PNO5_96_0.8 L8_PNO5_9
1c1
< 1 2.9
---
> 1 2.99
diff L8_PNO5_96_0.8_robin//K.xl L8_PNO5_96_0.8_robin_new//K.xl
1,3c1,3
< sdhfsfhjsd
< 1 2
< 2 4
---
> kljgklfjgkf
> 1 2.1
> 2 4.2
---------------------------------------------------------------------
In particular the file must contain:
stuff like  
1c1 
< something
---
> something else

floatdiff looks for numbers with a c in the middle and also for --- 
*/


#define STRLEN 65535

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/wait.h>
/* #include "wolfio.h" */   /* commented out because of next line */
#include "../src/utility/checkpoint/wolfio.c"  /* include some funcs directly instead of using makefile */


int main(int argc, char *argv[])
{
 FILE *difffile; 
 char difffilename[STRLEN];
 char *filename1;
 char *filename2;
 char currentdiff[STRLEN];
 char str[STRLEN];
 char str1[STRLEN];
 char str2[STRLEN];
 char substr[STRLEN];
 char substr1[STRLEN];
 char substr2[STRLEN];
 char pos1[STRLEN];
 char pos2[STRLEN];
 int begline1, endline1, begline2, endline2;
 long filepos1;
 long filepos2;
 int line, col;
 int n, p1,p2;
 double x1,x2;
 double reldiff, absdiff, reltol=0, abstol=0;
 int i, firstline, details=1;
 int exerr=0, exwarn=0;
 char *astr;

 if(argc<3 || argc>9)
 {  
  printf("Usage: floatdiff [options] file1 file2\n");
  printf("or:    diff file1 file2 > diff.out ; \n");
  printf("       floatdiff [options] --diff-out diff.out\n\n");
  printf("options: -rt rel.tolerance   to set relative tolerence\n");
  printf("         -at abs.tolerance   to set absolute tolerence\n");
  printf("         -d  detail-level    0 quiet, 1 show number, 2 show line\n");
  return -1;
 } 


 /* set the name of the diff output file to NULL */
 difffilename[0]=0;
 
 /* parse command line options, which start with - */
 for(i=1; (i<argc)&&(argv[i][0] == '-'); i++)
 {
  astr = argv[i];

  if( (strcmp(astr+1,"rt")==0) )
  {
    if(i>=argc-1) 
    {
      printf("no relative tolerence value specified, after -rt\n");
      return -1;
    } 
    reltol=atof(argv[i+1]);
    i++;
  }
  else if( (strcmp(astr+1,"at")==0) )
  {
    if(i>=argc-1) 
    {
      printf("no absolute tolerence value specified, after -at\n");
      return -1;
    } 
    abstol=atof(argv[i+1]);
    i++;
  }
  else if( astr[1]=='d' )
  {
    if(i>=argc-1) 
    {
      printf("no detail-level specified, after -d\n");
      return -1;
    } 
    /* set details flag, if 2 print more details on diffs */
    details=atoi(argv[i+1]);
    i++;
  }
  else if( (strcmp(astr+1,"-diff-out")==0) )
  {
    if(i>=argc-1) 
    {
      printf("no filename specified, after --diff-out\n");
      return -1;
    }
    snprintf(difffilename,STRLEN,"%s",argv[i+1]);
    i++;
  }
  else 
  {
    printf("floatdiff: unknown argument %s\n", astr);
    return -1;
  }
 }

 /* parse the remaining two arguments if we don't have a difffilename yet */
 if(difffilename[0]==0)
 {
   for (n=0;  (i < argc)&&(n<2);  i++, n++)
   {
     if(n==0) filename1= (char *) strdup(argv[i]);
     if(n==1) filename2= (char *) strdup(argv[i]);
   }
  if(n!=2)
  {
    printf("floatdiff: 2 file or directory names are needed as arguments!\n");
    return -1;
  }
  
  snprintf(difffilename,STRLEN,"%s",".last_tmp_diffout.floatdiff");
  snprintf(str,STRLEN,"rm -f %s", difffilename);
  /* printf("%s\n",str); */
  system(str);

  /* works only with bash: */  
  snprintf(str,STRLEN,"diff -r %s %s >& %s", 
           filename1, filename2, difffilename);
  /* works with all shells: */
  snprintf(str,STRLEN,"diff -r %s %s > %s 2>&1", 
           filename1, filename2, difffilename);
  /* printf("%s\n",str); */

  n=system(str);
  if( (WEXITSTATUS(n)>1)||(WEXITSTATUS(n)<0) )
  {
    printf("floatdiff: Trouble executing\n %s\n",str);
    printf(" -> Return value n=%d     WEXITSTATUS(n)=%d\n", n, WEXITSTATUS(n));
  }
 }
 else
  printf("reading diffs from %s \n",difffilename);

 /* printf("filename1=%s  filename2=%s\n",filename1,filename2);
    printf("difffilename[1]=%s \n",difffilename);                */

 
 /* open difffile */
 difffile=fopen(difffilename,"r");
 if(difffile==NULL)
 {
  printf("%s not found. \n",difffilename);
  return -1;
 }

 printf("floatdiff:  reltol=%e  abstol=%e  details=%d\n",
        reltol, abstol, details);

 /* now start comparing numbers */
 firstline=1;
 currentdiff[0]=0;
 while(fscanline(difffile, str)!=EOF)
 {
  sscanf(str,"%s",substr);

  /* check if a diff has been done between two files */
  if(strncmp(str,"diff ",5)==0) 
  {
   strcpy(currentdiff,str);
   /* find_before_after(str,substr,currentdiff,"diff"); */
   firstline=1;
   continue;
  }
  
  /* print any warnings diff could have given */
  if(isalpha(str[0])) 
  {
   if(details>0) printf("%s\n",str);
   exwarn++;
   continue;
  }
  
  /* check if the first char is a number, which could indicate a line number */
  if(isdigit(substr[0]))
  {
   /* check if lines were compared by diff, i.e. if c is contained in substr */
   n=extrstr_before_after(substr,pos1,pos2, 'c');
   if(n==EOF) 
   {
     if(firstline==1)
     {
       if(details>0) printf("\n%s\n",currentdiff);
       firstline=0;
     }
     if(details>0) 
       printf("WARNING (%s): One of the files contains additional lines!\n",
              substr);
     exerr++;
     continue;
   }
   
   /* extract first and last lines compared by diff in first file*/
   n=extrstr_before_after(pos1,substr1,substr2, ',');
   begline1=atoi(substr1); 
   if(n!=EOF) endline1=atoi(substr2);
   else       endline1=begline1;
   
   /* extract first and last lines compared by diff in second file*/
   n=extrstr_before_after(pos2,substr1,substr2, ',');
   begline2=atoi(substr1); 
   if(n!=EOF) endline2=atoi(substr2);
   else       endline2=begline2;
   
   if(endline1-begline1 != endline2-begline2)
   {
     if(firstline==1)
     {
      if(details>0) printf("\n%s\n",currentdiff);
      firstline=0;
     }
     if(details>0)
     {
      printf("WARNING: Files contain different number of lines! ");
      printf("All comparisons may fail!\n");
     }
     exerr++;
   }
   
   /* save starting position in first file */
   filepos1=ftell(difffile);
   /* save starting position in second file */
   while(fscanline(difffile, str)!=EOF && strncmp(str,"---",3)!=0 );
   filepos2=ftell(difffile);
   
   /* go through lines of first file */
   for(line=begline1; line<=endline1; line++)
   {
    fseek(difffile,filepos1,SEEK_SET);
    for(i=0;i<=line-begline1;i++)   fscanline(difffile, str1);

    fseek(difffile,filepos2,SEEK_SET);
    for(i=0;i<=line-begline1;i++)   fscanline(difffile, str2);
    
    /* set file pointer to where it was before, so it's on the right place 
       on break */
    fseek(difffile,filepos1,SEEK_SET);
    
    if( str1[0]!='<' || str2[0]!='>' ) break;
    
    /* go through columns of line in first file 
       and compare to columns in second file    */
    p1=0; p2=0; col=-1;
    while( p1!=EOF )
    {
     col++;
     p1=sscan_word_at_p(str1,p1,substr1); 
     p2=sscan_word_at_p(str2,p2,substr2);
     if(p1==EOF && p2==EOF) break;
     if(p1==EOF || p2==EOF)
     {
      if(firstline==1)
      {
       if(details>0) printf("\n%s\n",currentdiff);
       firstline=0;
      }
      if(details>0)
        printf("L%dC%d WARNING: Files contain different number of columns!\n",
               line,col);
      if(details>1) printf("%s\n%s\n",str1,str2);
      exerr++;
      break;
     }
     x1=atof(substr1);
     x2=atof(substr2);
     if(x1==0.0 && x2==0.0) continue;
     if(x1!=0.0) reldiff=fabs( (x2-x1)/x1 );
     else        reldiff=fabs( (x2-x1)/x2 );
     absdiff=fabs( x2-x1 );
     if(reldiff<=reltol) continue;
     if(absdiff<=abstol) continue;
     if(firstline==1)
     {
      if(details>0) printf("\n%s\n",currentdiff);
      if(details>0) 
        printf("%s%s\n","Pos.       < number1             > number2",  
               "            rel.diff    abs.diff");
      firstline=0;
     }
     if(details>0)
       printf("L%dC%d %.15e %.15e  r=%1.2e  a=%1.2e\n",
              line,col,x1,x2,reldiff,absdiff);  
     if(details>1) printf("%s\n%s\n",str1,str2);
     exerr++;
    } 
   }/* end for */
  } /* end if */
 }
 
 fclose(difffile);
 
 if(exerr>0) 
 {
  printf("\n=> %d relevant difference(s), and %d warning(s) found! \n",
         exerr, exwarn);
  return 2;
 }
 else if(exwarn>0) 
 {
  printf("\n=> %d warning(s)\n",exwarn);
  return 1;
 }
 else return 0;
}

