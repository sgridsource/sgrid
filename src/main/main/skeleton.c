/* skeleton.c */
/* Wolfgang Tichy, April 2005  &  Bernd Bruegmann 12/99 */

/* better name? fps function pointer skeleton */

#include "sgrid.h"
#include "main.h"

/* global var that tells us how often we have restarted sgrid */
extern int sgrid_restarts;


typedef struct tTODO {
  struct tTODO *next;
  int (*f)(tGrid *);
  char *name;
} tTodo;

tTodo *fps[NFUNCTIONS];



void AddFun(int step, int (*f)(tGrid *), char *name)
{
  tTodo *t;

  /* do not add func again if this is a restart of sgrid */
  if(sgrid_restarts>0) return;

  if (1) printf("  function  %s\n", name);

  if (!fps[step]) fps[step] = (tTodo *) calloc(1, sizeof(tTodo));

  for (t = fps[step]; t->next; t = t->next);
  t->next = (tTodo *) calloc(1, sizeof(tTodo));
  t->f = f;
  t->name = (char *) calloc(strlen(name)+1, sizeof(char));
  strcpy(t->name, name);
}

/* remove all functions from all func bins */
void remove_all_Funs(void)
{
  int step;

  for(step=0; step<NFUNCTIONS; step++)
  {
    tTodo *t, *p=NULL;

    for(t = fps[step]; t; )
    {
      //printf("step%d %p %s\n", step, t, t->name);
      p = t;         /* previous entry */
      t = t->next;   /* get next entry */
      free(p->name); /* free previous name */
      free(p);       /* free previous entry */
    }
  }
}


void RunFun(int step, tGrid *grid) 
{
  tTodo *t;

  if (!fps[step]) return;

  for (t = fps[step]; t->next; t = t->next) {
    (*(t->f))(grid);
  }
}



/* print functions in func bin step */
void PrintFun(int step)
{
  tTodo *t;

  if(!fps[step]) return;

  printf("func_T%02d:\n", step);
  for(t = fps[step]; t->next; t = t->next)
  {
    printf("          %s\n", t->name);
  }
}

/* print functions in all function bins */
void PrintFuncs(void)
{
  int step;

  prdivider(0);
  printf("Function bin skeleton: functions are called in this order:\n");
  for(step=0; step<NFUNCTIONS; step++) PrintFun(step);
  prdivider(0);
}
