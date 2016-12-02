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


void RunFun(int step, tGrid *grid) 
{
  tTodo *t;

  if (!fps[step]) return;

  for (t = fps[step]; t->next; t = t->next) {
    (*(t->f))(grid);
  }
}
