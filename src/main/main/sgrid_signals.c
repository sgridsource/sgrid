/* sgrid_signals.c.c */
/* Wolfgang Tichy, Dec 2010 */

#include "sgrid.h"
#include "main.h"
#include <signal.h>



/* handler function for SIGFPE */
/* To use it we need something like this in in e.g. main:
...
#include <fenv.h>  
#include <signal.h>

...
  struct sigaction sa;

  sa.sa_flags = SA_SIGINFO;
  sigemptyset(&sa.sa_mask);
  sa.sa_sigaction = sgrid_SIGFPE_handler;
  //sigaction(SIGSEGV, &sa, NULL);
  sigaction(SIGFPE, &sa, NULL);

  feenableexcept(FE_ALL_EXCEPT);
*/
static void sgrid_SIGFPE_handler(int sig, siginfo_t *si, void *unused)
{
  printf("sgrid_SIGFPE_handler: Got SIGFPE at address: 0x%lx\n",
         (long) si->si_addr);
  exit(EXIT_FAILURE);
}
