//#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <stdlib.h>

#define STAT_WAITING       1
#define STAT_WORKING       2
#define STAT_FINISHED      3
#define STAT_NO_MORE       4

struct t_params
{
	int    id;
	int    status;	
	int    iloop;
	int    ne_loc;
	int    ne;
	int    nt;
	double *fce;
	double *qce;
	double *uce;
	double *vce;
	void   (*emissivity)();
	void   (*element)();
};

static pthread_cond_t  thread_work = PTHREAD_COND_INITIALIZER;
static pthread_cond_t  main_work = PTHREAD_COND_INITIALIZER;
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

static pthread_t       *th;
static struct t_params *params;


void *thread(void *arg){

 struct t_params *my_params;

 my_params = (struct t_params*)arg;

 pthread_mutex_lock(&lock);
 while(my_params->status != STAT_NO_MORE){
  if(my_params->status == STAT_WORKING){
   pthread_mutex_unlock(&lock);
   my_params->element(my_params->iloop,my_params->ne_loc,my_params->ne,
            my_params->nt,my_params->fce,my_params->qce,my_params->uce,
            my_params->vce,my_params->emissivity);
   pthread_mutex_lock(&lock);
   my_params->status = STAT_FINISHED;
   pthread_cond_signal(&main_work);
  }else pthread_cond_wait(&thread_work, &lock);
 };
 pthread_mutex_unlock(&lock);
 pthread_exit(NULL);
}


void create_ide_threads(const int nthreads, const int nthreads_old, 
                        const int ne, const int ne_old,
                        const int nt, const int nt_old, 
                        const int polar, const int polar_old){

int   i;
void  *return_pointer;

// Let's free pointers if we change nthreads, ne, nt or polar
 if(nthreads_old > 1 && (nthreads != nthreads_old || ne != ne_old || 
                         nt != nt_old || (polar == 0 && polar_old == 1))){
  pthread_mutex_lock(&lock); 
  if(nthreads != nthreads_old || ne != ne_old || nt != nt_old) 
   for(i = 0; i < nthreads_old; i++) free((double *) params[i].fce);
  if(polar_old) for(i = 0; i < nthreads_old; i++){
    free((double *) params[i].qce);
    free((double *) params[i].uce);
    free((double *) params[i].vce);
  };
  pthread_mutex_unlock(&lock); 
  if(nthreads != nthreads_old){
   pthread_mutex_lock(&lock); 
   for(i = 0; i < nthreads_old; i++) params[i].status = STAT_NO_MORE;
   pthread_cond_broadcast(&thread_work);
   pthread_mutex_unlock(&lock); 
   for(i = 0; i < nthreads_old; i++) pthread_join(*(th + i), &return_pointer);
   free(th);
   free(params);
  };
 };
// Let's create pointers if we change nthreads, ne, nt or polar
 if(nthreads > 1){
  if(nthreads != nthreads_old){
   th = (pthread_t *)malloc(nthreads * sizeof(pthread_t));
   params = (struct t_params *)malloc(nthreads * sizeof(struct t_params));
   for(i = 0; i < nthreads; i++){
    params[i].id = i;
    params[i].status = STAT_WAITING;
    pthread_create(th + i, NULL, thread, (void *)(params + i));
   }; 
  };
  if(nthreads != nthreads_old || ne != ne_old || nt != nt_old)
   for(i = 0; i < nthreads; i++){
    params[i].ne = ne;
    params[i].nt = nt;
    params[i].fce = (double *)malloc(ne*nt*sizeof(double));
   };
  if(polar && (nthreads != nthreads_old || ne != ne_old || nt != nt_old || 
     polar_old == 0)){
   for(i = 0; i < nthreads; i++){
    params[i].qce = (double *)malloc(ne*nt*sizeof(double));
    params[i].uce = (double *)malloc(ne*nt*sizeof(double));
    params[i].vce = (double *)malloc(ne*nt*sizeof(double));
   };
  };
 };
}


void ide_threads(const int nthreads, const int nloops, const int polar, 
                 const int ne_loc, const int ne, const int nt, double *fc, 
                 double *qc, double *uc, double *vc, void (*emissivity)(), 
                 void (*element)()){
                  
 int       i, ie, it;
 int       iloop, finished;
 
 finished = 0;
 iloop = 1;
 pthread_mutex_lock(&lock);
 for(i = 0; i < nthreads; i++){
  params[i].status = STAT_WAITING;
  params[i].ne_loc = ne_loc;
  for(ie = 0; ie < ne; ie++) for(it = 0; it < nt; it++){
   *(params[i].fce+ie+ne*it)=0.;
   if(polar){
    *(params[i].qce+ie+ne*it)=0.;
    *(params[i].uce+ie+ne*it)=0.;
    *(params[i].vce+ie+ne*it)=0.;
   };
  };
  params[i].emissivity=emissivity;
  params[i].element=element;
 };
 while(finished < nloops){
  for(i = 0; (i < nthreads) && (iloop <= nloops); i++) 
   if(params[i].status != STAT_WORKING){
    params[i].status = STAT_WORKING;
    params[i].iloop = iloop;
    iloop++;
  };
  pthread_cond_broadcast(&thread_work);
  pthread_cond_wait(&main_work, &lock);
  for(i = 0; i < nthreads; i++) if(params[i].status == STAT_FINISHED){
   finished++;
   params[i].status = STAT_WAITING;
  };
 };
 for(i = 0; i < nthreads; i++){
  for(ie = 0; ie < ne; ie++) for(it = 0; it < nt; it++){
   *(fc+ie+ne*it) += *(params[i].fce+ie+ne*it);
   if(polar){
    *(qc+ie+ne*it) += *(params[i].qce+ie+ne*it);
    *(uc+ie+ne*it) += *(params[i].uce+ie+ne*it);
    *(vc+ie+ne*it) += *(params[i].vce+ie+ne*it);
   };
  };
 };
 pthread_mutex_unlock(&lock);
 return;
}
