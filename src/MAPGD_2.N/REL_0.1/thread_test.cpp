#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>

using namespace std;

#define NUM_THREADS 5

typedef struct str_thdata
{
    int no;
    char message[100];
    int ret;
} thdata;

void *wait(void *t)
{
	str_thdata *data;

	data = (thdata *) t;

	sleep(5);
	cout << "Sleeping " << data->no << " ... " << data->message << endl;
	data->ret=data->no;
	pthread_exit(NULL);
}

int main ()
{
	int rc;
	int i;
	int sum=0;
	pthread_t threads[NUM_THREADS];
	int list[]={1,2,3,4,5,6,7,8,9,10};
	int pos=0;
	pthread_attr_t attr;
	thdata data[NUM_THREADS];
	for( i=0; i < NUM_THREADS; i++ ){
		data[i].no=i;
		sprintf(data[i].message, " goodbye!");
 	}

	void *status;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
 
 for( i=0; i < NUM_THREADS; i++ ){
 cout << "main() : creating thread, " << i << endl;
 data[i].no=list[pos];
 rc = pthread_create(&threads[i], NULL, wait, (void *)&data[i] );
 pos+=1;
 if (rc){
 cout << "Error:unable to create thread," << rc << endl;
 exit(-1);
 }
 }
 
 pthread_attr_destroy(&attr);

 bool running=true;
 while(running){
 running=false;
 for( i=0; i < NUM_THREADS; i++ ){
 rc = pthread_join(threads[i], &status);
 if (rc){
 cout << "Error:unable to join," << rc << endl;
 exit(-1);
 }
 cout << "Main: completed thread id :" << i ;
 cout << " exiting with status :" << status << endl;
 cout << " return value :" << data[i].ret << endl;
 sum+=data[i].ret;
 if (pos<10){
 running=true;
 data[i].no=list[pos];
 rc = pthread_create(&threads[i], NULL, wait, (void *)&data[i] );
 pos+=1;
 pthread_attr_init(&attr);
 pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
 if (rc){
 cout << "Error:unable to create thread," << rc << endl;
 exit(-1);
 }
 }
 }
 }
 
 cout << "Main: program exiting : sum " << sum << endl;
 pthread_exit(NULL);
 }
