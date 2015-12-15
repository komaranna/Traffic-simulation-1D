#include <iostream>
#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#define N 1000
#define alpha_max 1
#define beta_max 1
#define alpha_min 0.001
#define beta_min 0.001
#define alpha_N 20
#define beta_N 20



using namespace std;

void exit(string message, int exitcode);
double *vector_double(long m);
int *vector_int(long m);
double **matrix_double(long m, long n);
int **matrix_int(long m, long n);
string convertInt2Str(int number);
void *simulate ( void *arg );



////Global variable
   double Tlength, Tavstart; //Time of simulation
   int *threadind;	//Thread index
   double *vector_alpha; //Vector containing values of alpha (incoming rate)
   double *vector_beta; //Vector containing values of beta (leaving rate)
   double *vector_phase; //Parameters defining the phase (J, current)
   double *density_avr; //Average density (rho)
   double **density_profile; //Stationary density profiles (rho_i)


int main(int argc, char *argv[]) {

	
   int nprocs=0; //Number of threads
   int i,j;	//Generic integers
   pthread_t *pth; //Thread pointers
   pthread_attr_t pth_attr;
   clock_t start, finish; // Time variables
   ofstream myfile; //Output file handle
   

if ( argc != 3 ) // 2 arguments are the default (nprocs, Tlength)
    {
        exit("Not appropriate amount of arguments", 1);
    }
    else 
    {
        nprocs = atoi(argv[1]); // Setting the number of threads to use
        Tlength = strtod(argv[2],NULL); // Setting the time (length) of simulation
	}
    
    //Allocating memory for variables, and setting the array to zero
    vector_alpha = vector_double(alpha_N);
    vector_beta = vector_double(beta_N);
    vector_phase = vector_double(alpha_N*beta_N);
    density_avr = vector_double(alpha_N*beta_N);
    density_profile = matrix_double(alpha_N*beta_N, N);

    //Pointer allocation for threads
    pth = ( pthread_t * ) malloc ( nprocs * sizeof ( pthread_t ) );
    
    //Start the clock!
    start = time(NULL);
    printf("Number of threads: %i \n",nprocs);
    printf("Time length: %f \n",Tlength);
    
    //set the attr to be joinable thread
    pthread_attr_init(&pth_attr);
    pthread_attr_setdetachstate(&pth_attr, PTHREAD_CREATE_JOINABLE);

    Tavstart=Tlength*0.75; //Setting the start of averaging
    
//    //Index array
    threadind = ( int * ) malloc ( alpha_N * beta_N * sizeof ( int ) );
    
    // Assigning the alpha and beta values to the vectors
    for(i=0; i<alpha_N; i++)
    {
    vector_alpha[i] = alpha_min + (double)i*((double)alpha_max - alpha_min)/((double)alpha_N);
    }
    for(i=0; i<beta_N; i++)
    {
    vector_beta[i] = beta_min + (double)i*((double)beta_max - beta_min)/((double)beta_N);
    }
    
    
    j=0;
    while(j < alpha_N*beta_N-((alpha_N*beta_N)%nprocs)){ //Will repeat until most of the graphs are ready
    //Starting threads
         for ( i = 0; i < nprocs; i++ ){	  
    	    threadind[j]=j;
   	        pthread_create(&pth[i],&pth_attr,simulate,&threadind[j]);
   	        j+=1;
         }
	
	//Waiting for threads
 	    for ( i = 0; i < nprocs; i++ ){
		    (void) pthread_join(pth[i], NULL);
	    }
	}
	

    while(j < alpha_N*beta_N){ //rest of the graphs
    //Starting thread	  
    	    threadind[j]=j;
   	        pthread_create(&pth[i],&pth_attr,simulate,&threadind[j]);
   	        j+=1;	
	//Waiting for thread
		    (void) pthread_join(pth[i], NULL);
	    
	}

	
	// Exporting data for each alpha and beta pair: alpha, beta, J (current), average density (rho)
	myfile.open("phase_diagram.txt");
    for ( i = 0; i < alpha_N; i++ )
    {
        for ( j = 0; j < beta_N; j++ )
        {
		myfile <<vector_alpha[i]<<"\t"<<vector_beta[j]<<"\t"<<vector_phase[alpha_N*j+i]<<"\t"<<density_avr[alpha_N*j+i]<<"\n";
        }
	}
    myfile.flush();
    myfile.close();
	
		
	//Exiting
	finish = time(NULL);
    printf("Total elapsed time: %i minutes %i seconds\n",(int)((long)(finish - start)/60),(int)((long)(finish - start) % 60));
    printf("Press any key to continue!");
    getchar();
    return 0;
    
}

//---------------------------------------------------------------------    
// Simulating a system for fix alpha and beta

void *simulate( void *arg){
  int threadID; //ID of the thread, handed down from main
  long i,k;	//Generic long integers for stepping
  double x, rand_alpha, rand_beta, u, deltat;	//Generic doubles
  double minimum; //Time of next step out of all points
  int minind; //index of the next step
  double t; //Time
  double alpha, beta; //values of alpha and beta
  int *systems; //Matrix containing the system (1 means there is a car in that position, 0 means no car)
  double *next_step; //Time of next step for each point
  
  ofstream myfile; //Output file handle
  string filename; //Filename string
  
//  //Assigning threadID from input
  threadID=*((int*)arg);
  
  alpha =  vector_alpha[threadID%alpha_N]; //Defining alpha for the calculations
  beta = vector_beta[(int)(threadID/alpha_N)]; //Defining beta for the calculations
    
      //Initializing random number generator
    srand(time(NULL)+threadID);
    
    //Throwing away the first 10000 numbers, as it is customary
    for (i = 0; i < 10000; i++) {
        u = rand();
    }
    
    u=(double)rand()/((double)RAND_MAX + 1); //Uniform random numbers
    x=-log(1-u); //Exponential random numbers, lambda=1
    rand_alpha=-log(1-u)/alpha; //lambda=alpha
    rand_beta=-log(1-u)/beta;    //lambda=beta

  
  //Allocating memory for variables, and setting the array to zero
  systems=vector_int(N+1); //zero array means the initial condition is an empty road
  next_step=vector_double(N+1);
  
  
  next_step[0]=0+rand_alpha; //Ensuring the first step is an incoming car
  for(i=1; i<N+1; i++)
  {
           next_step[i]=Tlength+2; //Excluding steps which are not possible (there is no car in that position)
  }

  t=0; //time
  while(t < Tlength) //simulation is done for Tlength seconds
  {
          //Finding the next event
          minimum=next_step[0];
          minind=0;
          for(i=1; i<N+1; i++) //for all points of the road
          {
           if(next_step[i] < minimum) //if the time of the examined event is smaller than the minimum, reassign the minimum
           {
            minimum=next_step[i];
            minind=i;
           }
          }
          deltat = minimum - t; //time length between the last and the new step
          t=minimum; //Setting time to the time of next event
          
          
          if(minind==N)//Last car is leaving
          {
           systems[N]=0;
           next_step[N]=Tlength+2; //Excluding this point from the possible steps
          }
          else //Any other car is moving
          {          
              if(systems[minind+1] == 0) //There is space for the step
              {
                     if((minind != 0) && (minind != N-1)) //Steps in between the boundaries
                     {
                      systems[minind]=0;
                      systems[minind+1]=1;
                      next_step[minind]=Tlength+2; //Excluding the new void point from the possible steps

					  u=(double)rand()/((double)RAND_MAX + 1); //Uniform random number
					  x=-log(1-u); //Exponential random number, lambda=1
                      next_step[minind+1]=t+x; //Defining a new step time for the car
                     }
                     else 
                     {
                          if(minind == 0) //Incoming car
                          {
                           systems[0]=1; // #0 is a fictional position, there is always a car there, waiting to come
                           systems[1]=1;
						   u=(double)rand()/((double)RAND_MAX + 1); //Uniform random number
						   rand_alpha=-log(1-u)/alpha; //lambda=alpha
                           next_step[0]=t+rand_alpha; //Defining a new time for the incoming car
                           
						   u=(double)rand()/((double)RAND_MAX + 1); //Uniform random number
						   x=-log(1-u); //Exponential random number, lambda=1
                           next_step[1]=t+x; //Defining a new step time for the car already inside
                          }
                          else //Car stepping to the last position
                          {
                              systems[N-1]=0;
                              systems[N]=1;
                              next_step[N-1]=Tlength+2; //Excluding the void position from the possible steps

							  u=(double)rand()/((double)RAND_MAX + 1); //Uniform random number
						      rand_beta=-log(1-u)/beta;    //lambda=beta
                              next_step[N]=t+rand_beta; //Defining the new step time for the last (leaving) car
                          }
                     }
              }
              else //There is no space for the step
              {
                   if(minind != 0) //Steps in between the boundaries
                     {
                      //No need to change the system
					  u=(double)rand()/((double)RAND_MAX + 1); //Uniform random number
					  x=-log(1-u); //Exponential random numbers, lambda=1
                      next_step[minind]=t+x; //Only redefining the time of the next step
                     }
                     else //Incoming car
                     {
                      //No need to change the system
					  u=(double)rand()/((double)RAND_MAX + 1); //Uniform random number
					  rand_alpha=-log(1-u)/alpha; //lambda=alpha
                      next_step[0]=t+rand_alpha; //Only redefining the time for the incoming car
                     }
                   
              }
               
               
           }
                      
           
   if(t>Tavstart) //the time reached the value where averaging will start
   {
    //time averaging the stationary density profiles (only adding their weighted values)
    for(k=0; k<N; k++)
    {
    density_profile[threadID][k] = density_profile[threadID][k] + (double)(systems[k+1])*deltat;
    }
   }




 }
 
   //Finishing the time average: denoting the weighted sum by the whole length of averaging
   for(k=0; k<N; k++)
   {
    density_profile[threadID][k] = density_profile[threadID][k]/(Tlength-Tavstart);
   }
  
  
    //Calculating the average density of the stationary profile and the J current
    density_avr[threadID] = density_profile[threadID][0]/N;
    for(k=0; k<N-1; k++)
    {
    vector_phase[threadID] = vector_phase[threadID] + (density_profile[threadID][k] * (1 - density_profile[threadID][k+1]))/(N-1);
    density_avr[threadID] = density_avr[threadID] + density_profile[threadID][k+1]/N;
    }
    
    //Exporting the stationary density profile: first row is alpha and beta, then position (i) and density (rho_i)
    filename = "densityprofile" + convertInt2Str(threadID) + ".txt";
    myfile.open(filename.c_str());
    myfile <<alpha<<"\t"<<beta<<"\n";
    for ( i = 0; i < N; i++ )
    {
        myfile <<i+1<<"\t"<<density_profile[threadID][i]<<"\n";
	}
    myfile.flush();
    myfile.close();
    
    
        	
     free(next_step);
	free(systems);

	printf("Thread %i ready\n",threadID);
    pthread_exit((void*) threadID);		
}
    

    




//--------------------------------------------------------------------------------
//General routines, useful for everything

//Writing exit messages
void exit(string message, int exitcode = 0) {
    cout << "Application exited with error: " << endl << message << endl;
    exit(exitcode);
}

//Converting integer to string
string convertInt2Str(int number) {
    stringstream ss;
    ss << number;
    return ss.str();
}

//Allocating vector of integers, size: m
int *vector_int(long m) {
    int *res = (int *) calloc(m, sizeof (int));
    return (res);
}

//Allocating matrix of doubles, size: m x n
double **matrix_double(long m, long n) {
    double **res = (double **) calloc(m, sizeof (double *));
    long i;
    for (i = 0; i < m; i++)
        res[i] = (double *) calloc(n, sizeof (double));
    return (res);
}

//Allocating matrix of integers, size: m x n
int **matrix_int(long m, long n) {
    int **res = (int **) calloc(m, sizeof (int *));
    long i;
    for (i = 0; i < m; i++)
        res[i] = (int *) calloc(n, sizeof (int));
    return (res);
}

//Allocating vector of doubles, size: m
double *vector_double(long m) {

    double *res = (double *) calloc(m, sizeof (double));
    return (res);
}



