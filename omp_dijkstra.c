/* Program     : Dijkstra's Algorithm               *
 * Description : To find the shortest path from     *
 *               source node to all other vertices. *
 * Usage       : Keep the input files as mentioned  *
 *               in report in the same directory    *
 *               where the binary is placed.        *
 *               ./binary_name #threads file_name   *
 * ------------------------------------------------ */              
  
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<omp.h>
#define INFINITY 9999
#define MAX_LINE_SIZE 10000
#define RAND_NUM 5
 
void omp_dijkstra(int **G, int th, double *eta, int *dist, int *pre);
void serial_dijkstra(int **G,int *dist, int *pre);
void printMatrix(void);
void printVector(int *arr,char *str);
void print_shortest_path(int *dist, int *pre);
void get_matrix_type1(char *file_name);
void get_matrix_type2(char *file_name);
void get_matrix_type3(char *file_name);
void usage(char *prog_name);

int rows, columns, src_node, dyn_size; 
double omp_time;
int **M;
int firstline =1;

int main(int argc, char* argv[])
{
    int i, choice, file_ch, loop;
    char ch;
    double ser_time;
    char *file_name[]={"soc-sign-bitcoinalpha.csv",
                       "artist_edges.csv",
                       "facebook_combined.txt",
                       "soc-sign-Slashdot081106.txt",
                       "soc-Slashdot0902.txt",
                       "file_data1.txt",
                       "file_data.txt"};
   
    if(argc != 3)
      usage(argv[0]);

    /* This is meant for comparing the file and calling relevant filehandling function */
    if((strncmp(argv[2],file_name[0],strlen(file_name[0])))==0)
       get_matrix_type2(argv[2]);
    else if((strncmp(argv[2],file_name[1],strlen(file_name[1])))==0)
       get_matrix_type2(argv[2]);
    else if((strncmp(argv[2],file_name[2],strlen(file_name[2])))==0)
       get_matrix_type3(argv[2]);
    else if((strncmp(argv[2],file_name[3],strlen(file_name[3])))==0)
       get_matrix_type3(argv[2]);
    else if((strncmp(argv[2],file_name[4],strlen(file_name[4])))==0)
       get_matrix_type3(argv[2]);
    else if((strncmp(argv[2],file_name[5],strlen(file_name[5])))==0)
       get_matrix_type1(argv[2]);
    else if((strncmp(argv[2],file_name[6],strlen(file_name[6])))==0)
       get_matrix_type1(argv[2]);
    else
       printf("File_name did not match\n");
   
   /* It will print the adjacency Matrix */ 
   //printMatrix();    
    
    /* Find the thread_value */ 
    int thread_count = strtol(argv[1],NULL,10);
    if(thread_count <=0)
      usage(argv[0]);
   
    /* Allocating memory for distance and predecessor */ 
    int *omp_distance = (int *)malloc(rows*sizeof(int));
    int *omp_pred = (int *)malloc(rows*sizeof(int));

    /* Allocating a vector to record time for each thread inside the parallel region */
    double *execution_time_array=(double *)malloc(thread_count*sizeof(double));
 
    /* intializing the execution time array */
    for(i=0;i<thread_count;i++)
    {
      execution_time_array[i]=0.0;
    }

    /* The Matrix is suppose to be square matrix */
    if(rows!=columns) /* Condition check */
    {
      printf("The Input Matrix data is not square matrix\nCannot Proceed.......\n");
      /* free Resources before exiting */
      free (execution_time_array);
      free (omp_distance);
      free (omp_pred);
      free (M); 
      return 0;
    }
    printf("The starting node is %d\n",src_node);

    omp_time = 0.0;
    omp_time = omp_get_wtime(); 
    omp_dijkstra(M, thread_count, execution_time_array, omp_distance,omp_pred);  /* function call to Dijkstra */
    omp_time = omp_get_wtime() - omp_time;

    /* print the path and distance of each node */
    //print_shortest_path(omp_distance,omp_pred);

    /* Print the execution time of each thread */  
    printf("\n\nExecution time of each thread\n");
    for(i=0;i<thread_count;i++)
    {
      //printf("Thread %d: %lf\n",i,(avg_time_array[i]/10));
      printf("Thread %d: %lf\n",i,execution_time_array[i]);
    }
   
    printf("overall: %lf\n",omp_time);

    /* allocating the memory for the serial code */
    int *ser_distance = (int *)malloc(rows*sizeof(int));
    int *ser_pred = (int *)malloc(rows*sizeof(int));
     
    ser_time = omp_get_wtime();
    serial_dijkstra(M,ser_distance,ser_pred);      
    ser_time = omp_get_wtime()-ser_time;

    /* print the path and distance of each node */
    //print_shortest_path(ser_distance,ser_pred);
      
    printf("Total Execution time of serial code:%lf\n",ser_time);
    printf("Comparing the output of serial and parallel code\n");

    /* Comparing the distances for both serial and parallel code */
    for(i=0;i<rows;i++)
    {
      if(omp_distance[i]!=ser_distance[i])
      {
        printf("Value is not matching at position <%d> omp_distance[i]: %d, ser_distance[i]: %d\n",i,omp_distance[i],ser_distance[i]);
      }
    }
     
    printf("Total Speed up is: %lf\n",ser_time/omp_time);
    printf("efficiency is: %lf\n",((ser_time/omp_time)/thread_count));
    free (ser_distance);
    free (ser_pred);
 
    free (execution_time_array);
    free (omp_distance);
    free (omp_pred);
    free (M); 
    return 0;
}

void get_matrix_type2(char*file_name)
{
  int i,j;
  int index_i, index_j,edge_weight;
  char line[100];
  FILE *Inputfile;
  int max=0,max1=0;
  Inputfile = fopen(file_name,"r");
  if(Inputfile==NULL)
  {
    printf("Error in opening file\n");
    exit(0);
  }
  while(1)
  {
    if (fgets(line, 100, Inputfile) == NULL)
       break;
    /* Read the first line as the src_node */
    if(firstline ==1)
    {
      if(sscanf(line,"%d",&src_node)!=1)
         printf("Error in reading first line\n");
      firstline = 0;
      continue;
    }
    if (sscanf(line,"%d,%d", &index_i, &index_j) != 2)
      printf("Error in reading file\n");
   /*  printf("<%d> <%d>\n",index_i, index_j);*/
      if(max<index_i)
         max=index_i;
      if(max1<index_j)
         max1=index_j;
  }
  //printf("The max value is %d %d\n",max,max1);
  if(max<max1)
     max=max1;
  printf("The max node value is %d\n",max);
  columns = max+1;
  rows = max+1; /* Populating the value of rows after reading the maximum value */
  fclose(Inputfile);

  /* Allocating the space for Adjacency Matrix */
  M = (int **)malloc(rows*sizeof(int *));
  if(M==NULL){
    printf("Out of Memory\n");
    exit(-1);
  }
  for(i=0;i<rows;i++)
  {
    M[i] = (int *)malloc(rows*sizeof(int));
    if(M[i]==NULL){
      printf("Out of Memory\n");
      exit(-1);
    }
    for(j=0;j<rows;j++)
    {
      M[i][j]=0;
    }
  }
  firstline =1;
  Inputfile = fopen(file_name,"r");
  if(Inputfile==NULL)
  {
    printf("Error in opening file\n");
    exit(0);
  }
  /*printf("Read the file\n");*/
  while(1)
  {
    if (fgets(line, 100, Inputfile) == NULL)
       break;
    /* Redundant line just to avoid the first line although we have already captured it in first fetch */
    if(firstline ==1)
    {
      if(sscanf(line,"%d",&src_node)!=1)
         printf("Error in reading first line\n");
      firstline = 0;
      continue;
    }
    if (sscanf(line,"%d,%d", &index_i, &index_j) != 2)
      printf("Error in reading file\n");
   /*   printf("<%d> <%d> <%d>\n",index_i, index_j, edge_weight);*/
      M[index_i][index_j]=rand() % RAND_NUM;
  }
  printf("Matrix function type2 completed successfully\n"); 
  fclose(Inputfile); 
}

void get_matrix_type3(char*file_name)
{
  int i,j;
  int index_i, index_j,edge_weight;
  char line[100];
  FILE *Inputfile;
  int max=0,max1=0;
  Inputfile = fopen(file_name,"r");
  if(Inputfile==NULL)
  {
    printf("Error in opening file\n");
    exit(0);
  }
  while(1)
  {
    if (fgets(line, 100, Inputfile) == NULL)
       break;
    /* Read the first line as the src_node */ 
    if(firstline ==1)
    {
      if(sscanf(line,"%d",&src_node)!=1)
         printf("Error in reading first line\n");
      firstline = 0;
      continue;
    }
    if (sscanf(line,"%d %d", &index_i, &index_j) != 2)
      printf("Error in reading file\n");
   /*  printf("<%d> <%d>\n",index_i, index_j);*/
      if(max<index_i)
         max=index_i;
      if(max1<index_j)
         max1=index_j;
  }
  //printf("The max value is %d %d\n",max,max1);
  if(max<max1)
     max=max1;
  printf("The max node value is %d\n",max);
  columns = max+1;
  rows = max+1; /* Populating the value of rows after reading the maximum value */
  fclose(Inputfile);

  /* Allocating the space for Adjacency Matrix */
  M = (int **)malloc(rows*sizeof(int *));
  if(M==NULL){
    printf("Out of Memory\n");
    exit(-1);
  }
  for(i=0;i<rows;i++)
  {
    M[i] = (int *)malloc(rows*sizeof(int));
    if(M[i]==NULL){
      printf("Out of Memory\n");
      exit(-1);
    }
    for(j=0;j<rows;j++)
    {
      M[i][j]=0;
    }
  }
  firstline =1;
  Inputfile = fopen(file_name,"r");
  if(Inputfile==NULL)
  {
    printf("Error in opening file\n");
    exit(0);
  }
  /*printf("Read the file\n");*/
  while(1)
  {
    if (fgets(line, 100, Inputfile) == NULL)
       break;
    /* Redundant line just to avoid the first line although we have already captured it in first fetch */
    if(firstline ==1)
    {
      if(sscanf(line,"%d",&src_node)!=1)
         printf("Error in reading first line\n");
      firstline = 0;
      continue;
    }
    if (sscanf(line,"%d %d", &index_i, &index_j) != 2 )
      printf("Error in reading file\n");
   /*   printf("<%d> <%d> <%d>\n",index_i, index_j, edge_weight);*/
      M[index_i][index_j]=rand() % RAND_NUM;
  }
  printf("Facebook Matrix function completed successfully\n");
  fclose(Inputfile);
}

void get_matrix_type1(char *file_name)
{
  int i,j,firstline=1,count;
  FILE *Inputfile;
  char help[MAX_LINE_SIZE], *token, ch;
  
  /* Open the file to read the number of lines so as to allocate the array dynamically */
  //Inputfile = fopen("file_data.txt","r");
  Inputfile = fopen(file_name,"r");
  if(Inputfile==NULL)
  {
    printf("Error in opening file\n");
    exit(0);
  }
  while(1)  /* Continuous loop that breaks at EOF and counts the newline encountered */
  {
    ch = fgetc(Inputfile);
    if(ch == '\n')
    {
      count++;
    }
    if(feof(Inputfile))
      break;
  }
  dyn_size = count-1;
  printf("The count is %d\n",dyn_size);
  fclose(Inputfile);   /* Close the file as count is done */
 
  /* Allocate the 2-D array depending on the count value */
  M = (int **)malloc(dyn_size*sizeof(int *));
  for(i=0;i<dyn_size;i++)
  {
    M[i] = (int *)malloc(dyn_size*sizeof(int));  //Allocates memory for column in each row
    for(j=0;j<dyn_size;j++)
    {
      M[i][j]=0;
    }
  }

  /* Read the file again to read the content */
  Inputfile = fopen(file_name,"r");
  //Inputfile = fopen("file_data.txt","r");
  if(Inputfile==NULL)
  {
    printf("Error in opening file\n");
    exit(0);
  }

  i=0;
  while(1)
  {
    if(firstline ==1)
    {
      fscanf(Inputfile,"%d",&src_node);
      firstline = 0;
    }
    fscanf(Inputfile,"%s",help);
    if(feof(Inputfile))
      break;
   // printf("Read line %d: <%s>\n",i,help);
    token = strtok(help, ",");
    j=0;
    while(token!=NULL)
    {
      M[i][j]=atoi(token);
      token = strtok(NULL,",");
      j++;
    }
    i++;
  }
  fclose(Inputfile);
  rows = i;
  columns = j;
  printf("Successfully loaded random matrix\n");
}

void usage(char *prog_name)
{
  fprintf(stderr, "-----USAGE: %s <thread_count> [filename]\n",prog_name);
  fprintf(stderr, " Dijkstra's Algorithm: To find the shortest path from\n a given source node to rest of the nodes\n");
  fprintf(stderr, "  ----<thread_count> should be positive\n");
  exit(0);
}

void printMatrix()   /* To print the MAtrix */
{
  int i, j;
  printf("The value of rows: <%d> and column: <%d>\n",rows,columns);
  for (i=0;i<rows;i++)
  {
    for(j=0;j<columns;j++)
    {
      printf("%ld\t",M[i][j]);
    }
    printf("\n");
  }
}

void printVector(int *arr,char *str)   /* To print the MAtrix */
{
  int i;
  printf("<%s> Vector\n",str);
  for (i=0;i<rows;i++)
  {
      printf("%d\t",arr[i]);
  }
  printf("\n");
}

void print_shortest_path(int *dist, int *pre)
{
  int i,j;

  for(i=0;i<rows;i++)
  {
    if(i!=src_node)
    {
       printf("\nDistance of node <%d> from node <%d> =%d",i,src_node,dist[i]);
       printf("\nPath=%d",i);
       
       j=i;
       do
       {
           j=pre[j];
           printf("<-%d",j);
       }while(j!=src_node);
    }
  }
}

void serial_dijkstra(int **G,int* distance, int* pred)
{

    int count,mindistance,nextnode,i,j;
    /* Allocated the memory for few of the local variables */
    int **cost = (int **)malloc(rows*sizeof(int *));
    int *visited = (int *)malloc(rows*sizeof(int));

    /* Initialized all the allocated variable */
    for(i=0;i<rows;i++)
    {
      cost[i] = (int *)malloc(rows*sizeof(int));  //Allocates memory for column in each row
      for(j=0;j<rows;j++)
      {
        cost[i][j]=0;
      }
    }
    /* Intitalize the allocated array */
    for(i=0;i<rows;i++)
    {
      distance[i]=0;
      pred[i]=0;
    }

    /* pred[] stores the predecessor of each node
 *     count gives the number of nodes seen so far
 *         create the cost matrix */
    for(i=0;i<rows;i++)
        for(j=0;j<rows;j++)
            if(G[i][j]==0)
                cost[i][j]=INFINITY;
            else
                cost[i][j]=G[i][j];   /* The matrix value will become the weight of the edge */

    /* initialize pred[],distance[] and visited[] with the required value */
    for(i=0;i<rows;i++)
    {
        distance[i]=cost[src_node][i];
        pred[i]=src_node;
        visited[i]=0;
    }

    distance[src_node]=0;
    visited[src_node]=1;
    /*count=1;*/

  /*  while(count<n-1)*/
    for(count=1;count<rows-1;count++)
    {
        mindistance=INFINITY;

        /* nextnode gives the node at minimum distance */
        for(i=0;i<rows;i++)
            if(distance[i]<mindistance&&!visited[i])
            {
                mindistance=distance[i];
                nextnode=i;
            }

            /* check if a better path exists through nextnode */
            visited[nextnode]=1;
            for(i=0;i<rows;i++)
                if(!visited[i])
                    if(mindistance+cost[nextnode][i]<distance[i])
                    {
                        distance[i]=mindistance+cost[nextnode][i];
                        pred[i]=nextnode;
                    }
      /* count++;*/
    }

    free (cost);
    free (visited);
}
 
void omp_dijkstra(int **G, int thread_count, double *eta, int *distance, int *pred)
{
 
    int count,mindistance,nextnode,i,j;
    int my_mindist,my_nextnode;

    /* Allocating the memory fot start_time_array and end_time_array */
    double *start_time_array=(double *)malloc(thread_count*sizeof(double));
    double *end_time_array=(double *)malloc(thread_count*sizeof(double));

    /* Allocated the memory for few of the local variables */
    int **cost = (int **)malloc(rows*sizeof(int *));
    int *visited = (int *)malloc(rows*sizeof(int));

    #pragma omp parallel num_threads(thread_count) default(none) private(count, i, j, my_mindist, my_nextnode) shared(G, src_node, rows, mindistance, visited, distance, nextnode, cost, pred,start_time_array, end_time_array, eta,thread_count)
    { 
      int my_rank = omp_get_thread_num();

      /* intializing the execution time array */
      #pragma omp for schedule(guided,16)
      for(i=0;i<thread_count;i++)
      {
        start_time_array[i]=0.0;
        end_time_array[i]=0.0;
      } 

      start_time_array[my_rank] = omp_get_wtime(); /* Record the start time of all the threads. */

      /* Initialized all the allocated variable */
      #pragma omp for schedule(guided,16)
      for(i=0;i<rows;i++)
      {
        cost[i] = (int *)malloc(rows*sizeof(int));  //Allocates memory for column in each row
        for(j=0;j<rows;j++)
        {
          cost[i][j]=0;
        }
      }

      /* Intitalize the allocated array */
      #pragma omp for schedule(guided,16)
      for(i=0;i<rows;i++)
      { 
        distance[i]=0;
        pred[i]=0;
      }
    
      /* pred[] stores the predecessor of each node
      count gives the number of nodes seen so far
      create the cost matrix */
      #pragma omp for collapse(2) schedule(guided,16)
      for(i=0;i<rows;i++)
          for(j=0;j<rows;j++)
              if(G[i][j]==0)
                  cost[i][j]=INFINITY;
              else
                  cost[i][j]=G[i][j];   /* The matrix value will become the weight of the edge */ 
    
      /* initialize pred[],distance[] and visited[] with the required value */
      #pragma omp for schedule(guided,16)
      //#pragma omp simd for schedule(guided,16)
      for(i=0;i<rows;i++)
      {
        distance[i]=cost[src_node][i];
        pred[i]=src_node;
        visited[i]=0;
      }
      /* printf("After intialization\n");
         printVector(distance,"distance"); 
         printVector(pred,"pred"); */ 
      #pragma omp single
      { 
        distance[src_node]=0;
        visited[src_node]=1;
      }
      //printVector(visited,"visited"); 
      
      /* Picked the pragma from this location */
      for(count=1;count<rows-1;count++)
      {
        #pragma omp single
        {
          mindistance=INFINITY;
         // printf("count=%d\n",count);
         // printVector(distance,"distance"); 
        }

        my_mindist = INFINITY;
        my_nextnode = -1;

        /* nextnode gives the node at minimum distance */
        //#pragma omp for schedule (static,2)
        //#pragma omp for schedule(dynamic)
        #pragma omp for schedule(guided,16)
        for(i=0;i<rows;i++)
        {  // printf("distance[%d]=%d, mindistance=%d, visited[%d]=%d, if<%d < %d>\n",i,distance[i],mindistance,i,visited[i],distance[i],mindistance);
            if(distance[i]<my_mindist&&!visited[i])
            {
                my_mindist=distance[i];
                my_nextnode=i;
            }
        }   
        
        /* We need this region to be parallel because all the thread will have its own minimum and we want minimum of all */ 
        #pragma omp critical
        {
          if(my_mindist < mindistance)
          {
            mindistance = my_mindist;
            nextnode = my_nextnode;
          }
        }
        /* We need to make sure that all the threads have updated the critical section and then only we need to move further */
        #pragma omp barrier

        /* This needs to be updated by only single thread as this is shared among all the threads */
        #pragma omp single
        {
          visited[nextnode]=1;
         // printf("thread: %d, nextnode=%d,\n",my_rank,nextnode);
         // printf("thread: %d, mindistance=%d\n",my_rank,mindistance);
         // printVector(visited,"visited"); 
        }
        
        /* Again we don't want any thread to move further without the updation of the visited node */
        //#pragma omp barrier

        /* check if a better path exists through nextnode */            
        //#pragma omp for schedule(static,2)
        //#pragma omp for schedule(dynamic)
        #pragma omp for schedule(guided,16)
        for(i=0;i<rows;i++)
        {   
            //printf("visited[%d]=%d, if<%d>\n",i,visited[i],!visited[i]);
            if(!visited[i])
            {
               // printf("mindistance=%d, cost[%d][%d]=%d, distance[%d]=%d,\nif<%d < %d>\n",mindistance,nextnode,i,cost[nextnode][i],i,distance[i],mindistance+cost[nextnode][i],distance[i]);
                if(mindistance+cost[nextnode][i]<distance[i])
                {
                    distance[i]=mindistance+cost[nextnode][i];
                    pred[i]=nextnode;
                   // printf("pred[%d]=%d, nextnode=%d\n",i,pred[i],nextnode);
                }
            }
        }
      
        end_time_array[my_rank] = omp_get_wtime(); /* Record the end time of all the threads. */
        /* We don't want any thread to move to next step before each thread find its better path and update the distance vector */
        //#pragma omp barrier
      }
      /* Once all the threads completed their work it is now time to calculate their execution time */
      #pragma omp single
      {
        for(i=0;i<thread_count;i++)
        {
          eta[i]=end_time_array[i] - start_time_array[i];
        }
      } 
    } 

    /* Cleanup */
    free (cost);
    free (visited);
    free (start_time_array);
    free (end_time_array);
}
