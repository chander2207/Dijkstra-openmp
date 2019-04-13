# Dijkstra-openmp
Dijkstra's Algorithm using openMP for parallelization

This program has 2 functions in it.
     -- serial_dijkstra: It is a normal execution of the code
     -- omp_dijkstra: This function uses the openMP commands to run Dijkstra into multiple threads to achieve better efficiency.

Usage:
      ./binary_file <thread_count> [filename]
     
      where filename could be:
      -- facebook_combined.txt : This dataset consists of 'circles' (or 'friend’s lists') from Facebook.
                                 Facebook data was collected from survey participants using this Facebook app.
                          Nodes: 4039
                          Edges: 88234
      -- soc-sign-Slashdot081106.txt : Slashdot is a technology-related news website known for its specific user community. 
                                       In 2002 Slashdot introduced the Slashdot Zoo feature which allows users to tag each 
                                       other as friends or foes. The network contains friend/foe links between the users of Slashdot.
                               Nodes : 82168
                               Edges : 948464 
      -- soc-sign-bitcoinalpha.csv : This is who-trusts-whom network of people who trade using Bitcoin on a platform called Bitcoin OTC.
                       Data format : SOURCE, TARGET, RATING, TIME
                             Nodes : 5,881
                             Edges : 35,592

The first line of the each file consist of the source node which can be changed in the file. This has been done to make sure 
that all the input should be provided to the program through file.

OMP clauses selected for parallelization
These lines have been taken from the code:
#pragma omp parallel num_threads(thread_count) default(none) private(count, i, j, my_mindist, my_nextnode) shared(G, src_node, rows, mindistance, visited, distance, nextnode, cost, pred,start_time_array, end_time_array, eta,thread_count)
 This part is used for the parallelization of the code. Rest all the code came under this section.

#pragma omp for schedule(guided,16)
  The schedule clause has been used to determine how loop iterations are assigned to threads. Here guided has been selected with chunk size 16, the iteration will be assigned to threads in the team of chunk, until no chunks remain to be assigned. This has been used majorly for all for loops.

#pragma omp for collapse(2) schedule(guided,16)
  The collapse was used to for two consecutive nested loop where it count both of them in single space and then divide the iterations. It performed faster. In code, the initialization of the cost matrix which contain the weight of all the edges. This was a big chunk was initialized quickly.

#pragma omp single
  This single statement was used to update few values that were supposed to be executed only once. Hence using this helped execute the code section by only single thread. In code, the “visited” node vector was marked by the single thread.

#pragma omp critical
  This clause was used to update the “mindistance” and “nextnode” variable once all the thread calculated their own values.

#pragma omp barrier
  This was used to restrict any thread moving further without updating the value of mindistance and the nextnode. This is because the further calculation was based on these values.

    
