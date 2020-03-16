
/* 
AUTHORS: Mikael Karlsson , Felix Magnell, Alexander Engström
 */

/**
Compile:  g++ TSP.cpp -o TSP
Run: ./TSP < input.bin > test.txt
*/
#include <iostream>
#include <vector>
#include <array>
#include <iterator>
#include <string>
#include <deque>
#include <numeric>
#include <math.h>       /* pow */
//#include<bits/stdc++.h>
using namespace std;
using std::cin;
using std::cout;
using std::cerr;
using std::array;
int N;
 
 
typedef pair<double, double> Pair;
vector<Pair> graph;
array<int, 1000> tour;
array<int, 1000> best_tour_ever;
array<int, 1000> new_route;
array<int, 1000> used;
int distanceMatrix[1000][1000];
 
//vector<int> tour(1000, -1);
//vector<int> new_route(1000, -1);
//vector<int> used(1000, 0);
vector<int> edges(1000, -1);
 
 
/* void printResult(){
 
  int currentNode = 0;
  cout << 0 << endl;
  currentNode = edges.at(0);
 
  for(int i = 0; i < N-1; i++){
      cout << currentNode << endl;
      currentNode = edges.at(currentNode);
  }
} */
 
 
void readInput() {
 
  // Läs antal
  cin >> N;
 
  //Läs in x och y (punkterna)
  for (int i = 0; i < N; ++i) {
    double x, y;
    cin >> x >> y;
    graph.push_back(make_pair(x, y));
  }
}
 
 
/* void writeGraph() {
 
  for (int i = 0; i < N; ++i) {
    cout << " " << graph[i].first << " " << graph[i].second << "\n";
  }
  // Var noggrann med att flusha utdata när grafen skrivits ut!
  cout.flush();
 
} */
 
 
 
 
double dist(double p1x, double p1y, double p2x, double p2y){
  return sqrt(pow((p1x - p2x),2) + pow((p1y - p2y),2));
}
 
double distance(int p1, int p2){
  return dist(graph.at(p1).first,graph.at(p1).second, graph.at(p2).first, graph.at(p2).second);
}


void buildDistanceMatrix(){
  for(int i = 0; i < N; i++ ){
    for(int j = 0; j < N; j++){
      distanceMatrix[i][j] = distance(i, j);
    }
  }
}


/* void writeDistanceMatrix(){
  for(int i = 0; i < N; i++ ){
    for(int j = 0; j < N; j++){
      cout << distanceMatrix[i][j]  << " ";
    }
    cout << "\n";
  }
  cout << "\n";
} */
 
void writeTour(){
  for (int i = 0; i < N; ++i) {
    cout << best_tour_ever[i] << "\n";
  }
  cout.flush();
}
 
 
void greedyTour(){
  tour[0] = 0;
  used[0] = 1;
  double best;
  for(int i = 1; i < N; ++i){
      best = -1;
      for(int j = 0; j < N; ++j){
        if(used[j] == 0 && (best == -1)){
          best = j;
        }else if(used[j] == 0  && (distance(tour[i-1], j) < distance(tour[i-1], best))){
            best = j;
        }
       
      }
        tour[i] = best;
        used[best] = 1;
  }
}

 
void twoOptSwap(int i, int k) {

        double d = 0;
        
       for(int j=0; j < i; j++){
         new_route[j] = tour[j];

       }
             
       int k2 = k;
       int j = i;

       while(k2 >= i){
         new_route[j] =  tour[k2];
         k2--;
         j++;

       }
       
       for(int j = k+1; j < N; j++){
          new_route[j] =  tour[j];
       }

   }
 
 
/* double calcTotDistance(std::array<int, 1000> &tour){
  double totDist=0;
  for(int i = 0; i < N-2; i++ ){
    totDist += distance(tour[i], tour[i+1]);
  }
 
  totDist += distance(tour[N-1], tour[0]);
  return totDist;
} */


double calcTotDistanceFast(std::array<int, 1000> &tour){
  double totDist=0;
  for(int i = 0; i < N-2; i++ ){
    //totDist += distance(tour[i], tour[i+1]);
    totDist += distanceMatrix[tour[i]][tour[i+1]];
  }
  //totDist += distance(tour[N-1], tour[0]);
  totDist += distanceMatrix[tour[N-1]][tour[0]];
  return totDist;
}


void random_swap(){
    int r1;
    int j;
    int k;

    for(int i = 0; i < 25; i++){
        r1 = rand() % (N - 4);
        j = tour[r1];
        k = tour[r1 + 1];
        tour[r1] = k;
        tour[r1 + 1] = j;
    }
}
 
int main(void) {
    int cc = 0;
    
  clock_t timer;
  timer = clock();
  double duration;

  std::ios::sync_with_stdio(false);
 
 
  bool first_iteration = true;
  cin.tie(0);
  duration = ( clock() - timer ) / (double) CLOCKS_PER_SEC;
  readInput();
  greedyTour();
  buildDistanceMatrix();
  best_tour_ever = tour;

  double dd = calcTotDistanceFast(tour);
  double old_dd = 999999999999999999;
 
  bool swapped = false;

  double new_distance = 0;
  double current_distance = 0;
  bool shouldSwap = false;
  int best_i = 0;
  int best_k = 1;
  

   while(duration < 1.995){
    


    cc += 1;
    first_iteration = false;


     for (int i = 0; i < N ; i++) {
       duration = (clock() - timer ) / (double) CLOCKS_PER_SEC;
         
          if(duration > 1.995) break;
           
           for (int k = i + 1; k < N-1; k++) {
             
                duration = (clock() - timer ) / (double) CLOCKS_PER_SEC;
                if(duration > 1.995) break;
                
                //om man gatt runt ett varv ?
                if(i == N-1){
                k = 0;
                }
                
                //if(i == k) continue;
                int i2 = i - 1;

                //om i ar 0 satt i2 till sista noden i listan  
                if(i2 == -1){
                i2 = N-1;
                }

                current_distance = distanceMatrix[tour[i]][tour[i2]] + distanceMatrix[tour[k]][tour[k+1]];
                new_distance = distanceMatrix[tour[i2]][tour[k]] + distanceMatrix[tour[i]][tour[k+1]];

                if (new_distance < current_distance) {
                    best_i = i;
                    best_k = k;
                    shouldSwap = true;   
                }
           


           if(shouldSwap == true){

               dd = dd - current_distance;
               dd = dd + new_distance;

              twoOptSwap(best_i, best_k);
              tour = new_route;
              shouldSwap = false;
              swapped = true;
           } 

     }  
               
           }   

            if(swapped == true && dd < old_dd ){
               best_tour_ever = tour;
               swapped = false;
               old_dd = dd;

               
           }else{
            
               //if(N > 9){
                    random_swap();
                    dd = calcTotDistanceFast(tour);
              // }
           
       }

 
    duration = (clock() - timer ) / (double) CLOCKS_PER_SEC;
 
   
   
  }

    //cout << "dd: " << dd << endl;
    writeTour();
 
  return 0;
}