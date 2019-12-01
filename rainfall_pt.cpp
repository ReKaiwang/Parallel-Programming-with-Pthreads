#include <iostream>
#include <fstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>
#include <pthread.h>
#include <algorithm>
using namespace std;

int num_threads;
int M;
double A;
int N;

vector<vector<vector<double >>> dir;
vector<vector<double>> update;
vector<vector<double>> cur;
vector<vector<double>> absorb;
double total_absorb = 0;
int timesteps = 0;
vector<vector<int>> D = {{-1,0}, {1, 0}, {0, -1}, {0, 1}};
//initialize the lock for total_absorb
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
//define the lock matrix for update matrix
pthread_mutex_t* locks;
//define the barrier for update matrix clearance
pthread_barrier_t mybarrier;
double calc_time(struct timespec start, struct timespec end) {
    double start_sec = (double)start.tv_sec+ (double)start.tv_nsec / 1000000000.0 ;
    double end_sec = (double)end.tv_sec+ (double)end.tv_nsec / 1000000000.0 ;

    if (end_sec < start_sec) {
        return 0;
    } else {
        return end_sec - start_sec;
    }
}

//parallel
void clear_update(int startrow, int endrow){
    for(int i = startrow; i <= endrow; ++i){
        for(int j = 0; j < N; ++j){
            update[i][j] = 0;
        }
    }
}
void update_pal(int i, int j){
    cur[i][j] += update[i][j];
}

void raining_pal(int i, int j){
    cur[i][j]  += 1.0;
}

void absorb_pal(int i, int j, double& cur_sum){
    double value = min(A, cur[i][j]);
    if(!value){return;}
//    pthread_mutex_lock(&lock);
    cur_sum += value;
//    pthread_mutex_unlock(&lock);
    absorb[i][j] += value;
    cur[i][j] -= value;
}
void trickle_pal_lock_free(int i, int j){
    if(!dir[i][j].empty()) {
        double trickle_num = min(double(1.0), cur[i][j]);
        for (int k : dir[i][j]) {
            int row = k / N;
            int col = k % N;
            update[row][col] += trickle_num / dir[i][j].size();
        }
        cur[i][j] -= trickle_num;
    }
}
void trickle_pal(int i, int j){
    if(!dir[i][j].empty()) {
        double trickle_num = min(double(1.0), cur[i][j]);
        for (int k : dir[i][j]) {
            int row = k / N;
            int col = k % N;
            pthread_mutex_lock(&locks[k]);
            update[row][col] += trickle_num / dir[i][j].size();
            pthread_mutex_unlock(&locks[k]);
        }
        cur[i][j] -= trickle_num;
    }
}
void* rainaction(void* argv){
    double cur_sum;
    int num  = *((int *)argv);
    int startrow = num * N / num_threads;
    int endrow = (num+1)* N / num_threads - 1;
    for(int r = 0; r < M; ++r) {
        cur_sum = 0.0;
        //clear the update matrix
        clear_update(startrow, endrow);
        pthread_barrier_wait(&mybarrier);
        //only lock startrow and endrow
        for (int i = startrow; i <= endrow; ++i) {
            if(i == startrow || i == endrow){
                for (int j = 0; j < N; ++j) {
                    raining_pal(i, j);
                    absorb_pal(i, j, cur_sum);
                    trickle_pal(i, j);
                }
            }
            else {
                for (int j = 0; j < N; ++j) {
                    raining_pal(i, j);
                    absorb_pal(i, j, cur_sum);
                    trickle_pal_lock_free(i, j);
                }
            }
        }
        pthread_mutex_lock(&lock);
        total_absorb += cur_sum;
        pthread_mutex_unlock(&lock);
        pthread_barrier_wait(&mybarrier);
        for (int i = startrow; i <= endrow; ++i) {
            for (int j = 0; j < N; ++j) {
                update_pal(i, j);
            }
        }
    }
}
void* unrainaction(void* argv){
    double cur_sum;
    double total_rain = M * N * N;
    int num  = *((int *)argv);
    int startrow = num * N / num_threads;
    int endrow = (num+1)* N / num_threads - 1;
    while(1) {
        cur_sum = 0.0;
        //clear the update matrix
        clear_update(startrow, endrow);
        pthread_barrier_wait(&mybarrier);
        //only lock startrow and endrow
        for (int i = startrow; i <= endrow; ++i) {
            if(i == startrow || i == endrow){
                for (int j = 0; j < N; ++j) {
                    absorb_pal(i, j, cur_sum);
                    trickle_pal(i, j);
                }
            }
            else {
                for (int j = 0; j < N; ++j) {
                    absorb_pal(i, j, cur_sum);
                    trickle_pal_lock_free(i, j);
                }
            }
        }
        pthread_mutex_lock(&lock);
        total_absorb += cur_sum;
        pthread_mutex_unlock(&lock);
        pthread_barrier_wait(&mybarrier);
        if(total_rain - total_absorb < 1e-5){break;}
        for (int i = startrow; i <= endrow; ++i) {
            for (int j = 0; j < N; ++j) {
                update_pal(i, j);
            }
        }
    }
}
void* unrainaction_plus(void* argv){
    double cur_sum;
    double total_rain = M * N * N;
    int num  = *((int *)argv);
    int startrow = num * N / num_threads;
    int endrow = (num+1)* N / num_threads - 1;
    while(1) {
        cur_sum = 0.0;
        //clear the update matrix
        clear_update(startrow, endrow);
        pthread_barrier_wait(&mybarrier);
        //only lock startrow and endrow
        for (int i = startrow; i <= endrow; ++i) {
            if(i == startrow || i == endrow){
                for (int j = 0; j < N; ++j) {
                    absorb_pal(i, j, cur_sum);
                    trickle_pal(i, j);
                }
            }
            else {
                for (int j = 0; j < N; ++j) {
                    absorb_pal(i, j, cur_sum);
                    trickle_pal_lock_free(i, j);
                }
            }
        }
        pthread_mutex_lock(&lock);
        total_absorb += cur_sum;
        pthread_mutex_unlock(&lock);
        pthread_barrier_wait(&mybarrier);
        ++timesteps;
        if(total_rain - total_absorb < 1e-5){break;}
        for (int i = startrow; i <= endrow; ++i) {
            for (int j = 0; j < N; ++j) {
                update_pal(i, j);
            }
        }
    }
}
int main(const int argc, char** argv){
    if(argc != 6){
        cout << "invalid argument" << endl;
        return 1;
    }
    num_threads = stoi(argv[1]);
    M = stoi(argv[2]);
    A = stod(argv[3]);
    N = stoi(argv[4]);
    string file_name = argv[5];
    num_threads = min(num_threads, N);
    //read from file
    //cout << M << " " << A << " " <<N <<  file_name << endl;
    vector<vector<int>> land(N, vector<int>(N));
    string line;
    ifstream myfile (file_name);
    if(myfile.is_open()){
        int row = 0;
        while(getline(myfile, line)){
            int col = 0;
            int prev = 0;
            while(prev < line.size() && line[prev] == ' '){
                ++prev;
            }
            int cur = prev;
            while(col < N){
                while(cur < line.size() && line[cur] != ' '){
                    ++cur;
                }
                if(cur == line.size()){
                    land[row][col] = stoi(line.substr(prev));
                }else {
                    land[row][col] = stoi(line.substr(prev, cur - prev));
                    while(cur < line.size() && line[cur] == ' '){
                        ++cur;
                    }
                    prev = cur;
                }
                ++col;
            }
            ++row;
        }

    }else{
        cout << "Unable to open the file" << endl;
        return 1;
    }
    //initialize vectors
    dir = vector<vector<vector<double >>>(N, vector<vector<double>>(N));
    update = vector<vector<double>>(N, vector<double>(N, 0));
    cur = vector<vector<double>>(N, vector<double>(N, 0));
    absorb = vector<vector<double>>(N, vector<double>(N, 0));
    //initialize rwlock matrix
    locks = (pthread_mutex_t*) malloc(N * N * sizeof(pthread_mutex_t));
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            locks[i*N+j] = PTHREAD_MUTEX_INITIALIZER;
        }
    }
    //initialize the barrier
    pthread_barrier_init(&mybarrier, NULL, num_threads);
    // calculate direction
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            int cur_lowest = land[i][j];
            //find lowest
            for(int k = 0; k < 4; ++k){
                int x = i + D[k][0];
                int y = j + D[k][1];
                if(x >= 0 && x < N && y >=0 && y < N && land[x][y] < cur_lowest){
                    cur_lowest = land[x][y];
                }
            }
            //cout << i << " " << j << ":";
            if(cur_lowest < land[i][j]){
                for(int k = 0; k < 4; ++k) {
                    int x = i + D[k][0];
                    int y = j + D[k][1];
                    if (x >= 0 && x < N && y >= 0 && y < N && land[x][y] == cur_lowest) {

                        //cout << "("<<x << " " << y << ") ";
                        dir[i][j].push_back(x * N + y);
                    }
                }
            }
            //cout << endl;

        }
    }
    vector<int> ID(num_threads, 0);
    for(int i = 0; i < num_threads; i++){
        ID[i] = i;
    }
    struct timespec start_time, end_time;
    double elapsed_s;
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    pthread_t * threads = (pthread_t *) malloc(num_threads * sizeof(pthread_t));
//    cout << "preparation done" << endl;
    //begin rain
    //create threads
    for(int i = 0; i < num_threads; i++){
        pthread_create(&threads[i], NULL, rainaction,(void *)(&ID[i]));
    }
    //wait for subjob
    for(int i = 0; i < num_threads; i++){
        pthread_join(threads[i], NULL);
    }
    timesteps += M;
    cout << "raining is over" << endl;

    //only absort
    //create threads
    pthread_create(&threads[0], NULL, unrainaction_plus,(void *)(&ID[0]));
    for(int i = 1; i < num_threads; i++){
        pthread_create(&threads[i], NULL, unrainaction,(void *)(&ID[i]));
    }
    //wait for subjob
    for(int i = 0; i < num_threads; i++){
        pthread_join(threads[i], NULL);
    }
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    //output result
    elapsed_s = calc_time(start_time, end_time);
    cout << "Rainfall simulation took " << timesteps <<" time steps to complete." << endl;
    cout << "Runtime = "<<elapsed_s << " seconds" << endl;
    cout<< "The following grid shows the number of raindrops absorbed at each point:" << endl;

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; ++j){
            cout << absorb[i][j];
            if(j != N - 1){
                cout << " ";
            }
        }
        cout << endl;
    }
    free(threads);
    free(locks);
    return 0;

}
