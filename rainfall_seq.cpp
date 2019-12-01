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
#include <exception>
using namespace std;

int M;
double A;
int N;

vector<vector<int>> D = {{-1,0}, {1, 0}, {0, -1}, {0, 1}};
double calc_time(struct timespec start, struct timespec end) {
    double start_sec = (double)start.tv_sec+ (double)start.tv_nsec / 1000000000.0 ;
    double end_sec = (double)end.tv_sec+ (double)end.tv_nsec / 1000000000.0 ;

    if (end_sec < start_sec) {
        return 0;
    } else {
        return end_sec - start_sec;
    }
}


//sequence
void update_seq(vector<vector<double>>& cur, vector<vector<double>>& update){
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            cur[i][j] += update[i][j];
        }
    }
}

void raining_seq(vector<vector<double >>& cur){
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            cur[i][j]  += 1.0;
        }
    }
}

void absorb_seq(vector<vector<double>>& cur, vector<vector<double>>& absorb, double& total_absorb){
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            total_absorb += min(A, cur[i][j]);
            absorb[i][j] += min(A, cur[i][j]);
            cur[i][j] -= min(A, cur[i][j]);
        }
    }
}

void trickle_seq(vector<vector<double>>& cur, vector<vector<double>>& update, vector<vector<vector<double>>>& dir){

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            update[i][j] = 0;
        }
    }

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            if(!dir[i][j].empty()){
                double trickle_num = min(double(1.0),cur[i][j]);
                update[i][j] -= trickle_num;
                for(int k : dir[i][j]){
                    int row = k / N;
                    int col = k % N;

                    update[row][col] += trickle_num / dir[i][j].size();
                }
                //cur[i][j] = 0;
            }
        }
    }
}

int main(const int argc, char** argv){
    if(argc != 6){
        cout << "invalid argument" << endl;
        return 1;
    }

    M = stoi(argv[2]);
    A = stod(argv[3]);
    N = stoi(argv[4]);
    string file_name = argv[5];
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

    vector<vector<vector<double>>> dir(N, vector<vector<double>>(N));
    vector<vector<double>> update(N, vector<double>(N, 0));
    vector<vector<double>> cur(N, vector<double>(N, 0));
    vector<vector<double>> absorb(N, vector<double>(N, 0));
    int timesteps = 0;
    double total_absorb = 0;
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
            if(cur_lowest < land[i][j]){
                for(int k = 0; k < 4; ++k) {
                    int x = i + D[k][0];
                    int y = j + D[k][1];
                    if (x >= 0 && x < N && y >= 0 && y < N && land[x][y] == cur_lowest) {

                        dir[i][j].push_back(x * N + y);
                    }
                }
            }

        }
    }

    struct timespec start_time, end_time;
    double elapsed_s;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    //begin rain
    for(int r = 0; r < M; ++r){
        raining_seq(cur);
        absorb_seq(cur, absorb, total_absorb);
        trickle_seq(cur, update, dir);
        update_seq(cur, update);
        ++timesteps;
    }

    //only absort

    while(M * N * N  - total_absorb > 1e-5){
        absorb_seq(cur, absorb, total_absorb);
        trickle_seq(cur, update, dir);
        update_seq(cur, update);
        //cout << timesteps <<endl;
        ++timesteps;
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
    return 0;

}