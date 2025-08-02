#include<bits/stdc++.h>
using namespace std;
set<int> part[1200000];
int main(){
    ifstream file("./dblp.txt.1000.pedges", ios::in);
    int u, v, w;
    int mx = -1;
    int cnt = 0;
    while(file >> u){
        file >> v;
        file >> w;
        mx = max(mx, max(u, v));
        ++cnt;
        part[u].insert(w);
        part[v].insert(w);
    }
    file.close();
    std::vector<int> iso[1205];
    for(int i = 0; i < mx; ++i){
        if(part[i].size() > 1) iso[*part[i].begin()].push_back(i);
        part[i].clear();
    }
    for(int i = 0; i < 1000; ++i){
        ofstream outp("./node_replica/" + to_string(i) + ".txt", ios::out);
        for(auto j = iso[i].begin(); j != iso[i].end(); ++j){
            outp << *j << ' ';
        }
        outp << '\n';
        outp.close();
    }
    return 0;
}