#include "network.h"
#include <cmath>

#define N 303
#define G 320.
#define V 105.
#define L 4
#define INF 67108864.
#define R 6356.752
#define SQ(x) ((x)*(x))
#define d2r(x) (x/180.*3.14159265358979323846264338327950288419716939937510)

#define LOOP(x,y) for(int x=0;x<y;++x)
#define LOOP1 LOOP(i,N)
#define LOOP2 LOOP1 LOOP(j,N)
#define LOOP3 LOOP(k,N) LOOP2
#define RELAX(x,y,u,v) if(y<x){x=y;u=v;}

using namespace std;

string sn, tn;
int s,t;
double a[L][N][N], d[N][N], g[N];
int A[L][N][N];

bool input(int argc, char** argv) {
    if (argc != 3) return false;
    sn = argv[1];
    tn = argv[2];
    
    s = t = 0;
    while (s < network.size() && network[s].name != sn) ++s;
    while (t < network.size() && network[t].name != tn) ++t;
    if (s >= network.size() || t >= network.size()) return false;
    
    return true;
}

string output(int l, int i, int j) {
    if (i == j) return "";
    if (l == 0) {
        if (A[l][i][j] < 0)
            return " " + to_string(d[i][j]*g[i]/V) + " " + network[j].name;
        return output(0, i, A[l][i][j]) + output(0, A[l][i][j], j);
    }
    if (l == 1)
        return " " + network[A[l][i][j]].name + " "
        + to_string((d[i][A[l][i][j]]+d[A[l][i][j]][j]-G) * g[A[l][i][j]] / V)
        + " " + network[j].name;
    if (l == 2)
        return output(1, i, A[l][i][j]) + output(0, A[l][i][j], j);
    if (A[l][i][j] == -1)
        return output(2, i, j) + " " + to_string(g[j]/V);
    if (A[l][i][j] == -2)
        return " " + network[j].name + " " + to_string(d[i][j] * g[j] / V);
    return output(3, i, A[l][i][j]) + output(3, A[l][i][j], j);
}

int main(int argc, char** argv)
{
    if (!input(argc, argv)) {
        cout << "Error: requires initial and final supercharger names" << endl;
        return -1;
    }
    
    LOOP1 g[i] = V / network[i].rate;
    LOOP2 d[i][j]
        = R*2.*asin(sqrt(SQ(sin(d2r((network[i].lat-network[j].lat)/2.)))
        + cos(d2r(network[i].lat))*cos(d2r(network[j].lat))
        * SQ(sin(d2r((network[i].lon-network[j].lon)/2.)))));
    LOOP2 if (d[i][j] > G) d[i][j] = INF;
    if (d[s][t] <= G) {
        cout << sn + " " + tn << endl;
        return 0;
    }
        
    LOOP(l,L) LOOP2 a[l][i][j] = INF;
    LOOP(l,L) LOOP2 A[l][i][j] = -1;
    
    LOOP2 a[0][i][j] = d[i][j] * (g[i]+1.);
    LOOP3 RELAX(a[0][i][j], a[0][i][k] + a[0][k][j], A[0][i][j], k);
    
    LOOP3 if (d[i][k] + d[k][j] >= G)
        RELAX(a[1][i][j], G + (d[i][k]+d[k][j]-G)*(g[k]+1.), A[1][i][j], k);
    
    LOOP3 RELAX(a[2][i][j], a[1][i][k] + a[0][k][j], A[2][i][j], k);
    
    LOOP2 a[3][i][j] = a[2][i][j] + g[j] * G;
    LOOP2 RELAX(a[3][i][j], d[i][j] * (1. + g[j]), A[3][i][j], -2);
    LOOP3 RELAX(a[3][i][j], a[3][i][k] + a[3][k][j], A[3][i][j], k);
    
    int m = INF, M = -1;
    LOOP1 RELAX(m, a[3][s][i] + a[2][i][t], M, i);
    if (M < 0)
        cout << "No path exists." << endl;
    else
        cout << sn + output(3, s, M) + output(2, M, t) << endl;
    
    return 0;
}
