#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <cassert>

#define GROUP_SIZE 2000
#define GENERATIONS 2000
#define ABNORMAL 0.6
using namespace std;

pair<double, double> pos[30];
double dis[30][30];             //all x to y distance
string id2name[30];
int n;


class Randomizer{
  //this is a randomizer which generates int in range[0, 2^64 - 1]
  unsigned long long num;
public:
  const unsigned long long MAX;
  Randomizer() : MAX(-1) {
    num = time(0);
  }
  unsigned long long next_random() {
    return num = num * (unsigned long long)25214903917 + 11;
  }
} randomizer;

struct Generation {

  // This is the class to store one generation of the specie
  // Encoding method:
  //    a sequence of city ids to represent the TSP route
  // Mating method:
  //    for FA[x1, x2, ..., xp, ..., xn] and
  //        FB[y1, y2, ..., yp, ..., yn]
  //    get SA[x1, x2, ..., xp, yt1, yt2, ..., yt(n-p)] where (t1 < t2 < t3 ...)  and
  //        SB[yt1, yt2, ..., yt(n-p), xp + 1, ..., xn] where (t1 < t2 < t3 ...)
  // Variation method:
  //    for A[x1, x2, ..., xa, xa+1, ..., xb, ..., xn]
  //    get B[x1, x2, ..., xb, xb-1, ..., xa, ..., xn]
  // Generating method:
  //    get 2 parent under probability density 1/cost, where cost is the length of the route
  //    generate two son
  //    random variation
  // Terminate condition:
  //    Finish loop

  static int best[30];
  static double bestv;
  int status[GROUP_SIZE][30];
  double pro[GROUP_SIZE];
  int cnt;
  void clear() {
    cnt = 0;
  }
  double evaluate(int p) {
    double res = 0;
    for (int i = 0; i < n; i++) {
      res += dis[status[p][i]][status[p][(i + 1 == n) ? 0 : i + 1]];
    }
    if (res < bestv) {
      bestv = res;
      for (int i = 0; i < n; i++)
        best[i] = status[p][i];
    }
    return res;
  }
  void flushpro() {
    for (int i = 0; i < GROUP_SIZE; i++)
      pro[i] = evaluate(i);
    for (int i = 0; i < GROUP_SIZE; i++)
      pro[i] = 1.0 / pro[i] + ((i > 0) ? pro[i - 1] : 0);
    for (int i = 0; i < GROUP_SIZE; i++)
      pro[i] /= pro[GROUP_SIZE - 1];
  }
  int* pick() {
    double x = (double)randomizer.next_random() / randomizer.MAX;
    int p = lower_bound(pro, pro + GROUP_SIZE, x) - pro;
    if (p == GROUP_SIZE) p--;
    return status[p];
  }
  void add(int *p1, int *p2) {
    int p = randomizer.next_random() % (n - 1);
    bool h[30];
    memset(h, 0, sizeof(h));
    for (int i = 0; i < p; i++) {
      h[ status[cnt][i] = p1[i] ] = 1;
    }
    for (int i = p, x = 0; i < n; i++) {
      while (h[p2[x]]) x++;
      h[status[cnt][i] = p2[x]] = 1;
    }
    cnt++;

    memset(h, 0, sizeof(h));
    for (int i = p; i < n; i++) {
      h[ status[cnt][i] = p1[i] ] = 1;
    }
    for (int i = 0, x = 0; i < p; i++) {
      while (h[p2[x]]) x++;
      h[status[cnt][i] = p2[x]] = 1;
    }
    cnt++;
  }
  void abnormal(int p) {
    int x1 = randomizer.next_random() % (n - 1);
    int x2 = randomizer.next_random() % (n - 1);
    if (x1 > x2) swap(x1, x2);
    double pre = evaluate(p);
    reverse(status[p] + x1, status[p] + x2);
    double now = evaluate(p);
    if (now > pre && (double)randomizer.next_random() / randomizer.MAX > exp((pre - now) * 5))
      reverse(status[p] + x1, status[p] + x2);
  }
};

double Generation::bestv;
int Generation::best[30];
Generation *g1, *g2;


double sqr(double x) {return x * x; }

double calcdis(pair<double, double> &a, pair<double, double> &b) {
  return sqrt(sqr(a.first - b.first) + sqr(a.second - b.second));
}

void init() {

  //read position information
  for (int i = 0; i < n; i++) {
    string name;
    double x, y;
    cin >> name >> x >> y;
    id2name[i] = name;
    pos[i] = make_pair(x, y);
  }

  //calc all distance
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      dis[i][j] = calcdis(pos[i], pos[j]);

  Generation::bestv = 1e8;
}

void Gen() {
  g1 = new Generation;
  g2 = new Generation;
  for (int i = 0; i < GROUP_SIZE; i++) {
    for (int j = 0; j < n; j++)
      g1->status[i][j] = j;
    random_shuffle(g1->status[i], g1->status[i] + n);
  }
  g1->cnt = GROUP_SIZE;
  g1->flushpro();

  for (int i = 0; i < GENERATIONS; i++) {
    g2->clear();
    for (int j = 0; j < GROUP_SIZE / 2; j++) {
      int *p1 = g1->pick();
      int *p2 = g1->pick();
      while (p1 == p2) {
        p2 = g1->pick();
      }
      g2->add(p1, p2);
    }
    g2->flushpro();
    for (int j = 0; j < GROUP_SIZE; j++)
      if ((double)randomizer.next_random() / randomizer.MAX < ABNORMAL)
        g2->abnormal(j);
    g2->flushpro();
    swap(g1, g2);
  }
}

int main(int argc, char **argv) {
  if (argc == 3) {
    freopen(argv[1], "r", stdin);
    freopen(argv[2], "w", stdout);
  }
  scanf("%d", &n);

  init();

  for (int i = 0; i < 10; i++) Gen();

  for (int i = 0; i < n; i++)
    printf("%s ", id2name[Generation::best[i]].c_str());
  printf("%.8lf\n", Generation::bestv);
}
