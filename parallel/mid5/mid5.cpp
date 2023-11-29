#include <bits/stdc++.h>
#include "mpi.h"

#define WORLD MPI_COMM_WORLD

using namespace std;

using ll = long long int;
typedef pair<double, int> CostTo;

struct cmp {
	bool operator() (const CostTo& a, const CostTo& b) const {
		return a.first > b.first;
	}
};

struct Point{
    int x, y;
    double operator * (Point q) { return x * q.x + y * q.y; }
    Point operator + (Point b) { return {x + b.x, y + b.y}; }
    Point operator - (Point b) { return {x - b.x, y - b.y}; }
    double operator ^ (Point q) { return x * q.y - y * q.x; }
    bool operator < (const Point& b) const { return x == b.x ? y < b.y : x < b.x; }
	double dis(Point& b) {
		return floor(sqrt(pow(x - b.x, 2) + pow(y - b.y, 2)) * 10000.0) / 10000.0 + numeric_limits<double>::epsilon();
	}
};

void convexHull(vector<Point>& ps, vector<Point>& res) {
    sort(ps.begin(), ps.end());
    res.clear();

    int n = 0;
    for(int i = 0; i < ps.size(); ++i) {
        while(n >= 2 && ((res[n - 1] - res[n - 2]) ^ (ps[i] - res[n - 2])) < 0) {
            n--;
            res.pop_back();
        }
        res.push_back(ps[i]);
        n++;
    }

    for(int i = ps.size() - 2, n2 = n + 1; i >= 0; --i) {

        while(n >= n2 && ((res[n - 1] - res[n - 2]) ^ (ps[i] - res[n - 2])) < 0) {
            n--;
            res.pop_back();
        }
        res.push_back(ps[i]);
        n++;
    }
    res.pop_back();
}

bool vis[20] = {0};
priority_queue<CostTo, vector<CostTo>, cmp> pq;

double prim(vector<Point>& ps) {
	while(!pq.empty()) {
		pq.pop();
	}
    pq.push({0, 0});
	memset(vis, 0, 20 * sizeof(bool));

    double total = 0;
    int cnt = 0;
	int n = ps.size();

    while(!pq.empty()) {
        int curr = pq.top().second; 

        if(vis[curr]) { 
            pq.pop();
            continue;
        }

        total += pq.top().first; pq.pop();
        vis[curr] = true;
        cnt++;

        for(int i = 0; i < n; ++i) {
			Point& nei = ps[i];
            if(!vis[i]) {
                pq.push({nei.dis(ps[curr]), i});
            }
        }

        if(cnt == n) {
            break;
        }
    }

    if(cnt == n) {
        return total;
    } else {
        return INT_MAX;
    }
}

int main(void) {

	int nums, id;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(WORLD, &nums);
	MPI_Comm_rank(WORLD, &id);

	vector<Point> ps;
	int n;
	if(id == 0) {
		string s;
		cin >> s;

		FILE* f = fopen(s.c_str(), "r");
		fscanf(f, "%d", &n);
		ps.reserve(n);

		int x, y;
		for(int i = 0; i < n; ++i) {
			fscanf(f, "%d %d", &x, &y);
			ps.push_back({x, y});
		}
	}

	vector<Point> ch;	
	vector<Point> internal;
	
	if(id == 0) {
		convexHull(ps, ch);
		set<Point> hull(ch.begin(), ch.end());
		for(int i = 0; i < n; ++i) {
			if(hull.find(ps[i]) == hull.end()) {
				internal.push_back(ps[i]);
			}
		}
	}

	int intlen = internal.size(), chlen = ch.size();
	int arr[2] = {intlen, chlen};

	MPI_Bcast(arr, 2, MPI_INT, 0, WORLD);
	intlen = arr[0], chlen = arr[1];
	
	if(id != 0) {
		ch.resize(chlen);
		internal.resize(intlen);
	}

	MPI_Bcast(ch.data(), chlen * sizeof(Point), MPI_BYTE, 0, WORLD);
	MPI_Bcast(internal.data(), intlen * sizeof(Point), MPI_BYTE, 0, WORLD);

	// calc
	ll finalState = (1 << intlen) - 1;
	vector<Point> p(ch.begin(), ch.end());
	double minn = INT_MAX;

	for(ll i = id; i <= finalState; i += nums) {
		// add points to p
		for(int j = 0; j < intlen; ++j) {
			if(i & (1 << j)) {
				p.push_back(internal[j]);
			}
		}
		
		// do prim algor
		minn = min(minn, prim(p));

		for(int j = p.size(); j > chlen; --j) {
			p.pop_back();
		}
	}

	double ans;
	MPI_Reduce(&minn, &ans, 1, MPI_DOUBLE, MPI_MIN, 0, WORLD);

	if(id == 0) {
		ans *= 10000;
		ans = floor(ans);
		ans /= 10000.0;
		ans += numeric_limits<double>::epsilon();
		cout << fixed << setprecision(4) << ans;
	}

	MPI_Finalize();
	return 0;
}
