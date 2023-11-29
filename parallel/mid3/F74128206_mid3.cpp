#include <bits/stdc++.h>
#include <mpi.h>

#define WORLD MPI_COMM_WORLD
using ll = long long int;

using namespace std;
void mergeSort(vector<int>& a, vector<int>& b, vector<int>& res) {
	res.clear();
	res.reserve(a.size() + b.size());

	int i = 0, j = 0;
	while(i < a.size() || j < b.size()) {
		if(i < a.size() && j < b.size()) {
			if(a[i] > b[j]) {
				res.push_back(b[j++]);
			} else {
				res.push_back(a[i++]);
			}
		} else if(i < a.size()) {
			res.push_back(a[i++]);
		} else if(j < b.size()) {
			res.push_back(b[j++]);
		}
	}
}
int main(void) {

	int nums, id;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(WORLD, &nums);
	MPI_Comm_rank(WORLD, &id);
	
	int n;
	vector<int> arr;
	if(id == 0) {
		string fn;
		cin >> fn;

		FILE* f = fopen(fn.c_str(), "r");
		fscanf(f, "%d", &n);
		arr.resize(n);

		for(int i = 0; i < n; ++i) {
			fscanf(f, "%d", &arr[i]);
		}
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, WORLD);

	int epp = n / nums;
	int remain = n % nums;
	int from = id * epp + min(id, remain);
	int size = epp + (id < remain ? 1: 0);

	int cnt[nums] = {0};
	int disp[nums] = {0};
	for(int i = 0; i < nums; ++i) {
		cnt[i] = epp + (i < remain ? 1 : 0);
		if(i > 0) {
			disp[i] = disp[i - 1] + cnt[i - 1];
		}
	}

	vector<int> a(size);
	MPI_Scatterv(arr.data(), cnt, disp, MPI_INT, a.data(), size, MPI_INT, 0, WORLD);

	sort(a.begin(), a.end());

	vector<int> b(n);
	vector<int> res(n);

	int step = 1;
	while(step < nums) {
		for(int i = 0; i < nums; i+=2*step) {
			if(i + step < nums) {
//				cout << i << " " << i + step << "\n";
				if(id == i + step) {
					MPI_Send(a.data(), a.size(), MPI_INT, i, i, WORLD);
				} else if(id == i) {
					MPI_Status s;
					b.resize(n);
					MPI_Recv(b.data(), n, MPI_INT, i + step, i, WORLD, &s);
					int len;
					MPI_Get_count(&s, MPI_INT, &len);

					b.resize(len);

					mergeSort(a, b, res);
					a = res;
				}
			}
		}
		step *= 2;
	}

	if(id == 0) {

		for(int i = 0; i < a.size(); ++i) {
			cout << a[i] << " ";
		}
	}

	MPI_Finalize();
	return 0;
}
