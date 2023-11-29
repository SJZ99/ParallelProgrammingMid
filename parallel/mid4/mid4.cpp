#include <bits/stdc++.h>
#include <mpi.h>

#define WORLD MPI_COMM_WORLD

using namespace std;
using ll = long long int;

int nums, id;

inline void calc(int* a, int* b, int* k, int& n, int& m, int& d1, int& d2) {
	int hd1 = (d1 - 1) >> 1;
	int hd2 = (d2 - 1) >> 1;

	register int index = 0;
	for(register int i = 0; i < n; ++i) {
		for(register int j = 0; j < m; ++j) {
			index = i * m + j;
			b[index] = 0;
			for(register int l = -hd1; l <= hd1; ++l) {
				for(register int o = -hd2; o <= hd2; ++o) {
					int c = (j + o + m) % m;
					b[index] += k[(hd1 + l) * d2 + (hd2 +o)] * *(a + (i + l) * m + c);
				}
			}
			b[index] = b[index] / (d1 * d2);

		}
	}
}

int main(void) {

	cin.tie(0) -> sync_with_stdio(false);

	MPI_Init(NULL, NULL);
	MPI_Comm_size(WORLD, &nums);
	MPI_Comm_rank(WORLD, &id);

	int t, n, m, d1, d2;
	int *fullA, *k;

	if(id == 0) {
		string fileName;
		cin >> fileName;

		FILE* file = fopen(fileName.c_str(), "r");
		fscanf(file, "%d", &t);
		fscanf(file, "%d %d", &n, &m);

		fullA = (int*) calloc(n * m, sizeof(int));
		for(int i = 0, end = n * m; i < end; ++i) {
			fscanf(file, "%d", fullA + i);
		}
		
		fscanf(file, "%d %d", &d1, &d2);
		k = (int*) calloc(d1 * d2, sizeof(int));
		for(int i = 0, end = d1 * d2; i < end; ++i) {
			fscanf(file, "%d", k + i);
		}

		// make fullA contain head and tail
		int *temp = (int*) calloc((n + d1 - 1) * m, sizeof(int));
		int hk = (d1 - 1) / 2;

		memcpy(temp + hk * m, fullA, n * m * sizeof(int));
		free(fullA);
		fullA = temp;
		memcpy(fullA, fullA + n * m, hk * m * sizeof(int));
		memcpy(fullA + (n + hk) * m, fullA + hk * m, hk * m * sizeof(int));
	}

	int arr[5] = {t, n, m, d1, d2};
	MPI_Bcast(arr, 5, MPI_INT, 0, WORLD);

	t = arr[0], n = arr[1], m = arr[2], d1 = arr[3], d2 = arr[4];

	if(id != 0) k = (int*) calloc(d1 * d2, sizeof(int));
	MPI_Bcast(k, d1 * d2, MPI_INT, 0, WORLD);

	int rowPerProcess = n / nums;
	int remain = n % nums;
	int size = rowPerProcess + (id < remain ? 1 : 0);
	int from = size * id + min(id, remain);
	int halfK = (d1 - 1) / 2;

	int *a1, *a, *a2;
	a1 = (int*) calloc((size + d1 - 1) * m, sizeof(int));
	a = a1 + (halfK * m);
	a2 = a + size * m;

	int *b1, *b, *b2;
	b1 = (int*) calloc((size + d1 - 1) * m, sizeof(int));
	b = b1 + (halfK * m);
	b2 = b + size * m;


	// modify to distribute fullA
	int cnt[nums] = {0};
	int disp[nums] = {0};
	for(int i = 0; i < nums; ++i) {
		cnt[i] = rowPerProcess + (i < remain ? 1 : 0) + d1 - 1;
		cnt[i] *= m;
		if(i > 0) {
			disp[i] = disp[i - 1] + cnt[i - 1] - (d1 - 1) * m;
		}
	}

	MPI_Scatterv(fullA, cnt, disp, MPI_INT,
                 a1, (d1 - 1 + size) * m, MPI_INT,
                 0, WORLD);

	int cnt2[nums] = {0};
	int disp2[nums] = {0};
	for(int i = 0; i < nums; ++i) {
		cnt2[i] = rowPerProcess + (i < remain ? 1 : 0);
		cnt2[i] *= m;
		
		if(i > 0) {
			disp2[i] = disp2[i - 1] + cnt2[i - 1];
		}

	}

	int *temp;
	for(int i = 1; i <= t; ++i) {
		calc(a, b, k, size, m, d1, d2);

		
		MPI_Gatherv(b, size * m, MPI_INT, (fullA + halfK * m), cnt2, disp2, MPI_INT, 0, WORLD);

		if(i != t) {
			if(id == 0) {
				memcpy(fullA, fullA + n * m, halfK * m * sizeof(int));
				memcpy((fullA + ((n + halfK) * m)), (fullA + halfK * m), (halfK * m * sizeof(int)));
			}
			MPI_Scatterv(fullA, cnt, disp, MPI_INT, a1, (d1 - 1 + size) * m, MPI_INT, 0, WORLD);
		}
	}
	
	if(id == 0) {
		for(int i = 0, end = n * m; i < end; ++i) {
			cout << fullA[i + halfK * m] << " ";
		}
	}
	
	MPI_Finalize();
	return 0;
}

