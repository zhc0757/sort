//#define TESTING

//#define SERIAL

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<algorithm>
#include<chrono>
#include<functional>
#include<list>
#include<mutex>
#include<forward_list>
#include<thread>
#include<omp.h>
#include<climits>
#include<atomic>
#include<numeric>
#include<vector>
using namespace std;
using namespace chrono;
namespace nsSerial
{
	/*
	Initialize integer array d[N] with seed
	*/
	void sort_gen(int *d, int N, int seed)
	{
		srand(seed);
		for (int i = 0; i < N; i++) {
			d[i] = rand();
		}
	}
	int cmp(const void*a, const void*b)
	{
		return (*(int*)a) - (*(int*)b);
	}
	const unsigned BUCKET_LEN = ((unsigned)RAND_MAX) + 1;
	int bucketSort(int n, int seed)
	{
		unsigned*bucket = new unsigned[BUCKET_LEN];
		int y = -1;
		unsigned i, m = n >> 1, cnt = 0, n2 = n;
		fill_n(bucket, BUCKET_LEN, 0u);
		srand(seed);
		for (i = 0; i < n2; ++i) {
			++bucket[rand()];
		}
		/*for (i = 0; i < BUCKET_LEN; ++i) {
			cnt += bucket[i];
			if (cnt >= m) {
				y = (int)i;
				break;
			}
		}*/
		//considering d[0]
		for (i = 0; cnt <= m; ++i)cnt += bucket[i];
		y = i - 1;
		delete[]bucket;
		return y;
	}
}
namespace nsParallelRadix
{
	typedef list<int>* plist;
	const int THREAD_NUM = 16,
		RADIX_LEN = 8,
		RADIX = 1 << RADIX_LEN,
		R_MAX = RADIX - 1,
		T0 = sizeof(int) * 8 / RADIX_LEN - 1,
		RD0 = T0*RADIX_LEN,
		R0 = R_MAX << RD0;
	template<typename T>void sortGen(T from, T to, int seed)
	{
		srand(seed);
		//for (T i = from; i != to; ++i)*i = rand();
		generate(from, to, rand());
	}
	void bucketRadixSort(list<int>&b0)
	{
		if (!b0.empty()) {
			int buckets[RADIX] = { 0 }, base = b0.back()&(~R_MAX), i, j, t;
			list<int>::iterator it;
			for (it = b0.begin(); it != b0.end(); ++it)++buckets[(*it)&R_MAX];
			b0.clear();
			for (i = 0; i < RADIX; ++i) {
				t = base + i;
				for (j = 0; j < buckets[j]; ++j)b0.push_back(t);
			}
		}
	}
	void subRadixSort(list<int>&b0, int rd)
	{
		list<int>buckets[RADIX];
		list<int>::iterator it;
		int r = R_MAX << rd, rd_next = rd - RADIX_LEN, i;
		for (it = b0.begin(); it != b0.end(); it++)buckets[((*it)&r) >> rd].push_back(*it);
		b0.clear();
		if (rd_next > 0)for (i = 0; i < RADIX; ++i)subRadixSort(buckets[i], rd_next);
		else for (i = 0; i < RADIX; ++i)bucketRadixSort(buckets[i]);
		for (i = 0; i < RADIX; ++i) {
			/*for (it = buckets[i].begin(); it != buckets[i].end(); ++it)
				b0.push_back(*it);*/
			b0.splice(b0.end(), buckets[i]);
		}
	}
	const int NOT_R_MAX = ~R_MAX;
	void bucketRadixSort(int*d, int n)
	{
		if (n) {
			int bucket[RADIX] = { 0 }, base = (*d)&NOT_R_MAX, i, j, k = 0, t;
			for (i = 0; i < n; ++i)
				++bucket[d[i] & R_MAX];
			for (i = 0; i < RADIX; ++i) {
				t = base + i;
				for (j = 0; j < bucket[i]; ++j) {
					d[k] = t;
					k++;
				}
			}
		}
	}
	void subRadixSort(int*&d, int n, int rd)
	{
		if (n) {
			int*bucket[RADIX];
			int r = R_MAX << rd, rd_next = rd - RADIX_LEN, i, t;
			int len[RADIX] = { 0 };
			int*p;
			for (i = 0; i < n; ++i)++len[(d[i] & r) >> rd];
			for (i = 0; i < RADIX; ++i)if (len[i])bucket[i] = new int[len[i]];
			fill_n(len, RADIX, 0);
			for (i = 0; i < n; ++i) {
				t = (d[i] & r) >> rd;
				bucket[t][len[t]++] = d[i];
			}
			delete[]d;
			if (rd_next/* > RADIX_LEN*/)for (i = 0; i < RADIX; ++i)subRadixSort(bucket[i], len[i], rd_next);
			else for (i = 0; i < RADIX; ++i)if (len[i])bucketRadixSort(bucket[i], len[i]);
			d = new int[n];
			p = d;
			for (i = 0; i < RADIX; ++i) {
				if (len[i]) {
					p = copy_n(bucket[i], len[i], p);
					delete[]bucket[i];
				}
			}
		}
	}
	const int L = RADIX >> 1, RD_N = RD0 - RADIX_LEN;
	void radixSort(int*&d, int n)
	{
		int i, j = 0;
		/*list<int>buckets[L];
		list<int>::iterator it;*/
		mutex mtx[L];
		int len[L] = { 0 };
		int*bucket[L];
		int*p;
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < n; ++i) {
			int t = (d[i] & R0) >> RD0;
			mtx[t].lock();
			++len[t];
			mtx[t].unlock();
		}
		for (i = 0; i < L; ++i)if (len[i])bucket[i] = new int[len[i]];
		fill_n(len, L, 0);
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < n; ++i) {
			int t = (d[i] & R0) >> RD0;
			mtx[t].lock();
			//buckets[t].push_back(d[i]);
			bucket[t][len[t]++] = d[i];
			mtx[t].unlock();
		}
		delete[]d;
//#pragma omp parallel for num_threads(THREAD_NUM)
//		for (i = 0; i < L; ++i)subRadixSort(buckets[i], RD_N);
//		d = new int[n];
//		for (i = 0; i < L; ++i) {
//			for (it = buckets[i].begin(); it != buckets[i].end(); ++it) {
//				d[j] = *it;
//				++j;
//			}
//			buckets[i].clear();
//		}
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < L; ++i)subRadixSort(bucket[i], len[i], RD_N);
		d = new int[n];
		p = d;
		for (i = 0; i < L; ++i) {
			if (len[i]) {
				p = copy_n(bucket[i], len[i], p);
				delete[]bucket[i];
			}
		}
	}
	/*void radixSort(list<int>&d)
	{
		const int L = RADIX >> 1, RD_N = RD0 - RADIX_LEN;
		int i;
		list<int>buckets[L];
		list<int>::iterator it;
		mutex mtx[L];
#pragma omp parallel for num_threads(THREAD_NUM)
		for (it = d.begin(); it != d.end(); ++it) {
			int t = (*it & R0) >> RD0;
			mtx[t].lock();
			buckets[t].push_back(*it);
			mtx[t].unlock();
		}
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < L; ++i)subRadixSort(buckets[i], RD_N);
		for (i = 0; i < L; ++i) {
			d.splice(d.end(), buckets[i]);
		}
	}*/
	const double P = 1.1, P_RD8 = 1.5;
	const unsigned Q = nsSerial::BUCKET_LEN / THREAD_NUM,
		R = nsSerial::BUCKET_LEN%THREAD_NUM,
		NR = nsSerial::BUCKET_LEN - R,
		Q2 = Q - 1,
		NR2 = NR + 1;
	int*p;
#ifdef TESTING
	atomic_int tooSmallTime[4];
#endif
	void bucketRadixSort2(vector<int>&d)
	{
		if (d.empty())d.shrink_to_fit();
		else {
			int*begin = p + (d.back()&NOT_R_MAX);
			int i;
			for (i = 0; i < d.size(); ++i)begin[d[i] & R_MAX]++;
			d.clear();
			d.shrink_to_fit();
		}
	}
	void subRadixSort2(vector<int>&d, int rd)
	{
		if (d.empty())d.shrink_to_fit();
		else {
			int r = R_MAX << rd, rd_next = rd - RADIX_LEN, i;
			int bucketSize0 = (int)((d.size() / RADIX) * (rd_next ? P : P_RD8));
			vector<int>bucket[RADIX];
			for (i = 0; i < RADIX; ++i)bucket[i].reserve(bucketSize0);
			for (i = 0; i < d.size(); ++i)bucket[(d[i] & r) >> rd].push_back(d[i]);
			d.clear();
			d.shrink_to_fit();
			if (rd_next > 0)
				for (i = 0; i < RADIX; ++i) {
#ifdef TESTING
					//if (bucket[i].size() > bucketSize0)cout << "Too small: " << rd << endl;
					if (bucket[i].size() > bucketSize0)tooSmallTime[rd / RADIX_LEN]++;
#endif
					subRadixSort2(bucket[i], rd_next);
				}
			else for (i = 0; i < RADIX; ++i) {
#ifdef TESTING
				//if (bucket[i].size() > bucketSize0)cout << "Too small: " << rd << endl;
				if (bucket[i].size() > bucketSize0)tooSmallTime[rd / RADIX_LEN]++;
#endif
				bucketRadixSort2(bucket[i]);
			}
		}
	}
	int radixSort2(int n, int seed)
	{
		int bucketSize0 = (int)((n / L)*P);
		vector<int>bucket[L];
		int i, t;
#ifdef TESTING
		for (auto& ele : tooSmallTime)ele = 0;
		auto t0 = high_resolution_clock::now();
#endif
		for (i = 0; i < L; ++i)bucket[i].reserve(bucketSize0);
		srand(seed);
		for (i = 0; i < n; ++i) {
			t = rand();
			bucket[t >> RD0].push_back(t);
		}
		p = new int[nsSerial::BUCKET_LEN];
		fill_n(p, nsSerial::BUCKET_LEN, 0);
#ifdef TESTING
		cout << duration_cast<duration<double>>(high_resolution_clock::now() - t0).count() << endl;
#endif
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < L; ++i) {
#ifdef TESTING
			//if (bucket[i].size() > bucketSize0)cout << "Too small: " << RD0 << endl;
			if (bucket[i].size() > bucketSize0)tooSmallTime[RD0 / RADIX_LEN]++;
#endif
			subRadixSort2(bucket[i], RD_N);
		}
#ifdef TESTING
		cout << duration_cast<duration<double>>(high_resolution_clock::now() - t0).count() << endl;
#endif
		int s[THREAD_NUM];
#pragma omp parallel num_threads(THREAD_NUM)
		{
			int thNum = omp_get_thread_num();
			unsigned tmp = thNum*Q;
			s[thNum] = accumulate(p + tmp, p + (tmp + Q), 0);
		}
		int m = n >> 1, y = -1;
		unsigned j;
		t = 0;
		for (i = 0; i < THREAD_NUM; ++i) {
			t += s[i];
			if (t > m) {
				m -= t - s[i];
				t = 0;
				for (j = i*Q; t <= m; ++j)t += p[j];
				y = (int)j - 1;
				goto end;
			}
		}
		m -= t;
		t = 0;
		for (j = NR; t <= m; ++j)t += p[j];
		y = (int)j - 1;
	end:
#ifdef TESTING
		cout << duration_cast<duration<double>>(high_resolution_clock::now() - t0).count() << endl;
		for (auto& ele : tooSmallTime)cout << ele << ' ';
		cout << endl;
#endif
		delete[]p;
		return y;
	}
}
namespace nsParallelMerge
{
	void mergeSort(int*d, int n, int threadCnt = 1)
	{
		if (threadCnt >= nsParallelRadix::THREAD_NUM)
			sort(d, d + n);
		else {
			int thCnt_next = threadCnt << 1, n2 = n >> 1, n3 = n - n2;
			int*d2 = d + n2;
			thread th1(mergeSort, d, n2, thCnt_next), th2(mergeSort, d2, n3, thCnt_next);
			th1.join();
			th2.join();
			int*p = new int[n];
			merge(d, d + n2, d2, d2 + n3, p);
			copy_n(p, n, d);
			delete[]p;
		}
	}
	const int L = nsParallelRadix::RADIX >> 1, THREAD_NUM = 12;
	const int L2 = L - 1;
	union ULen {
		static const int CACHE_LINE_LEN = 64 * 2;
		atomic_int l;
		char padding[CACHE_LINE_LEN];
	};
	void radixMergeSort(int*&d, int n)
	{
		ULen len[L];
		int len2[L];
		int*bucket[L];
		int i;
		for (i = 0; i < L; ++i)len[i].l = 0;
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < n; ++i)++len[d[i] >> nsParallelRadix::RD0].l;
		for (i = 0; i < L; ++i)if (len[i].l)bucket[i] = new int[len[i].l];
		for (i = 0; i < L; ++i)len[i].l = 0;
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < n; ++i) {
			int t = d[i] >> nsParallelRadix::RD0;
			bucket[t][len[t].l++] = d[i];
		}
		delete[]d;
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < L; ++i)
			if (len[i].l)
				sort(bucket[i], bucket[i] + len[i].l);
		d = new int[n];
		len2[0] = 0;
		for (i = 0; i < L2; ++i)len2[i + 1] = len[i].l + len2[i];
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < L; ++i) {
			if (len[i].l) {
				copy_n(bucket[i], (int)len[i].l, d + len2[i]);
				delete[]bucket[i];
			}
		}
	}
	void radixMergeSort2(int*&d, int n)
	{
#define CLEAR_LEN for (i = 0; i < L; ++i)len[i].l = 0;
		using namespace chrono;
		ULen len[L];
		int*bucket[L + 1];
		int i;
		CLEAR_LEN;
		//for (i = 0; i < L; ++i)len[L] = 0;
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < n; ++i)
			len[d[i] >> nsParallelRadix::RD0].l++;
#ifdef TESTING
		auto t0 = high_resolution_clock::now();
#endif
		bucket[0] = new int[n];
#ifdef TESTING
		cout << duration_cast<duration<double>>(high_resolution_clock::now() - t0).count() << endl;
#endif
		for (i = 0; i < L; ++i) {
			bucket[i + 1] = bucket[i] + len[i].l;
		}
		CLEAR_LEN;
		//for (i = 0; i < L; ++i)len[L] = 0;
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < n; ++i) {
			int t = d[i] >> nsParallelRadix::RD0;
			bucket[t][len[t].l++] = d[i];
		}
		delete[]d;
		d = bucket[0];
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < L; ++i)
			if (len[i].l)
				sort(bucket[i], bucket[i + 1]);
	}
	void radixMergeSort3(int*&d, int n)
	{
		int len[L] = { 0 };
		int*bucket[L + 1];
		int i, t;
		for (i = 0; i < n; ++i)len[d[i] >> nsParallelRadix::RD0]++;
		bucket[0] = new int[n];
		for (i = 0; i < L; ++i)bucket[i + 1] = bucket[i] + len[i];
		fill_n(len, L, 0);
		for (i = 0; i < n; ++i) {
			t = d[i] >> nsParallelRadix::RD0;
			bucket[t][len[t]++] = d[i];
		}
		delete[]d;
		d = bucket[0];
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < L; ++i)
			if (len[i])
				sort(bucket[i], bucket[i + 1]);
	}
	int radixMergeSort4(int n, int seed)
	{
		const double P = 1.1;
		int bucketSize0 = (int)((n / L) * P);
		vector<int>bucket[L];
		int i, t, m = n >> 1, s = 0;
		for (i = 0; i < L; ++i)bucket[i].reserve(bucketSize0);
		srand(seed);
		for (i = 0; i < n; ++i) {
			t = rand();
			bucket[t >> nsParallelRadix::RD0].push_back(t);
		}
#pragma omp parallel for num_threads(THREAD_NUM)
		for (i = 0; i < L; ++i)sort(bucket[i].begin(), bucket[i].end());
		for (i = 0; s <= m; ++i)s += (int)bucket[i].size();
		t = i - 1;
		return bucket[t][bucket[t].size() - (s - m)];
	}
}
namespace nsParaBucket
{
	const int THREAD_CNT = 12;
	const unsigned BUCKET_LEN =
#ifdef TESTING
		INT_MAX + 1u
#else
		nsSerial::BUCKET_LEN
#endif
		;
	const unsigned Q = BUCKET_LEN / THREAD_CNT,
		R = BUCKET_LEN%THREAD_CNT,
		NR = BUCKET_LEN - R,
		Q2 = Q - 1,
		NR2 = NR + 1;
	//const int I_NR = NR, I_Q = Q;
	void subFill(int*bucket)
	{
		fill_n(bucket + (omp_get_thread_num()*Q), Q, 0);
	}
	void getSum(int*bucket)
	{
		unsigned i = Q*omp_get_thread_num(), last = i + Q;
		for (++i; i < last; ++i)bucket[i] += bucket[i - 1];
	}
	int paraBucketSort(int n, int seed)
	{
		using namespace chrono;
		int*bucket = new int[BUCKET_LEN];
		int y = -1, m = n >> 1, cnt = 0, k;
		/*const lldiv_t divt = lldiv(nsSerial::BUCKET_LEN, THREAD_CNT);
		const int nr = nsSerial::BUCKET_LEN - divt.rem;*/
		unsigned i;
#ifdef TESTING
		auto t0 = high_resolution_clock::now();
#endif
		//fill_n(bucket, BUCKET_LEN, 0);
#pragma omp parallel num_threads(THREAD_CNT)
		subFill(bucket);
		fill_n(bucket + NR, R, 0);
		srand(seed);
		for (k = 0; k < n; ++k)++bucket[rand()];
#ifdef TESTING
		cout << duration_cast<duration<double>>(high_resolution_clock::now() - t0).count() << endl;
#endif
//#pragma omp parallel for num_threads(THREAD_CNT)
//		for (k = 0; k < I_NR; k += I_Q) {
//			//printf("%d\n", k);
//			unsigned j, last = k + Q;
//			for (j = k + 1u; j < last; ++j)
//				bucket[j] += bucket[j - 1];
//		}
#pragma omp parallel num_threads(THREAD_CNT)
		getSum(bucket);
		for (i = NR2; i < BUCKET_LEN; ++i)bucket[i] += bucket[i - 1];
		for (i = Q2; i < NR; i += Q) {
			cnt += bucket[i];
			if (cnt > m) {
				unsigned j;
				m -= cnt - bucket[i];
				//cnt = 0;
				for (j = i - Q2; /*cnt*/bucket[j] <= m; ++j)/*cnt += bucket[j]*/;
				//y = j - 1;
				y = (int)j;
				goto end;
			}
		}
		m -= cnt;
		//cnt = 0;
		for (i = NR; /*cnt*/bucket[i] <= m; ++i)/*cnt += bucket[i]*/;
		//y = i - 1;
		y = (int)i;
	end:
		delete[]bucket;
		return y;
	}
	int paraBucketSort2(int n, int seed)
	{
		using namespace chrono;
		int*bucket = new int[BUCKET_LEN];
		int y = -1, m = n >> 1, cnt = 0, k;
		unsigned i;
		thread*th[THREAD_CNT];
		int s[THREAD_CNT];
		fill_n(bucket, BUCKET_LEN, 0);
		srand(seed);
#ifdef TESTING
		auto t0 = high_resolution_clock::now();
#endif
		for (k = 0; k < n; ++k)bucket[rand()]++;
#ifdef TESTING
		cout << duration_cast<duration<double>>(high_resolution_clock::now() - t0).count() << endl;
#endif
		for (k = 0; k < THREAD_CNT; ++k) {
			th[k] = new thread([&](int thNum) {
				unsigned t = thNum*Q;
				s[thNum] = accumulate(bucket + t, bucket + (t + Q), 0);
			}, k);
		}
		for (k = 0; k < THREAD_CNT; ++k) {
			th[k]->join();
			delete th[k];
		}
		for (k = 0; k < THREAD_CNT; ++k) {
			cnt += s[k];
			if (cnt > m) {
				m -= cnt - s[k];
				cnt = 0;
				for (i = k*Q; cnt <= m; ++i)cnt += bucket[i];
				y = (int)i - 1;
				goto end;
			}
		}
		m -= cnt;
		cnt = 0;
		for (i = NR; cnt <= m; ++i)cnt += bucket[i];
		y = (int)i - 1;
	end:
		delete[]bucket;
		return y;
	}
}
int calculate(int n, int seed)
{
	//return nsSerial::bucketSort(n, seed);
	//return nsParaBucket::paraBucketSort2(n, seed);
	//return nsParallelMerge::radixMergeSort4(n, seed);
	//return nsParallelRadix::radixSort2(n, seed);
	int mid;
	int*d = new int[n];
	nsSerial::sort_gen(d, n, seed);
#ifdef SERIAL
	//qsort(d, n, sizeof(int), nsSerial::cmp);
	sort(d, d + n);
#endif
	//nsParallelRadix::radixSort(d, n);
	nsParallelMerge::mergeSort(d, n);
	//nsParallelMerge::radixMergeSort3(d, n);
	mid = d[n >> 1];
	delete[]d;
	return mid;
}
/*
input: 1000000 1
output: 16386
input: 1000000000 1
output: 16384
*/
int main(int argc, char*argv[])
{
	using namespace chrono;
	int n, seed;
	switch (argc) {
	case 3:
		n = atoi(argv[1]);
		seed = atoi(argv[2]);
		break;
	default:
		cin >> n >> seed;
	}
	auto t0 = high_resolution_clock::now();
	printf("%d", calculate(n, seed));
	printf(" %lf", duration_cast<duration<double>>(high_resolution_clock::now() - t0).count());
	printf("\n");
	system("pause");
	return 0;
}
