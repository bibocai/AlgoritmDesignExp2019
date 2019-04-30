// Exp_quicksort.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include <chrono>   
#include <iostream>
using namespace std;
#include <vector>
#include<random>
using namespace chrono;
#include<algorithm>
#include<tuple>
void printVector(const std::vector<int>& L) {
	int len = L.size();
	for (int i = 0; i < len; i++) {
		cout << L[i] << ' ';
		
	}
	cout << endl;
}


std::random_device rd;
std::default_random_engine e = std::default_random_engine(rd());
inline int Random(int p, int r) {
	 //生成随机数 随机数中可以包含p和r
	uniform_int_distribution<int> u(p,r);

	return u(e);

}

void _exchange(std::vector<int> &L, int p, int r) {
	int temp = L[p]; L[p] = L[r]; L[r] = temp;

}


std::tuple<int, int> newPartition(std::vector<int>& A, int p, int r) {
	int i = Random(p, r);
	_exchange(A, r, i);
	int x = A[r];
	//    cout<<"choosen ele"<<x<<endl;
	i = p - 1;
	int k = p - 1;
	for (int j = p; j <= r - 1; j++) {
		if (A[j] < x) {
			i++;
			_exchange(A, j, i);
			k++;
			if (i < k) {
				_exchange(A, j, k);
			}
		}
		else if (A[j] == x)
		{
			k++;
			_exchange(A, j, k);
		}	
	}
	_exchange(A, k + 1, r);

	return { i,k + 1 };
}

void newQuickSort(std::vector<int>& L, int p, int r) {
	if (p < r) {
		std::tuple<int,int> rn = newPartition(L, p, r);
		int rn1, rn2;
		std::tie(rn1, rn2) = rn;
		newQuickSort(L, p, rn1);
		newQuickSort(L, rn2+1, r);
	}
}

int Modify_Rand_Partition(std::vector<int>& L, int p, int r) {
	int i = Random(p, r);
	int allSameFlag = 1;
	//	cout << "generate random "<<i<< endl ;
	_exchange(L, i, r);
	i = p - 1;
	for (int j = p; j <= r - 1; j++) {
		if (L[j] <= L[r]) {
			if (L[j] < L[r]) allSameFlag = 0;//只要有一次不是等于 则说明不是全等
			i++;
			_exchange(L, i, j);
		}
		else {
			allSameFlag = 0;
		}
	}
	_exchange(L, i + 1, r);
	if (allSameFlag == 1)
		return -1;
	return i + 1;
}

void Modify_QuickSort(std::vector<int>& L, int p, int r) {
	if (p < r) {
		int q = Modify_Rand_Partition(L, p, r);
		if (q == -1) return;
		Modify_QuickSort(L, p, q - 1);
		Modify_QuickSort(L, q + 1, r);
	}
}
int Rand_Partition(std::vector<int> &L, int p, int r){
	int i = Random(p, r);
//	cout << "generate random "<<i<< endl ;
	_exchange(L, i, r);
	i = p - 1;
	for (int j = p; j <= r - 1;j++) {
		if (L[j] <= L[r]) {
			i++;
			_exchange(L, i, j);
		}
	}
	_exchange(L, i + 1, r);
	return i + 1;
}

void QuickSort(std::vector<int>& L, int p, int r) {
	if (p < r) {
		int q = Rand_Partition(L, p, r);
		QuickSort(L, p, q - 1);
		QuickSort(L, q + 1, r);
	}
}
std::vector<int> gen_randVec(int length) {
	std::vector<int> vec(length);
	std::random_device rd;
	std::default_random_engine e = std::default_random_engine(rd());
	uniform_int_distribution<int> u; //生成随机数 随机数中可以包含p和r
	cout << "minimum " << u.min() << endl;
	cout << "maximum" << u.max() << endl;
	for (int i = 0; i < length; i++) {
		vec[i] = u(e);
	}
	return vec;
}
std::vector<int> genMixVec(int length, float percentOfOne) {
	std::vector<int> vec(length, 10000);
	std::random_device rd;
	std::default_random_engine e = std::default_random_engine(rd());
	uniform_int_distribution<int> u; //生成随机数 随机数中可以包含p和r
	for (int i = int(length * percentOfOne); i < length; i++) {
		vec[i] = u(e);
	}
	return vec;
}

int main()
{
	std::vector<double> pct = { 0.5,0.6,0.7,0.8,0.9,1 };
	for (int i = 0; i < pct.size(); i++) {
		std::vector<int> L = genMixVec(1000000, pct[i]);
//		std::vector<int> L = gen_randVec(100);
		cout << "数据规模" << L.size() << endl;
		cout << "输入vec中1的比例" << pct[i]<<endl;
		auto start = system_clock::now();
//		QuickSort(L, 0, L.size() - 1);
		std::sort(L.begin(), L.end());
//		Modify_QuickSort(L, 0, L.size() - 1);
//		newQuickSort(L, 0, L.size() - 1);
//		cout << endl;
		auto end = system_clock::now();
		auto duration = duration_cast<microseconds>(end - start);

		cout << "花费了"
			<< double(duration.count()) * microseconds::period::num / microseconds::period::den
			<< "秒" << endl;
//	printVector(L);

	}
	//std:vector<int> L = gen_randVec(1000000);


	

	




}	

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
