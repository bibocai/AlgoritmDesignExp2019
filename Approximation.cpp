	// Approximation.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
	//
	#include <chrono>   

	#include <iostream>
	#include <string>
	#include <unordered_set>
	#include <stdlib.h>
	#include<vector>
	#include<iterator>
	#include <map>
	using namespace std;
	using namespace chrono;

	#define TotalEleNumOfF 1000
	#define MaxNumOfX 10000

	int SetCoverLP(int* ia, int* ja, int numMatrixValueNotZero, double* ar, double* x, int numX);//ia,ja,ar描述参数矩阵 从1开始 numMatrixValueNotZero描述参数矩阵中的有效参数个数 

	void genSetFamily(unordered_set<int>& X, vector<unordered_set<int> >& F,int n) {
	
		for (int i = 0; i < n; i++) {
			X.insert(i);
		}  //从0到n-1
		//cout << "************************输出X中的所有元素*****************" << endl;
		//copy(X.begin(), X.end(), ostream_iterator<int>(cout, " "));
		
	
		for (int i = 0; i < n; i++) {
			unordered_set<int> temp;
			int tempSize = rand() % n;
			for(int j=0;j<tempSize;j++){
				temp.insert(rand() % n);
			}
			F.push_back(temp);
		}
		//cout << endl;
		//cout << "*************************输出F中的所有集合:******************" << endl;
		//int Fcout = 0;
		//for (auto it = F.begin(); it != F.end(); it++) {
		//	cout << "集合" << Fcout++ << "中的元素为" << " ";
		//	copy((*it).begin(), (*it).end(), ostream_iterator<int>(cout, " "));
		//	cout << endl;
		//}
		//for (auto it = X.begin(); it != X.end(); ++it)
		//	std::cout << " " << *it;
	}

	int numIntersection(unordered_set<int>& x, unordered_set<int>& y) {
		int i = 0;
		for (auto it = x.begin(); it != x.end(); ++it) {
			if (y.find(*it) != y.end()) {//如果在y中找到x中的元素
				i++;
			}
		}
		return i;
	}
	vector<int> GreedySetCover(unordered_set<int>& X,vector<unordered_set<int>> & F) {
		unordered_set<int> U(X);
		vector<int> C;
	//	vector<bool> mask(F.size(),0); //以免有不中指的风险 
	//	bool correct = 1;
		while (!U.empty()) {
			int maxNum = numIntersection(F[0], U); int maxIdx = 0;
			int newMax = numIntersection(F[maxIdx], U);

			for (int idx = 1; idx < F.size(); idx++) {
				newMax = numIntersection(F[idx], U);
				if (newMax  > maxNum) { maxNum = newMax; maxIdx = idx; } 
			}//第3步

			for (auto it = F[maxIdx].begin(); it != F[maxIdx].end(); ++it) {
				unordered_set<int>::iterator itU = U.find(*it);
				if (itU != U.end()) { U.erase(itU); }
			} //第4步
	//		mask[maxIdx] = 1;
			C.push_back(maxIdx);
		}
//		cout << endl << "选中的集合的索引是" << endl;
		copy(C.begin(), C.end(), ostream_iterator<int>(cout, " "));

		return C;

	}


	vector<int>  SetCoverRounding(unordered_set<int>& X,vector<unordered_set<int>> & F, double* x){ //x是输出的解 x的大小得够大
		//生成线性规划相关矩阵 初始化都是0;

		//int ia[1 + TotalEleNumOfF], ja[1 + TotalEleNumOfF];
		//double ar[1 + TotalEleNumOfF]; //这三个函数都得索引从1开始
		//double x[MaxNumOfX];
		int numX = X.size();
		int numMatrixValueNotZero = 0;
		for (auto it = F.begin(); it != F.end(); it++) {
			numMatrixValueNotZero += (*it).size();
		}
		int* ia = new int[1 + numMatrixValueNotZero];
		int* ja = new int[1 + numMatrixValueNotZero];
		double* ar = new double[1 + numMatrixValueNotZero];
		//统计每个元素在集族中的最大频率
		map<int,int> counter;
		for(auto it=X.begin();it!=X.end();it++){
			counter[ *it]=0;
		}
		map<int, int> num2idx;
		int invertIdx = 0;
		for (auto it = X.begin(); it != X.end(); it++) {
			num2idx[*it] = invertIdx++;
		}

		//在遍历过程中 初始化 参数矩阵ia,ja,ar 和 counter
		int matrixIdx = 0;
		for(auto it=F.begin();it!=F.end();++it){
			//横轴代表点 纵轴代表集合
			int mj = it - F.begin();
			for(auto setIt = (*it).begin();setIt!=(*it).end();++setIt){
				ia[1 + matrixIdx] = 1 + num2idx[ *setIt ]; ja[1 + matrixIdx] = 1+ mj; ar[1 + matrixIdx] = 1.0;
				 matrixIdx++;
				counter[ *setIt ]++;
			}
		
		}
		//输出数据倒排索引
		//cout << "输出倒排索引即 std::map num2idx的值为（从元素的值映射到元素在unorder_set中的位置）" << endl;
		//for (auto it = num2idx.begin(); it != num2idx.end(); ++it) {
		//	cout << it->first << "->" << it->second << "    ";
		//}
		//cout << endl;
		//输出参数矩阵
		//cout << "输出参数矩阵" << endl;
		//float(*mx)[MaxNumOfX];
		//mx = new float[numX][MaxNumOfX];
		//for (int i = 0; i < numX; i++) {
		//	for (int j = 0; j < numX; j++) {
		//	 mx[i][j]=0;
		//	}
		//}
		//for (int i = 0; i < numMatrixValueNotZero; i++) {
		//	int mi = ia[1 + i]; int mj = ja[1 + i]; float r = ar[1 + i];
		//	mx[mi-1][mj-1] = r;
		//}
		//for (int i = 0; i < numX; i++) {
		//	for (int j = 0; j < numX; j++) {
		//		cout << mx[i][j] << " ";			
		//	}
		//	cout << endl;
		//}
		//delete[] mx;
		//cout << "输出ia，ja，ar中的值" << endl;
		//for (int i = 0; i < numMatrixValueNotZero; i++) {
		//	cout <<"("<< *(ia + i + 1) <<", "<< *(ja + i + 1) <<") "<< *(ar + i + 1) << endl;
		//}



		int f=0;
		for(auto it=counter.begin();it!=counter.end();++it){
			if( it->second > f) f=it->second;
		}
		cout << "元素的出现的最大频率是 " << f << endl;
		//线性规划
		SetCoverLP(ia, ja, numMatrixValueNotZero, ar, x, numX);
		delete[] ia, ja, ar;
		//输出最后结果

		std::vector<int> out;
		cout << "prime solution" << endl;
		for (int i = 0; i < numX; i++) {
			cout << x[i] << ' ';
			
		}
		for (int i = 0; i < numX; i++) {
			if (x[i] > 1 / float(f))  out.push_back(i);
		}
		cout << endl << "选中的集合的索引是" << endl;
		copy(out.begin(), out.end(), ostream_iterator<int>(cout, " "));

		return out;
	}

	int main()
	{

		unordered_set<int> X;
		vector<unordered_set<int>> F;
		genSetFamily(X,F,100);
		auto start = system_clock::now();

		std::vector<int> C = GreedySetCover(X ,F);

		auto end = system_clock::now();
		auto duration = duration_cast<microseconds>(end - start);

		cout << "花费了"
			<< double(duration.count()) * microseconds::period::num / microseconds::period::den
			<< "秒" << endl;

		double* x = new double[5000];

		start = system_clock::now();

		std::vector<int> out = SetCoverRounding( X, F, x);  //x是输出的解 x的大小得够大
		delete[] x;


		end = system_clock::now();		duration = duration_cast<microseconds>(end - start);		cout << "花费了"
			<< double(duration.count()) * microseconds::period::num / microseconds::period::den
			<< "秒" << endl;
	}


	// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
	// 调试程序: F5 或调试 >“开始调试”菜单

	//对于unorderset 迭代每次的输出是否是一致的的？
