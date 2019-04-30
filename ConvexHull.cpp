
#include <iostream>
#include <vector>
using namespace std;
#include<random>
#include <chrono>   
#include <deque>
#include<algorithm>
#include <chrono>   
#include<iterator>
using namespace chrono;

class Vector2d
{
public:
	double x_;
	double y_;

public:
	Vector2d(double x, double y) :x_(x), y_(y) {}
	Vector2d() {
		std::random_device rd;
		std::default_random_engine e = std::default_random_engine(rd());
		uniform_real_distribution<double> u(0.0, 100.0);
		x_ = u(e);
		y_ = u(e);
	}
	double x() const { return x_; }
	double y() const { return y_; }
	//二维向量叉乘, 叉乘的结果其实是向量，方向垂直于两个向量组成的平面，这里我们只需要其大小和方向
	double CrossProduct(const Vector2d vec)
	{
		return x_ * vec.y_ - y_ * vec.x_;
	}

	//二维向量点积
	double DotProduct(const Vector2d vec)
	{
		return x_ * vec.x_ + y_ * vec.y_;
	}

	//二维向量减法
	Vector2d Minus(const Vector2d vec) const
	{
		return Vector2d(x_ - vec.x_, y_ - vec.y_);
	}

	//判断点M,N是否在直线AB的同一侧
	static bool IsPointAtSameSideOfLine(const Vector2d & pointM, const Vector2d & pointN,
		const Vector2d & pointA, const Vector2d & pointB)
	{
		Vector2d AB = pointB.Minus(pointA);
		Vector2d AM = pointM.Minus(pointA);
		Vector2d AN = pointN.Minus(pointA);

		//等于0时表示某个点在直线上
		return AB.CrossProduct(AM)* AB.CrossProduct(AN) >= 0;
	}
};

void ShowPoints(const std::vector <Vector2d>& Points) {
	cout << "生成的点中X的坐标序列是" << endl;

	for (int i = 0; i < Points.size(); i++) {
		if (i == Points.size() - 1) cout << Points[i].x() << endl;
		else	    cout << Points[i].x() << ",";
	}
	cout << endl;
	cout << "生成序列中Y的坐标序列是" << endl;
	for (int i = 0; i < Points.size(); i++) {
		if (i == Points.size() - 1) cout << Points[i].y() << endl;
		else	    cout << Points[i].y() << ",";
	}
	cout << endl;
}


void ShowPoint(const Vector2d & p)
{
	cout <<  "Point : (" << p.x() << ", " << p.y() << ")" << endl;
}
//三个条件均成立时认为成立
bool IsPointInTriInner(const Vector2d& PointI, const Vector2d& PointJ, const Vector2d& PointK, const Vector2d& PointX) {
	return Vector2d::IsPointAtSameSideOfLine(PointK, PointX,PointI, PointJ ) && \
		Vector2d::IsPointAtSameSideOfLine(PointJ, PointX,PointI, PointK) && \
		Vector2d::IsPointAtSameSideOfLine(PointI, PointX,PointK, PointJ);
}

std::vector<int> BruteForceConvexHull(const std::vector<Vector2d>& Points) {
	//ShowPoints(Points);
	std::vector<bool> mark(Points.size(), 0); //如果确定是某个三角形内部的点 则置为1
	for (int i = 0; i < Points.size(); i++) {
		if (mark[i]) continue;
		for (int j = 0; j < Points.size(); j++) {
			if (j == i || mark[j]) continue;
			for (int k = 0; k < Points.size(); k++) {
				if (k == i || k == j || mark[k]) continue;
				for (int x = 0; x < Points.size(); x++) {
					if (x == i || x == j || x == k || mark[x]) continue;
					if (IsPointInTriInner(Points[i], Points[j], Points[k], Points[x])) {

						mark[x] = 1;
					}

				}
			}
		}
	}

	std::vector<int> convexHullIndex;
	for (int i = 0; i < Points.size(); i++) {
		if (!mark[i]) { convexHullIndex.push_back(i); }
	}

	return convexHullIndex;
}

void testBruteForce() {
	auto start = system_clock::now();
	auto end = system_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	vector<int> setSize;
	for (int i =500 ; i <501 ; i++) {
		setSize.push_back(1000 + 1000 * i);
	}
	cout << "花费时间" << endl;
	for (int i = 0; i < setSize.size(); i++) {
		std::vector<Vector2d> Points(setSize[i]);
		start = system_clock::now();
		std::vector<int> CH = BruteForceConvexHull(Points);
		end = system_clock::now();
		duration = duration_cast<microseconds>(end - start);
		cout << double(duration.count()) * microseconds::period::num / microseconds::period::den << " ";
	}
}

void test() {
	std::vector<Vector2d> x(10);
	cout << "测试生成的Point点是否随机" << endl;
//	for (int i = 0; i < 10; i++) {
//		ShowPoint(x[i]);
//  }

	cout << "测试凸包算法是否有效" << endl;
	std::vector<Vector2d> Points(30);
	ShowPoints(Points);
	std::vector<bool> mark(30, 0); //如果确定是某个三角形内部的点 则置为1
	
	///*设定固定的x,y供测试用*/
	//std::vector<double> xx = { 98.0,42.0,6,34 };
	//std::vector<double> yy = { 70,60.0,81.0,8 };
	//for (int i = 0; i < 4; i++) {
	//	Points[i].x_ = xx[i];
	//	Points[i].y_ = yy[i];
	//}
	///*测试代码完毕*/

	for (int i = 0; i < Points.size(); i++) {
		if (mark[i]) continue;
		for (int j = 0; j < Points.size(); j++) {
			if (j == i || mark[j]) continue;
			for (int k = 0; k < Points.size(); k++) {
				if (k == i || k == j || mark[k]) continue;
				for (int x = 0; x < Points.size(); x++) {
					if (x == i || x == j || x == k || mark[x]) continue;
					if (IsPointInTriInner(Points[i], Points[j], Points[k], Points[x])) {
						
						mark[x] = 1;
					}
						
				}
			}
		}
	}
	
	std::vector<Vector2d> convexHull;
	for (int i = 0; i < Points.size(); i++) {
		if (!mark[i]) { cout << i << ","; convexHull.push_back(Vector2d(Points[i])); }
	}
//	ShowPoints(convexHull);
	cout << endl;


	
}



struct Less
{
	Vector2d s0;

	Less(Vector2d& p0) { s0 = p0; }
	bool operator()(const Vector2d& s1, const Vector2d& s2)
	{
		Vector2d AB = s1.Minus(s0);
		Vector2d AC = s2.Minus(s0);
		return AB.CrossProduct(AC)>0 ; //从小到大排序
	}
};


void CoordinateCouterClockSort(std::vector<Vector2d>& Points, Vector2d& p0) {

	std::sort(Points.begin(), Points.end(), Less(p0)); //极角排序
//先不管去重的操作了。。	
}
bool IsMakeNonLeftRun(Vector2d& top, Vector2d& nextToTop, Vector2d& newPoint) {
	Vector2d AB = top.Minus(nextToTop);
	Vector2d BC = newPoint.Minus(top);
	return !(AB.CrossProduct(BC) > 0);
}

std::deque<Vector2d> GrahamScanConvexHull(std::vector<Vector2d>& Points) {
	Vector2d p0;
	int minYIdx = 0;
	double minY = Points[0].y(); //y的范围是从0到100
	for (int i = 1; i < Points.size(); i++) {
		if (Points[i].y() < minY)
		{
			minY = Points[i].y(); minYIdx = i;
		}
		else if (Points[i].y() == minY) {
			if (Points[i].x() < Points[minYIdx].x()) minYIdx = i;

		}
	} //寻找极点
	p0 = Points[minYIdx];
	Points.erase(Points.begin() + minYIdx);
	CoordinateCouterClockSort(Points, p0); //极角排序
	//先不管去重的问题 去重可以用vector的去重算法来写
	Less cmp(p0);//极角比较大小比较器
	
	std::deque<Vector2d> out;
	out.push_front(p0);
	out.push_front(Points[0]);
	out.push_front(Points[1]);
	for (int i = 2; i < Points.size(); i++) {

		while (IsMakeNonLeftRun(out.front(), *(out.begin()+1), Points[i])) {
			out.pop_front();
		}
		out.push_front(Points[i]);
	}
	return out;
}

void showdequePoints(std::deque<Vector2d> input) {//用的副本
	std::vector<Vector2d> Points;
	while (!input.empty()) {
		Points.push_back(input.front());
		input.pop_front();
	}
	ShowPoints(Points);
}
void testGraham() {


	/*设定固定的x,y供测试用*/
	//std::vector<double> xx = { 98.0,42.0,6,34 };
	//std::vector<double> yy = { 70,60.0,81.0,8 };
	//for (int i = 0; i < 4; i++) {
	//	Points[i].x_ = xx[i];
	//	Points[i].y_ = yy[i];
	//}
	/*测试代码完毕*/
	//ShowPoints(Points);
	auto start = system_clock::now();
	auto end = system_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	vector<int> setSize;
	for (int i = 0; i < 50; i++) {
		setSize.push_back(1000 + 1000 * i);
	}
	cout << "花费时间" << endl;
	for (int i = 0; i < setSize.size(); i++) {
		std::vector<Vector2d> Points(setSize[i]);
		start = system_clock::now();
		std::deque<Vector2d> x = GrahamScanConvexHull(Points);
		end = system_clock::now();
		duration = duration_cast<microseconds>(end - start);
		cout << double(duration.count()) * microseconds::period::num / microseconds::period::den << "\t";
	}
	cout << endl;

//	showdequePoints(x);
}

void lambdaCoutDouble(double i) {
	cout << i << ',';
}

void _exchange(std::vector<double>& L, int p, int r) {
	double temp = L[p]; L[p] = L[r]; L[r] = temp;

}

void _exchange(std::vector<Vector2d>& L, int p, int r) {
	Vector2d temp;
	temp = L[p]; L[p] = L[r]; L[r] = temp;
}

int RandomPartition(std::vector<Vector2d> &L, int p, int r) {
	//std::random_device rd;
	std::default_random_engine e = std::default_random_engine();
	uniform_int_distribution<int> u(p, r);

	int i = u(e);
	_exchange(L, i, r);
	i = p - 1;
	for (int j = p; j <= r - 1; j++) {
		if (L[j].x() <= L[r].x()) {
			i++;
			_exchange(L, i, j);
		}

	}
	_exchange(L, i + 1, r);
	
	return i+1;

}


//第i大的数 而非索引
Vector2d& RandomSelect(std::vector<Vector2d> &A, int p, int r, int i) {
	if (p == r) return A[p];
	int q = RandomPartition(A, p, r);
	int k = q - p + 1;
	if (i == k) return A[q];
	else if (i < k) {
		return RandomSelect(A, p, q - 1, i);
	}
	else {
		return RandomSelect(A, q + 1, r, i - k);
	}

}

int minYPointIdx(std::vector<Vector2d>& Points) { //最下最左点

	int minYIdx = 0;
	double minY = Points[0].y(); //y的范围是从0到100
	for (int i = 1; i < Points.size(); i++) {
		if (Points[i].y() < minY)
		{
			minY = Points[i].y(); minYIdx = i;
		}
		else if (Points[i].y() == minY) {
			if (Points[i].x() < Points[minYIdx].x()) minYIdx = i;

		}
	} //寻找极点
	return minYIdx;
}


//std::vector<Vector2d> DivideConquerConvexHull(std::vector<Vector2d>& Points,int i,int j) {
//
//	if ((j - i + 1) <= 3) {
//		std::vector<Vector2d> CH(Points.begin() + i, Points.begin() + j + 1); 
//		return CH;
//	}
//	
//	RandomSelect(Points, i, j, ((j-i+1) + 1) / 2 ); //下中位数
//	
//	std::vector<Vector2d> CHl = DivideConquerConvexHull(Points, i, i+ ((j - i + 1) + 1) / 2 - 1);
//	std::vector<Vector2d> CHr = DivideConquerConvexHull(Points, i+ ((j - i + 1) + 1) / 2, j);
//
//	//从CH1和CH2中找出纵坐标最小的点Q
//	int idxL = minYPointIdx(CHl);
//	int idxR = minYPointIdx(CHr);
//	
//	bool PolarInLeft = ( CHl[idxL].y() <= CHr[idxR].y() );
//	Vector2d Polar = PolarInLeft ? CHl[idxL] : CHr[idxR];
//
//
//
//	std::vector<Vector2d> PolarSideHull;
//	std::vector<Vector2d> notPolarSideHull;
//	int PolarIdx;
//	if (PolarInLeft) {
//		PolarSideHull.assign(CHl.begin(), CHl.end());
//		notPolarSideHull.assign(CHr.begin(), CHr.end());
//		PolarIdx = idxL;
//	}
//	else {
//		PolarSideHull.assign(CHr.begin(), CHr.end());
//		notPolarSideHull.assign(CHl.begin(), CHl.end());
//		PolarIdx = idxR;
//	}
//
//
//	Less AngleLess(PolarSideHull[PolarIdx]);
//	/*******************************从PolarSideHull中得到序列S1***************************/
//	int size = PolarSideHull.size();
//	std::vector<Vector2d> S1(size-1);
//	for (int i = 0; i < PolarSideHull.size() - 1; i++) {
//		S1[i] = PolarSideHull[(PolarIdx + i + 1) % size];
//	} // 不含polar
//	if (size == 3) {
//		if (!AngleLess(S1[0], S1[1])) _exchange(S1, 0, 1);
//	}
//
//
//
//	/*******************************从PolarSideHull中得到序列S2***************************/
//	//从notPolarSideHull中得到逆时针序列S2和顺时针序列S3
//		//找到notPolarSideHull中的最大极角和最小极角点。
//	size = notPolarSideHull.size();
//	int maxAngleIdx, minAngleIdx,start_idx;
//	if (size == 3) {
//		if (!AngleLess(notPolarSideHull[0], notPolarSideHull[1])) _exchange(notPolarSideHull, 0, 1);
//		if (!AngleLess(notPolarSideHull[1], notPolarSideHull[2])) _exchange(notPolarSideHull, 1, 2);
//		if (!AngleLess(notPolarSideHull[0], notPolarSideHull[2])) _exchange(notPolarSideHull, 0, 2);
//	}
//
//	if (size % 2 != 0) {
//		maxAngleIdx = minAngleIdx = 0;
//		start_idx = 1;
//	}
//	else {
//		start_idx = 2;
//		minAngleIdx = 1; maxAngleIdx = 0;
//		if (AngleLess(notPolarSideHull[0], notPolarSideHull[1])) {
//				minAngleIdx = 0; maxAngleIdx = 1;
//		}
//	}
//	//每次取一对数比较大小，复杂度3[n/2] 这部分有毛病！！
//	for (int i = start_idx; i < size - 1; i++,i++) {
//		if (AngleLess(notPolarSideHull[i], notPolarSideHull[i+1])) {
//			if (AngleLess(notPolarSideHull[maxAngleIdx], notPolarSideHull[i+1])) { maxAngleIdx = i+1; }
//			if (AngleLess(notPolarSideHull[i], notPolarSideHull[minAngleIdx])) { minAngleIdx = i; }
//		}
//		else {
//			if ( AngleLess(notPolarSideHull[maxAngleIdx], notPolarSideHull[i])) { maxAngleIdx = i; }
//			if ( AngleLess(notPolarSideHull[i+1], notPolarSideHull[minAngleIdx])) { minAngleIdx = i+1; }
//		}
//	}
//	
//	   //由maxIdx和minIdx生成S2和S3
//	//由最小极角到最大极角的逆时针序列 
//	std::vector<Vector2d> S2;
//	for (int i = 0; (minAngleIdx + i) % size != maxAngleIdx; i++ ) {
//		S2.push_back(notPolarSideHull[(minAngleIdx + i) % size]);
//	} //S2 = Qr[y],Qr[y+1],……Qr[z-1]
//	std::vector<Vector2d> S3;
//	for (int i = 1; (minAngleIdx - i+size) % size != (maxAngleIdx-1+size) %size; i++) {
//		S3.push_back(notPolarSideHull[(minAngleIdx - i+size) % size]);
//	}//S3 = Qr[y-1],Qr[y-2],……Qr[z] 
//
//	//cout << "the true size of S3 " << S3.size() << endl;
//	//cout << "the correct size of S3 " << size - S2.size() << endl;
//	//对上述三个数组进行归并排序
//
//	std::vector<Vector2d> S23(S2.size() + S3.size());
//	merge(S3.begin(),S3.end(), S2.begin(),S2.end(),S23.begin(), AngleLess);
//	std::vector<Vector2d> S123(S23.size() + S1.size());
//	merge(S23.begin(), S23.end(), S1.begin(), S1.end(), S123.begin(), AngleLess);
//
//	//graham 排序
//	std::deque<Vector2d> out;
//	out.push_front(PolarSideHull[PolarIdx]);
//	out.push_front(S123[0]);
//	out.push_front(S123[1]);
//	for (int i = 2; i < S123.size(); i++) {
//
//		while (IsMakeNonLeftRun(out.front(), *(out.begin() + 1), S123[i])) {
//			out.pop_front();
//		}
//		out.push_front(S123[i]);
//	}
//
//	std::vector<Vector2d> vecCH(out.size());
//	std::deque<Vector2d>::reverse_iterator rit = out.rbegin();
//	for (int i = 0; i < out.size(); i++) {
//		vecCH[i] = *(rit++);
//	}
//	return vecCH;
//}
Vector2d* Conquer(vector<Vector2d>& PolarSideHull, vector<Vector2d>& notPolarSideHull, int PolarIdx,Less& AngleLess) {
	/*******************************从PolarSideHull中得到序列S1***************************/
	int size = PolarSideHull.size();
	std::vector<Vector2d> S1(size - 1);
	for (int i = 0; i < PolarSideHull.size() - 1; i++) {
		S1[i] = PolarSideHull[(PolarIdx + i + 1) % size];
	} // 不含polar
	if (size == 3) {
		if (!AngleLess(S1[0], S1[1])) _exchange(S1, 0, 1);
	}



	/*******************************从PolarSideHull中得到序列S2***************************/
	//从notPolarSideHull中得到逆时针序列S2和顺时针序列S3
		//找到notPolarSideHull中的最大极角和最小极角点。
	size = notPolarSideHull.size();
	int maxAngleIdx, minAngleIdx, start_idx;
	if (size == 3) {
		if (!AngleLess(notPolarSideHull[0], notPolarSideHull[1])) _exchange(notPolarSideHull, 0, 1);
		if (!AngleLess(notPolarSideHull[1], notPolarSideHull[2])) _exchange(notPolarSideHull, 1, 2);
		if (!AngleLess(notPolarSideHull[0], notPolarSideHull[2])) _exchange(notPolarSideHull, 0, 2);
	}

	if (size % 2 != 0) {
		maxAngleIdx = minAngleIdx = 0;
		start_idx = 1;
	}
	else {
		start_idx = 2;
		minAngleIdx = 1; maxAngleIdx = 0;
		if (AngleLess(notPolarSideHull[0], notPolarSideHull[1])) {
			minAngleIdx = 0; maxAngleIdx = 1;
		}
	}
	//每次取一对数比较大小，复杂度3[n/2] 这部分有毛病！！
	for (int i = start_idx; i < size - 1; i++, i++) {
		if (AngleLess(notPolarSideHull[i], notPolarSideHull[i + 1])) {
			if (AngleLess(notPolarSideHull[maxAngleIdx], notPolarSideHull[i + 1])) { maxAngleIdx = i + 1; }
			if (AngleLess(notPolarSideHull[i], notPolarSideHull[minAngleIdx])) { minAngleIdx = i; }
		}
		else {
			if (AngleLess(notPolarSideHull[maxAngleIdx], notPolarSideHull[i])) { maxAngleIdx = i; }
			if (AngleLess(notPolarSideHull[i + 1], notPolarSideHull[minAngleIdx])) { minAngleIdx = i + 1; }
		}
	}

	//由maxIdx和minIdx生成S2和S3
 //由最小极角到最大极角的逆时针序列 
	std::vector<Vector2d> S2;
	for (int i = 0; (minAngleIdx + i) % size != maxAngleIdx; i++) {
		S2.push_back(notPolarSideHull[(minAngleIdx + i) % size]);
	} //S2 = Qr[y],Qr[y+1],……Qr[z-1]
	std::vector<Vector2d> S3;
	for (int i = 1; (minAngleIdx - i + size) % size != (maxAngleIdx - 1 + size) % size; i++) {
		S3.push_back(notPolarSideHull[(minAngleIdx - i + size) % size]);
	}//S3 = Qr[y-1],Qr[y-2],……Qr[z] 

	//cout << "the true size of S3 " << S3.size() << endl;
	//cout << "the correct size of S3 " << size - S2.size() << endl;
	//对上述三个数组进行归并排序

	std::vector<Vector2d> S23(S2.size() + S3.size());
	merge(S3.begin(), S3.end(), S2.begin(), S2.end(), S23.begin(), AngleLess);
	Vector2d * S123 = new Vector2d[S23.size() + S1.size()];
	merge(S23.begin(), S23.end(), S1.begin(), S1.end(), S123, AngleLess);

	return S123;
}
std::vector<Vector2d> DivideConquerConvexHull(std::vector<Vector2d> & Points, int i, int j) {

	if ((j - i + 1) <= 3) {
		std::vector<Vector2d> CH(Points.begin() + i, Points.begin() + j + 1);
		return CH;
	}

	RandomSelect(Points, i, j, ((j - i + 1) + 1) / 2); //下中位数

	std::vector<Vector2d> CHl = DivideConquerConvexHull(Points, i, i + ((j - i + 1) + 1) / 2 - 1);
	std::vector<Vector2d> CHr = DivideConquerConvexHull(Points, i + ((j - i + 1) + 1) / 2, j);

	//从CH1和CH2中找出纵坐标最小的点Q
	int idxL = minYPointIdx(CHl);
	int idxR = minYPointIdx(CHr);

	bool PolarInLeft = (CHl[idxL].y() <= CHr[idxR].y());
	Vector2d Polar = PolarInLeft ? CHl[idxL] : CHr[idxR];


	Vector2d * S123;

	Less AngleLess(Polar);


	if (PolarInLeft) {
		S123 = Conquer(CHl, CHr, idxL, AngleLess);

	}
	else {
		S123 = Conquer(CHr, CHl, idxR, AngleLess);

	}







	//graham 排序
	std::deque<Vector2d> out;
	out.push_front(Polar);
	out.push_front(S123[0]);
	out.push_front(S123[1]);
	for (int i = 2; i < sizeof(S123) / sizeof(S123[0]); i++) {

		while (IsMakeNonLeftRun(out.front(), *(out.begin() + 1), S123[i])) {
			out.pop_front();
		}
		out.push_front(S123[i]);
	}

	std::vector<Vector2d> vecCH(out.size());
	std::deque<Vector2d>::reverse_iterator rit = out.rbegin();
	for (int i = 0; i < out.size(); i++) {
		vecCH[i] = *(rit++);
	}
	delete[] S123;
	return vecCH;
}
void testDivideConquerCH() {

	//std::vector<double> xx = { 16.1686,8.14013,86.4206,22.7636,18.5437,28.5542,1.54808,80.4782,47.948,8.13792 };
	//std::vector<double> yy = { 68.7322,93.8056,11.7102,79.7337,2.68421,70.001,90.4233,99.1033,16.7415,96.827 };
	//for (int i = 0; i < xx.size(); i++) {
	//	Points[i].x_ = xx[i];
	//	Points[i].y_ = yy[i];
	//}


	auto start = system_clock::now();
	auto end = system_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	vector<int> setSize;
	for (int i = 500; i < 501; i++) {
		setSize.push_back(1000 + 1000 * i);
	}
	cout << "花费时间" << endl;
	for (int i = 0; i < setSize.size(); i++) {
		std::vector<Vector2d> Points(setSize[i]);
		start = system_clock::now();
		std::vector<Vector2d> CH = DivideConquerConvexHull(Points, 0, Points.size() - 1);
		end = system_clock::now();
		duration = duration_cast<microseconds>(end - start);
		cout << double(duration.count()) * microseconds::period::num / microseconds::period::den << " ";
	}
	//ShowPoints(CH);
}
int main()
{
	//std::vector<int> dataSetSize = { 1000,2000,3000,4000 };


	//for (int i = 0; i < dataSetSize.size(); i++) {
	//	std::vector<Vector2d> dataSet(dataSetSize[i]);		
	//}

//	test();

//	testGraham();

//	std::vector<Vector2d> Points(31);
//	RandomSelect(Points, 0, Points.size()-1,(Points.size()+1)/2 ) ;

	//cout << "下面是graham的表演时间" << endl;

	//testGraham();
	cout << "下面是分治的表演时间" << endl;
	testDivideConquerCH();
	cout << "下面是brute的表演时间" << endl;

	testBruteForce();
	return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单
