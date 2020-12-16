// C++ program to convert a decimal 
// number to binary number 

#include <iostream>
//#include <algorithm>    // std::reverse
#include <iomanip>      // std::setprecision
#include<vector>
using namespace std; 

// function to convert decimal to binary 
vector<int> decToBinary(int n) 
{ 
	// array to store binary number 
	vector<int> binaryNum, binaryNumReversed;

	for (int i=0; i<8; i++){
		// storing remainder in binary array 
		binaryNumReversed.push_back(n % 2);
		n = n / 2; 
	}

	for (int i=7; i>=0; i--)
		binaryNum.push_back(binaryNumReversed[i]);
	
	return binaryNum;
} 

// Driver program to test above function 
int main() 
{ 
	int n;
	//int* binaryNum;

	cout << "Enter a number: ";
	cin >> n;

	vector<int> binaryNum = decToBinary(n);

	//binaryNum.insert( binaryNum.end(), binaryNum.begin(), binaryNum.end() );

	// printing binary array in reverse order 
	for (int j = 0; j < 8; j++) 
		cout << binaryNum[j];

	cout << endl;

	int sec = 0;
	for (int j = 0; j < 7; j++) 
		sec = 2 * sec + binaryNum[j];


	int msec = 0;
	for (int j = 7; j < 8; j++) // change it to 20
		msec = 2 * msec + binaryNum[j];

	double time = sec + 125e-6 * msec;
	setprecision(9);

	cout << sec << " " << 125e-6 * msec << " " << time << " " << 125e-6 << endl;
	printf("%3.6f\n",100*time);


	return 0; 
} 
