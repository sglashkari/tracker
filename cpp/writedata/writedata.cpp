
// basic file operations
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main (int argc, char** argv) {
	
	string directory = argv[1];
    string data_filename = directory + "/tracking.dat";
    data_filename.erase(data_filename.find("//"),1);

    /*
	ofstream data_file1;
	data_file1.open (data_filename, ios::out | ios::binary);

	vector<int> gpio;
  	for(int i = 0; i < 5; i++){
  		double time = (double) (4-i)*1.1;
        vector<int> gpio{ 1*i, 2*i, 3*i, 4*i};
        data_file1.write((char*) &time, sizeof(time));
        data_file1.write((char*) &gpio[0], gpio.size() * sizeof(int));
        float x = 0.252 * i;
        float y = 0.329 * i;
        data_file1.write((char*) &x, sizeof(x));
        data_file1.write((char*) &y, sizeof(y));
  		cout << gpio.size() << " " << time << ", " << gpio[0] << ", " << gpio[1]<< ", " << gpio[2]<< ", " << gpio[3] << endl;
  	}
	
  	data_file1.close();
	*/

  	ifstream data_file;
  	data_file.open (data_filename, ios::in | ios::binary);

  	double time;
  	int pin;
  	
  	while (true){

  		data_file.read((char*) &time, sizeof(double));
  		if (data_file.eof()) break;
  		cout << "t = " << time;
  		for(int j = 0; j < 4; j++){
  			data_file.read((char*) &pin, sizeof(int));
  			cout << ", " << pin;
  		}
  		float x , y;
  		data_file.read((char*) &x, sizeof(x));
  		data_file.read((char*) &y, sizeof(float));

  		cout << ", x = " << x << ", y = " << y << endl;
  	}

  	data_file.close();


  	return 0;
}