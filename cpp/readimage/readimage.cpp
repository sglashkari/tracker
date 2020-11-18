#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <vector>
#include <iostream>
#include <string>
#include <tuple>
#include <ctime> // time calculation
#include <sstream>  // std::ostringstream
#include <iomanip> // std::setw
#include <fstream> // file operations

using namespace cv;
using namespace std;
//! [includes]

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

// function to calculate timestmp and gpio
tuple<double, vector<int>> readimageinfo (Mat image)
{
	int intensity;
	vector<int> binaryNum, binaryNum32;
    
    for (int i=0; i<4; i++){
    	intensity = image.at<uchar>(0, i);
    	binaryNum = decToBinary(intensity);
    	binaryNum32.insert( binaryNum32.end(), binaryNum.begin(), binaryNum.end());
    }

	int sec = 0;
	for (int j = 0; j < 7; j++)
		sec = 2 * sec + binaryNum32[j];

	int msec = 0;
	for (int j = 7; j < 20; j++)
		msec = 2 * msec + binaryNum32[j];

	double time = sec + 125e-6 * msec;

	// GPIO
    int intensity4 = image.at<uchar>(0, 4);
	vector<int> gpio = decToBinary(intensity4);

	return {time, gpio};
}

int main(int argc, char** argv)
{

    if ( argc != 2 && argc != 3)
    {
        printf("usage: readimage <Image_Path> [Option]\n");
        return -1;
    }

    if ( argc ==3 && strcmp(argv[2],"-r") != 0 && strcmp(argv[2],"--raw") != 0)
    {
        std::cerr << "Usage: " << argv[0] << " " << argv[1] << " [Option]"
            << "Options:\n"
              << "\t-r,--raw RAW File\tSpecify the image file type"
              << std::endl;
        return 1;
    }    


    string path = argv[1], filename;
    int row = 1440, col = 708;

    int imageCnt = 0;
    clock_t t_start = clock();
    while(true){ 

        //string s_number = to_string(imageCnt++);
        //string str_number = string(4 - s_number.length(), '0') + s_number;

        Mat image;

        if (argc ==3 && (strcmp(argv[2],"-r") == 0 || strcmp(argv[2],"--raw") == 0)){
            filename = path + "frame-" + to_string(imageCnt) + ".raw";

            ifstream rawimagefile;
            rawimagefile.open (filename, ios::in | ios::binary);

            if (rawimagefile.fail()){
                if (imageCnt == 0){
                    cout << "Could not read the image: " << filename << endl;
                    return -1;
                } else {
                    clock_t t_end = clock();
                    double time_taken = double(t_end - t_start) / double(CLOCKS_PER_SEC);
                    cout << "Time taken by program is : " << fixed  << time_taken << setprecision(5);
                    cout << " sec for " << imageCnt << " frames.\n" << fixed  << setprecision(0) << time_taken/imageCnt*1e6;
                    cout << " microseconds per frame." << endl; 
                    return 0;
                }
            }

            char pixels[row*col];
            //pixels = new char [size];
            rawimagefile.read (pixels, row*col);

            
            //fread(pixels,row*col,1,f);
            rawimagefile.close();

            image = Mat(row, col, CV_8UC1, pixels);

        } else {
            filename = path + "frame-" + to_string(imageCnt) + ".pgm";
            image = imread(filename, IMREAD_GRAYSCALE);
        }

        cout << filename << endl;
        
        if (!image.data ){
            clock_t t_end = clock();
            double time_taken = double(t_end - t_start) / double(CLOCKS_PER_SEC);
            cout << "Time taken by program is : " << fixed  << time_taken << setprecision(5);
            cout << " sec for " << imageCnt << " frames.\n" << fixed  << setprecision(0) << time_taken/imageCnt*1e6;
            cout << " microseconds per frame." << endl; 
            return 0;
        }

        //! [empty]
        if(image.empty()){
            cout << "Could not read the image: " << argv[1] << endl;
            return 1;
        }

        auto [time, gpio] = readimageinfo (image);

        // printing out the values
        printf("%3.6f\n",time);
        for (int j = 0; j < 4; j++) 
    		cout << gpio[j];

        imageCnt++;
    	cout << endl;
    }

    return 0;
}