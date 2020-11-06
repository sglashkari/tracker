#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include<vector>
#include <iostream>
#include <string>
#include <tuple>

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
    if ( argc != 2 )
    {
        printf("usage: DisplayImage.out <Image_Path>\n");
        return -1;
    }
    Mat image = imread( argv[1], IMREAD_GRAYSCALE); //IMREAD_COLOR);
    
    auto [time, gpio] = readimageinfo (image);


    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }


    //! [empty]
    if(image.empty())
    {
        std::cout << "Could not read the image: " << argv[1] << std::endl;
        return 1;
    }
    //! [empty]

    // printing out the values
    cout << __DATE__ << " " << __TIME__ << endl;
    printf("%3.6f\n",time);
    for (int j = 0; j < 4; j++) 
		cout << gpio[j];

	cout << endl;

    //! [imshow]
    //namedWindow("Display Image", WINDOW_AUTOSIZE );
    imshow("Display window: by Shahin", image);
    
    int k = waitKey(0); // Wait for a keystroke in the window
    //! [imshow]

    //! [imsave]

    string str = argv[1];
    if(k == 's')
    {
        imwrite(str.append(".png"), image);
    }
    //! [imsave]

    return 0;
}