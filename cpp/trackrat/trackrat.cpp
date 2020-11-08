#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <vector>
#include <iostream>
#include <string>
#include <tuple>
#include <chrono> // time, chrono
#include <sstream>  // std::ostringstream
#include <iomanip> // std::setw

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

    int imageCnt = 0;
    auto start = chrono::high_resolution_clock::now();
    while(true){ 

        string filename = argv[1];

        filename = filename + to_string(imageCnt++) + ".pgm";
        Mat image = imread(filename, IMREAD_GRAYSCALE), image_bin;
        
        if (!image.data ){
            //printf("No image data \n");
            auto finish = chrono::high_resolution_clock::now();
            double duration = chrono::duration_cast<chrono::nanoseconds>(finish-start).count()*1e-9;
            cout << "Time taken by program is : " << fixed  << duration << setprecision(6);
            cout << " sec for " << --imageCnt << " frame(s).\n" << fixed  << setprecision(0) << duration/imageCnt*1e6;
            cout << " microseconds per frame." << endl;
            return -1;
        }

        //! [empty]
        if(image.empty()){
            std::cout << "Could not read the image: " << argv[1] << std::endl;
            return 1;
        }

        auto [time, gpio] = readimageinfo (image);

        cout << filename << endl;
        // printing out the values
        printf("time = %3.6f, gpio = ",time);
        for (int j = 0; j < 4; j++) 
    		cout << gpio[j];

    	cout << endl;

        // make binary image (see trackled.cpp)
        int threshold_value = 200;
        int threshold_type = 1;
        int const max_binary_value = 255;

        threshold(image, image_bin, threshold_value, max_binary_value, threshold_type );

        // rat detection (see blob.cpp)
            // Setup SimpleBlobDetector parameters.
        SimpleBlobDetector::Params params;

        // Change thresholds
        params.minThreshold = 10;
        params.maxThreshold = 200;

        // Filter by Area.
        params.filterByArea = true;
        params.minArea = 100;
        params.maxArea = 300;

        // Filter by Circularity
        params.filterByCircularity = true;
        params.minCircularity = 0.1;

        // Filter by Convexity
        params.filterByConvexity = false;

        // Filter by Inertia
        params.filterByInertia = true;
        params.minInertiaRatio = 0.01;


        // Storage for blobs
        vector<KeyPoint> keypoints;

        // Set up detector with params
        Ptr<SimpleBlobDetector> detector = SimpleBlobDetector::create(params);   

        // Detect blobs
        detector->detect( image_bin, keypoints);

        // Draw detected blobs as red circles.
        // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures
        // the size of the circle corresponds to the size of blob

        //Mat im_with_keypoints;
        //drawKeypoints( image_bin, keypoints, im_with_keypoints, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );

        // Show blobs
        //imshow("keypoints", im_with_keypoints );
        //int k = waitKey(0);

        for (int i=0;i<keypoints.size();i++)
            cout << "x = " << keypoints[i].pt.x << ", y = " << keypoints[i].pt.y << ", r = " << keypoints[i].size << endl;

    }

    return 0;
}