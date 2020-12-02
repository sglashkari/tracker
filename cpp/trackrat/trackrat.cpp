#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <vector>
#include <iostream>
#include <fstream> // file operations
#include <string>
#include <tuple>
#include <chrono> // time, chrono
#include <sstream>  // std::ostringstream
#include <iomanip> // std::setw
#include <algorithm>    // std::max

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
        printf("usage: trackrat <Image_Path> [Option]\n");
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

    string directory = argv[1], filename;
    ofstream data_file;
    string data_filename = directory + "tracking.dat";
    data_file.open (data_filename, ios::out | ios::binary);

    int col = 1600, row = 580;
    //int col = 1536, row = 740;
    //int col = 1440, row = 708;
    char pixels[row*col];
    ifstream rawimagefile;

    Mat image = Mat(row, col, CV_8U, pixels), image_bin, image_circle, image_cropped, image_RGB;
    float x = -1, y = -1, x1 = 0, x2 = col, y1 = 0, y2 = row;
    int flag = 0;
    
    int fps = 240;
    VideoWriter video1(directory + "video.avi", CV_FOURCC('X','2','6','4'),fps, Size(col,row)); // 'M','J','P','G' // CV_FOURCC('F','F','V','1') // CV_FOURCC('X','2','6','4')
    VideoWriter video2(directory + "tracked-video.avi", CV_FOURCC('X','2','6','4'),fps, Size(col,row));

    int imageCnt = 0;
    auto start = chrono::high_resolution_clock::now(); // timer
    while(true){

        if (argc ==3 && (strcmp(argv[2],"-r") == 0 || strcmp(argv[2],"--raw") == 0)){
            filename = directory + "frame-" + to_string(imageCnt) + ".raw";

            ifstream rawimagefile;
            rawimagefile.open (filename, ios::in | ios::binary);

            if (rawimagefile.fail()){
                if (imageCnt == 0){
                    cout << "Could not read the image: " << filename << endl;
                    return -1;
                } else {
                    auto finish = chrono::high_resolution_clock::now();
                    double duration = chrono::duration_cast<chrono::nanoseconds>(finish-start).count()*1e-9;
                    cout << "Time taken by program is : " << fixed  << duration << setprecision(6);
                    cout << " sec for " << imageCnt << " frame(s).\n" << fixed  << setprecision(3) << duration/imageCnt*1e3;
                    cout << " milliseconds per frame." << endl;
                    data_file.close();
                    cout << "Tracking data saved in " << data_filename << endl;
                    return 0;
                }
            }

            rawimagefile.seekg (0, rawimagefile.end);
            int length = rawimagefile.tellg();
            if (length != row*col){
                cout << "Image size (" << length << ") is not equal to row x col (" << row*col << ")." << endl;
                return -1;
            }
            rawimagefile.seekg (0, rawimagefile.beg);

            rawimagefile.read (pixels, row*col);
            rawimagefile.close();
            image = Mat(row, col, CV_8U, pixels); 

        } else {
            filename = directory + "frame-" + to_string(imageCnt) + ".pgm";
            image = imread(filename, IMREAD_GRAYSCALE);

            if (!image.data ){
                if (imageCnt == 0){
                    cout << "No image data" << endl;
                    return -1;
                } else {
                    auto finish = chrono::high_resolution_clock::now();
                    double duration = chrono::duration_cast<chrono::nanoseconds>(finish-start).count()*1e-9;
                    cout << "Time taken by program is : " << fixed  << duration << setprecision(6);
                    cout << " sec for " << imageCnt << " frame(s).\n" << fixed  << setprecision(3) << duration/imageCnt*1e3;
                    cout << " milliseconds per frame." << endl;
                    data_file.close();
                    cout << "Tracking data saved in " << data_filename << endl;
                    return 0;
                }
            }

            if(image.empty()){
                std::cout << "Could not read the image: " << filename << std::endl;
                return 1;
            }
        } 

        auto [time, gpio] = readimageinfo (image);

        cout << filename << endl;
        // printing out the values
        printf("time = %3.6f, gpio = ",time);
        for (int j = 0; j < 4; j++) 
    		cout << gpio[j];

    	cout << endl;

        //crop image to the region of interest (ROI)
        if (flag == 1){
            x1 = (x>30) ? x-30 : 0;
            y1 = (y>30) ? y-30 : 0; 
            x2 = (x+30> col) ? col-x1 : 60; 
            y2 = (y+30> row) ? row-y1 : 60;

            Rect roi(x1, y1, x2, y2);
            image_cropped = image(roi);
        } else {
            x1 = 0;
            y1 = 0;
            x2 = col;
            y2 = row;
            image_cropped = image;
        }

        // make binary image (see trackled.cpp)
        int threshold_value = 200;
        int threshold_type = 1;
        int const max_binary_value = 255;

        threshold(image_cropped, image_bin, threshold_value, max_binary_value, threshold_type );


        // rat detection (see blob.cpp)
            // Setup SimpleBlobDetector parameters.
        SimpleBlobDetector::Params params;

        // Change thresholds
        params.minThreshold = 10;
        params.maxThreshold = 200;

        // Filter by Area.
        params.filterByArea = true;
        params.minArea = 150;
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
        if (countNonZero(Scalar::all(255) - image_bin) > 0){ // check if there is any blob
            detector->detect( image_bin, keypoints);

            for (int i=0;i<keypoints.size();i++)
                cout << "x = " << keypoints[i].pt.x << ", y = " << keypoints[i].pt.y << ", r = " << keypoints[i].size << endl;

            flag = (int) (keypoints.size()==1); // successful (1) or not (0,2,3,..)
        } else {
            flag = 0;
        }


        // write data to file 
        if (flag == 1){
            keypoints[0].pt.x = x1 + keypoints[0].pt.x;
            keypoints[0].pt.y = y1 + keypoints[0].pt.y;
            x = (float) keypoints[0].pt.x;
            y = (float) keypoints[0].pt.y;
        } else {
            x = -1;
            y = -1;
        }

        data_file.write((char*) &imageCnt, sizeof(int));
        data_file.write((char*) &flag, sizeof(int));
        data_file.write((char*) &time, sizeof(double));
        data_file.write((char*) &gpio[0], 4 * sizeof(int));
        data_file.write((char*) &x, sizeof(float));
        data_file.write((char*) &y, sizeof(float));

        // Draw detected blobs as red circles.
        // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures
        // the size of the circle corresponds to the size of blob
        drawKeypoints(image, keypoints, image_circle, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );

        cvtColor(image, image_RGB, CV_GRAY2RGB);

        video1.write(image_RGB);
        video2.write(image_circle);

        // Show blobs
        //imshow("keypoints", image_circle );
        //int k = waitKey(1);

        imageCnt++;
    }

    video1.release();
    video2.release();

    return 0;
}