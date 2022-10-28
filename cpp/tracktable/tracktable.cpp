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
#include <algorithm>    // std::max, std::sort

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

// function to calculate offset X and Y
vector<int> readRoiOffset (Mat image){
    // ROI Position
    int intensity;
    vector<int> binaryNum, binaryNum32, RoiOffset;

    binaryNum32.clear();
    for (int i=8; i<12; i++){
        intensity = image.at<uchar>(0, i);
        binaryNum = decToBinary(intensity);
        binaryNum32.insert( binaryNum32.end(), binaryNum.begin(), binaryNum.end());
    }

    int OffsetX = 0;
    for (int j = 0; j < 16; j++)
        OffsetX = 2 * OffsetX + binaryNum32[j];

    int OffsetY = 0;
    for (int j = 16; j < 32; j++)
        OffsetY = 2 * OffsetY + binaryNum32[j];

    cout << "Offset X = " << OffsetX << ", Offset Y = " << OffsetY << endl;
    RoiOffset.push_back(OffsetX);
    RoiOffset.push_back(OffsetY);

    return RoiOffset;
}

// function to calculate timestmp and gpio
tuple<double, vector<int>> readimageinfo (Mat image)
{
	int intensity;
	vector<int> binaryNum, binaryNum32;

    // Time
    for (int i=0; i<4; i++){
    	intensity = image.at<uchar>(0, i);
    	binaryNum = decToBinary(intensity);
    	binaryNum32.insert( binaryNum32.end(), binaryNum.begin(), binaryNum.end());
    }

    int nSecond = 0;
    for (int j = 0; j < 7; j++)
      nSecond = 2 * nSecond + binaryNum32[j];

  int nCycleCount = 0;
  for (int j = 7; j < 20; j++)
      nCycleCount = 2 * nCycleCount + binaryNum32[j];

    int nCycleOffset = 0; // https://www.flir.com/support-center/iis/machine-vision/knowledge-base/imaging-products-timestamping-and-different-timestamp-mechanisms/
    for (int j = 20; j < 32; j++)
        nCycleOffset = 2 * nCycleOffset + binaryNum32[j];

    double time = (double) nSecond + 1.25e-4 * (double) nCycleCount + 4.069e-8 * (double) nCycleOffset;

	// GPIO
    int intensity4 = image.at<uchar>(0, 4);
    vector<int> gpio = decToBinary(intensity4);



    return {time, gpio};
}

int main(int argc, char** argv)
{

    if ( argc < 2 || argc > 4)
    {
        cerr << "Usage: trackrat <Frames_Directory_Path> [File Type] [Size]\ntrackrat ./ -r 1200x600" << endl;
        return 1;
    }

    if ( argc >=3 && strcmp(argv[2],"-p") != 0 && strcmp(argv[2],"-r") != 0 && strcmp(argv[2],"--raw") != 0)
    {
        cerr << "Usage: " << argv[0] << " " << argv[1] << " [File Type] [Size]\n"
        << "File Type:\n\t-r,--raw RAW File\tSpecify the image file type \n"
        << "\t(Default): pgm file type\n"
        << "Size:\n\tcolumn x row (e.g. 1200x600)\n"
        << endl;
        return 1;
    }

    string directory = argv[1], filename;
    ofstream data_file, log_file;
    string data_filename = directory + "tracking.dat", log_filename = directory + "log.txt";
    data_file.open (data_filename, ios::out | ios::binary);

    //int col = 1600, row = 580;
    //int col = 1536, row = 740;
    //int col = 1440, row = 708;
    int col = 1200, row = 700;
    //int col = 1920, row = 600;

    if (argc ==4){
        string str =argv[3];
        size_t pos = str.find("x");
        col = stoi(str.substr (0,pos));
        row = stoi(str.substr (pos+1));
    }

    char pixels[row*col];
    ifstream rawimagefile;

    Mat image = Mat(row, col, CV_8U, pixels), image_bin, image_circle, image_cropped;
    float xr = -1, xl = -1, xrc = col;
    int flag = 0, failure = 0;
    
    int fps = 30;
    VideoWriter video1(directory + "video.avi", CV_FOURCC('X','2','6','4'),fps, Size(col,row),false); // 'M','J','P','G' // CV_FOURCC('F','F','V','1') // CV_FOURCC('X','2','6','4') // false for grayscale (isColor)
    VideoWriter video2(directory + "tracked-video.avi", CV_FOURCC('X','2','6','4'),fps, Size(col,row));

    int imageCnt = 0;
    auto start = chrono::high_resolution_clock::now(); // timer
    while(true){
        if (argc >=3 && (strcmp(argv[2],"-r") == 0 || strcmp(argv[2],"--raw") == 0)){
            filename = directory + "frame-" + to_string(imageCnt) + ".raw";

            ifstream rawimagefile;
            rawimagefile.open (filename, ios::in | ios::binary);

            if (rawimagefile.fail()){
                cerr << "Could not read the image: " << filename << endl;
                if (imageCnt == 0){
                    return 1;
                } else break;
            }

            rawimagefile.seekg (0, rawimagefile.end);
            int length = rawimagefile.tellg();
            if (length != row*col){
                cerr << "Image size (" << length << ") is not equal to col x row (" << col << " x " << row << " = " << col*row << ")." << endl;
                if (imageCnt == 0){
                    return 1;
                } else {
                    cerr << "Could not read the image: " << filename << endl;
                    cerr << "\nProgram aborted!\n" << endl;
                    break;
                }
            }
            rawimagefile.seekg (0, rawimagefile.beg);

            rawimagefile.read (pixels, row*col);
            rawimagefile.close();
            image = Mat(row, col, CV_8U, pixels); 

        } else {
            filename = directory + "frame-" + to_string(imageCnt) + ".pgm";
            try {
                image = imread(filename, IMREAD_GRAYSCALE);
            } catch(...) {
                cerr << "Could not read the image: " << filename << endl;
                if (imageCnt == 0){
                    return 1;
                } else break;
            }

            if (!image.data ){
                if (imageCnt == 0){
                    cerr << "No image data" << endl;
                    return 1;
                } else break;
            }

            if(image.empty()){
                cerr << "Could not read the image: " << filename << endl;
                return 1;
            }
        } 

        auto [time, gpio] = readimageinfo (image);
        if (imageCnt == 0){
            vector<int> RoiOffset = readRoiOffset(image);
        }

        cout << filename << endl;
        // printing out the values
        printf("time = %3.6f, gpio = ",time);
        for (int j = 0; j < 4; j++) 
          cout << gpio[j];

      cout << endl;

        //crop image to the region of interest (ROI): dynamic size
      Rect roi(0, 0, col, 300);
      image_cropped = image(roi);

        // make binary image (see trackled.cpp)
      int threshold_value = 50;
      int threshold_type = 0;
      int const max_binary_value = 255;

      threshold(image_cropped, image_bin, threshold_value, max_binary_value, threshold_type );


        // rat detection (see blob.cpp)
        // Setup SimpleBlobDetector parameters.
      SimpleBlobDetector::Params params;

        // Change thresholds
      params.minThreshold = 30;
      params.maxThreshold = 100;

        // Filter by Area.
      params.filterByArea = true;
      params.minArea = 30;
      params.maxArea = 70;

        // Filter by Circularity
      params.filterByCircularity = true;
      params.minCircularity = 0.75;

        // Filter by Convexity
      params.filterByConvexity = true;
      params.minConvexity = 0.6;

        // Filter by Inertia
      params.filterByInertia = true;
      params.minInertiaRatio = 0.5;


        // Storage for blobs
      vector<KeyPoint> keypoints;

        // Set up detector with params
      Ptr<SimpleBlobDetector> detector = SimpleBlobDetector::create(params);   

      int left_nail_idx=0, right_nail_idx=1; // 1 & 2 for most of the days for 980, 0 & 1 for 1024, 1055
        // Detect blobs
        if (countNonZero(Scalar::all(255) - image_bin) > 0){ // check if there is any blob
            detector->detect( image_bin, keypoints);
            std::sort(keypoints.begin(),keypoints.end(),    // sorting
                [](const KeyPoint &first, const KeyPoint &second){
                    return first.pt.x < second.pt.x;
                }); 


            for (int i=0;i<keypoints.size();i++)
                cout << "x = " << keypoints[i].pt.x << ", y = " << keypoints[i].pt.y << ", r = " << keypoints[i].size << endl;
            cout << "right_nail_idx: " << right_nail_idx << " left_nail_idx: " << left_nail_idx<< endl;
            
            for (int i=right_nail_idx+1;i<keypoints.size();i++){
                //if (keypoints[right_nail_idx].pt.y > (keypoints[0].pt.y + 10) || keypoints[right_nail_idx].pt.y < (keypoints[0].pt.y - 30) ) //|| keypoints[right_nail_idx].pt.x < (xr - 50))  // 2021-12-09 detecting the largest marker
                if (keypoints[right_nail_idx].pt.y < (keypoints[left_nail_idx].pt.y - 20) || keypoints[right_nail_idx].pt.y > (keypoints[left_nail_idx].pt.y + 20))
                    right_nail_idx = i;
            }

            if (keypoints.size()<2)
                flag = 0;
            //else if (keypoints[right_nail_idx].pt.x > (xrc + 50))
            //    flag = 0;
            else if (keypoints[right_nail_idx].pt.y < (keypoints[left_nail_idx].pt.y - 20) || keypoints[right_nail_idx].pt.y > (keypoints[left_nail_idx].pt.y + 20))
                flag = 0;
            else
                flag = (int) keypoints.size(); // updated Dec 13, 2021 previously keypoints.size()== 1 successful (1) or not (0,2,3,..)

            if (keypoints.size()>=2 && keypoints[left_nail_idx].pt.x > col/2){
                flag = 0;
            }


        } else {
            flag = 0;
        }

        drawKeypoints(image, keypoints, image_circle, Scalar(0,255,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );

        // write data to file 
        if (!flag){
            failure++;
        } else {
            for (int i=0;i<keypoints.size();i++){
                if (i!=right_nail_idx && i!=left_nail_idx){
                    keypoints[i].size = 0;
                }
            }
            xr = (float) keypoints[right_nail_idx].pt.x;
            xl = (float) keypoints[left_nail_idx].pt.x;

            if (xr > col/2) 
                xrc = xr;
            else
                xrc = col;
            

            // Draw detected blobs as red circles.
            // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures
            // the size of the circle corresponds to the size of blob
            drawKeypoints(image_circle, keypoints, image_circle, Scalar(50,205,50), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );

            keypoints[left_nail_idx].size = 0;
            drawKeypoints(image_circle, keypoints, image_circle, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
        }
        

        cout << "right = " << xr << ", left = " << xl << endl;
        data_file.write((char*) &imageCnt, sizeof(int));
        data_file.write((char*) &flag, sizeof(int));
        data_file.write((char*) &time, sizeof(double));
        data_file.write((char*) &gpio[0], 4 * sizeof(int));
        data_file.write((char*) &xr, sizeof(float));
        data_file.write((char*) &xl, sizeof(float));

        video1.write(image);
        video2.write(image_circle);
        // Show blobs

        if (imageCnt % 10 == 0){
            //imshow("binary image", image_bin);
            //imshow("cropped image", image_cropped);
            imshow("keypoints", image_circle );
            int k = waitKey(1);
            if (k == 27) {
                cout << "\n User aborted the program!\n" << endl;
                return 0; // Esc key pressed
            }
        }
        imageCnt++;
    }

    video1.release();
    video2.release();

    log_file.open(log_filename);
    //vector<int> RoiOffset = readRoiOffset(image);
    //log_file << "ROI Offset = (" << *(RoiOffset.begin()) << "," << *(RoiOffset.begin()+1) << ")"<< endl;
    log_file << "Frame Size = " << col << "x" << row << endl;
    double success = 100 * (1 - (double) failure / imageCnt);

    auto finish = chrono::high_resolution_clock::now();
    double duration = chrono::duration_cast<chrono::nanoseconds>(finish-start).count()*1e-9;
    cout << endl;
    cout << "Time taken by program is : " << fixed  << duration << setprecision(6);
    cout << " sec for " << imageCnt << " frame(s).\n" << fixed  << setprecision(3) << duration/imageCnt*1e3;
    cout << " milliseconds per frame." << endl;
    cout << "Success rate was " << success << "%" << endl;
    log_file << "Success rate was " << success << "%" << endl;
    log_file << "Processing speed: " << setprecision(3) << duration/imageCnt*1e3 << " milliseconds per frame." << endl;
    log_file.close();
    data_file.close();
    cout << "Tracking data saved in " << data_filename << endl;

    return 0;
}