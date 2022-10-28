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
        << "Size:\n\tcolumn x row (e.g. 2048x400)\n"
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
    int col = 2048, row = 400;
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
    float x = -1, y = -1, x1 = 0, x2 = col, y1 = 0, y2 = row, xc = round(col/2), yc = round(row/2), center_x, center_y;
    int no_markers = 0, flag, failure = 0, d = 0, d2 = 0;
    double last_time;
    bool certain = false;
    
    int fps = 30;
    VideoWriter video1(directory + "video.avi", CV_FOURCC('X','2','6','4'),fps, Size(col,row),false); // 'M','J','P','G' // CV_FOURCC('F','F','V','1') // CV_FOURCC('X','2','6','4') // false for grayscale (isColor)
    VideoWriter video2(directory + "tracked-video.avi", CV_FOURCC('X','2','6','4'),fps, Size(col,row));

    int imageCnt = 0;
    auto start = chrono::high_resolution_clock::now(); // timer
    while(true){
        //if (imageCnt == 168111) continue;
        if (argc >=3 && (strcmp(argv[2],"-r") == 0 || strcmp(argv[2],"--raw") == 0)){
            filename = directory + "frame-" + to_string(imageCnt) + ".raw";

            ifstream rawimagefile;
            rawimagefile.open (filename, ios::in | ios::binary);

            if (rawimagefile.fail()){
                if (imageCnt == 0){
                    cerr << "Could not read the image: " << filename << endl;
                    return 1;
                } else {
                    cout << "Finished!\n" << endl;
                    break;
                }
            }

            rawimagefile.seekg (0, rawimagefile.end);
            int length = rawimagefile.tellg();
            if (length != row*col){
                cerr << "Image size (" << length << ") is not equal to col x row (" << col << " x " << row << " = " << col*row << ")." << endl;
                if (imageCnt == 0){
                    return 1;
                } else {
                    cerr << "Could not read the image: " << filename << endl;
                    cout << "Finished!" << endl;
                    break;
                }
            }
            rawimagefile.seekg (0, rawimagefile.beg);

            rawimagefile.read (pixels, row*col);
            rawimagefile.close();
            image = Mat(row, col, CV_8U, pixels); 

        } else {
            filename = directory + "frame-" + to_string(imageCnt) + ".pgm";
            try{
                image = imread(filename, IMREAD_GRAYSCALE);
            }
            catch(cv::Exception& e){
                cout << "exception caught: " << e.what() << endl;
                return 1;
            }

            if (!image.data ){
                if (imageCnt == 0){
                    cerr << "No image data" << endl;
                    return -1;
                } else {
                    cout << "Finished!" << endl;
                    filename = directory + "frame-" + to_string(imageCnt-1) + ".pgm";
                    image = imread(filename, IMREAD_GRAYSCALE);
                    break;
                }
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
        //cout << "x1 = " << x1 << ", x2 = " << x2 << endl;
        if (no_markers>0){ // previously || x2 < col

            int dt = (int) (1000 * (time - last_time + 128)) % 128000; // time difference (in msec), time resets every 128 seconds
            d = 2 * dt; // size of the search window, max speed of rat is 2 pixels per milliseconds (~ 3 m/s)
            //cout << "d = " << d  << ", dt = " << dt << " milliseconds." << endl;
            if (certain)
                d = (d>=16) ? d : 16;
            else
                d = (d>=75) ? d : 75;
            cout << "dt = " << time - last_time << ", d = " << d << endl;
            d2 = d;

            x1 = (xc>d) ? xc-d : 0;
            y1 = (yc>d2) ? yc-d2 : 0; 
            x2 = (x1+2*d > col) ? col-x1 : 2*d; 
            y2 = (y1+2*d2 > row) ? row-y1 : 2*d2;

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
        int threshold_value = 60;
        int threshold_type = 1;
        int const max_binary_value = 255;

        threshold(image_cropped, image_bin, threshold_value, max_binary_value, threshold_type );


        // rat detection (see blob.cpp)
            // Setup SimpleBlobDetector parameters.
        SimpleBlobDetector::Params params;

        // Change thresholds
        params.minThreshold = 20;
        params.maxThreshold = 255;

        // Filter by Area.
        params.filterByArea = true;
        params.minArea = 20; // 10: 3mm markers
        params.maxArea = 200;

        // Filter by Circularity
        params.filterByCircularity = true;
        params.minCircularity = 0.2;

        // Filter by Convexity
        params.filterByConvexity = true;
        params.minConvexity = 0.2;

        // Filter by Inertia
        params.filterByInertia = true;
        params.minInertiaRatio = 0.1;


        // Storage for blobs
        vector<KeyPoint> keypoints;

        // Set up detector with params
        Ptr<SimpleBlobDetector> detector = SimpleBlobDetector::create(params);   

        int largest_idx=0;
        no_markers = 0;
        // Detect blobs

        detector->detect( image_bin, keypoints);
        if (countNonZero(Scalar::all(255) - image_bin) > 6) {// check if there is any blob (at least 6 pixels)
            no_markers = (int) keypoints.size(); // updated Dec 13, 2021 previously keypoints.size()== 1 successful (1) or not (0,2,3,..)
        } else {
            no_markers = -1;
            certain = false;
        }
        //cout << "index: " << largest_idx << ", no_markers : "<< no_markers << ", certain : "<< certain << endl;

        if (no_markers > 0){
            for (int i=0;i<no_markers;i++){
                if (keypoints[largest_idx].size<keypoints[i].size) // 2021-12-09 detecting the largest marker
                    largest_idx = i;
                keypoints[i].pt.x = x1 + keypoints[i].pt.x;
                keypoints[i].pt.y = y1 + keypoints[i].pt.y;
                cout << "x = " << keypoints[i].pt.x << ", y = " << keypoints[i].pt.y << ", r = " << keypoints[i].size << endl;
            }
            //cout << "-1! " << (keypoints[largest_idx].size < 9.5 && keypoints[largest_idx].pt.x >  800 && keypoints[largest_idx].pt.x < 1200 ) << endl;
            if (keypoints[largest_idx].size < 9 && keypoints[largest_idx].pt.x >  750 && keypoints[largest_idx].pt.x < 1300 ){ // Jan 2, 2022 size (especially in ditch)
                if (no_markers == 3){
                    center_x = (keypoints[0].pt.x + keypoints[1].pt.x + keypoints[2].pt.x)/3;
                    center_y = (keypoints[0].pt.y + keypoints[1].pt.y + keypoints[2].pt.y)/3;
                    double max_dist = 0;
                    for (int i=0;i<no_markers;i++){
                        //double dist = sqrt(4*pow(center_x - keypoints[i].pt.x,2) + pow(center_y - keypoints[i].pt.y,2));
                        double dist = abs(center_x - keypoints[i].pt.x);
                        cout << "dist = " << dist << endl;
                        if (max_dist<dist){
                            max_dist = dist;
                            largest_idx = i;
                        }
                    }
                /*} else if (keypoints.size() == 2 && d < 100) { // closer to the previous one (2 markers)
                    cout << "check here" << endl;
                    //waitKey(0);
                    if (abs(keypoints[0].pt.x - (x1+d)) < abs(keypoints[1].pt.x - (x1+d)))
                        largest_idx = 0;
                    else
                        largest_idx = 1;
                    if (keypoints[largest_idx].size<keypoints[1-largest_idx].size) // check if the largest idx is the larger of the 2 markers
                        no_markers = 0;
                */
                    // if no_markers == 3 certainty is the same as before
                    //cout << "0!" << endl;
                } else if (no_markers == 2){
                    certain = false;
                }
                    // if no_markers == 1 certainty is the same as before
            } else {
                certain = true;
                //cout << "2!" << endl;
            }
        } else{
            certain = false;
            //cout << "3!" << endl;
        }
        // write data to file
        
        drawKeypoints(image, keypoints, image_circle, Scalar(0,255,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
        if (!certain){
            center_x = 0;
            center_y = 0;
            for (int i=0;i<no_markers;i++){
                if (i!=largest_idx){
                    keypoints[i].size = 0;
                }
            }
            if (no_markers==2 || no_markers < 1){
                x = -1;
                y = -1;
            } else {
                x = (float) keypoints[largest_idx].pt.x;
                y = (float) keypoints[largest_idx].pt.y;
                xc = round(x);
                yc = round(y);
                last_time = time;
            }
            failure++;
            drawKeypoints(image_circle, keypoints, image_circle, Scalar(255,0,0), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
            //cout << "index: " << largest_idx << ", no_markers : "<< no_markers << ", certain : "<< certain << endl;

        } else {
            
            for (int i=0;i<no_markers;i++){
                if (i!=largest_idx){
                    keypoints[i].size = 0;
                } else {
                    x = (float) keypoints[largest_idx].pt.x;
                    y = (float) keypoints[largest_idx].pt.y;
                    xc = round(x);
                    yc = round(y);
                }
            }
            
            // Draw detected blobs as red circles.
            // DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures
            // the size of the circle corresponds to the size of blob
            drawKeypoints(image_circle, keypoints, image_circle, Scalar(0,0,255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS );
            last_time = time;
        }
        
        flag = (int) certain;

        data_file.write((char*) &imageCnt, sizeof(int));
        data_file.write((char*) &flag, sizeof(int));
        data_file.write((char*) &time, sizeof(double));
        data_file.write((char*) &gpio[0], 4 * sizeof(int));
        data_file.write((char*) &x, sizeof(float));
        data_file.write((char*) &y, sizeof(float));

        video1.write(image);
        video2.write(image_circle);

        // Show blobs
        /*if (!certain && no_markers>0) {
            //imshow("binary image", image_bin);
            imshow("keypoints", image_circle );
            waitKey(1);
        }
        */
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
    vector<int> RoiOffset = readRoiOffset(image);

    log_file << "ROI Offset = (" << *(RoiOffset.begin()) << "," << *(RoiOffset.begin()+1) << ")"<< endl;
    log_file << "Frame Size = " << col << "x" << row << endl;
    double success = 100 * (1 - (double) failure / imageCnt);

    auto finish = chrono::high_resolution_clock::now();
    double duration = chrono::duration_cast<chrono::nanoseconds>(finish-start).count()*1e-9;
    cout << "Time taken by program is " << fixed  << duration << setprecision(6);
    cout << " sec for " << imageCnt << " frame(s).\n" << fixed  << setprecision(3) << duration/imageCnt*1e3;
    cout << " milliseconds per frame." << endl;
    cout << "Success rate: " << success << "%" << endl;
    log_file << "Success rate: " << success << "%" << endl;
    log_file << "Processing speed: " << setprecision(3) << duration/imageCnt*1e3 << " milliseconds per frame." << endl;
    log_file.close();
    data_file.close();
    cout << "Tracking data saved in " << data_filename << endl;

    return 0;
}