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

int main(int argc, char** argv){

    if ( argc != 2 && argc != 3){
        printf("usage: readimage <Image_Path> [Option]\n");
        return -1;
    }

    if ( argc ==3 && strcmp(argv[2],"-r") != 0 && strcmp(argv[2],"--raw") != 0 && strcmp(argv[2],"-v") != 0){
        std::cerr << "Usage: " << argv[0] << " " << argv[1] << " [Option]"
        << "Options:\n"
        << "\t-r,--raw RAW File\tSpecify the image file type as raw\n"
        << "\t-v\tSpecify the file type as video\n"
        << std::endl;
        return 1;
    }

    // Create a VideoCapture object and open the input file
    // If the input is the web camera, pass 0 instead of the video file name
    VideoCapture cap(argv[1]); 

    // Check if camera opened successfully
    if(!cap.isOpened() && argc == 3 && strcmp(argv[2],"-v") == 0){
        cout << "Error opening video stream or file" << endl;
        return -1;
    }

    //int col = 1536, row = 740;
    int col = 1440, row = 708;
    char pixels[row*col];
    string directory = argv[1], filename;
    ifstream rawimagefile;

    Mat image = Mat(row, col, CV_8UC1, pixels), imageRGB;

    size_t found = directory.find_last_of("/\\");
    directory = directory.substr(0,found+1);
    std::cout << " path: " << directory << '\n';

    //string directory = text.substr(0, directory.find("what "));
    string videoname = directory + "output.avi";

    int fps = 1;
    VideoWriter video(videoname,CV_FOURCC('M','J','P','G'),fps, Size(col,row));
    cout << "Video file created!" << endl;

    int imageCnt = 0;
    clock_t t_start = clock();

    while(true){ 

        if ((argc ==3) && (strcmp(argv[2],"-r") == 0 || strcmp(argv[2],"--raw") == 0)){
            filename = directory + "frame-" + to_string(imageCnt) + ".raw";
            rawimagefile.open (filename, ios::in | ios::binary);
            cout << filename << endl;
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

            rawimagefile.seekg (0, rawimagefile.end);
            int length = rawimagefile.tellg();
            if (length != row*col){
                cout << "Image size (" << length << ") is not equal to row x col (" << row*col << ")." << endl;
                return -1;
            }
            rawimagefile.seekg (0, rawimagefile.beg);

            rawimagefile.read (pixels, row*col);
            rawimagefile.close();
            image = Mat(row, col, CV_8UC1, pixels);   


        } else {
            if (argc ==3 && strcmp(argv[2],"-v") == 0) {
                // Capture frame-by-frame
                cap >> image;
            } else {
                filename = directory + "frame-" + to_string(imageCnt) + ".pgm";
                image = imread(filename, IMREAD_GRAYSCALE);
            }

            if (!image.data ){
                if (imageCnt == 0){
                    cout << "No image data" << endl;
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

            if(image.empty()){
                cout << "Could not read the image: " << filename << endl;
                return 1;
            }

        }

        auto [time, gpio] = readimageinfo (image);

        // printing out the values
        printf("%3.6f\n",time);
        for (int j = 0; j < 4; j++) 
          cout << gpio[j];

        cout << endl;

        if (strcmp(argv[2],"-v") != 0){
            //cvtColor(imageBGR, image, CV_GRAY2BGR);
            cvtColor(image, imageRGB, CV_GRAY2RGB);
        } else {
            imageRGB = image;
        }

        video.write(imageRGB);

        imshow("Display window: by Shahin", image);
        int k = waitKey(0); // Wait for a keystroke in the window
        imageCnt++;
    }

    // When everything done, release the video capture object
    cap.release();
    video.release();

    // Closes all the frames
    destroyAllWindows();

    return 0;
}