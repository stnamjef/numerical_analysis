#include <iostream>
#include "opencv2/opencv.hpp"
using namespace std;

int main()
{
	cv::Mat in = cv::imread("sunflower.bmp");
	cv::Mat mask(3, 3, CV_32F, cv::Scalar(0));
	mask.at<float>(0, 0) = -1.0;
	mask.at<float>(2, 2) = 1.0;

	int pad = 4;

	cv::Mat out(in.size(), CV_32FC1);
	for (int c = 0; c < 3; c++) {
		for (int i = 0; i < out.rows; i++) {
			for (int j = 0; j < out.cols; j++) {
				out.at<float>(i, j) = 0;
				for (int y = 0; y < mask.rows; y++) {
					for (int x = 0; x < mask.cols; x++) {
						int ii = i + y - pad;
						int jj = j + x - pad;
						if (ii >= 0 && ii < in.rows &&
							jj >= 0 && jj < in.cols) {
							out.at<float>(i, j) += in.at<cv::Vec3b>(ii, jj)[c] * mask.at<float>(y, x);
						}
					}
				}
			}
		}
	}

	out.convertTo(out, CV_8U, 1, 128);
	cv::imwrite("emboss_opencv.bmp", out);

	return 0;
}