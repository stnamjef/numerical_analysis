#include <iostream>
#include <fstream>
using namespace std;

int main()
{
	// sprint 1
	float fs, dt, f, A, t;

	fs = 44100.F;
	dt = 1.F / fs;
	f = 440.F;
	A = 10000.F;

	int n = (int)fs * 1;
	short* data = new short[n];

	for (int i = 0; i < n; i++) {
		t = i * dt;
		data[i] = (short)(A * sin(2 * 3.141592 * f * t));
	}

	// sprint2

	ifstream src("Beatles.wav", ios::binary | ios::in);

	WaveHeader head;
	src.read((char*)&head, sizeof(head));
	src.close();

	ofstream tmp("mywave_multi_channels.wav", ios::binary | ios::out);
	tmp.write((char*)&head, sizeof(head));
	tmp.write((char*)data, sizeof(short) * n);
	tmp.close();

	delete[] data;

	return 0;
}