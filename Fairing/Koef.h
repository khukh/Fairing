#pragma once
class DampingKoef{
public:
	DampingKoef();
	
	~DampingKoef();
	double linInterp(double x1, double x2, double y1, double y2, double x);
	double MrFr(double mach, double alpha);
	double MpFr(double mach, double alpha);
	double MyFr(double mach, double alpha);

	double MrFp(double mach, double alpha);
	double MpFp(double mach, double alpha);
	double MyFp(double mach, double alpha);

	double MrFy(double mach, double alpha);
	double MpFy(double mach, double alpha);
	double MyFy(double mach, double alpha);

	double MrB(double mach, double beta);
	double MpB(double mach, double beta);
	double MyB(double mach, double beta);

private:
	double ArMrFr[5][10], ArMpFr[5][10], ArMyFr[5][10];
	double ArMrFp[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, 0, 0, 0.1412, -0.2625, -0.0187, -0.075, -0.6, 0.2471, 0.0},
		{0.98, 0.0, -0.0187, 0.4059, -0.4235, 0, -0.3, -2.8125, 0.0529, 0},
		{1.2, 0.0, -0.1687, -1.0412, -0.3, 0, 0.0529, -0.9882, 0.7412, 0},
		{2, 0.0, 0.5471, -0.075, -0.1875, 0, 0.1412, -3.0176, 0.4588, 0}
	};
	double ArMpFp[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -5.187, -3.882, -0.562, -2.625, -6.294, -8.75, -1.625, -9.812, -5.187},
		{0.98, -9, -5.25, 1.438, -1.125, -10, -12.375, 1.938, -10.706, -9},
		{1.2, -10.706, -4.937, -6.176, -2.625, -10.471, -9.875, -4.187, -10.824, -10.706},
		{2, -8.25, -5.812, -8.5, -8.937, -5.812, -14, -13.353, -11.125, -8.25}
	};
	double ArMyFp[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -0.0, -0.0159, -0.1263, -0.4737, -0.0474, -0.0474, -0.4896, 0.3474, 0},
		{0.98, -0.0, -0.0789, 0.2211, -1.1367, -0.0159, -0.6159, -3.6159, -0.0474, 0.0},
		{1.2, -0.0, -0.2211, 0.1578, 1.3896, 0.0474, 0, -0.1737, -0.2211, -0.0},
		{2, -0.0, 0.3315, -0.0315, 0.8685, 0.0315, -0.0474, -1.1526, 1.5, -0.0}

	};

	double ArMrFy[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, 1.686, -0.232, 0.371, -0.087, -1.536, -6.265, -0.319, 4.629, 1.686},
		{0.98, 1.257, -1.217, -0.261, 1.171, -0.783, -2.171, -0.29, -3.6, 1.257},
		{1.2, 1.229, -0.377, -0.029, 0.571, -0.464, -1.275, -0.435, 0.6, 1.229},
		{2, 1.486, 0.114, -0.29, 0.143, -0.58, -0.638, 0.457, 0.829, 1.486}
	};
	double ArMpFy[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -0.516, -0.548, -0.548, -0.6915, -0.468, -0.484, -0.452, -0.3405, -0.516},
		{0.98, -0.5, -0.532, -0.564, 0.8885, -0.516, -0.165, -0.532, -0.2285, -0.5},
		{1.2, -0.5, -0.819, -0.149, -2.516, -0.532, -0.436, -0.4045, -0.532, -0.5},
		{2, -0.516, -0.947, -0.787, -2.8225, -0.516, -0.452, 2.3245, 0.1385, -0.516}
	};
	double ArMyFy[5][10];

	double ArMrFb[4][10] = {
		{0, -90, -60, -30, 0, 30, 60, 90},
		{0.6, 1.333, 2.604, 0.125, -0.208, 0.125, 2.625, 1.354},
		{0.98, -0.354, 0.542, -0.354, -0.292, -0.333, 0.542, -0.396},
		{1.2, -0.792, -0.375, -0.729, -0.375, -0.75, -0.396, -0.792}
	};
	double ArMpFb[4][10] = {
		{0, -90, -60, -30, 0, 30, 60, 90},
		{0.6, -1.84, 1.85, -2.519, -3.897, -2.519, 1.85, -1.783},
		{0.98, -2.745, -0.028, -4.935, -8.355, -4.935, 0, -2.774},
		{1.2, -3.28, -2.972, -4.991, -8.888, -4.935, -2.943, -3.308}
	};
	double ArMyFb[4][10] = {
		{0, -90, -60, -30, 0, 30, 60, 90},
		{0.6, 2.477, 0.209, -3.244, -3.244, -3.349, 0.174, 2.512},
		{0.98, 0.663, -3.593, -4.744, -4.605, -4.779, -3.628, 0.663},
		{1.2, 0.628, -2.442, -4.535, -4.5, -4.535, -2.477, 0.523},
	};
};
