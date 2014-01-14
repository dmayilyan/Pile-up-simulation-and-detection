#define LD_MAX    1.0 	//  Landau curve surrounding box
#define AV_POINTS 7		//  Should be odd
#define PEAK_WIDTH 50	//  Peak width
#define av_p 47			//  Points by which the signal will be averaged

const int GATE =  1024;
double peak_av = 0;

Double_t tau_BGO,Npe,firstPETime;


Double_t get_landau_distr(Double_t start_point, Int_t event_count, Double_t tau, Double_t Npe, TH1F* m)
{
	Double_t t, y, x;
	Int_t count = 0;
	Int_t event_no = 0;
	// Big number to avoid infinite loop
	while ((event_no < event_count) && (count < 1000000))
	{
		t = gRandom->Rndm()*(GATE-1);	// Random point in the gate
		y = gRandom->Rndm()*LD_MAX; 	// Random point in the Box
		x = start_point + t;
		// Checking if the point is under the curve
		if (y <= Npe*t/tau*exp(-Npe*t/tau))
		{
			event_no++;
			if ((x < GATE) && (x >= 0))
				m->Fill(x);				// Filling histogram
		}
		count++;
	}
}

void get_statistics(int count)
{
	double T, peak_min,peak_max = 0;
	double temp_max;
	peak_min = 1.0e+6;
	peak_max = -1.0e+6;
	peak_av = 0.0;
	double local_sum = 0.0;
	int real_count = 0;
	double mean_array[GATE] = {0.,};
	double std_dev = 0.0;

	TH1F* h=new TH1F("stat data","Temporary histogram",GATE,0,GATE);

	for (int j = 0; j < count; ++j)
	{
		h->Reset();
		temp_max = -1.0e+6;

		T = gRandom->Rndm()*GATE;

		if (T > 900)
			continue;

		get_landau_distr(T, 2000, tau_BGO, Npe, h);

		// average the data

		for (int q = 0; q < GATE-1; q++)
			mean_array[q] = 0;

		for (int s = av_p/2; s <= GATE-1-av_p/2; s++)
		{
			for (int ss = s-av_p/2; ss <= s+av_p/2; ss++)
				mean_array[s] += h->GetBinContent(ss);
			mean_array[s] /= av_p;
			// cout << mean_array[s] << endl;
		}

		for (int kk = 0; kk <= GATE-1; ++kk)
			if (mean_array[kk] > temp_max)
				temp_max = mean_array[kk];
	
		real_count++;

		peak_av += temp_max;
		local_sum += temp_max;
		std_dev += (temp_max - local_sum / real_count) * (temp_max - local_sum / real_count);
		if (temp_max > peak_max)
			peak_max = temp_max;
		if (temp_max < peak_min)
			peak_min = temp_max;
	}
/* // for debugging
	h->Reset();
	for (int j = 0; j <= GATE-1; j++)
		h->SetBinContent(j,mean_array[j]);

	h->Draw();
*/
	peak_av /= real_count;
	std_dev = sqrt(std_dev / real_count);
	cout << peak_min << "\t" << peak_max << "\t" << peak_av << "\t" << std_dev << endl;
}


void PhotoElec(Double_t Energy)
{
	Int_t Nevents=1000/*00*/; // number of events

	Double_t tdc,tdc_s;
	
	// Typical number of photo electrons for our BGOs at 511 keV is 100.
	// (this can be calculated with the resolution of the crystals
	// that is approx sigma = 50 keV)
	Npe=Energy/511.*100.;
	// Decay time of BGO in ns (300)
	// We artificialy increase this number as we need longer tails
	tau_BGO=35*300.; 

	// Creating canvas
	c = new TCanvas("c","Event");

	TH1F* k=new TH1F("k","Number of photons I",GATE,0,GATE);
	TH1F* l=new TH1F("l","Number of photons II",GATE,0,GATE);
	TH1F* m=new TH1F("Hist data","Number of photons I+II",GATE,0,GATE);
	
	// Getting highest peak average for double peak exclusion
	cout << "Getting peak statistics..." << endl;
	get_statistics(1000);

	double sumj = 0.;
	int pileup_count = 0;
	for(Int_t ev_num = 0; ev_num < Nevents; ev_num++)
	{
		// Getting start points for two Landaus.
		double T1 = gRandom->Rndm() * GATE;
		double T2 = gRandom->Rndm() * GATE;

		// Neglecting curves that start over 900, as they won't give pile-up
		if ((T1 > 900) || (T2 > 900))
			continue;

		get_landau_distr(T1, 2000, tau_BGO, Npe, m);
		// Start time, number of points, tau_BGO, Number of phot_el, histogram

		get_landau_distr(T2, 2000, tau_BGO, Npe, m);

		std::stringstream out;       //Converting integer to string
		out << ev_num+1;

		// Making filenames for histograms
		string filename = "SimPics/Histo_";
		filename += out.str() + ".png";

		double values[4] = {0,};
		double slope[GATE] = {0,};
		slope[0] = 1;  //to avoid having a peak on 0 point
		int peak_count = 0;
		int peak_start = 0;
		double mean_array[GATE] = {0.,};

		// Boundary values difinition
		for (int q = 0; q < GATE-1; q++)
		{
			mean_array[q] = 0;
			slope[q] = 0.;
		}

		// Averaging
		for (int s = av_p/2; s <= GATE-1-av_p/2; s++)
		{
			mean_array[s] = m->GetBinContent(s);
/*			for (int ss = s-av_p/2; ss <= s+av_p/2; ss++)
				mean_array[s] += m->GetBinContent(ss);
			mean_array[s] /= av_p;
*/			// cout << mean_array[s] << endl; // For debugging
		}

		// Slope calculation
		for (int j = AV_POINTS/2; j <= GATE-1 - AV_POINTS/2; j++)
		{
			for (int jjj = 0; jjj <= 3; jjj++)
			{
				values[jjj] = 0;
			}

			for (int jj = j - AV_POINTS/2; jj <= j + AV_POINTS/2; jj++)
			{
				values[0] += jj;					//Calculating sum(x)
				values[1] += mean_array[jj];		//Calculating sum(y)
				values[2] += jj * mean_array[jj];	//Calculating sum(x*y)
				values[3] += jj * jj;				//Calculating sum(x^2)

			}
			slope[j] = (values[2] - values[0] * values[1] / AV_POINTS)/(values[3] - values[0] * values[0]/AV_POINTS);
		}

		int j0 = 0;		// Highest peak bin
		int j1 = 0;		// Second peak bin
		int j_min = 0;	// Min point bin between peaks
		int CDI = 0;	// Channel Data Index(Bin)
		double peak0 = 0.;
		double peak1 = 0.;
		for (int j = AV_POINTS/2; j <= GATE-1 - AV_POINTS/2; j++)
		{
			if (slope[j-1] >= 0.0 && slope[j] < 0.0)
			{
				CDI = 0;
				if (mean_array[j] >= mean_array[j-1])
					CDI = j;
				else
					CDI = j-1;

				if (mean_array[CDI] >= peak0)  //  Getting highest peak
				{
					peak0 = mean_array[CDI];
					j0 = CDI;
				}
			}
		}

		// cout << "!!! " << peak_av * 1.5 << "\t" << peak0 << endl;
		// if (peak0 > peak_av*1.7)  // Exclusion of double peak
		// {
		// 	// cout << ev_num+1<< "\t" << peak_av * 1.7 << "\t" << j0 << "\t" << peak0 << endl;
		// 	m->Reset();
		// 	continue;
		// }

		int i = 0, j = 0, j_min_i;
		double max_diff = 0.0, max_diff_i;

		j1 = -1;
		j_min = -1;
		// Try to find the second peak left to the main peak
		if (j0 > PEAK_WIDTH)
		{
			max_diff_i = -1e+6;
			for (i = 1; i < j0 - PEAK_WIDTH; ++i)
			{
				max_diff_i = -1e+6;

				if (slope[i-1] >= 0.0 && slope[i] < 0.0)  // some local peak
				{
					for (j = i; j < j0; ++j)  // find the largest difference
					{
						if (mean_array[i] - mean_array[j] > max_diff_i)
						{
							j_min_i = j;
							max_diff_i = mean_array[i] - mean_array[j];
						}
					}
				}

				if (max_diff_i > max_diff) // This difference is more, exchange
				{
					max_diff = max_diff_i;
					j1 = i;
					j_min = j_min_i;
				}
			}
		}

		// Try to find the second peak right to the main peak
		if (j0 < GATE - PEAK_WIDTH)
		{
			max_diff_i = -1e+6;
			for (i = GATE; i > j0 + PEAK_WIDTH; --i)
			{
				max_diff_i = -1e+6;

				if (slope[i-1] >= 0.0 && slope[i] < 0.0)  // Some local peak
				{
					// Find the largest difference
					for (j = i; j > j0 + PEAK_WIDTH; --j)
					{
						if (mean_array[i] - mean_array[j] > max_diff_i)
						{
							j_min_i = j;
							max_diff_i = mean_array[i] - mean_array[j];
						}
					}
				}

				// This difference is more exchange
				if (max_diff_i > max_diff)
				{
					max_diff = max_diff_i;
					j1 = i;
					j_min = j_min_i;
				}
			}
		}

/*
		if (j1 != -1)
		{
			peak1 = mean_array[j1];
			cout << ev_num + 1 << "\t" << j0 << "\t" << peak0 << "\t" << j1;
			cout << "\t" << peak1 << "\t" << j_min << "\t";
			cout << mean_array[j_min] << "\t"  << max_diff << "\t" << endl;
		}
*/
		// cout << j0 << "\t" << peak0 << "\t" << j1 << "\t" << peak1 << endl;  // For debugging


		if ((peak0 > (peak1 * 1.3) || peak1 > (peak0 * 1.3)))
		{
			pileup_count++;
			sumj += abs(j0 - j1);
		}

		m->Reset();
		for (int j = 0; j <= GATE-1; j++)
			m->SetBinContent(j,mean_array[j]);

		m->Draw();
		m->GetXaxis()->SetTitle("Time (ns)");
		m->GetXaxis()->SetTitleOffset(1.2);
		m->GetYaxis()->SetTitle("Signal simulation");
		m->GetYaxis()->SetTitleOffset(0.8);
		c->SaveAs(filename.c_str()); // Saving histograms to disk
		k->Reset();
		l->Reset();
		m->Reset();

	}

	cout << "Pile-up count: " << pileup_count << endl;
	if (pileup_count != 0)
		cout << (double)sumj / pileup_count << endl; // Average time
}
