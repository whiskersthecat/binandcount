// give two categorized datasets from happypenguin
// randomly subsample molecules in binsizes and count
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sys/stat.h>
#include <iomanip>
#include <climits>
#include <random>
#include <cmath>

using namespace std;

// copied from happypenguin2.0.2.cpp
#define NUM_CATEGORIES 9
#define NUM_DATASETS 2
static string CATEGORY_NAMES[] = {"Genotype1", "Genotype2", "LowQualityCrossover", "MediumQualityCrossover", "HighQualityCrossover","ChompBadMapping", "ChompNoSNPs", "ChompBadGT", "ChompMultimapRec"};

int binsize = 0;
int num_bins;
int MAX_BINSIZE = 100000;
int FACTOR = 2;
string directories[NUM_DATASETS];

vector<vector<int>> binned_lengths[NUM_DATASETS];
vector<int> sampling_size;

ofstream logs[NUM_CATEGORIES];

// formatting
#define SPACING 15

int NUM_EXPERIMENTS = 3;

// vector<int> experiment_counts[NUM_EXPERIMENTS][NUM_CATEGORIES][NUM_DATASETS];

vector<vector<vector<vector<int>>>> experiment_counts;

// experiment totals are in the last column
vector<double> experiment_stats[2][NUM_CATEGORIES][NUM_DATASETS];
string stat_names[2] = {"Average", "STDev"};

// vector<vector<int>> experiment_totals;

random_device rd;
mt19937 PRNG(rd());

// [dataset][length_bin_number] -> category e.g. 0, 0, 0, 5, 2, 3, 8, ...

void log(int category, string text, int width) {
	// write to output file(s)
	if(category == -1) {
		for(int i = 0; i < NUM_CATEGORIES; ++i) {
			if(width != -1)
				logs[i] << left << setw(width) << text;
			else
				logs[i] << text;
		}
	}
	else
		if(width != -1)
			logs[category] << left << setw(width) << text;
		else
			logs[category] << text;
}

void populateLengths(int dataset_num) {
	for (int category = 0; category < NUM_CATEGORIES; ++category) {
		string ifile_path = directories[dataset_num] + CATEGORY_NAMES[category] + "_lengths.tsv";
		ifstream ifile(ifile_path);
		if(!ifile) {
			cout << "[error] Trouble opening file " << ifile_path << endl;
			exit(1);
		}
		int len;

		// cout << "[parse] populating from file " << ifile_path << endl;

		while(ifile >> len) {
			int bin = 0;
			bin = len / binsize; // will floor to the closest bin
			if (bin > num_bins) {
				bin = num_bins;
			}
			binned_lengths[dataset_num][bin].push_back(category);
		}
	}
	return;
}

void countTotalMolecules() {

	cout << "[info] counting reads" << endl;

	vector <int> smallest(num_bins, INT_MAX);

	for (int j = 0; j < NUM_DATASETS; ++j) {
		int total  = 0;
		log(-1, "Dataset_" + to_string(j + 1) + "_total", 15);
		for (int i = 0; i < num_bins; ++i) {
			int count = binned_lengths[j][i].size();
			// cout << i << '\t' << count << '\t' << smallest[i] << endl;
			if (count < smallest[i]) {
				
				smallest[i] = count;
			}
			total += count;
			log(-1, "\t" + to_string(count), SPACING);
		}
		log(-1, "\t" + to_string(total) + "\n", -1);
	}

	log(-1, "Sampling_size", 15);
	int total = 0;
	for (int i = 0; i < num_bins; ++i) {
		sampling_size[i] = smallest[i] / FACTOR;
		total += sampling_size[i];
		log(-1, "\t" + to_string(sampling_size[i]), SPACING);
	}
	log(-1, "\t" + to_string(total) + "\n", -1);
	
}

void runExperiment(int experiment_num) {

	cout << "[info] starting experiment number " << experiment_num + 1 << endl;
	
	// stores the counts for this experiment
	// vector<int> experiment_counts[NUM_CATEGORIES][NUM_DATASETS];

	for (int category = 0; category < NUM_CATEGORIES; ++category) {
		for(int j = 0; j < NUM_DATASETS; ++j) {
			for (int i = 0; i < num_bins + 1; ++i) {
				experiment_counts[experiment_num][category][j].push_back(0);
			}
		}
	}

	for (int j = 0; j < NUM_DATASETS; ++j) {
		for (int i = 0; i < num_bins; ++i) {
			// shuffle the vector
			shuffle(binned_lengths[j][i].begin(), binned_lengths[j][i].end(), PRNG);
			for(int n = 0; n < sampling_size[i]; ++n) {
				// draw and count the top n reads
				int read_category = binned_lengths[j][i][n];
				experiment_counts[experiment_num][read_category][j][i] += 1;
				// increment the total
				experiment_counts[experiment_num][read_category][j][num_bins] += 1;
			}
		}
	}

}

void calculateStatistics() {

	cout << "[info] calculating average and standard deviations" << endl;

	for (int category = 0; category < NUM_CATEGORIES; ++category) {
		for(int j = 0; j < NUM_DATASETS; ++j) {
			for (int i = 0; i < num_bins + 1; ++i) {
				for(int num_stats = 0; num_stats < 2; ++num_stats) {
					experiment_stats[num_stats][category][j].push_back(0);
				}
			}
		}
	}
	
	for (int category = 0; category < NUM_CATEGORIES; ++category) {
		for (int i = 0; i < num_bins + 1; ++i) {
			for(int j = 0; j < NUM_DATASETS; ++j) {
				int total = 0;
				for(int experiment_num = 0; experiment_num < NUM_EXPERIMENTS; ++experiment_num) {
					total += experiment_counts[experiment_num][category][j][i];
				}
				double mean = (double) total / (double) NUM_EXPERIMENTS;

				double std = 0;
				for(int experiment_num = 0; experiment_num < NUM_EXPERIMENTS; ++experiment_num) {
					std += pow(experiment_counts[experiment_num][category][j][i] - mean, 2);
				}
				std /= (double) NUM_EXPERIMENTS;
				std = pow(std, 0.5);
				experiment_stats[0][category][j][i] = mean;
				experiment_stats[1][category][j][i] = std;
			}
		}
	}

}

void logResults() {
	cout << "[info] writing output to files" << endl;

	for(int experiment_num = 0; experiment_num < NUM_EXPERIMENTS + 2; ++experiment_num) {

		if(experiment_num < NUM_EXPERIMENTS) {
			log(-1 , "Experiment_" + to_string(experiment_num + 1), 15);
		}
		else {
			if(experiment_num == NUM_EXPERIMENTS)
				log(-1, "\n", -1);
			log(-1 , stat_names[experiment_num - NUM_EXPERIMENTS], 15);
		}
		for (int category = 0; category < NUM_CATEGORIES; ++category) {

			for (int i = 0; i < num_bins + 1; ++i) {
				string results_str = "";
				for(int j = 0; j < NUM_DATASETS; ++j) {
					if(experiment_num < NUM_EXPERIMENTS)
						results_str += to_string(experiment_counts[experiment_num][category][j][i]);
					else
						results_str += to_string(experiment_stats[experiment_num - NUM_EXPERIMENTS][category][j][i]).substr(0, 4);
					if (j != NUM_DATASETS - 1) {
						results_str += ",";
					}
				}
				log(category, "\t" + (results_str), SPACING);
			}

			log(category, "\n", -1);
		}
	}

}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "Usage: ./binandcount binsize num_experiments foldername1 foldername2" << endl;
	}
	NUM_EXPERIMENTS = atoi(argv[2]);
	directories[0] = string(argv[3]) + "/misc/read_lengths/";
	directories[1] = string(argv[4]) + "/misc/read_lengths/";
	binsize = atoi(argv[1]);
	cout << "[info] max bin catches sequence more than " << MAX_BINSIZE << endl;
	cout << "[info] binsize set to " << binsize << endl;
	num_bins = (MAX_BINSIZE / binsize) + 1;
	cout << "[info] total number of bins " << num_bins << endl;

	for (int i = 0; i < num_bins; ++i) {
		for(int j = 0; j < NUM_DATASETS; ++j) {
			binned_lengths[j].push_back(vector<int>());
		}
		sampling_size.push_back(-1);
	}

	string output_dir = "binandcount_output_" + string(argv[1]) + "_" + argv[2] + "_" + argv[3] + "_" + argv[4];
	cout << "[info] creating output files in directory " << output_dir << endl;
	mkdir((output_dir).c_str(), 0777);
	for(int i = 0; i < NUM_CATEGORIES; ++i) {
		logs[i].open(output_dir + "/" + CATEGORY_NAMES[i]);
		// log(i, CATEGORY_NAMES[i]);
	}

	log(-1 , "Bin_Ranges", 15);
	for (int i = 0; i < num_bins; ++i) {
		if (i == num_bins - 1) {
			log(-1, "\t" + to_string(binsize * i) + "+", SPACING);
			log(-1, "\tTotal\n", -1);
		}
		else
			log(-1, "\t" + to_string(binsize * i), SPACING);
	}

	// read in lengths into memory
	populateLengths(0);
	populateLengths(1);

	// get sample size and distribution
	countTotalMolecules();

	// run experiments
	experiment_counts.resize(NUM_EXPERIMENTS);
	for (int i = 0; i < NUM_EXPERIMENTS; ++i) {
		experiment_counts[i].resize(NUM_CATEGORIES);
		for (int j = 0; j < NUM_CATEGORIES; ++j) {
			experiment_counts[i][j].resize(NUM_DATASETS);
		}
	}

	log(-1, "\n", -1);

	cout << "[info] running " << NUM_EXPERIMENTS << " total experiments" << endl;
	for (int i = 0; i < NUM_EXPERIMENTS; ++i) {
		runExperiment(i);
	}

	// calculate average and standard deviation
	calculateStatistics();

	// write experiment results to output files
	logResults();

	return 0;
}