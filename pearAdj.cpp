//////////////////////////////////////////////////////////////////////////////
// This was drated using ChatGPT, but required editing/debugging to work
// Author: ChatGPT and me, Jared Streich
// Version: 0.1.0
// email: streich.jared@gmail.com or at ORNL ju0@ornl.gov
// This should return an unordered pairwise matrix of pearson values.
// Date: 2023-02-19
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
////////////////////////////// Inlude Libraries //////////////////////////////
//////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Start Script ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



using namespace std;

int main(int argc, char* argv[]) {

    ///// Check if filename is provided
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <data_file>" << endl;
        return 1;
    }

    ////// Open file for reading
    ifstream input_file(argv[1]);
    if (!input_file.is_open()) {
        cerr << "Failed to open file: " << argv[1] << endl;
        return 1;
    }

    ////// Read data into vector
    vector<vector<double> > data;
    vector<string> sample_names;
    string line;
    while (getline(input_file, line)) {
        vector<double> row;
        string sample_name;
        double value;
        istringstream iss(line);
        iss >> sample_name;
        while (iss >> value) {
            row.push_back(value);
        }
        data.push_back(row);
        sample_names.push_back(sample_name);
    }
    input_file.close();

    ///// Calculate pairwise Pearson correlation coefficient
    int n = data.size();
    vector<vector<double> > corr(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            double x_sum = 0, y_sum = 0, xy_sum = 0, x2_sum = 0, y2_sum = 0;
            int k;
            for (k = 0; k < data[i].size() && k < data[j].size(); k++) {
                x_sum += data[i][k];
                y_sum += data[j][k];
                xy_sum += data[i][k] * data[j][k];
                x2_sum += data[i][k] * data[i][k];
                y2_sum += data[j][k] * data[j][k];
            }
            double denom = sqrt(x2_sum - x_sum * x_sum / k) * sqrt(y2_sum - y_sum * y_sum / k);
            corr[i][j] = denom == 0 ? 0 : (xy_sum - x_sum * y_sum / k) / denom;
            corr[j][i] = corr[i][j];
        }
    }

    ///// Write output to file
    ofstream output_file("pearson_corr_output.txt");
    output_file << "Sample Names\t";
    for (string sample_name:sample_names) {
        output_file << sample_name << "\t";
    }
    output_file << endl;
    for (int i = 0; i < n; i++) {
        output_file << sample_names[i] << "\t";
        for (int j = 0; j < n; j++) {
            output_file << corr[i][j] << "\t";
        }
        output_file << endl;
    }
    output_file.close();

    return 0;
}
