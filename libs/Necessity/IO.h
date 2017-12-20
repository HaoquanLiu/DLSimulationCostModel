#ifndef IO_H_
#define IO_H_

#include <string>
#include <vector>

using namespace std;

// read a file as a string
string ReadFileAsString(string);

// read data from the input csv file
vector<vector<double> > ReadCSV(const string);

// read data from the input csv file
vector<vector<unsigned int> > ReadCSVTemp(const string);

// write data to csv file
void WriteToCSV(string, vector<vector<double> >&);

// write data to csv file
void WriteToCSV(string, vector<vector<unsigned int> >&);

// write data to csv file
void WriteToCSV(string,vector<double>&);

// write data to csv file
void WriteToCSV(string, vector<unsigned int>&);

#endif /* IO_H_ */
