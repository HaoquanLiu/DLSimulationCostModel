#include <IO.h>
#include <fstream>
#include <sstream>
#include <iomanip>

// read a file as a string
string ReadFileAsString(string fname)
{
  ifstream in(fname.c_str(), ios::in);
  string contents;
  in.seekg(0, ios::end);
  contents.resize(in.tellg());
  in.seekg(0, ios::beg);
  in.read(&contents[0], contents.size());
  in.close();
  return (contents);
}

// read data from the input csv file
vector<vector<double> > ReadCSV(const string fname)
{
  std::stringstream stringstreamwhole(ReadFileAsString(fname));
  vector<vector<double> > data;
  while(true)
    {
      string stringline;
      if(!getline(stringstreamwhole,stringline))
	break;
      stringstream stringstreamline(stringline);
      vector<double> temp_data;
      while(stringstreamline)
	{
	  string stringval;
	  if(!getline(stringstreamline,stringval,','))
	    break;
	  stringstream stringstreamval(stringval);
	  double temp;
	  stringstreamval >> std::setprecision(15) >> temp;
	  temp_data.push_back(temp);
	}
      data.push_back(temp_data);
    }

  return data;
}

// read data from the input csv file
vector<vector<unsigned int> > ReadCSVTemp(const string fname)
{
  std::stringstream stringstreamwhole(ReadFileAsString(fname));
  vector<vector<unsigned int> > data;
  while(true)
    {
      string stringline;
      if(!getline(stringstreamwhole,stringline))
	break;
      stringstream stringstreamline(stringline);
      vector<unsigned int> temp_data;
      while(stringstreamline)
	{
	  string stringval;
	  if(!getline(stringstreamline,stringval,','))
	    break;
	  stringstream stringstreamval(stringval);
	  unsigned int temp;
	  stringstreamval >> temp;
	  temp_data.push_back(temp);
	}
      data.push_back(temp_data);
    }

  return data;
}


void WriteToCSV(string fname, vector<vector<double> >& data)
{
  ofstream out_stream;
  out_stream.open(fname);
  for(unsigned int i = 0; i < data.size(); i++)
    {
      for(unsigned int j = 0; j < data[i].size(); j++)
	{
	  if(j != data[i].size()-1)
	    out_stream << std::setprecision(15) << data[i][j] << ",";
	  else
	    out_stream << std::setprecision(15) << data[i][j] << endl;
	}
    }
  out_stream.close();
}

void WriteToCSV(string fname, vector<vector<unsigned int> >& data)
{
  ofstream out_stream;
  out_stream.open(fname);
  for(unsigned int i = 0; i < data.size(); i++)
    {
      for(unsigned int j = 0; j < data[i].size(); j++)
	{
	  if(j != data[i].size()-1)
	    out_stream << data[i][j] << ",";
	  else
	    out_stream << data[i][j] << endl;
	}
    }
  out_stream.close();
}

void WriteToCSV(string fname, vector<double>& data)
{
  ofstream out_stream;
  out_stream.open(fname);
  for(unsigned int i = 0; i < data.size(); i++)
    {
      out_stream << data[i] << endl;
    }
  out_stream.close();
}

void WriteToCSV(string fname, vector<unsigned int>& data)
{
  ofstream out_stream;
  out_stream.open(fname);
  for(unsigned int i = 0; i < data.size(); i++)
    {
      out_stream << data[i] << endl;
    }
  out_stream.close();
}
