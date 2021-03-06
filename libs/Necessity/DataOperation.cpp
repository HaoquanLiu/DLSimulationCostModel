#include <DataOperation.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <exception>

// check whether vector v2 is clockwise/anticloskwise from vector v1 (2D vectors)
// v[0]: x coordinate, v[1]: y coordinate
int SwingDir(vector<double> v1, vector<double> v2)
{
  if(v2[0]*v1[1] > v2[1]*v1[0])
    return anticlockwise;
  else
    return clockwise;
}

// safe version of acos in case of rounding error
double SafeCos (double x)
{
  if (x < -1.0) x = -1.0 ;
  else if (x > 1.0) x = 1.0 ;
  return acos (x) ;
}

// prevent negative zeros when computing extraction of two values
double SafeZero(double val)
{
  if(fabs(val) <= 0.0001)
    val = 0.0;
  return val;
}

// prevent negative zeros when computing extraction of two values
vector<vector<double> > SafeZero(vector<vector<double> >& val)
{
  for(unsigned int i = 0; i < val.size(); i++)
    {
      for(unsigned int j = 0; j < val[i].size(); j++)
	{
	  val[i][j] = SafeZero(val[i][j]);
	}
    }

  return val;
}

double Dot(vector<double>& vec_1, vector<double>& vec_2)
{
  if(vec_1.size() != vec_2.size())
    {
      cout << "Dimensions of the two vectors do not match." << endl;
      terminate();
    }
  double result = 0.0;
  for(unsigned int i = 0; i < vec_1.size(); i++)
    {
      result += vec_1[i] * vec_2[i];
    }
  return result;
}

double Norm(vector<double>& vec)
{
  if(vec.empty())
    {
      cout << "Empty vector." << endl;
      terminate();
    }
  double result = 0.0;
  for(unsigned int i = 0; i < vec.size(); i++)
    {
      result += vec[i] * vec[i];
    }
  return sqrt(result);
}

//method to check if a Coordinate is located in a polygon
bool InPolygon(vector<vector<double> >& polygon_vert,double test_point_x, double test_point_y)
{
  //this method uses the ray tracing algorithm to determine if the point is in the polygon
  unsigned int nPoints=polygon_vert.size();
  unsigned int j=-999;
  unsigned int i=-999;
  bool locatedInPolygon=false;
  for(i=0;i<(nPoints);i++)
    {
      //repeat loop for all sets of points
      if(i==(nPoints-1))
	//if i is the last vertex, let j be the first vertex
	j= 0;
      else
	//for all-else, let j=(i+1)th vertex
	j=i+1;

      double vertY_i= polygon_vert[i][1];
      double vertX_i= polygon_vert[i][0];
      double vertY_j= polygon_vert[j][1];
      double vertX_j= polygon_vert[j][0];

      // following statement checks if testPoint.Y is below Y-coord of i-th vertex
      bool belowLowY=vertY_i>test_point_y;
      // following statement checks if testPoint.Y is below Y-coord of i+1-th vertex
      bool belowHighY=vertY_j>test_point_y;

      /* following statement is true if testPoint.Y satisfies either (only one is possible)
	 -->(i).Y < testPoint.Y < (i+1).Y        OR
	 -->(i).Y > testPoint.Y > (i+1).Y

	 (Note)
	 Both of the conditions indicate that a point is located within the edges of the Y-th coordinate
	 of the (i)-th and the (i+1)- th vertices of the polygon. If neither of the above
	 conditions is satisfied, then it is assured that a semi-infinite horizontal line draw
	 to the right from the testpoint will NOT cross the line that connects vertices i and i+1
	 of the polygon
      */
      bool withinYsEdges= (belowLowY != belowHighY);

      if( withinYsEdges)
	{
	  // this is the slope of the line that connects vertices i and i+1 of the polygon
	  double slopeOfLine   = ( vertX_j-vertX_i )/ (vertY_j-vertY_i) ;

	  // this looks up the x-coord of a point lying on the above line, given its y-coord
	  double pointOnLine   = ( slopeOfLine* (test_point_y - vertY_i) )+vertX_i;

	  //checks to see if x-coord of testPoint is smaller than the point on the line with the same y-coord
	  bool isLeftToLine= test_point_x < pointOnLine;

	  if(isLeftToLine)
	    {
	      //this statement changes true to false (and vice-versa)
	      locatedInPolygon= !locatedInPolygon;
	    }//end if (isLeftToLine)
	}//end if (withinYsEdges
    }

  return locatedInPolygon;
}

double SumVector(const vector<double>& vec)
{
  return accumulate(vec.begin(),vec.end(),0.0);
}

double MinVector(vector<double>& vec)
{
  return *min_element(vec.begin(),vec.end());
}

double MaxVector(vector<double>& vec)
{
  return *max_element(vec.begin(),vec.end());
}

vector<double> SubtractVector(const vector<double>& v1,const  vector<double>& v2)
{
  if(v1.size() != v2.size())
    {
      cout << "Dimensions of these two vectors do not match." << endl;
      terminate();
    }
  vector<double> result(v1.size(),0.0);
  for(unsigned int i = 0; i < v1.size(); i++)
    {
      result[i] = v1[i] - v2[i];
    }
  return result;
}

vector<unsigned int> SubtractVector(const vector<unsigned int>& v1,const  vector<unsigned int>& v2)
{
  if(v1.size() != v2.size())
    {
      cout << "Dimensions of these two vectors do not match." << endl;
      terminate();
    }
  vector<unsigned int> result(v1.size(),0.0);
  for(unsigned int i = 0; i < v1.size(); i++)
    {
      result[i] = v1[i] - v2[i];
    }
  return result;
}

vector<double> SumVector(const vector<double>& v1,const  vector<vector<double> >& v2)
{
  if(v1.size() != v2.front().size())
    {
      cout << "Dimensions of these two vectors do not match." << endl;
      terminate();
    }
  vector<double> result(v1.size(),0.0);
  for(unsigned int i = 0; i < v1.size(); i++)
    {
      result[i] += v1[i];
      for(unsigned int j = 0; j < v2.size(); j++)
	result[i] += v2[j][i];
    }
  return result;
}

vector<double> SumVector(const vector<double>& v1,const vector<double>& v2)
{
  if(v1.size() != v2.size())
    {
      cout << "Dimensions of these two vectors do not match." << endl;
      terminate();
    }
  vector<double> result(v1.size(),0.0);
  for(unsigned int i = 0; i < v1.size(); i++)
    {
      result[i] += v1[i]+v2[i];
    }
  return result;
}


vector<unsigned int> SumVector(const vector<unsigned int>& v1,const  vector<vector<unsigned int> >& v2)
{
  if(v1.size() != v2.front().size())
    {
      cout << "Dimensions of these two vectors do not match." << endl;
      terminate();
    }
  vector<unsigned int> result(v1.size(),0.0);
  for(unsigned int i = 0; i < v1.size(); i++)
    {
      result[i] += v1[i];
      for(unsigned int j = 0; j < v2.size(); j++)
	result[i] += v2[j][i];
    }
  return result;
}

void PolyTran(vector<vector<double> >& poly, vector<double>& direction, double length)
{
  // compute the unit vector along the translation direction
  vector<double> unit_vector(2,0.0);
  unit_vector[0] = direction[0]/sqrt(pow(direction[0],2)+pow(direction[1],2));
  unit_vector[1] = direction[1]/sqrt(pow(direction[0],2)+pow(direction[1],2));

  for(unsigned int i = 0; i < poly.size(); i++)
    {
      for(unsigned int j = 0; j < 2; j++)
	{
	  poly[i][j] += length*unit_vector[j];
	}
    }
}

vector<unsigned int> IndicesCoal(vector<double>& h_target,double coal_height)
{
  vector<unsigned int> indices_coal;
  for(unsigned int i = 0; i < h_target.size(); i++)
    {
      if(h_target[i] <= coal_height)
	indices_coal.push_back(i);
    }
  return indices_coal;
}

double ExposeRatio(vector<double>& h_current,vector<double>& h_target,vector<unsigned int>& indices_coal)
{
  double expose_ratio = 0.0;
  double count = 0.0;
  for(unsigned int i = 0; i < indices_coal.size(); i++)
    {
      if(SafeZero(h_current[indices_coal[i]] - h_target[indices_coal[i]]) == 0.0)
	count++;
    }
  expose_ratio = count/(double)(indices_coal.size());
  return expose_ratio;
}
