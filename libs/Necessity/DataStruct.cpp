#include <DataOperation.h>
#include <IO.h>
#include <DataStruct.h>
#include <cmath>
#include <algorithm>
#include <fstream>

point_cloud::point_cloud(double x_cellsize_, double y_cellsize_)
{
  x_cellsize = x_cellsize_;
  y_cellsize = y_cellsize_;
  created_from_polygon_extraction = false;
}

point_cloud::point_cloud(const string xyz_file_name, double x_cellsize_, double y_cellsize_)
{
  x_cellsize = x_cellsize_;
  y_cellsize = y_cellsize_;
  points = ReadCSV(xyz_file_name);
  for(unsigned int i = 0; i < points.size(); i++)
    {
      if(points[i].size() <= 2)
	{
	  cout << "Dimensions of the input data are not sufficient." << endl;
	  terminate();
	}
      if(points[i].size() >= 4)
	points[i].resize(3);
    }
  created_from_polygon_extraction = false;
  polygon.assign(4,vector<double>(2,0.0));
  polygon[0][0] = MinX();
  polygon[0][1] = MinY();
  polygon[1][0] = MinX();
  polygon[1][1] = MaxY();
  polygon[2][0] = MaxX();
  polygon[2][1] = MaxY();
  polygon[3][0] = MaxX();
  polygon[3][1] = MinY();
}

point_cloud::point_cloud(const string xyz_file_name)
{
  x_cellsize = 0.0;
  y_cellsize = 0.0;
  points = ReadCSV(xyz_file_name);
  for(unsigned int i = 0; i < points.size(); i++)
    {
      if(points[i].size() <= 2)
	{
	  cout << "Dimensions of the input data are not sufficient." << endl;
	  exit(0);
	}
      if(points[i].size() >= 4)
	points[i].resize(3);
    }
  created_from_polygon_extraction = false;
  polygon.assign(4,vector<double>(2,0.0));
  polygon[0][0] = MinX();
  polygon[0][1] = MinY();
  polygon[1][0] = MinX();
  polygon[1][1] = MaxY();
  polygon[2][0] = MaxX();
  polygon[2][1] = MaxY();
  polygon[3][0] = MaxX();
  polygon[3][1] = MinY();
}


point_cloud::~point_cloud(){}

void point_cloud::ResetPolygon()
{
  polygon[0][0] = MinX();
  polygon[0][1] = MinY();
  polygon[1][0] = MinX();
  polygon[1][1] = MaxY();
  polygon[2][0] = MaxX();
  polygon[2][1] = MaxY();
  polygon[3][0] = MaxX();
  polygon[3][1] = MinY();
  created_from_polygon_extraction = false;
}

void point_cloud::Append(vector<double>& point)
{
  points.push_back(point);
  ResetPolygon();
}

vector<double> point_cloud::XCoordinates()
{
  if(points.empty())
    {
      cout << "Point cloud is empty." << endl;
      terminate();
    }
  vector<double> x_coordinates(points.size(),0.0);
  for(unsigned int i = 0; i < points.size(); i++)
    x_coordinates[i] = points[i][0];

  return x_coordinates;
}

vector<double> point_cloud::YCoordinates()
{
  if(points.empty())
    {
      cout << "Point cloud is empty." << endl;
      terminate();
    }
  vector<double> y_coordinates(points.size(),0.0);
  for(unsigned int i = 0; i < points.size(); i++)
    y_coordinates[i] = points[i][1];

  return y_coordinates;
}
 
vector<double> point_cloud::ZCoordinates() 
{
  if(points.empty())
    {
      cout << "Point cloud is empty." << endl;
      terminate();
    }
  vector<double> z_coordinates(points.size(),0.0);
  for(unsigned int i = 0; i < points.size(); i++)
    z_coordinates[i] = points[i][2];

  return z_coordinates;
}

void point_cloud::AssignXCoordinates(const vector<double>& new_z)
{
  if(new_z.size() != points.size())
    {
      cout << "Dimensions do not match." << endl;
      terminate();
    }
  for(unsigned int i = 0; i < new_z.size(); i++)
    points[i][0] = new_z[i];
}

void point_cloud::AssignXCoordinates(unsigned int index, double new_z)
{
    points[index][0] = new_z;
}

void point_cloud::AssignYCoordinates(const vector<double>& new_z)
{
  if(new_z.size() != points.size())
    {
      cout << "Dimensions do not match." << endl;
      terminate();
    }
  for(unsigned int i = 0; i < new_z.size(); i++)
    points[i][1] = new_z[i];
}

void point_cloud::AssignYCoordinates(unsigned int index, double new_z)
{
    points[index][1] = new_z;
}

void point_cloud::AssignZCoordinates(const vector<double>& new_z)
{
  if(new_z.size() != points.size())
    {
      cout << "Dimensions do not match." << endl;
      terminate();
    }
  for(unsigned int i = 0; i < new_z.size(); i++)
    points[i][2] = new_z[i];
}

void point_cloud::AssignZCoordinates(unsigned int index, double new_z)
{
  points[index][2] = new_z;
}

vector<double> point_cloud::Centroid2D()
{
  vector<double> x_coordinates = XCoordinates();
  vector<double> y_coordinates = YCoordinates();
  vector<double> centroid(2,0.0);
  centroid[0] = SumVector(x_coordinates) / points.size();
  centroid[1] = SumVector(y_coordinates) / points.size();
  return centroid;
}

vector<double> point_cloud::Centroid3D()
{
  vector<double> x_coordinates = XCoordinates();
  vector<double> y_coordinates = YCoordinates();
  vector<double> z_coordinates = ZCoordinates();
  vector<double> centroid(3,0.0);
  centroid[0] = SumVector(x_coordinates) / points.size();
  centroid[1] = SumVector(y_coordinates) / points.size();
  centroid[2] = SumVector(z_coordinates) / points.size();
  return centroid;
}

double point_cloud::MinX()
{
  vector<double> x_coordinates = XCoordinates();
  return *min_element(x_coordinates.begin(),x_coordinates.end());
}

double point_cloud::MaxX()
{
  vector<double> x_coordinates = XCoordinates();
  return *max_element(x_coordinates.begin(),x_coordinates.end());
}

double point_cloud::MinY()
{
  vector<double> y_coordinates = YCoordinates();
  return *min_element(y_coordinates.begin(),y_coordinates.end());
}

double point_cloud::MaxY()
{
  vector<double> y_coordinates = YCoordinates();
  return *max_element(y_coordinates.begin(),y_coordinates.end());
}

double point_cloud::MinZ()
{
  vector<double> z_coordinates = ZCoordinates();
  return *min_element(z_coordinates.begin(),z_coordinates.end());
}

double point_cloud::MaxZ()
{
  vector<double> z_coordinates = ZCoordinates();
  return *max_element(z_coordinates.begin(),z_coordinates.end());
}

void point_cloud::UpdateRegion(point_cloud& another_point_cloud)
{
  double min_x = another_point_cloud.MinX();
  double max_x = another_point_cloud.MaxX();
  double min_y = another_point_cloud.MinY();
  double max_y = another_point_cloud.MaxY();
  for(unsigned int i = 0; i < points.size(); i++)
    {
      if(points[i][0] <= max_x && points[i][0] >= min_x && points[i][1] <= max_y && points[i][1] >= min_y)
	{
	  for(unsigned int j = 0; j < another_point_cloud.points.size(); j++)
	    {
	      if(SafeZero(another_point_cloud.points[j][0]-points[i][0]) == 0.0 && SafeZero(another_point_cloud.points[j][1]-points[i][1]) == 0.0)
		{
		  points[i][2] = another_point_cloud.points[j][2];
		  break;
		}
	    }
	}
    }

}

point_cloud point_cloud::ExtractPointCloud(vector<vector<double> > specified_polygon)
{
  vector<double> x_polygon(specified_polygon.size(),0.0);
  vector<double> y_polygon(specified_polygon.size(),0.0);
  for(unsigned int i = 0; i < specified_polygon.size(); i++)
    {
      x_polygon[i] = specified_polygon[i][0];
      y_polygon[i] = specified_polygon[i][1];
    }
  double min_x = MinVector(x_polygon);
  double min_y = MinVector(y_polygon);
  double max_x = MaxVector(x_polygon);
  double max_y = MaxVector(y_polygon);
  point_cloud result(this->x_cellsize,this->y_cellsize);
  for(unsigned int i = 0; i < points.size(); i++)
    {
      // if it is possible that the element is within the polygon
      if(points[i][0] <= max_x && points[i][0] >= min_x && points[i][1] <= max_y && points[i][1] >= min_y)
	{
	  if(InPolygon(specified_polygon,points[i][0],points[i][1]))
	    {
	      result.points.push_back(points[i]);
	    }	
	}
    }
  result.created_from_polygon_extraction = true;
  result.polygon = specified_polygon;
  return result;
}

unsigned int point_cloud::IndexAtCoordinates(const double x,const double y)
{
  // if((x < this->MinX() || x > this->MaxX() || y < this->MinY() || y > this->MaxY()) &&
  //    !InPolygon(this->polygon,x,y))
  //   {
  //     cout << "Coordinates of the point is " << x << ", " << y << endl;
  //     cout << "The specified point is outside the possible range of the point cloud." << endl;
  //     terminate();
  //   }
  unsigned int index = 9999;
  double min_dist_metric = 1000000;
  for(unsigned int i = 0; i < points.size(); i++)
    {
      if(fabs(points[i][0]-x) <= x_cellsize && fabs(points[i][1]-y) <= y_cellsize) // trick to accelerate to computation
	{
	  double dist_metric = pow(x-points[i][0],2) + pow(y-points[i][1],2);
	  if(dist_metric < min_dist_metric)
	    {
	      index = i;
	      min_dist_metric = dist_metric;
	    }
	}
    }

  return index;
}

vector<unsigned int> point_cloud::IndicesAtCoordinates(const vector<vector<double> >& xys)
{
  vector<unsigned int> indices(xys.size(),9999);
  for(unsigned int i = 0; i < indices.size(); i++)
    {
      if(xys[i].size() <= 1)
	{
	  cout << "Please provide at least x,y coordinates." << endl;
	  terminate();
	}
      indices[i] = IndexAtCoordinates(xys[i][0],xys[i][1]);
    }
    return indices;
}

void point_cloud::ApplyDiggingActionsInPlace(vector<double>& ug)
{
  
  if(ug.size() != points.size())
    {
      cout << "Number of points does not match." << endl;
      terminate();
    }
  for(unsigned int j = 0; j < ug.size(); j ++)
    points[j][2] -= ug[j];
}

point_cloud point_cloud::ApplyDiggingActions(vector<double>& ug)
{
  point_cloud new_point_cloud(*this);
  new_point_cloud.ApplyDiggingActionsInPlace(ug);
  return new_point_cloud;
}

void point_cloud::ApplyDiggingActionsInPlace(vector<vector<double> >& ug_seq)
{  
  for(unsigned int i = 0; i < ug_seq.size(); i++)
    ApplyDiggingActionsInPlace(ug_seq[i]);
}

point_cloud point_cloud::ApplyDiggingActions(vector<vector<double> >& ug_seq)
{
  point_cloud new_point_cloud(*this);
  new_point_cloud.ApplyDiggingActionsInPlace(ug_seq);
  return new_point_cloud;
}

void point_cloud::ApplyDumpingActionsInPlace(vector<vector<double> >& up_seq)
{
  for(unsigned int i = 0; i < up_seq.size(); i++)
    {
      if(up_seq[i].size() != points.size())
	{
	  cout << "Number of points does not match." << endl;
	  terminate();
	}
      for(unsigned int j = 0; j < up_seq[i].size(); j ++)
	points[j][2] += up_seq[i][j];
    }
}

point_cloud point_cloud::ApplyDumpingActions(vector<vector<double> >& up_seq)
{
  point_cloud new_point_cloud(*this);
  new_point_cloud.ApplyDumpingActionsInPlace(up_seq);
  return new_point_cloud;
}

void point_cloud::ApplyDumpingActionsInPlace(vector<double>& up)
{
  for(unsigned int i = 0; i < up.size(); i++)
    {
      if(up.size() != points.size())
	{
	  cout << "Number of points does not match." << endl;
	  terminate();
	}
	points[i][2] += up[i];
    }
}

point_cloud point_cloud::ApplyDumpingActions(vector<double>& up)
{
  point_cloud new_point_cloud(*this);
  new_point_cloud.ApplyDumpingActionsInPlace(up);
  return new_point_cloud;
}

void point_cloud::ApplyDumpingActionsInPlace(vector<vector<vector<double> > >& up_seq)
{
  for(unsigned int i = 0; i < up_seq.size(); i++)
    {
      ApplyDumpingActionsInPlace(up_seq[i]);
    }
}

point_cloud point_cloud::ApplyDumpingActions(vector<vector<vector<double> > >& up_seq)
{
  point_cloud new_point_cloud(*this);
  for(unsigned int i = 0; i < up_seq.size(); i++)
    {
      new_point_cloud.ApplyDumpingActionsInPlace(up_seq[i]);
    }
  return new_point_cloud;
}

void point_cloud::WriteToXYZ(string fname)
{
  WriteToCSV(fname,points);
}

optdata optdata::operator()(unsigned int ind)
{
  optdata result;
  result.d_fg.push_back(this->d_fg[ind]);
  result.d_rp.push_back(this->d_rp[ind]);
  result.t_mdr.push_back(this->t_mdr[ind]);
  result.beta_1.push_back(this->beta_1[ind]);
  result.sa_1.push_back(this->sa_1[ind]);
  result.sa_2.push_back(this->sa_2[ind]);
  result.sa_1_nohoist.push_back(this->sa_1_nohoist[ind]);
  result.sa_2_nohoist.push_back(this->sa_2_nohoist[ind]);
  result.dumping_points.push_back(this->dumping_points[ind]);
  result.d_dumping_points_to_elements.push_back(this->d_dumping_points_to_elements[ind]);
  result.preceding_elem.push_back(this->preceding_elem[ind]);
  result.M_dg = this->M_dg;
  result.spoil_shell = this->spoil_shell;
  result.M_dp = this->M_dp;
  result.max_height_dumping_points.push_back(this->max_height_dumping_points[ind]);
  result.d_dragline.push_back(this->d_dragline[ind]);
  result.alpha_1.push_back(this->alpha_1[ind]);
  result.floor_dumping_points.push_back(this->floor_dumping_points[ind]);
  result.indices_dppt.push_back(this->indices_dppt[ind]);
  return result;
}

const optdata optdata::operator()(unsigned int ind) const
{
  optdata result;
  result.d_fg.push_back(this->d_fg[ind]);
  result.d_rp.push_back(this->d_rp[ind]);
  result.t_mdr.push_back(this->t_mdr[ind]);
  result.beta_1.push_back(this->beta_1[ind]);
  result.sa_1.push_back(this->sa_1[ind]);
  result.sa_2.push_back(this->sa_2[ind]);
  result.sa_1_nohoist.push_back(this->sa_1_nohoist[ind]);
  result.sa_2_nohoist.push_back(this->sa_2_nohoist[ind]);
  result.dumping_points.push_back(this->dumping_points[ind]);
  result.d_dumping_points_to_elements.push_back(this->d_dumping_points_to_elements[ind]);
  result.preceding_elem.push_back(this->preceding_elem[ind]);
  result.M_dg = this->M_dg;
  result.spoil_shell = this->spoil_shell;
  result.M_dp = this->M_dp;
  result.max_height_dumping_points.push_back(this->max_height_dumping_points[ind]);
  result.d_dragline.push_back(this->d_dragline[ind]);
  result.alpha_1.push_back(this->alpha_1[ind]);
  result.floor_dumping_points.push_back(this->floor_dumping_points[ind]);
  result.indices_dppt.push_back(this->indices_dppt[ind]);
  return result;
}

optdata optdata::operator()(vector<unsigned int> indices)
{
  optdata result;
  for(unsigned int i = 0; i < indices.size(); i++)
    {
      result.d_fg.push_back(this->d_fg[indices[i]]);
      result.d_rp.push_back(this->d_rp[indices[i]]);
      result.t_mdr.push_back(this->t_mdr[indices[i]]);
      result.beta_1.push_back(this->beta_1[indices[i]]);
      result.sa_1.push_back(this->sa_1[indices[i]]);
      result.sa_2.push_back(this->sa_2[indices[i]]);
      result.sa_1_nohoist.push_back(this->sa_1_nohoist[indices[i]]);
      result.sa_2_nohoist.push_back(this->sa_2_nohoist[indices[i]]);
      result.dumping_points.push_back(this->dumping_points[indices[i]]);
      result.d_dumping_points_to_elements.push_back(this->d_dumping_points_to_elements[indices[i]]);
      result.preceding_elem.push_back(this->preceding_elem[indices[i]]);
      result.M_dg = this->M_dg;
      result.spoil_shell = this->spoil_shell;
      result.M_dp = this->M_dp;
      result.max_height_dumping_points.push_back(this->max_height_dumping_points[indices[i]]);
      result.d_dragline.push_back(this->d_dragline[indices[i]]);
      result.alpha_1.push_back(this->alpha_1[indices[i]]);
      result.floor_dumping_points.push_back(this->floor_dumping_points[indices[i]]);
      result.indices_dppt.push_back(this->indices_dppt[indices[i]]);
    }
  return result;
}

const optdata optdata::operator()(vector<unsigned int> indices) const
{
  optdata result;
  for(unsigned int i = 0; i < indices.size(); i++)
    {
      result.d_fg.push_back(this->d_fg[indices[i]]);
      result.d_rp.push_back(this->d_rp[indices[i]]);
      result.t_mdr.push_back(this->t_mdr[indices[i]]);
      result.beta_1.push_back(this->beta_1[indices[i]]);
      result.sa_1.push_back(this->sa_1[indices[i]]);
      result.sa_2.push_back(this->sa_2[indices[i]]);
      result.sa_1_nohoist.push_back(this->sa_1_nohoist[indices[i]]);
      result.sa_2_nohoist.push_back(this->sa_2_nohoist[indices[i]]);
      result.dumping_points.push_back(this->dumping_points[indices[i]]);
      result.d_dumping_points_to_elements.push_back(this->d_dumping_points_to_elements[indices[i]]);
      result.preceding_elem.push_back(this->preceding_elem[indices[i]]);
      result.M_dg = this->M_dg;
      result.spoil_shell = this->spoil_shell;
      result.M_dp = this->M_dp;
      result.max_height_dumping_points.push_back(this->max_height_dumping_points[indices[i]]);
      result.d_dragline.push_back(this->d_dragline[indices[i]]);
      result.alpha_1.push_back(this->alpha_1[indices[i]]);
      result.floor_dumping_points.push_back(this->floor_dumping_points[indices[i]]);
      result.indices_dppt.push_back(this->indices_dppt[indices[i]]);
    }
  return result;
}

double sol_digging::UgSum()
{
  double result = 0.0;
  for(unsigned int i = 0; i < ug_seq.size(); i++)
    result += SumVector(ug_seq[i]);
  return result;
}

double sol_dumping::VdpSum()
{
  double result = 0.0;
  for(unsigned int i = 0; i < V_dp_seq.size(); i++)
    for(unsigned int j = 0; j < V_dp_seq[i].size(); j++)
      result += V_dp_seq[i][j];
  return result;
}

double sol_digging_dumping::UgSum()
{
  double result = 0.0;
  for(unsigned int i = 0; i < ug_seq.size(); i++)
    result += SumVector(ug_seq[i]);
  return result;
}

double sol_digging_dumping::VdpSum()
{
  double result = 0.0;
  for(unsigned int i = 0; i < V_dp_seq.size(); i++)
    for(unsigned int j = 0; j < V_dp_seq[i].size(); j++)
      result += V_dp_seq[i][j];
  return result;
}
