#ifndef DATASTRUCT_H_
#define DATASTRUCT_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <exception>

using namespace std;

// class that stores dragline parameters
class dl_para
{
 public:
  double r; // operating radius
  double r_f; // fairlead radius
  double h_f; // fairlead height
  double mdr; // minimum disengaging radius
  double V_ca; // bucket capacity
  double alpha_b; // boom angle (unit in rad)
  double w_buc; // dragline bucket width
  double max_dumping_height;
  double max_swing_speed; // degrees per second
  double max_hoist_speed; // metres per second
  double max_swing_accel; // degrees per second square
  double max_hoist_accel; // metres per second square
  double walking_speed; // metres per minute

  dl_para(double v1,double v2,double v3,double v4,double v5,double v6,double v7, double v8,
	  double v9, double v10, double v11, double v12, double v13)
    {
      r = v1;
      r_f = v2;
      h_f = v3;
      mdr = v4;
      V_ca = v5;
      alpha_b = v6*180/3.1415925;
      w_buc = v7;
      max_dumping_height = v8;
      max_swing_speed = v9;
      max_hoist_speed = v10;
      max_swing_accel = v11;
      max_hoist_accel = v12;
      walking_speed = v13;
    }
 
 private:
  dl_para(); // avoid empty initialization
};

// class that stores the parameters of the block to be excavated
class blk_para
{
 public:
  double swell_factor; 
  double angle_of_repose; // unit in rad
  vector<double> start_point; // starting point of a vector indicating the excavation direction of this block, it is also the spoil shell design line start point
  vector<double> end_point; // end point of a vector indicating the excavation direction of this block, it is also the spoil shell design line end point
  blk_para(double v1,double v2,double v3,double v4, double v5, double v6)
    {
      swell_factor = v1;
      angle_of_repose = v2*3.1415925/180.0;
      start_point.push_back(v3);
      start_point.push_back(v4);
      end_point.push_back(v5);
      end_point.push_back(v6);
    }

 private:
  blk_para(); // avoid empty initialization
};

// class related to the 3D terrain data points
class point_cloud
{
 public:
  point_cloud(const string, double, double); // initialize the point_cloud from a xyz file from a regular grid (or part of a regular grid)
  point_cloud(const string); // initialize the point_cloud from a xyz file (can be random points)
  point_cloud(unsigned int); // initialize the point_cloud with a specified number of points;
  point_cloud(const point_cloud& another_point_cloud) // copy constructor
    {
      x_cellsize = another_point_cloud.x_cellsize;
      y_cellsize = another_point_cloud.x_cellsize;
      points = another_point_cloud.points;
      created_from_polygon_extraction = another_point_cloud.created_from_polygon_extraction;
      polygon = another_point_cloud.polygon;
    }

  ~point_cloud(); // destructor

  vector<double> &operator()(unsigned int ind) // access the point cloud by indices
  {
    if(!points.empty())
      {
	if(ind > points.size())
	  {
	    cout << "Index outside the bound." << endl;
	    terminate();
	  }
      }
    else
      {
	cout << "Point cloud is empty." << endl;
	terminate();
      }
    return points[ind];
  }

  const vector<double> &operator()(unsigned int ind) const // access the point cloud by indices
  {
    if(!points.empty())
      {
	if(ind > points.size())
	  {
	    cout << "Index outside the bound." << endl;
	    terminate();
	  }
      }
    else
      {
	cout << "Point cloud is empty." << endl;
	terminate();
      }
    return points[ind];
  }

  // assignment operator?

  unsigned int NumberOfPoints()
  {
    return points.size();
  }

  void ResetPolygon(); // reset the polygon using the extremum of the x,y coordinates

  void Append(vector<double>&); // append a point to the points

  double X(unsigned int index) // return the x coordinates of the point specified by index
  {
    return points[index][0];
  }

  double Y(unsigned int index) // return the y coordinates of the point specified by index
  {
    return points[index][1];
  }

  double Z(unsigned int index) // return the z coordinates of the point specified by index
  {
    return points[index][2];
  }

  vector<double> XCoordinates(); // return the x coordinates of the points

  vector<double> YCoordinates(); // return the y coordinates of the points

  vector<double> ZCoordinates(); // return the z coordinates of the points

  void AssignXCoordinates(unsigned int, double); // assign a new value of z coordinate for a single point specified by index

  void AssignXCoordinates(const vector<double>&); // assign a new vector of z coordinates

  void AssignYCoordinates(unsigned int, double); // assign a new value of z coordinate for a single point specified by index

  void AssignYCoordinates(const vector<double>&); // assign a new vector of z coordinates

  void AssignZCoordinates(unsigned int, double); // assign a new value of z coordinate for a single point specified by index

  void AssignZCoordinates(const vector<double>&); // assign a new vector of z coordinates

  vector<double> Centroid2D(); // return the centroid of the point cloud on the x-y plane 

  vector<double> Centroid3D(); // return the centroid of the point cloud

  double MinX(); // return the minimum x coordinates in points

  double MaxX(); // return the maximum x coordinates in points

  double MinY(); // return the minimum y coordinates in points

  double MaxY(); // return the maximum y coordinates in points

  double MinZ(); // return the minimum z coordinates in points

  double MaxZ(); // return the maximum z coordinates in points

  void UpdateRegion(point_cloud&); // update the terrain of the current point_cloud using another point_cloud

  point_cloud ExtractPointCloud(vector<vector<double> >); // extract the points within a specified polygon and form another point_cloud

  unsigned int IndexAtCoordinates(const double,const double); // return the index of a point nearest to the specified xy coordinates

  vector<unsigned int> IndicesAtCoordinates(const vector<vector<double> >&); // return the indices of points nearest to a set of xy coordinates

  void ApplyDiggingActionsInPlace(vector<double>&); // apply digging actions on the terrain data in points by manipulating the original data directly

  point_cloud ApplyDiggingActions(vector<double>&); // apply digging actions on the terrain data in points by manipulating the copied of the original data

  void ApplyDiggingActionsInPlace(vector<vector<double> >&); // apply digging actions on the terrain data in points by manipulating the original data directly

  point_cloud ApplyDiggingActions(vector<vector<double> >&); // apply digging actions on the terrain data in points by manipulating the copied of the original data

  void ApplyDumpingActionsInPlace(vector<vector<vector<double> > >&); // apply dumping actions on the terrain data in points by manipulating the original data directly

  point_cloud ApplyDumpingActions(vector<vector<vector<double> > >&); // apply dumping actions on the terrain data in points by manipulating the copied of the original data

  void ApplyDumpingActionsInPlace(vector<vector<double> >&); // apply dumping actions on the terrain data in points by manipulating the original data directly

  point_cloud ApplyDumpingActions(vector<vector<double> >&); // apply dumping actions on the terrain data in points by manipulating the copied of the original data

  void ApplyDumpingActionsInPlace(vector<double>&); // apply dumping actions on the terrain data in points by manipulating the original data directly

  point_cloud ApplyDumpingActions(vector<double>&); // apply dumping actions on the terrain data in points by manipulating the copied of the original data

  void WriteToXYZ(string); // write the terrain data to xyz file (all points)

  double x_cellsize;
  double y_cellsize;
  bool created_from_polygon_extraction;
  vector<vector<double> > polygon;

 private:
  point_cloud(); 
  point_cloud(double,double); // initialize the point_cloud with specified x_cellsize and y_cellsize
  vector<vector<double> > points;
};

// struct that stores the optimization data computed for all possible dragline positions
/* struct optdata_all */
/* { */
/*   vector<vector<double> > d_fg_all; // horizontal distances from the fairlead to the digging elements at all possible dragline positions */
/*   vector<vector<double> > d_fp_all; //  horizontal distances from the fairlead to the dumping elements at all possible dragline positions */
/*   vector<vector<double> > t_mdr_all; // computed terrain heights used in the dragline planning and optimization */
/*   vector<vector<unsigned int> > beta_1_all; // indices of digging elements available for digging determined by r, r_f and block boundary crest */
/*   vector<vector<double> > sa_1_all; // swing angles from the digging elements to the reference boom position */
/*   vector<vector<double> > sa_2_all; // swing angle from the reference boom position to the dumping points */
/*   vector<vector<double> > dumping_points_all; // coordinates of the feasible dumping points */
/*   vector<vector<double> > d_dumping_points_to_elements_all; // horizontal distances from the dumping points to the dumping elements */
/*   vector<double> spoil_shell; // maximum heights of the dumping elements */
/*   vector<vector<vector<unsigned int> > > preceding_elem_all; // indices of the preceding elements that are used to address drag rope bending constraint */
/*   vector<double> M_dg; // maximum digging depths of the digging elements */
/*   double M_dp; // total dumping volume of the block */
/*   vector<vector<double> > max_height_dumping_points_all; // maximum heights of the feasible dumping points */
/*   vector<vector<double> > d_dragline_all; // horizontal distances between different dragline positions */
/* }; */

// struct that stores the optimization data computed for the given sequence of dragline positions
struct optdata
{
  vector<vector<double> > d_fg; // horizontal distances from the fairlead to the digging elements
  vector<vector<double> > d_rp; //  horizontal distances from the rotation centre to the dumping elements
  vector<vector<double> > t_mdr; // computed terrain heights used in the dragline planning and optimization
  vector<vector<unsigned int> > beta_1; // indices of digging elements available for digging determined by r, r_f and block boundary crest
  vector<vector<double> > sa_1; // swing angles from the digging elements to the reference boom position
  vector<vector<double> > sa_2; // swing angle from the reference boom position to the dumping points
  vector<vector<double> > sa_1_nohoist; // swing angles from the digging elements to the reference boom position
  vector<vector<double> > sa_2_nohoist; // swing angle from the reference boom position to the dumping points
  vector<vector<vector<double> > >dumping_points; // coordinates of the feasible dumping points
  vector<vector<vector<double> > >d_dumping_points_to_elements; // horizontal distances from the dumping points to the dumping elements
  vector<double> spoil_shell; // maximum heights of the dumping elements
  vector<vector<vector<unsigned int> > > preceding_elem; // indices of the preceding elements that are used to address drag rope bending constraint
  vector<double> M_dg; // maximum digging depths of the digging elements
  double M_dp; // total dumping volume of the block
  vector<vector<double> > max_height_dumping_points; // maximum heights of the feasible dumping points
  vector<vector<double> > d_dragline; // horizontal distances between different dragline positions
  vector<vector<vector<unsigned int> > > alpha_1; // indices of dumping elements that are possibly affected by dumping at the dumping points
  vector<vector<double> > floor_dumping_points;
  vector<vector<unsigned int> > indices_dppt;
  vector<vector<double> > ref_points; // reference points used to compute the reference boom positions 
  
  optdata operator()(unsigned int);
  const optdata operator()(unsigned int) const; 
  optdata operator()(vector<unsigned int>);
  const optdata operator()(vector<unsigned int>) const; 
};

// struct that stores the inscreased heights of dumping elements and the dumping volumes for the considered dragline positions
struct sol_dumping
{
  vector<vector<vector<double> > > up_seq;
  vector<vector<double> > V_dp_seq;
  vector<vector<double> > height_dumping_points;
  double VdpSum();
};

// struct that stores the decreased heights of digging elements for the considered dragline positions
struct sol_digging
{
  vector<vector<double> > ug_seq;
  vector<vector<double> > sa_end_of_fill_points; // swing angles of the centroids of the end-of-fill points
  vector<vector<double> > height_end_of_fill_points; 
  double UgSum();
};

// struct that stores the results from maximizing the digging volume from ONE dragline position
// considering the available spoil room
struct sol_digging_dumping
{
  vector<vector<double> > ug_seq;
  vector<vector<double> > sa_end_of_fill_points; // swing angles of the centroids of the end-of-fill points
  vector<vector<double> > height_end_of_fill_points; 

  vector<vector<double> > V_dp_seq;
  vector<vector<vector<double> > > up_seq;
  vector<vector<double> > height_dumping_points;

  double UgSum();
  double VdpSum();
};

// struct that stores the costs of the solution
struct cost_sol
{
  double tot_op_time;
  vector<double> swing_time;
  vector<double> move_time;
};

#endif /* DATASTRUCT_H_ */
