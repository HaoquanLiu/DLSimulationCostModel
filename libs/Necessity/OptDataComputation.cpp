#include <OptDataComputation.h>
#include <DataOperation.h>
#include <DataStruct.h>
#include <IO.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <omp.h>

// compute the optimization data for one dragline position, specified in indices
optdata ComputeOptData(unsigned int index,point_cloud& ini_elements,
		       point_cloud& dl_candidate_pos,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,
		       dl_para& dragline_para,blk_para& block_para)
{
  optdata result;
  unsigned int num_dg = dg_elements_ini.NumberOfPoints();
  unsigned int num_dl = dl_candidate_pos.NumberOfPoints();
  unsigned int num_dp = dp_elements.NumberOfPoints();
  
  // compute d_fg, t_mdr and d_dragline
  vector<double> tmp_1(num_dg,0.0);
  vector<double> tmp_2(num_dg,0.0);
  vector<unsigned int> tmp_3;
  for(unsigned int i = 0; i < num_dg; i++)
    {
      tmp_1[i] = fabs(sqrt(pow(fabs(dl_candidate_pos.Y(index)-dg_elements_ini.Y(i)),2)+
			  pow(fabs(dl_candidate_pos.X(index)-dg_elements_ini.X(i)),2))-dragline_para.r_f);
      if(sqrt(pow(tmp_1[i],2)+pow(dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(i),2)) > dragline_para.mdr)
	tmp_2[i] = dg_elements_ini.Z(i);
      else
	tmp_2[i] = dl_candidate_pos.Z(index)+dragline_para.h_f-sqrt(pow(dragline_para.mdr,2)-pow(tmp_1[i],2));

      if(sqrt(pow(abs(dg_elements_ini.Y(i)-dl_candidate_pos.Y(index)),2)+pow(abs(dg_elements_ini.X(i)-dl_candidate_pos.X(index)),2)) > dragline_para.r_f &&
	 tmp_1[i]+dragline_para.r_f <= dragline_para.r)
	tmp_3.push_back(i);
    }
  result.d_fg.push_back(tmp_1); // d_fg //
  result.t_mdr.push_back(tmp_2); // t_mdr //
  result.beta_1.push_back(tmp_3); // beta_1 //

  vector<double> tmp_4(num_dl,0.0);
  for(unsigned int i = 0; i < num_dl; i++)
    tmp_4[i] = sqrt(pow(fabs(dl_candidate_pos.Y(index)-dl_candidate_pos.Y(i)),2)+pow(fabs(dl_candidate_pos.X(index)-dl_candidate_pos.X(i)),2));
  result.d_dragline.push_back(tmp_4); // d_dragline //

  // compute M_dg
  for(unsigned int i = 0; i < num_dg; i++)
    result.M_dg.push_back(dg_elements_ini.Z(i) - dg_elements_target.Z(i));

  // compute the coordinates of the dumping points and the reference point used to compute the reference boom position
  // iterating all sides of the dumping polygon to compute the intersection points between each side and the dumping arc
  vector<vector<double> > intersection_points;
  vector<double> tmp_ref_point(2,0.0);
  for(unsigned int i = 0; i < dp_elements.polygon.size(); i++)
    {
      // we solve a quadratic equation for each side
      if(i != dp_elements.polygon.size()-1)
	{
	  double a = pow(dp_elements.polygon[i][0]-dp_elements.polygon[i+1][0],2) + pow(dp_elements.polygon[i][1]-dp_elements.polygon[i+1][1],2);
	  double b = 2 * ((dp_elements.polygon[i+1][0]-dl_candidate_pos.X(index))*(dp_elements.polygon[i][0]-dp_elements.polygon[i+1][0])+
			  (dp_elements.polygon[i+1][1]-dl_candidate_pos.Y(index))*(dp_elements.polygon[i][1]-dp_elements.polygon[i+1][1]));
	  double c = pow(dp_elements.polygon[i+1][0]-dl_candidate_pos.X(index),2) + pow(dp_elements.polygon[i+1][1]-dl_candidate_pos.Y(index),2) - pow(dragline_para.r,2);

	  if(pow(b,2)-4*a*c < 0)
	    continue;
	  else if(pow(b,2)-4*a*c == 0)
	    {
	      vector<double> tmp_intersection(2,0.0);
	      double t = -b/2/a;
	      if(t >= 0 && t <= 1)
		{
		  tmp_intersection[0] = t*dp_elements.polygon[i][0]+(1-t)*dp_elements.polygon[i+1][0];
		  tmp_intersection[1] = t*dp_elements.polygon[i][1]+(1-t)*dp_elements.polygon[i+1][1];
		  intersection_points.push_back(tmp_intersection);
		}
	    }
	  else
	    {
	      vector<double> tmp_intersection(2,0.0);
	      double t = (-b+sqrt(pow(b,2)-4*a*c))/2/a;
	      if(t >= 0 && t <= 1)
		{
		  tmp_intersection[0] = t*dp_elements.polygon[i][0]+(1-t)*dp_elements.polygon[i+1][0];
		  tmp_intersection[1] = t*dp_elements.polygon[i][1]+(1-t)*dp_elements.polygon[i+1][1];
		  intersection_points.push_back(tmp_intersection);
		}
	      t = (-b-sqrt(pow(b,2)-4*a*c))/2/a;
	      if(t >= 0 && t <= 1)
		{
		  tmp_intersection[0] = t*dp_elements.polygon[i][0]+(1-t)*dp_elements.polygon[i+1][0];
		  tmp_intersection[1] = t*dp_elements.polygon[i][1]+(1-t)*dp_elements.polygon[i+1][1];
		  intersection_points.push_back(tmp_intersection);
		}
	    }
	}
      else
	{
	  double a = pow(dp_elements.polygon[i][0]-dp_elements.polygon[0][0],2) + pow(dp_elements.polygon[i][1]-dp_elements.polygon[0][1],2);
	  double b = 2 * ((dp_elements.polygon[0][0]-dl_candidate_pos.X(index))*(dp_elements.polygon[i][0]-dp_elements.polygon[0][0])+
			  (dp_elements.polygon[0][1]-dl_candidate_pos.Y(index))*(dp_elements.polygon[i][1]-dp_elements.polygon[0][1]));
	  double c = pow(dp_elements.polygon[0][0]-dl_candidate_pos.X(index),2) + pow(dp_elements.polygon[0][1]-dl_candidate_pos.Y(index),2) - pow(dragline_para.r,2);

	  if(pow(b,2)-4*a*c < 0)
	    continue;
	  else if(pow(b,2)-4*a*c == 0)
	    {
	      vector<double> tmp_intersection(2,0.0);
	      double t = -b/2/a;
	      if(t >= 0 && t <= 1)
		{
		  tmp_intersection[0] = t*dp_elements.polygon[i][0]+(1-t)*dp_elements.polygon[0][0];
		  tmp_intersection[1] = t*dp_elements.polygon[i][1]+(1-t)*dp_elements.polygon[0][1];
		  intersection_points.push_back(tmp_intersection);
		}
	    }
	  else
	    {
	      vector<double> tmp_intersection(2,0.0);
	      double t = (-b+sqrt(pow(b,2)-4*a*c))/2/a;
	      if(t >= 0 && t <= 1)
		{
		  tmp_intersection[0] = t*dp_elements.polygon[i][0]+(1-t)*dp_elements.polygon[0][0];
		  tmp_intersection[1] = t*dp_elements.polygon[i][1]+(1-t)*dp_elements.polygon[0][1];
		  intersection_points.push_back(tmp_intersection);
		}
	      t = (-b-sqrt(pow(b,2)-4*a*c))/2/a;
	      if(t >= 0 && t <= 1)
		{
		  tmp_intersection[0] = t*dp_elements.polygon[i][0]+(1-t)*dp_elements.polygon[0][0];
		  tmp_intersection[1] = t*dp_elements.polygon[i][1]+(1-t)*dp_elements.polygon[0][1];
		  intersection_points.push_back(tmp_intersection);
		}
	    }
	}
    }
  
  vector<double> sa_intersection(intersection_points.size(),0.0);
  vector<double> reverse_strip_excavation_direction(2,0.0);
  reverse_strip_excavation_direction[0] = block_para.start_point[0] - block_para.end_point[0];
  reverse_strip_excavation_direction[1] = block_para.start_point[1] - block_para.end_point[1];
  for(unsigned int i = 0; i < intersection_points.size(); i++)
    {
      vector<double> dl2intersection(2,0.0);  
      dl2intersection[0] = intersection_points[i][0] - dl_candidate_pos.X(index);
      dl2intersection[1] = intersection_points[i][1] - dl_candidate_pos.Y(index);
      sa_intersection[i] = SafeCos(Dot(reverse_strip_excavation_direction,dl2intersection)/(Norm(reverse_strip_excavation_direction)*Norm(dl2intersection)));      
    }
  unsigned int min_sa_ind = min_element(sa_intersection.begin(),sa_intersection.end())-sa_intersection.begin();
  tmp_ref_point[0] = intersection_points[min_sa_ind][0];
  tmp_ref_point[1] = intersection_points[min_sa_ind][1];
  result.ref_points.push_back(tmp_ref_point);
  double angle_covered = *max_element(sa_intersection.begin(),sa_intersection.end()) - *min_element(sa_intersection.begin(),sa_intersection.end());
  double delta_angle = 10; // angular distance between adjacent dumping points
  double angle_sum = delta_angle;
  vector<vector<double> > tmp_dumping_points;
  vector<unsigned int> tmp_indices_dppt;
  while(angle_sum <= angle_covered*180/3.1415925)
    {
      double x_dpp = (intersection_points[min_sa_ind][0]-dl_candidate_pos.X(index))*cos(angle_sum*3.1415925/180)+
	(intersection_points[min_sa_ind][1]-dl_candidate_pos.Y(index))*sin(angle_sum*3.1415925/180)+dl_candidate_pos.X(index);
      double y_dpp = -(intersection_points[min_sa_ind][0]-dl_candidate_pos.X(index))*sin(angle_sum*3.1415925/180)+
	(intersection_points[min_sa_ind][1]-dl_candidate_pos.Y(index))*cos(angle_sum*3.1415925/180)+dl_candidate_pos.Y(index);
      vector<double> dpp2dl(2,0.0);
      dpp2dl[0] = x_dpp - dl_candidate_pos.X(index);
      dpp2dl[1] = y_dpp - dl_candidate_pos.Y(index);
      if(Dot(dpp2dl,reverse_strip_excavation_direction) >= 0 && InPolygon(dp_elements.polygon,x_dpp,y_dpp)) // constrain the angle between the connection between the dragline and the dumping point and the excavation direction
	// if(InPolygon(dp_elements.polygon,x_dpp,y_dpp)) // constrain the angle between the connection between the dragline and the dumping point and the excavation direction 
	{
	  vector<double> tmp_point(2,0.0);
	  tmp_point[0] = x_dpp;
	  tmp_point[1] = y_dpp;
	  tmp_indices_dppt.push_back(dp_elements.IndexAtCoordinates(x_dpp,y_dpp));
	  tmp_dumping_points.push_back(tmp_point);  
      	}
      else
      	break;
      angle_sum += delta_angle;
    }
  result.dumping_points.push_back(tmp_dumping_points);
  result.indices_dppt.push_back(tmp_indices_dppt);

  // compute swing angles
  vector<double> centroid_dg = dg_elements_ini.Centroid2D();
  vector<double> centroid_dp = dp_elements.Centroid2D(); // temporarily disabled due to the definition of dumping area
  vector<double> dl2dgcentroid(2,0.0);
  vector<double> dl2dpcentroid(2,0.0);
  vector<double> ref_boom_vec(2,0.0);

  vector<double> dl2dg_vec(2,0.0);
  ref_boom_vec[0] = tmp_ref_point[0] - dl_candidate_pos.X(index);
  ref_boom_vec[1] = tmp_ref_point[1] - dl_candidate_pos.Y(index);
  dl2dgcentroid[0] = centroid_dg[0] - dl_candidate_pos.X(index);
  dl2dgcentroid[1] = centroid_dg[1] - dl_candidate_pos.Y(index);
  dl2dpcentroid[0] = centroid_dp[0] - dl_candidate_pos.X(index);
  dl2dpcentroid[1] = centroid_dp[1] - dl_candidate_pos.Y(index);
  int direction_1 = SwingDir(dl2dgcentroid,dl2dpcentroid); // swing direction from the digging area to the dumping area

  for(unsigned int i = 0; i < num_dg; i++)
    {
      dl2dg_vec[0] = dg_elements_ini.X(i) - dl_candidate_pos.X(index);
      dl2dg_vec[1] = dg_elements_ini.Y(i) - dl_candidate_pos.Y(index);
      int direction_2 = SwingDir(dl2dg_vec,ref_boom_vec);
      tmp_1[i] = (direction_1*direction_2) * SafeCos(Dot(ref_boom_vec,dl2dg_vec)/(Norm(ref_boom_vec)*Norm(dl2dg_vec)))*180/3.1415925;
    }
  result.sa_1_nohoist.push_back(tmp_1);

  // test calculating the equilavent swing angle based on the coincident limit graph relationship
  double swhoi_limit_a = 0.37;
  double swhoi_limit_b = -0.8;
  // check hoist-dependency. The hoist height is assumed to be half of the absolute difference between the element's initial height and its target height
  for(unsigned int i = 0; i < num_dg; i++)
    if(fabs(tmp_1[i])*swhoi_limit_a+swhoi_limit_b < dl_candidate_pos.Z(index)-(dg_elements_ini.Z(i)+dg_elements_target.Z(i))/2.0) 
      tmp_1[i] = (dl_candidate_pos.Z(index)-(dg_elements_ini.Z(i)+dg_elements_target.Z(i))/2.0-swhoi_limit_b)/swhoi_limit_a;
  result.sa_1.push_back(tmp_1);

  // compute the preceding elements
  vector<double> dl2dg_vec_1(2,0.0);
  vector<double> dl2dg_vec_2(2,0.0);

  double delta_r = dg_elements_ini.x_cellsize;
  vector<vector<unsigned int> > tmp_5(num_dg);
  // ofstream ofs;
  // ofs.open("/tmp/test.xyz");
  
  vector<vector<double> > dg_dl_polygon(4,vector<double>(2,0.0));
  double min_x = 0.0;
  double min_y = 0.0;
  double max_x = 0.0;
  double max_y = 0.0;
  if(dg_elements_ini.MinX() <= dl_candidate_pos.MinX())
    min_x = dg_elements_ini.MinX();
  else
    min_x = dl_candidate_pos.MinX();
  if(dg_elements_ini.MinY() <= dl_candidate_pos.MinY())
    min_y = dg_elements_ini.MinY();
  else
    min_y = dl_candidate_pos.MinY();
  if(dg_elements_ini.MaxX() >= dl_candidate_pos.MaxX())
    max_x = dg_elements_ini.MaxX();
  else
    max_x = dl_candidate_pos.MaxX();
  if(dg_elements_ini.MaxY() >= dl_candidate_pos.MaxY())
    max_y = dg_elements_ini.MaxY();
  else
    max_y = dl_candidate_pos.MaxY();
  dg_dl_polygon[0][0] = min_x;
  dg_dl_polygon[0][1] = min_y;
  dg_dl_polygon[1][0] = min_x;
  dg_dl_polygon[1][1] = max_y;
  dg_dl_polygon[2][0] = max_x;
  dg_dl_polygon[2][1] = max_y;
  dg_dl_polygon[3][0] = max_x;
  dg_dl_polygon[3][1] = min_y;
  point_cloud dg_dl_elements = ini_elements.ExtractPointCloud(dg_dl_polygon);// a point_cloud that is used to accelerate the computation of preceding elements
  
  for(unsigned int i = 0; i < num_dg; i++)
    {
      vector<unsigned int>::iterator it = find(result.beta_1[0].begin(),result.beta_1[0].end(),i);
      if(it != result.beta_1[0].end())
	{
	  bool only_affected_by_crest = false;
	  bool preceding_computation_continue = false;
	  // I think this is a better way to iterate the points along the line than the way in DLAuto
	  vector<double> unit_vector(2,0.0);
	  unit_vector[0] = (dg_elements_ini.X(i) - dl_candidate_pos.X(index))/sqrt(pow(dg_elements_ini.X(i) - dl_candidate_pos.X(index),2)+pow(dg_elements_ini.Y(i) - dl_candidate_pos.Y(index),2));
	  unit_vector[1] = (dg_elements_ini.Y(i) - dl_candidate_pos.Y(index))/sqrt(pow(dg_elements_ini.X(i) - dl_candidate_pos.X(index),2)+pow(dg_elements_ini.Y(i) - dl_candidate_pos.Y(index),2));
	  double min_grad = 999; // keep track of the minimum gradient between the fairlead and elements along line connecting the fairlead and the digging element
	  for(double r = 0; r <= result.d_fg[0][i]; r += delta_r)
	    {
	      double x = dl_candidate_pos.X(index) + (r + dragline_para.r_f) * unit_vector[0];
	      double y = dl_candidate_pos.Y(index) + (r + dragline_para.r_f) * unit_vector[1];
	      
	      bool in_poly = InPolygon(dg_elements_ini.polygon,x,y);
	      if(!in_poly)
	      	{
	      	  unsigned int ind_in_dg_dl = dg_dl_elements.IndexAtCoordinates(x,y);
		  
	      	  if(ind_in_dg_dl != 9999)
	      	    {
	      	      //ofs << x << "," << y << "," << dg_dl_elements.Z(ind_in_dg_dl)+dragline_para.h_f << endl;
	      	      if(x != dg_dl_elements.X(ind_in_dg_dl) || y != dg_dl_elements.Y(ind_in_dg_dl))
	      		{
	      		  double tmp_grad = (dl_candidate_pos.Z(index)+dragline_para.h_f-dg_dl_elements.Z(ind_in_dg_dl))/
	      		    sqrt(pow(x-dg_dl_elements.X(ind_in_dg_dl),2)+pow(y-dg_dl_elements.Y(ind_in_dg_dl),2));
	      		  if(tmp_grad < min_grad)
	      		    min_grad = tmp_grad;
	      		}
	      	    }
	      	}
	      
	      if(r >= dragline_para.mdr && !in_poly) // if this happens, it means theorectically the maximum digging depth of the current element is determined by the crest 
		{
		  only_affected_by_crest = true;
		}

	      if(in_poly)
		{
		  unsigned int ind_preceding = dg_elements_ini.IndexAtCoordinates(x,y);
		  if(only_affected_by_crest && !preceding_computation_continue)
		    {
		      if(ind_preceding != 9999 && result.d_fg[0][ind_preceding] <= result.d_fg[0][i])
			{
			  if((dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(i))/result.d_fg[0][i] >=
			     (dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(ind_preceding))/result.d_fg[0][ind_preceding] ||
			     (dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(i))/result.d_fg[0][i] >= min_grad)
			    {
			      result.beta_1[0].erase(it); // definition of beta_1 has been expanded
			      break;
			    }
			  else
			    preceding_computation_continue = true;
			}
		      else if(ind_preceding == 9999)
			continue;
		      else
			break;
		    }
		  else
		    {
		      if((dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(i))/result.d_fg[0][i] >= min_grad)
			{
			  result.beta_1[0].erase(it); // definition of beta_1 has been expanded
			  break;
			}

		      if(ind_preceding != 9999 && result.d_fg[0][ind_preceding] <= result.d_fg[0][i])
			{
			  tmp_5[i].push_back(ind_preceding);
			}		   
		      else if(ind_preceding == 9999)
			continue;
		      else
			break;
		    }
		}
	    }
	  // make sure it is its own preceding element
	  if(find(tmp_5[i].begin(),tmp_5[i].end(),i) == tmp_5[i].end())
	    tmp_5[i].push_back(i);
	}
    }
  result.preceding_elem.push_back(tmp_5);
  //ofs.close();

  // maximum dumping volume of material in the dumping area
  result.M_dp = accumulate(result.M_dg.begin(),result.M_dg.end(),0.0)*block_para.swell_factor*dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize;

  // compute d_rp
  vector<double> tmp_6(num_dp,0.0);
  for(unsigned int i = 0; i < num_dp; i++)
      tmp_6[i] = fabs(sqrt(pow(fabs(dl_candidate_pos.Y(index)-dp_elements.Y(i)),2)+pow(fabs(dl_candidate_pos.X(index)-dp_elements.X(i)),2)));
  result.d_rp.push_back(tmp_6);

  // compute spoil_shell
  vector<double> tmp_spoil_shell(num_dp,0.0);
  vector<double> vector_design_line(2,0.0);
  vector_design_line[0] = block_para.end_point[0]-block_para.start_point[0];
  vector_design_line[1] = block_para.end_point[1]-block_para.start_point[1];
  vector<double> unit_vector_design_line(2,0.0);
  double design_vector_length = sqrt(pow(vector_design_line[0],2)+pow(vector_design_line[1],2));
  unit_vector_design_line[0] =  vector_design_line[0]/ design_vector_length;
  unit_vector_design_line[1] =  vector_design_line[1]/ design_vector_length;
  //currently commented out because using the Geom library does not need the spoil shell, the material flow is contained by max_height_dumping_points------------------------------
  // for(unsigned int i = 0; i < num_dp; i++) 
  //   {
  //     vector<double> tmp_vec(2,0.0);
  //     tmp_vec[0] = dp_elements.X(i)-block_para.start_point[0];
  //     tmp_vec[1] = dp_elements.Y(i)-block_para.start_point[1];
  //     if(Norm(tmp_vec) != 0)
  // 	{
  // 	  double length_perpend_point = Norm(tmp_vec)*Dot(tmp_vec,vector_design_line)/(Norm(tmp_vec)*Norm(vector_design_line));
  // 	  double x_perpend = block_para.start_point[0] + length_perpend_point*unit_vector_design_line[0];
  // 	  double y_perpend = block_para.start_point[1] + length_perpend_point*unit_vector_design_line[1];
  // 	  double h_perpend = ini_elements.Z(ini_elements.IndexAtCoordinates(x_perpend,y_perpend)); // is this value nearly a constant? 0.0 in 3ns27???
  // 	  double horizontal_dist = sqrt(pow(Norm(tmp_vec),2)-pow(length_perpend_point,2));
  // 	  if(tan(block_para.angle_of_repose) * horizontal_dist + h_perpend < dp_elements.Z(i)) 
  // 	    tmp_spoil_shell[i] = dp_elements.Z(i);
  // 	  else
  // 	    tmp_spoil_shell[i] = tan(block_para.angle_of_repose) * horizontal_dist + h_perpend;
  // 	}
  //     else
  // 	tmp_spoil_shell[i] = dp_elements.Z(i);
  //   }
  //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  result.spoil_shell = tmp_spoil_shell;

  // compute the max_height_dumping_points (how to ensure this is consistent with the spoil shell when using it for simulation??)
  // only consider the angle of repose from the spoil toe design line ///////////////////////////////////////
  vector<double> tmp_max_height(tmp_dumping_points.size(),0.0);
  vector<double> tmp_floor(tmp_dumping_points.size(),0.0);
  for(unsigned int i = 0; i < tmp_dumping_points.size(); i++)
    {
      // vector<double> max_heights_dppt_from_side_lines(dp_elements.polygon.size(),0.0);
      // vector<double> h_perpend_from_side_lines(dp_elements.polygon.size(),0.0);
      vector<double> max_heights_dppt_from_side_lines(1,0.0);////////////////////////
      vector<double> h_perpend_from_side_lines(1,0.0);///////////////////////////
      for(unsigned int j = 0; j < 1; j++) // need to consider h_perpend for all sides of the dumping polygon//////////////////////////
	{
	  vector<double> start_point(2,0.0);
	  vector<double> end_point(2,0.0);
	  if(j < dp_elements.polygon.size()-1)
	    {
	      // start_point = dp_elements.polygon[j];
	      // end_point = dp_elements.polygon[j+1];
	      start_point = block_para.start_point;
	      end_point = block_para.end_point;
	    }
	  else
	    {
	      // start_point = dp_elements.polygon[j];
	      // end_point = dp_elements.polygon.front();
	      start_point = block_para.start_point;
	      end_point = block_para.end_point;
	    }
	  vector<double> vector_side_line(2,0.0);
	  vector<double> unit_vector_side_line(2,0.0);
	  vector_side_line[0] = end_point[0] - start_point[0];
	  vector_side_line[1] = end_point[1] - start_point[1];
	  double side_vector_length = sqrt(pow(vector_side_line[0],2)+pow(vector_side_line[1],2));
	  unit_vector_side_line[0] =  vector_side_line[0]/ side_vector_length;
	  unit_vector_side_line[1] =  vector_side_line[1]/ side_vector_length;

	  vector<double> tmp_vec(2,0.0);
	  tmp_vec[0] = tmp_dumping_points[i][0]-start_point[0];
	  tmp_vec[1] = tmp_dumping_points[i][1]-start_point[1];      
	  double length_perpend_point = Norm(tmp_vec)*Dot(tmp_vec,vector_side_line)/(Norm(tmp_vec)*Norm(vector_side_line));
	  double x_perpend = start_point[0] + length_perpend_point*unit_vector_side_line[0];
	  double y_perpend = start_point[1] + length_perpend_point*unit_vector_side_line[1];
	  h_perpend_from_side_lines[j] = ini_elements.Z(ini_elements.IndexAtCoordinates(x_perpend,y_perpend));
	  double horizontal_dist = sqrt(pow(Norm(tmp_vec),2)-pow(length_perpend_point,2));
	  if(tan(block_para.angle_of_repose) * horizontal_dist + h_perpend_from_side_lines[j] <= 
	     dragline_para.max_dumping_height + dl_candidate_pos.Z(index)) // cannot exceed the maximum dumping height of the dragline
	    max_heights_dppt_from_side_lines[j] = tan(block_para.angle_of_repose) * horizontal_dist + h_perpend_from_side_lines[j];
	  else
	    max_heights_dppt_from_side_lines[j] = dragline_para.max_dumping_height + dl_candidate_pos.Z(index);
	}
      tmp_max_height[i] = *min_element(max_heights_dppt_from_side_lines.begin(),max_heights_dppt_from_side_lines.end());
      double h_ini = dp_elements.Z(dp_elements.IndexAtCoordinates(tmp_dumping_points[i][0],tmp_dumping_points[i][1]));
      if(h_perpend_from_side_lines[distance(max_heights_dppt_from_side_lines.begin(),find(max_heights_dppt_from_side_lines.begin(),max_heights_dppt_from_side_lines.end(),tmp_max_height[i]))] >= h_ini)
	tmp_floor[i] = h_ini;
      else
	tmp_floor[i] = 
	  h_perpend_from_side_lines[distance(max_heights_dppt_from_side_lines.begin(),find(max_heights_dppt_from_side_lines.begin(),max_heights_dppt_from_side_lines.end(),tmp_max_height[i]))];
    }
  result.max_height_dumping_points.push_back(tmp_max_height);
  result.floor_dumping_points.push_back(tmp_floor);

  // compute sa_2
  vector<double> dl2dumpingpoint(2,0.0);
  vector<double> tmp_sa_2(tmp_dumping_points.size(),0.0);
  for(unsigned int i = 0; i < tmp_dumping_points.size(); i++)
    {
      dl2dumpingpoint[0] = tmp_dumping_points[i][0] - dl_candidate_pos.X(index);
      dl2dumpingpoint[1] = tmp_dumping_points[i][1] - dl_candidate_pos.Y(index);
      int direction_2 = SwingDir(ref_boom_vec,dl2dumpingpoint);
      tmp_sa_2[i] = (direction_1*direction_2) * SafeCos(Dot(ref_boom_vec,dl2dumpingpoint)/(Norm(ref_boom_vec)*Norm(dl2dumpingpoint)))*180/3.1415925;
    }
  result.sa_2_nohoist.push_back(tmp_sa_2);
  // test calculating the equilavent swing angle based on the coincident limit graph relationship
  // the dragline has to hoist to avoid the spoil pile from dumping actions from previous dumping points at this dragline position (assume the bucket stays straight below the boom tip)
  double hoist_height = 0.0; // hoist height required 
  double equivalent_swing_angle_hoist = 0.0;
  unsigned int ind_dppt_hoist_height = 0;
  vector<double> tmp_sa_2_copy(tmp_sa_2);
  for(unsigned int i = 0; i < tmp_dumping_points.size(); i++)
    {
      // check whether there is need for hoisting first for this dumping point
      if(dl_candidate_pos.Z(index) < tmp_max_height[i])
	{
	  if(tmp_max_height[i]-dl_candidate_pos.Z(index) > hoist_height)
	    {
	      hoist_height = tmp_max_height[i]-dl_candidate_pos.Z(index);
	      // hoist_height = (dl_candidate_pos.Z(index)+dragline_para.max_dumping_height+dp_elements.Z(dp_elements.IndexAtCoordinates(tmp_dumping_points[i][0],tmp_dumping_points[i][1])))/2;
	      ind_dppt_hoist_height = i;
	      if(fabs(tmp_sa_2[i])*swhoi_limit_a+swhoi_limit_b < hoist_height) // hoist-dependent in dumping
		tmp_sa_2[i] = (hoist_height-swhoi_limit_b)/swhoi_limit_a;
	      equivalent_swing_angle_hoist = tmp_sa_2[i];
	    }
	  else
	    {
	      tmp_sa_2[i] = equivalent_swing_angle_hoist + tmp_sa_2_copy[i] - tmp_sa_2_copy[ind_dppt_hoist_height];
	    }
	}
    }
  result.sa_2.push_back(tmp_sa_2);

  // compute d_dumping_points_to_elements
  vector<vector<double> > tmp_d_dppt_to_elem(tmp_dumping_points.size(),vector<double>(num_dp,0.0));
  for(unsigned int i = 0; i < tmp_dumping_points.size(); i++)
    {
      for(unsigned int j = 0; j < num_dp; j++)
	{
	  tmp_d_dppt_to_elem[i][j] = sqrt(pow(tmp_dumping_points[i][0]-dp_elements.X(j),2) + pow(tmp_dumping_points[i][1]-dp_elements.Y(j),2));
	}
    }
  result.d_dumping_points_to_elements.push_back(tmp_d_dppt_to_elem);

  // compute alpha_1
  vector<vector<unsigned int> > tmp_alpha_1(tmp_dumping_points.size(),vector<unsigned int>());
  for(unsigned int i = 0; i < tmp_dumping_points.size(); i++)
    {
      for(unsigned int j = 0; j < num_dp; j++)
	{
	  double h_post_max_dumping = result.max_height_dumping_points.front()[i] - tan(block_para.angle_of_repose) * result.d_dumping_points_to_elements.front()[i][j];
	  if(h_post_max_dumping >= dp_elements.Z(j))
	    {
	      tmp_alpha_1[i].push_back(j);
	    }
	}
    }
  result.alpha_1.push_back(tmp_alpha_1);

  // dumping points pruning-----------------------------------this may not be true in reality, but seems to make sense 
  // double min_sa_dg = *min_element(result.sa_1_nohoist.front().begin(),result.sa_1_nohoist.front().end());
  // double min_sa_dp = *min_element(result.sa_2_nohoist.front().begin(),result.sa_2_nohoist.front().end());
  // if(min_sa_dp < 0.0 && min_sa_dg < 0.0)
  //   {
  //     unsigned int feasible_dppt_ind = 0;
  //     while(result.sa_2_nohoist.front()[feasible_dppt_ind] < 0.0 && result.sa_2_nohoist.front()[feasible_dppt_ind] <= fabs(min_sa_dg))
  // 	feasible_dppt_ind++;
  //     if(feasible_dppt_ind != 0)
  // 	{
  // 	  result.sa_2.front().erase(result.sa_2.front().begin(),result.sa_2.front().begin()+feasible_dppt_ind-1);
  // 	  result.sa_2_nohoist.front().erase(result.sa_2_nohoist.front().begin(),result.sa_2_nohoist.front().begin()+feasible_dppt_ind-1);
  // 	  result.dumping_points.front().erase(result.dumping_points.front().begin(),result.dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	  result.d_dumping_points_to_elements.front().erase(result.d_dumping_points_to_elements.front().begin(),result.d_dumping_points_to_elements.front().begin()+feasible_dppt_ind-1);
  // 	  result.max_height_dumping_points.front().erase(result.max_height_dumping_points.front().begin(),result.max_height_dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	  result.alpha_1.front().erase(result.alpha_1.front().begin(),result.alpha_1.front().begin()+feasible_dppt_ind-1);
  // 	  result.floor_dumping_points.front().erase(result.floor_dumping_points.front().begin(),result.floor_dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	}
  //   }
  // else if(min_sa_dp > 0.0 && min_sa_dg < 0.0)
  //   {
  //     unsigned int feasible_dppt_ind = 0;
  //     while(result.sa_2_nohoist.front()[feasible_dppt_ind] <= fabs(min_sa_dg))
  // 	feasible_dppt_ind++;
  //     if(feasible_dppt_ind != 0)
  // 	{
  // 	  result.sa_2.front().erase(result.sa_2.front().begin(),result.sa_2.front().begin()+feasible_dppt_ind-1);
  // 	  result.sa_2_nohoist.front().erase(result.sa_2_nohoist.front().begin(),result.sa_2_nohoist.front().begin()+feasible_dppt_ind-1);
  // 	  result.dumping_points.front().erase(result.dumping_points.front().begin(),result.dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	  result.d_dumping_points_to_elements.front().erase(result.d_dumping_points_to_elements.front().begin(),result.d_dumping_points_to_elements.front().begin()+feasible_dppt_ind-1);
  // 	  result.max_height_dumping_points.front().erase(result.max_height_dumping_points.front().begin(),result.max_height_dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	  result.alpha_1.front().erase(result.alpha_1.front().begin(),result.alpha_1.front().begin()+feasible_dppt_ind-1);
  // 	  result.floor_dumping_points.front().erase(result.floor_dumping_points.front().begin(),result.floor_dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	}
  //   }
  // else if(min_sa_dp < 0.0 && min_sa_dg > 0.0)
  //   {
  //     unsigned int feasible_dppt_ind = 0;
  //     while(result.sa_2_nohoist.front()[feasible_dppt_ind] <= -min_sa_dg)
  // 	feasible_dppt_ind++;
  //     if(feasible_dppt_ind != 0)
  // 	{
  // 	  result.sa_2.front().erase(result.sa_2.front().begin(),result.sa_2.front().begin()+feasible_dppt_ind-1);
  // 	  result.sa_2_nohoist.front().erase(result.sa_2_nohoist.front().begin(),result.sa_2_nohoist.front().begin()+feasible_dppt_ind-1);
  // 	  result.dumping_points.front().erase(result.dumping_points.front().begin(),result.dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	  result.d_dumping_points_to_elements.front().erase(result.d_dumping_points_to_elements.front().begin(),result.d_dumping_points_to_elements.front().begin()+feasible_dppt_ind-1);
  // 	  result.max_height_dumping_points.front().erase(result.max_height_dumping_points.front().begin(),result.max_height_dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	  result.alpha_1.front().erase(result.alpha_1.front().begin(),result.alpha_1.front().begin()+feasible_dppt_ind-1);
  // 	  result.floor_dumping_points.front().erase(result.floor_dumping_points.front().begin(),result.floor_dumping_points.front().begin()+feasible_dppt_ind-1);
  // 	}
  //   }

  return result;
}

// compute the optimization data for one dragline position, specified in coordinates
optdata ComputeOptData(double x_dl,double y_dl,point_cloud& ini_elements,
		       point_cloud& dl_candidate_pos,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,
		       dl_para& dragline_para,blk_para& block_para)
{
  unsigned int index_dl = dl_candidate_pos.IndexAtCoordinates(x_dl,y_dl);
  return ComputeOptData(index_dl,ini_elements,dl_candidate_pos,dg_elements_ini,dg_elements_target,dp_elements,dragline_para,block_para); 
}

// compute the optimization data for a set of dragline positions, specified in indices
optdata ComputeOptData(vector<unsigned int>& indices,point_cloud& ini_elements,
		       point_cloud& dl_candidate_pos,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,
		       dl_para& dragline_para,blk_para& block_para)
{
  optdata result;
  result.d_fg = vector<vector<double> >(indices.size(),vector<double>()); 
  result.d_rp = vector<vector<double> >(indices.size(),vector<double>());  
  result.t_mdr = vector<vector<double> >(indices.size(),vector<double>()); 
  result.beta_1 = vector<vector<unsigned int> >(indices.size(),vector<unsigned int>()); 
  result.sa_1 = vector<vector<double> >(indices.size(),vector<double>()); 
  result.sa_1_nohoist = vector<vector<double> >(indices.size(),vector<double>()); 
  result.preceding_elem = vector<vector<vector<unsigned int> > >(indices.size(),vector<vector<unsigned int> >()); 
  result.d_dragline = vector<vector<double> >(indices.size(),vector<double>()); 
  result.dumping_points = vector<vector<vector<double> > >(indices.size(),vector<vector<double> >()); 
  result.max_height_dumping_points = vector<vector<double> >(indices.size(),vector<double>()); 
  result.sa_2 = vector<vector<double> >(indices.size(),vector<double>()); 
  result.sa_2_nohoist = vector<vector<double> >(indices.size(),vector<double>()); 
  result.d_dumping_points_to_elements = vector<vector<vector<double> > >(indices.size(),vector<vector<double> >()); 
  result.alpha_1 = vector<vector<vector<unsigned int> > >(indices.size(),vector<vector<unsigned int> >());
  result.floor_dumping_points  = vector<vector<double> >(indices.size(),vector<double>());
  result.indices_dppt = vector<vector<unsigned int> >(indices.size(),vector<unsigned int>());
  result.ref_points = vector<vector<double> >(indices.size(),vector<double>());
  
  #pragma omp parallel for
  for(unsigned int i = 0; i < indices.size(); i++)
    {
      optdata tmp_result = ComputeOptData(indices[i],ini_elements,dl_candidate_pos,dg_elements_ini,dg_elements_target,dp_elements,dragline_para,block_para); 
      result.d_fg[i] = (tmp_result.d_fg.front()); 
      result.d_rp[i] = (tmp_result.d_rp.front()); 
      result.t_mdr[i] = (tmp_result.t_mdr.front());
      result.beta_1[i] = (tmp_result.beta_1.front());
      result.sa_1[i] = (tmp_result.sa_1.front());
      result.sa_1_nohoist[i] = (tmp_result.sa_1_nohoist.front());
      result.preceding_elem[i] = (tmp_result.preceding_elem.front());
      result.M_dg = tmp_result.M_dg;
      result.M_dp = tmp_result.M_dp;
      result.d_dragline[i] = (tmp_result.d_dragline.front());
      result.dumping_points[i] = (tmp_result.dumping_points.front());
      result.spoil_shell = tmp_result.spoil_shell;
      result.max_height_dumping_points[i] = (tmp_result.max_height_dumping_points.front());
      result.sa_2[i] = (tmp_result.sa_2.front());
      result.sa_2_nohoist[i] = (tmp_result.sa_2_nohoist.front());
      result.d_dumping_points_to_elements[i] = (tmp_result.d_dumping_points_to_elements.front());
      result.alpha_1[i] = (tmp_result.alpha_1.front());
      result.floor_dumping_points[i] = (tmp_result.floor_dumping_points.front());
      result.indices_dppt[i] = tmp_result.indices_dppt.front();
      result.ref_points[i] = (tmp_result.ref_points.front());
    }
  return result;
}

// compute the optimization data for a set of dragline positions, specified in coordinates
optdata ComputeOptData(vector<vector<double> >& xy_dl,point_cloud& ini_elements,
		       point_cloud& dl_candidate_pos,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,
		       dl_para& dragline_para,blk_para& block_para)
{
  vector<unsigned int> indices(xy_dl.size(),0);
  for(unsigned int i = 0; i < xy_dl.size(); i++)
    indices[i] = dl_candidate_pos.IndexAtCoordinates(xy_dl[i][0],xy_dl[i][1]);
  return ComputeOptData(indices,ini_elements,dl_candidate_pos,dg_elements_ini,dg_elements_target,dp_elements,dragline_para,block_para);
}

// compute the optimization data for ALL candidate dragline positions
optdata ComputeOptData(point_cloud& ini_elements,point_cloud& dl_candidate_pos,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,
		       dl_para& dragline_para,blk_para& block_para)
{
  optdata result;
  result.d_fg = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>()); 
  result.d_rp = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>());  
  result.t_mdr = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>()); 
  result.beta_1 = vector<vector<unsigned int> >(dl_candidate_pos.NumberOfPoints(),vector<unsigned int>()); 
  result.sa_1 = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>()); 
  result.sa_1_nohoist = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>()); 
  result.preceding_elem = vector<vector<vector<unsigned int> > >(dl_candidate_pos.NumberOfPoints(),vector<vector<unsigned int> >()); 
  result.d_dragline = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>()); 
  result.dumping_points = vector<vector<vector<double> > >(dl_candidate_pos.NumberOfPoints(),vector<vector<double> >()); 
  result.max_height_dumping_points = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>()); 
  result.sa_2 = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>()); 
  result.sa_2_nohoist = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>()); 
  result.d_dumping_points_to_elements = vector<vector<vector<double> > >(dl_candidate_pos.NumberOfPoints(),vector<vector<double> >()); 
  result.alpha_1 = vector<vector<vector<unsigned int> > >(dl_candidate_pos.NumberOfPoints(),vector<vector<unsigned int> >());
  result.floor_dumping_points  = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>());
  result.indices_dppt = vector<vector<unsigned int> >(dl_candidate_pos.NumberOfPoints(),vector<unsigned int>());  
  result.ref_points = vector<vector<double> >(dl_candidate_pos.NumberOfPoints(),vector<double>());  

  #pragma omp parallel for
  for(unsigned int i = 0; i < dl_candidate_pos.NumberOfPoints(); i++)
    {
      optdata tmp_result = ComputeOptData(i,ini_elements,dl_candidate_pos,dg_elements_ini,dg_elements_target,dp_elements,dragline_para,block_para); 
      result.d_fg[i] = (tmp_result.d_fg.front()); 
      result.d_rp[i] = (tmp_result.d_rp.front()); 
      result.t_mdr[i] = (tmp_result.t_mdr.front());
      result.beta_1[i] = (tmp_result.beta_1.front());
      result.sa_1[i] = (tmp_result.sa_1.front());
      result.sa_1_nohoist[i] = (tmp_result.sa_1_nohoist.front());
      result.preceding_elem[i] = (tmp_result.preceding_elem.front());
      result.M_dg = tmp_result.M_dg;
      result.M_dp = tmp_result.M_dp;
      result.d_dragline[i] = (tmp_result.d_dragline.front());
      result.dumping_points[i] = (tmp_result.dumping_points.front());
      result.spoil_shell = tmp_result.spoil_shell;
      result.max_height_dumping_points[i] = (tmp_result.max_height_dumping_points.front());
      result.sa_2[i] = (tmp_result.sa_2.front());
      result.sa_2_nohoist[i] = (tmp_result.sa_2_nohoist.front());
      result.d_dumping_points_to_elements[i] = (tmp_result.d_dumping_points_to_elements.front());
      result.alpha_1[i] = (tmp_result.alpha_1.front());
      result.floor_dumping_points[i] = (tmp_result.floor_dumping_points.front());
      result.indices_dppt[i] = tmp_result.indices_dppt.front();
      result.ref_points[i] = (tmp_result.ref_points.front()); 
    }
  return result;
}
