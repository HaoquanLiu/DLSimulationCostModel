#include "GeomDigging.h"
#include "GeomDumping.h"
#include <OptDataComputation.h>
#include <DataOperation.h>
#include <cmath>
#include <algorithm>
#include <sys/time.h>

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

sol_digging GeomMaxDigging(unsigned int index,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
			   dl_para& dragline_para,blk_para& block_para)
{
  sol_digging sol;
  optdata optdata_ = optdata_all(index);

  // some useful variables that are computed beforehand
  unsigned int num_dg = dg_elements_ini.NumberOfPoints(); // number of digging elements
  unsigned int num_feasible_dg = optdata_.beta_1.front().size(); // number of feasible digging elements indicated by beta_1
  vector<unsigned int> ind_feasible_dg = optdata_.beta_1.front();

  vector<double> tmp_ug(num_dg,0.0);
  vector<double> tmp_sa_end_of_fill_points(num_dg,0.0);
  vector<double> tmp_height_end_of_fill_points(num_dg,0.0);

  // geometrically compute the maximum volume of material that can be removed from the block
  //#pragma omp parallel for
  for(unsigned int i = 0; i < num_feasible_dg; i++)
    {
      // compute the distance from the element to the fairlead
      double d_rf = sqrt(pow((dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(ind_feasible_dg[i])),2)+pow(optdata_.d_fg.front()[ind_feasible_dg[i]],2));
      if(d_rf > dragline_para.mdr) // ignore the element that is within the MDR
	{
	  // find the indices of the preceding elements
	  vector<unsigned int> ind_preceding = optdata_.preceding_elem.front()[ind_feasible_dg[i]];
	  // compute the gradient between the fairlead and the preceding elements
	  vector<double> grad_mdr;
	  vector<double> grad_outside_mdr;
	  bool mdr_affect = false;
	  //#pragma omp parallel for
	  for(unsigned int j = 0; j < ind_preceding.size(); j++)
	    {
	      if(sqrt(pow((dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(ind_preceding[j])),2)+pow(optdata_.d_fg.front()[ind_preceding[j]],2)) <= dragline_para.mdr) // if the preceding element is within the MDR
		{
		  grad_mdr.push_back((dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(ind_preceding[j]))/optdata_.d_fg.front()[ind_preceding[j]]);
		  mdr_affect = true;
		}
	      else
		{
		  if(!mdr_affect)
		    grad_outside_mdr.push_back((dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(ind_preceding[j]))/optdata_.d_fg.front()[ind_preceding[j]]);
		}
	    }
	  
	  if(mdr_affect) // if there are preceding elements within the MDR, the maximum digging depth is computed using the minimum gradient between
	    // the fairlead and the preceding elements within the MDR
	    {
	      double h_elem_post = dl_candidate_pos.Z(index)+dragline_para.h_f - *min_element(grad_mdr.begin(),grad_mdr.end())*optdata_.d_fg.front()[ind_feasible_dg[i]];
	      if(h_elem_post < dg_elements_ini.Z(ind_feasible_dg[i]))
		{
		  if(h_elem_post > dg_elements_target.Z(ind_feasible_dg[i]))
		    tmp_ug[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_feasible_dg[i]) - h_elem_post;
		  else
		    tmp_ug[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_feasible_dg[i]) - dg_elements_target.Z(ind_feasible_dg[i]);
		}
	    }
	  else // otherwise, the maximum digging depth is computed using the gradient between the fairlead and the element closest to the dragline
	    {
	      double h_elem_post = dl_candidate_pos.Z(index)+dragline_para.h_f - grad_outside_mdr.front()*optdata_.d_fg.front()[ind_feasible_dg[i]];
	      if(h_elem_post < dg_elements_ini.Z(ind_feasible_dg[i]))
		{
		  if(h_elem_post > dg_elements_target.Z(ind_feasible_dg[i]))
		    tmp_ug[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_feasible_dg[i]) - h_elem_post;
		  else
		    tmp_ug[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_feasible_dg[i]) - dg_elements_target.Z(ind_feasible_dg[i]);
		}
	    }
	}
	
    }

  // iteratively correct the digging depth of elements to avoid drag rope bending
  unsigned int max_iter_correction = 1;
  for(unsigned int iter = 0; iter < max_iter_correction; iter++)
    {
      ////#pragma omp parallel for
      for(unsigned int i = 0; i < num_feasible_dg; i++)
	{
	  if(tmp_ug[ind_feasible_dg[i]] != 0.0)
	    {
	      // find the indices of the preceding elements
	      vector<unsigned int> ind_preceding = optdata_.preceding_elem.front()[ind_feasible_dg[i]];
	      // compute the gradient between the fairlead and the preceding elements
	      vector<double> grad_preceding(ind_preceding.size(),0.0);
	      for(unsigned int j = 0; j < ind_preceding.size(); j++)
		{
		  grad_preceding[j] = (dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(ind_preceding[j])+tmp_ug[ind_preceding[j]])/optdata_.d_fg.front()[ind_preceding[j]];
		}

	      if((dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(ind_feasible_dg[i])+tmp_ug[ind_feasible_dg[i]])/optdata_.d_fg.front()[ind_feasible_dg[i]] > 
		 *min_element(grad_preceding.begin(),grad_preceding.end())) // if the digging depth on this element violates the constraint,correct it using the minimum gradient between the fairlead and the preceding elements
		{
		  double h_elem_post = -*min_element(grad_preceding.begin(),grad_preceding.end())*optdata_.d_fg.front()[ind_feasible_dg[i]]+dl_candidate_pos.Z(index)+dragline_para.h_f;
		  if(h_elem_post < dg_elements_ini.Z(ind_feasible_dg[i]))
		    {
		      if(h_elem_post > dg_elements_target.Z(ind_feasible_dg[i]))
			tmp_ug[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_feasible_dg[i])-h_elem_post;
		      else
			tmp_ug[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_feasible_dg[i]) - dg_elements_target.Z(ind_feasible_dg[i]);
		    }
		  else
		    tmp_ug[ind_feasible_dg[i]] = 0;
		}

	    }
	}
    }

  // ofstream ofs_test;
  // ofs_test.open("/tmp/end_of_fill_points.csv");
  for(unsigned int i = 0; i < num_feasible_dg; i++)
    {
      if(tmp_ug[ind_feasible_dg[i]] != 0.0)
	{
	  // find the indices of the preceding elements
	  vector<unsigned int> ind_preceding = optdata_.preceding_elem.front()[ind_feasible_dg[i]];
	  for(unsigned int j = 0; j < ind_preceding.size(); j++)
	    {
	      if(tmp_ug[ind_preceding[j]] != 0.0)
		{
		  if(j > 0)
		    {
		      tmp_sa_end_of_fill_points[ind_feasible_dg[i]] = optdata_.sa_1_nohoist.front()[ind_preceding[j-1]];
		      tmp_height_end_of_fill_points[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_preceding[j-1]);
		      // ofs_test << dg_elements_ini.X(ind_preceding[j-1]) << "," << dg_elements_ini.Y(ind_preceding[j-1]) << "," << dg_elements_ini.Z(ind_preceding[j-1]) << endl;
		    }
		  else
		    {
		      tmp_sa_end_of_fill_points[ind_feasible_dg[i]] = optdata_.sa_1_nohoist.front()[ind_preceding[j]];
		      tmp_height_end_of_fill_points[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_preceding[j]);
		      // ofs_test << dg_elements_ini.X(ind_preceding[j]) << "," << dg_elements_ini.Y(ind_preceding[j]) << "," << dg_elements_ini.Z(ind_preceding[j]) << endl;
		    }
		  break;
		}
	    }
	}
    }
  // ofs_test.close();
  sol.ug_seq.push_back(tmp_ug);
  sol.sa_end_of_fill_points.push_back(tmp_sa_end_of_fill_points);
  sol.height_end_of_fill_points.push_back(tmp_height_end_of_fill_points);
  return sol;
}


sol_digging GeomMaxDigging(double x,double y,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
			   dl_para& dragline_para,blk_para& block_para)
{
  unsigned int index = dl_candidate_pos.IndexAtCoordinates(x,y);
  return GeomMaxDigging(index,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para);
}


sol_digging GeomMaxDigging(vector<unsigned int>& indices,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
			   dl_para& dragline_para,blk_para& block_para)
{
  sol_digging sol;
  sol.ug_seq.assign(indices.size(),vector<double>(dg_elements_ini.NumberOfPoints(),0.0));
  point_cloud tmp_dg_elements_ini(dg_elements_ini); // currently the influence of dumping is not considered
  for(unsigned int i = 0; i < indices.size(); i++)
    {
      sol_digging tmp_sol = GeomMaxDigging(indices[i],optdata_all,tmp_dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para);
      sol.ug_seq[i] = tmp_sol.ug_seq.front();      
      sol.sa_end_of_fill_points[i] = tmp_sol.sa_end_of_fill_points.front();      
      sol.height_end_of_fill_points[i] = tmp_sol.height_end_of_fill_points.front();      
      tmp_dg_elements_ini.ApplyDiggingActionsInPlace(tmp_sol.ug_seq.front());
    }
  return sol;
}

sol_digging GeomMaxDigging(vector<vector<double> >& xy_dl,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
			   dl_para& dragline_para,blk_para& block_para)
{
  vector<unsigned int> indices(xy_dl.size(),0);
  for(unsigned int i = 0; i < xy_dl.size(); i++)
    indices[i] = dl_candidate_pos.IndexAtCoordinates(xy_dl[i][0],xy_dl[i][1]);
  return GeomMaxDigging(indices,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para);
}

sol_digging GeomDigGivenVolume(unsigned int index,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
			       dl_para& dragline_para,blk_para& block_para,double V_given,sol_digging& sol_max)
{
  sol_digging sol;
  optdata optdata_ = optdata_all(index);

  // some useful variables that are computed beforehand
  unsigned int num_dg = dg_elements_ini.NumberOfPoints(); // number of digging elements
  unsigned int num_feasible_dg = optdata_.beta_1.front().size(); // number of feasible digging elements indicated by beta_1
  vector<unsigned int> ind_feasible_dg = optdata_.beta_1.front();

  vector<double> tmp_ug(num_dg,0.0);
  vector<double> tmp_sa_end_of_fill_points(num_dg,0.0);
  vector<double> tmp_height_end_of_fill_points(num_dg,0.0);

  double lb = 0.0;
  double ub = 1.0;
  double current_V_dg = sol_max.UgSum()*dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize;
  double current_val = 1.0;
  double tol = dragline_para.V_ca;
    
  vector<double> angle_max(num_feasible_dg,0.0);
  for(unsigned int i = 0; i < num_feasible_dg; i++)
    {
      if(SafeZero(sol_max.ug_seq.front()[ind_feasible_dg[i]]) > 0.0)
	{
	  angle_max[i] = atan((dl_candidate_pos.Z(index)+dragline_para.h_f-dg_elements_ini.Z(ind_feasible_dg[i])+sol_max.ug_seq.front()[ind_feasible_dg[i]])/optdata_.d_fg.front()[ind_feasible_dg[i]]);
	}
    }

  while(fabs(current_V_dg-V_given) >= tol)
    {
      if(current_V_dg >= V_given)
	{
	  ub = current_val;
	  current_val = (current_val+lb)/2;
	}
      else
	{
	  lb = current_val;
	  current_val = (current_val+ub)/2;
	}

      current_V_dg = 0.0;
      for(unsigned int i = 0; i < num_feasible_dg; i++)
	{
	  if(SafeZero(sol_max.ug_seq.front()[ind_feasible_dg[i]]) > 0.0)
	    {
	      double h_post_digging = dl_candidate_pos.Z(index)+dragline_para.h_f-optdata_.d_fg.front()[ind_feasible_dg[i]]*tan(current_val*angle_max[i]);
	      if(h_post_digging < dg_elements_ini.Z(ind_feasible_dg[i]))
		{
		  if(h_post_digging < dg_elements_target.Z(ind_feasible_dg[i]))
		    tmp_ug[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_feasible_dg[i]) - dg_elements_target.Z(ind_feasible_dg[i]);
		  else
		    tmp_ug[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_feasible_dg[i]) - h_post_digging;

		  current_V_dg += tmp_ug[ind_feasible_dg[i]];
		}
	    }
	}
      current_V_dg *= dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize;
    }

  for(unsigned int i = 0; i < num_feasible_dg; i++)
    {
      if(tmp_ug[ind_feasible_dg[i]] != 0.0)
	{
	  // find the indices of the preceding elements
	  vector<unsigned int> ind_preceding = optdata_.preceding_elem.front()[ind_feasible_dg[i]];
	  for(unsigned int j = 0; j < ind_preceding.size(); j++)
	    {
	      if(tmp_ug[ind_preceding[j]] != 0.0)
		{
		  if(j > 0)
		    {
		      tmp_sa_end_of_fill_points[ind_feasible_dg[i]] = optdata_.sa_1_nohoist.front()[ind_preceding[j-1]];
		      tmp_height_end_of_fill_points[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_preceding[j-1]);
		    }
		  else
		    {
		      tmp_sa_end_of_fill_points[ind_feasible_dg[i]] = optdata_.sa_1_nohoist.front()[ind_preceding[j]];
		      tmp_height_end_of_fill_points[ind_feasible_dg[i]] = dg_elements_ini.Z(ind_preceding[j]);
		    }
		  break;
		}
	    }
	}
    }
  sol.sa_end_of_fill_points.push_back(tmp_sa_end_of_fill_points);
  sol.height_end_of_fill_points.push_back(tmp_height_end_of_fill_points);
  sol.ug_seq.push_back(tmp_ug);
  return sol;
}

sol_digging GeomDigGivenVolume(double x, double y,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
			       dl_para& dragline_para,blk_para& block_para,double V_given,sol_digging& sol_max)
{
  unsigned int index = dl_candidate_pos.IndexAtCoordinates(x,y);
  return GeomDigGivenVolume(index,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para,V_given,sol_max); 
}

sol_digging_dumping GeomMaxDiggingWithSpoil(unsigned int index,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
					    dl_para& dragline_para,blk_para& block_para)
{
  sol_digging_dumping sol;
  // compute the maximum digging volume with MILP and not considering dumping
  sol_digging sol_max_digging = GeomMaxDigging(index,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para);
  // compute the maximum dumping volume allowed with Geom
  sol_dumping sol_max_dumping = GeomMaxDumping(index,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para);
      
  if(sol_max_digging.UgSum()*dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize*block_para.swell_factor <= sol_max_dumping.VdpSum()) // sufficient spoil room
    {
      //      cout << "Sufficient spoil room." << endl;
      sol_dumping sol_actual_dumping = GeomDumpGivenVolume(index,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para,
							   sol_max_digging.UgSum()*dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize*block_para.swell_factor,sol_max_dumping);

      sol.ug_seq = sol_max_digging.ug_seq;
      sol.sa_end_of_fill_points = sol_max_digging.sa_end_of_fill_points;
      sol.height_end_of_fill_points = sol_max_digging.height_end_of_fill_points;
      sol.up_seq = sol_actual_dumping.up_seq;
      sol.V_dp_seq = sol_actual_dumping.V_dp_seq;
      sol.height_dumping_points = sol_actual_dumping.height_dumping_points;
    }
  else // not sufficient spoil room, compute the actual digging actions based on the given digging volume
    {
      //      cout << "Inufficient spoil room." << endl;
      sol.up_seq = sol_max_dumping.up_seq;
      sol.V_dp_seq = sol_max_dumping.V_dp_seq;
      sol.height_dumping_points = sol_max_dumping.height_dumping_points;
      sol_digging sol_actual_digging = GeomDigGivenVolume(index,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para,
							  sol_max_dumping.VdpSum()/block_para.swell_factor,sol_max_digging);
      sol.ug_seq = sol_actual_digging.ug_seq;
      sol.sa_end_of_fill_points = sol_actual_digging.sa_end_of_fill_points;
      sol.height_end_of_fill_points = sol_actual_digging.height_end_of_fill_points;
    }
  return sol;
}

sol_digging_dumping GeomMaxDiggingWithSpoil(double x_dl, double y_dl,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos, dl_para& dragline_para,blk_para& block_para)
{
  unsigned int index = dl_candidate_pos.IndexAtCoordinates(x_dl,y_dl);
  return GeomMaxDiggingWithSpoil(index,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para);
}

sol_digging_dumping GeomMaxDiggingWithSpoil(vector<unsigned int>& indices,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
					    dl_para& dragline_para,blk_para& block_para)
{
  sol_digging_dumping sol;
  point_cloud tmp_dg_elements_ini(dg_elements_ini); 
  point_cloud tmp_dp_elements(dp_elements); 
  for(unsigned int i = 0; i < indices.size(); i++)
    {
      sol_digging_dumping tmp_sol = GeomMaxDiggingWithSpoil(indices[i],optdata_all,tmp_dg_elements_ini,dg_elements_target,tmp_dp_elements,dl_candidate_pos,dragline_para,block_para);
      sol.ug_seq.push_back(tmp_sol.ug_seq.front());
      sol.sa_end_of_fill_points.push_back(tmp_sol.sa_end_of_fill_points.front());
      sol.height_end_of_fill_points.push_back(tmp_sol.height_end_of_fill_points.front());
      sol.up_seq.push_back(tmp_sol.up_seq.front());
      sol.V_dp_seq.push_back(tmp_sol.V_dp_seq.front());
      sol.height_dumping_points.push_back(tmp_sol.height_dumping_points.front());
      tmp_dg_elements_ini.ApplyDiggingActionsInPlace(tmp_sol.ug_seq.front());
      tmp_dp_elements.ApplyDumpingActionsInPlace(tmp_sol.up_seq.front());
    }
  return sol;
}

sol_digging_dumping GeomMaxDiggingWithSpoil(vector<vector<double> >& xy_dl,optdata& optdata_all,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,point_cloud& dl_candidate_pos,
					    dl_para& dragline_para,blk_para& block_para)
{
  vector<unsigned int> indices(xy_dl.size(),0);
  for(unsigned int i = 0; i < indices.size(); i++)
      indices[i] = dl_candidate_pos.IndexAtCoordinates(xy_dl[i][0],xy_dl[i][1]);

  return GeomMaxDiggingWithSpoil(indices,optdata_all,dg_elements_ini,dg_elements_target,dp_elements,dl_candidate_pos,dragline_para,block_para);
}
