#include <Cost.h>
#include <OptDataComputation.h>
#include <cmath>
#include <numeric>
#include <algorithm>

cost_sol Cost(sol_digging_dumping& sol,vector<unsigned int>& indices,optdata& optdata_all,
	      point_cloud& dl_candidate_pos,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,
	      dl_para& dragline_para,blk_para& block_para)
{
  cost_sol cost;
  cost.move_time = vector<double>(indices.size(),0.0);
  cost.swing_time = vector<double>(indices.size(),0.0);
  optdata optdata_ = optdata_all(indices);

  double k_sw = 1/dragline_para.max_swing_speed;
  double b_sw = dragline_para.max_swing_speed/dragline_para.max_swing_accel;
  double k_hoist = 1/dragline_para.max_hoist_speed;
  double b_hoist = dragline_para.max_hoist_speed/dragline_para.max_hoist_accel;
  double k_equiv = dragline_para.max_hoist_speed/dragline_para.max_swing_speed;
  double b_equiv = dragline_para.max_hoist_speed*dragline_para.max_swing_speed/dragline_para.max_swing_accel-pow(dragline_para.max_hoist_speed,2)/dragline_para.max_hoist_accel;
  double move_time_per_metre = 60/dragline_para.walking_speed;

  for(unsigned int i = 0; i < indices.size(); i++)
    {
      // swinging part in digging
      for(unsigned int j = 0; j < dg_elements_ini.NumberOfPoints(); j++)
  	{
	  if(sol.ug_seq[i][j] != 0.0)
  	    {
  	      if(k_sw*fabs(sol.sa_end_of_fill_points[i][j])+b_sw >= k_hoist*(dl_candidate_pos.Z(indices[i])-sol.height_end_of_fill_points[i][j])+b_hoist) // swing-dependent
		{
		  double equivalent_sa = (dl_candidate_pos.Z(indices[i])-sol.height_end_of_fill_points[i][j]-b_equiv)/k_equiv;
		  if(sol.sa_end_of_fill_points[i][j] < 0.0)
		    cost.swing_time[i] += (k_sw*(equivalent_sa+sol.sa_end_of_fill_points[i][j])+b_sw)*
		      dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize*sol.ug_seq[i][j]*block_para.swell_factor/dragline_para.V_ca;
		  else
		    cost.swing_time[i] += (k_sw*sol.sa_end_of_fill_points[i][j]+b_sw)*dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize*sol.ug_seq[i][j]*block_para.swell_factor/dragline_para.V_ca;

		}
  	      else // hoist-dependent
		{
		  double equivalent_sa = (dl_candidate_pos.Z(indices[i])-sol.height_end_of_fill_points[i][j]-b_equiv)/k_equiv;
		  if(sol.sa_end_of_fill_points[i][j] < 0.0)
		    cost.swing_time[i] += (k_sw*(equivalent_sa+sol.sa_end_of_fill_points[i][j])+b_sw)*
		      dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize*sol.ug_seq[i][j]*block_para.swell_factor/dragline_para.V_ca;
		  else
		    cost.swing_time[i] += (k_sw*equivalent_sa+b_sw)*dg_elements_ini.x_cellsize*dg_elements_ini.y_cellsize*sol.ug_seq[i][j]*block_para.swell_factor/dragline_para.V_ca;
		}
  	    }
  	}
      
      // swinging part in dumping
      vector<double> signed_equivalent_sa(optdata_.dumping_points[i].size(),0.0);
      double h_bucket = dl_candidate_pos.Z(indices[i]); // initial bucket height
      for(unsigned int j = 0; j < optdata_.dumping_points[i].size(); j++) // figure out the cycle dependency for each dumping point first and compute their signed equivalent swing angles
	{
	  if(sol.height_dumping_points[i][j] > h_bucket) // if it is necessary to hoist
	    {
	      if(j >= 1)
		{
		  if(k_hoist*(sol.height_dumping_points[i][j] - h_bucket)+b_hoist <= k_sw*(optdata_.sa_2_nohoist[i][j] - optdata_.sa_2_nohoist[i][j-1])+b_sw) // swing-dependent
		    {
		      signed_equivalent_sa[j] = signed_equivalent_sa[j-1] + optdata_.sa_2_nohoist[i][j] - optdata_.sa_2_nohoist[i][j-1];
		      h_bucket += k_equiv*(optdata_.sa_2_nohoist[i][j] - optdata_.sa_2_nohoist[i][j-1])+b_equiv; 
		    }
		  else // hoist-dependent
		    {
		      double delta_equivalent_sa = (sol.height_dumping_points[i][j]-h_bucket-b_equiv)/k_equiv;
		      signed_equivalent_sa[j] = signed_equivalent_sa[j-1] + delta_equivalent_sa;
		      h_bucket = sol.height_dumping_points[i][j];
		    }
		}
	      else
		{
		  if(k_hoist*(sol.height_dumping_points[i][j] - h_bucket)+b_hoist <= k_sw*(optdata_.sa_2_nohoist[i][j])+b_sw) // swing-dependent
		    {
		      signed_equivalent_sa[j] = optdata_.sa_2_nohoist[i][j];
		      h_bucket += k_equiv*(optdata_.sa_2_nohoist[i][j])+b_equiv; 
		    }
		  else
		    {
		      double delta_equivalent_sa = (sol.height_dumping_points[i][j]-h_bucket-b_equiv)/k_equiv;
		      if(optdata_.sa_2_nohoist[i][j] >= 0.0) // this seems unnecessary if using the boom entering point to the dumping area as the reference point
			signed_equivalent_sa[j] = delta_equivalent_sa;
		      else
			signed_equivalent_sa[j] = delta_equivalent_sa + optdata_.sa_2_nohoist[i][j];
		      h_bucket = sol.height_dumping_points[i][j];
		    }
		}
	    }
	  else
	    {
	      if(j >= 1)
		{
		  signed_equivalent_sa[j] = signed_equivalent_sa[j-1] + optdata_.sa_2_nohoist[i][j] - optdata_.sa_2_nohoist[i][j-1];
		  h_bucket += k_equiv*(optdata_.sa_2_nohoist[i][j] - optdata_.sa_2_nohoist[i][j-1])+b_equiv; 
		}
	      else
		{
		  signed_equivalent_sa[j] = optdata_.sa_2_nohoist[i][j];
		  h_bucket = sol.height_dumping_points[i][j];
		}
	    }
	}

      for(unsigned int j = 0; j < optdata_.dumping_points[i].size(); j++)
	  cost.swing_time[i] += (k_sw*signed_equivalent_sa[j]+b_sw)*(sol.V_dp_seq[i][j]/dragline_para.V_ca);

      if(i != 0)
	  cost.move_time[i] = optdata_.d_dragline[i][indices[i-1]]*move_time_per_metre;
	else
	  cost.move_time[i] = 0;
    }
  cost.tot_op_time = 2*std::accumulate(cost.swing_time.begin(),cost.swing_time.end(),0.0)+std::accumulate(cost.move_time.begin(),cost.move_time.end(),0.0)+600*(indices.size()-1);
  return cost;
}

cost_sol Cost(sol_digging_dumping& sol,vector<vector<double> >& xy_dl,optdata& optdata_all,
	      point_cloud& dl_candidate_pos,point_cloud& dg_elements_ini,point_cloud& dg_elements_target,point_cloud& dp_elements,
	      dl_para& dragline_para,blk_para& block_para)
{
  vector<unsigned int> indices = dl_candidate_pos.IndicesAtCoordinates(xy_dl);
  return Cost(sol,indices,optdata_all,dl_candidate_pos,dg_elements_ini,dg_elements_target,dp_elements,dragline_para,block_para);
}
