#ifndef OPTDATACOMPUTATION_H_
#define OPTDATACOMPUTATION_H_

#include "DataStruct.h"

// compute the optimization data for one dragline position, specified in indices
optdata ComputeOptData(unsigned int,point_cloud&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

// compute the optimization data for one dragline position, specified in coordinates
optdata ComputeOptData(double,double,point_cloud&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

// compute the optimization data for a set of dragline positions, specified in indices
optdata ComputeOptData(vector<unsigned int>&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

// compute the optimization data for a set of dragline positions, specified in coordinates
optdata ComputeOptData(vector<vector<double> >&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

// compute the optimization data for all candidate dragline positions
optdata ComputeOptData(point_cloud&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

#endif
