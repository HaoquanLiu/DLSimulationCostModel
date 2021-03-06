#ifndef COST_H_
#define COST_H_

#include <DataStruct.h>

// compute the operational cost considering both digging and dumping with equivalent swing angle online, along with the positioning and preparation cost
cost_sol Cost(sol_digging_dumping&,vector<unsigned int>&,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
cost_sol Cost(sol_digging_dumping&,vector<vector<double> >&,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

#endif 
