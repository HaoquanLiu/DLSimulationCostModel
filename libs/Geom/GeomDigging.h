#ifndef GEOMDIGGING_H_
#define GEOMDIGGING_H_

#include <DataStruct.h>

// compute the maximum accessible volume of material for digging (currently not considering dumping and layer excavation constraint)
// the optmization data for ALL candidate dragline positions are provided for the following four functions
sol_digging GeomMaxDigging(unsigned int,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_digging GeomMaxDigging(double x,double y,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_digging GeomMaxDigging(vector<unsigned int>&,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_digging GeomMaxDigging(vector<vector<double> >&,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

// compute the digging actions which finish the excavation of a specified amount of material. These functions achieve this by decreasing the digging angles 
// compared to the max digging actions
sol_digging GeomDigGivenVolume(unsigned int,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&,double,sol_digging&);
sol_digging GeomDigGivenVolume(double,double,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&,double,sol_digging&);

// compute the maximum volume of material available for digging at one dragline position considering the available spoil room (front to back spoiling strategy)
sol_digging_dumping GeomMaxDiggingWithSpoil(unsigned int,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_digging_dumping GeomMaxDiggingWithSpoil(double,double,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_digging_dumping GeomMaxDiggingWithSpoil(vector<unsigned int>&,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_digging_dumping GeomMaxDiggingWithSpoil(vector<vector<double> >&,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

#endif
