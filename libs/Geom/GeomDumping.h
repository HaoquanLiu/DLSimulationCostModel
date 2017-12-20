#ifndef GEOMMAXDUMPING_H_
#define GEOMMAXDUMPING_H_

#include <DataStruct.h>
// compute the maximum volume of material that can be dumped with front to back spoiling strategy
// the optmization data for ALL candidate dragline positions are provided for the following four functions
sol_dumping GeomMaxDumping(unsigned int,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_dumping GeomMaxDumping(double,double,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_dumping GeomMaxDumping(vector<unsigned int>&,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);
sol_dumping GeomMaxDumping(vector<vector<double> >&,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&);

sol_dumping GeomDumpGivenVolume(unsigned int,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&, double,sol_dumping&);
sol_dumping GeomDumpGivenVolume(double,double,optdata&,point_cloud&,point_cloud&,point_cloud&,point_cloud&,dl_para&,blk_para&, double,sol_dumping&);
#endif
