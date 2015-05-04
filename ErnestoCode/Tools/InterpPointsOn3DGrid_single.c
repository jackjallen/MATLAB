/*
  o = InterpPointsOn3DGrid( IM , X , Y , Z , [xyz1;xyz2...] )
                by default mode is 'linear'
 
  o = InterpPointsOn3DGrid( IM , X , Y , Z , [xyz1;xyz2...] , mode )
          mode = 'nearest'
                 '*nearest' (assuming regular grid)
                 'linear'
                 '*linear'  (assuming regular grid)

  o = InterpPointsOn3DGrid( IM , X , Y , Z , [xyz1;xyz2...] , mode , extrapolation_val )
          extrapolation_val is by default NaN

*/

#define   real       float
#define   mxREAL_CLASS     mxSINGLE_CLASS

#include "InterpPointsOn3DGrid.c"
