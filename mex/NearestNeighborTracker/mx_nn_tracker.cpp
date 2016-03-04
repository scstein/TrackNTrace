// mex
#include "mex.h"

// stdlib
#include <iostream> // cout
#include <vector>
#include <limits>

// our headers
#include "mexUtil.h"
#include "NearestNeighborTracker.h"

using namespace std;


// -- Type definitions -- //
typedef Array1D_t<double> Array1D;

// --  Global Variables  -- //
mxArray** plhs;  // Access to left hand side is global, so we can output to it in functions


// WARNING: The functions treat inherently integer parameters (like the frame number) as integers,
// although they come from a MATLAB double matrix. This might introduce problems, but as a 64-bit
// double captures all integer numbers up to 2^52 exact it's probably ok in practice..

/*
% Performs gap-closing nearest neighbor tracking of 2d point data recorded over time.
% Gap closing connects endpoints of tracks to the next startpoint of a
% track popping up within the specified distance.

% SYNTAX [ tracks ] = nn_tracker_cpp( Localizations, min_track_length,max_linking_distance,max_distance_gap,max_frame_gap,min_track_length_afterGapClosing,verbose )
%
% Input:
%     Localizations - (4+n)xN matrix of points, where N is the number of points 
%                     and the rows correspond to (frame,x,y,z + n additional data). 
%                     Data must be sorted with increasing frame number. 
%     min_track_length - Tracks with less then min_track_length connected
%                        positions are removed before gap closing. | default: 0
%     max_linking_distance - Points farther away then this distance will
%                            not be linked. | default: inf
%     max_distance_gap - Max distance to close gaps over | default: 0
%     max_frame_gap - Max distance in time to connect gaps | default: 0
%     min_track_length_afterGapClosing - Tracks shorter than this length
%                                        are rejected after gap closing | default: 0
%     verbose: If true, information is printed to the console | default: false
%
%  Parameters can be left empty [] to use their default values.
%
% Output:
%     tracks - (5+n)xP matrix of points, where P is the number of points in
%              the tracks and the rows correspond to (trackID, frame,x,y,z + n 
%              additional data from input). 
%
%  Example: to get the data for track 1 and plot its y-x trajectory:
%     track1data = tracks(tracks(1,:)==1, :);
%     plot(track1data(3,:), track1data(4,:));
%  To plot the y-t over time movement
%     plot(track1data(2,:), track1data(4,:));

% Make sure logicals are passed as correct datatype
%
% C++ implementation is written using the nanoflann library. Licensing see below.
*/
void mexFunction(int nlhs, mxArray* plhs_arg[], int nrhs, const mxArray* prhs[])
{	    

    // Access to left hand side is global, so we can output to it in functions
    plhs = plhs_arg;

	 //-- Check for proper number of arguments.
 	if(nrhs<1) {
    	mexErrMsgTxt("One input required.");
  	} 
    if(nlhs!=1) {
    	mexErrMsgTxt("One output required.");
  	}    
    
	// -- Input Parsing --
    PointCloud<double> input_pts( prhs[0] );
    if(input_pts.nRows < 4)
    {
        mexErrMsgTxt("Dimension mismatch: Input matrix must be 3xN, where N is the number of points and columns are Frame,x,y,z!");
    }    
    
    // Default values for non given arguments
    int    min_track_length = 0;
    double max_linking_dist = numeric_limits<double>::infinity();
    
    double max_distance_gap = 0;
	int    max_frame_gap    = 0;   
	int    min_track_length_afterGapClosing = 0;  
	
	bool verbose = false;

	// -- Input parsing -- //
    if(nrhs>1)
    {
        Array1D inputParam(prhs[1]);
        if(!inputParam.isEmpty())
            min_track_length = inputParam[0];
    }
    
    if(nrhs>2)
    {
        Array1D inputParam(prhs[2]);
        if(!inputParam.isEmpty())
        {
            max_linking_dist = inputParam[0];
        }        
    }
    
    if(nrhs>3)
    {
        Array1D inputParam(prhs[3]);
        if(!inputParam.isEmpty())
        {
            max_distance_gap = inputParam[0];
        }
    }

	if(nrhs>4)
    {
        Array1D inputParam(prhs[4]);
        if(!inputParam.isEmpty())
        {
            max_frame_gap = inputParam[0];
        }
    }

	if(nrhs>5)
    {
        Array1D inputParam(prhs[5]);
        if(!inputParam.isEmpty())
        {
            min_track_length_afterGapClosing = inputParam[0];
        }
    }

	if(nrhs>6)
    {
        Array1D inputParam(prhs[6]);
        if(!inputParam.isEmpty())
        {
            verbose = mxIsLogicalScalarTrue(prhs[6]);
        }
    }
    
	bool useGapClosing = ( (max_distance_gap > 0) && (max_frame_gap > 0) );

    // If user wants output, redirect std::cout to MATLAB console
    // MUST be undone by "std::cout.rdbuf(outbuf);" at end of program
    mexStream mout;
    nullStream nullOut;
    std::streambuf *outbuf;
	if(verbose)
	{
	  outbuf = std::cout.rdbuf(&mout);  //Replace the std stream with the 'matlab' stream
	}
    else
    {
      outbuf = std::cout.rdbuf(&nullOut);   
    }

    
    // -- Tracking -- //
    NearestNeighborTracker NNtracker(input_pts);    
    
    NNtracker.frameToframeLinking(max_linking_dist);
    NNtracker.buildTracks(min_track_length);
	if(useGapClosing)
	{
		NNtracker.gapClosing(max_distance_gap,max_frame_gap);    		
	}
    NNtracker.removeTracksShorterThan(min_track_length_afterGapClosing);

    vector< vector<int> > all_tracks = NNtracker.getTrackingData();
    
    
    // -- Copy Result to Matlab output -- //
    int nr_data_points = 0;
    for(int iTrack=0; iTrack<all_tracks.size(); ++iTrack) nr_data_points += all_tracks[iTrack].size();
    int nr_params = input_pts.nRows + 1; // one extra for the track ID
            
     mwSize dim_out[2] = { nr_params, nr_data_points };
     plhs[0] = mxCreateNumericArray( 2, dim_out , mxDOUBLE_CLASS, mxREAL);
     
     Array2D_t<double> mx_track_data ( plhs[0] );    // Map access to output     

	 int point_cnt = 0;
	 // For all tracks
	 for(int iTrack=0; iTrack<all_tracks.size(); ++iTrack)
	 {   // For all points within the track
		 for(int iP = 0; iP<all_tracks[iTrack].size(); ++iP)
		 {
			 int pointID = all_tracks[iTrack][iP];
			 // Write the track ID in the first row ...
			 mx_track_data(0,point_cnt) = iTrack+1; // We use MATLABs starting at 1 convention for the trackID
			 // ... then copy all parameters from the input file
			 for(int iParam = 0; iParam<input_pts.nRows; ++iParam)
			 {
				 mx_track_data(iParam+1,point_cnt) = input_pts(iParam,pointID);
			 }
			 point_cnt++;
		 }
	 }

     
     // restore standard stream
    std::cout.rdbuf(outbuf);
}






