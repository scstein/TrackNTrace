#pragma once

// std lib
#include <cmath>
#include <utility>
#include <algorithm>

// kd-tree library
#include "nanoflann.hpp"

using namespace std;
using namespace nanoflann;


// We need to define this class to work with the nanoflann kd-tree
template <typename T>
struct PointCloud : public Array2D_t<T>
{    
    using Array2D_t<T>::nRows;
    using Array2D_t<T>::nCols;
    using Array2D_t<T>::ptr;
    
     // Use constructor of base class
    PointCloud (const mxArray* mexPtr)
    :     Array2D_t<T>( mexPtr ),
             nr_points(nCols), // Note: nCols and nRows are set by the Array2D constructor
             nr_dimensions(nRows)
    {}
	
	// Use constructor of base class
    PointCloud (T* ptr_arg, int num_rows, int num_cols)
    :     Array2D_t<T>( ptr_arg, num_rows, num_cols ),
             nr_points(nCols), // Note: nCols and nRows are set by the Array2D constructor
             nr_dimensions(nRows)
    {}
    
    // Use constructor of base class. This one allocates a new array in memory!
    PointCloud (int num_rows, int num_cols)
    :     Array2D_t<T>( num_rows, num_cols ),
             nr_points(nCols), // Note: nCols and nRows are set by the Array2D constructor
             nr_dimensions(nRows)
    {}
    
    // Need a constant access function for the kd-tree
    inline T const_at(int iRow, int iCol) const
    {
        return *(ptr + iCol*nRows + iRow);
    }
    
	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return nr_points; }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t size) const
	{
		const T d0 = p1[0] - this->const_at(1, idx_p2);
		const T d1 = p1[1] - this->const_at(2, idx_p2);
		const T d2 = p1[2] - this->const_at(3, idx_p2);
		return d0*d0 + d1*d1 + d2*d2;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		return this->const_at(dim+1, idx);	//Note: dim+1, as our first row contains the frame number, then x,y,... etc.
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const { return false; }
    
    // Returns the distance between the two points the farthest apart
    double maxDistanceInCloud()
    {
        // Find maximum distance of all points
        double xmin = const_at(1,0);
        double xmax = const_at(1,0);
        double ymin = const_at(2,0);
        double ymax = const_at(2,0);
		double zmin = const_at(3,0);
        double zmax = const_at(3,0);
        
        for(unsigned int iP = 1; iP<nr_points; ++iP)
        {
            xmin = (xmin<const_at(1,iP)) ? xmin : const_at(1,iP);
            ymin = (ymin<const_at(2,iP)) ? ymin : const_at(2,iP);
			zmin = (zmin<const_at(3,iP)) ? zmin : const_at(3,iP);
            
            xmax = (xmax>const_at(1,iP)) ? xmax : const_at(1,iP);
            ymax = (ymax>const_at(2,iP)) ? ymax : const_at(2,iP);
			zmax = (zmax>const_at(3,iP)) ? zmax : const_at(3,iP);
        }
        return std::sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) + (zmax-zmin)*(zmax-zmin));
    }
    
    // Attributes
    const unsigned int nr_points;
    const unsigned int nr_dimensions;
};


class NearestNeighborTracker
{
public:
    NearestNeighborTracker(PointCloud<double>& input_pointcloud)
    : input_pts(input_pointcloud)
    {
        nr_frames = findNrFrames();
    }  
   
    // Go through the data frame by frame and perform nearest neighbor search to the next frame.
    // Point pairs farther then "max_linking_dist" appart are neglected.
    void frameToframeLinking(double max_linking_dist);
    
    // Must be called after frameToframeLinking. Builds the "all_tracks" data from the frame2frame
    // data, which contains for each track the global indices of all particle positions.
    // tracks below "min_track_length" are not saved.
    void buildTracks(int min_track_length);
    
    
    // Performs gap closing. Make sure frameToframeLinking and buildTracks is executed first!
    void gapClosing(double max_distance_gap, int max_frame_gap);
    
    
    void removeTracksShorterThan(int min_track_length)
    {
		if(min_track_length <= 1) // every trace has at least one point by nature
			return;

		cout << " Removing tracks shorter than " << min_track_length << " frames (";

            all_tracks.erase(std::remove_if(all_tracks.begin(), 
                                            all_tracks.end(),
                              [&](vector<int>& vec){return vec.size()<min_track_length;}),
                                      all_tracks.end());
          
		cout << all_tracks.size() << " tracks remaining)" << endl;
    }
    
    
    // Gives access to the trackingData
    vector< vector<int> >& getTrackingData()
    {
      return all_tracks;   
    }
    
private:
    // determine the frame span from the input data (maxframe-minframe+1)
    int findNrFrames();
    
    
    // Initialize the tracker. 
    // a) This counts the number of points per frame -> frame_pt_count
    // b) Counts the cumulative number of points before each frame -> global_index_offset
    void trackerInit();
    
    
// Nearest neighbor search connecting points in the two input clouds.
// Input:
//    pts_frame_1 - Source point cloud
//    pts_frame_2 - Target point cloud
//    max_linking_distance - Matches with a higher distance are neglected.
// Output:
//    out_match_IDs - Contains the matching from source to target: out_match_IDs[src_point] == tgt_match_point. If no match exists the value is -1.
void nn_search (PointCloud<double>& pts_frame_1,
                  PointCloud<double>& pts_frame_2,
                  double max_linking_dist,
                  vector<int>& out_match_IDs // Output: For every point in frame 1, ID of nearest neighbor in frame 2
                  );

// Radius search connecting endpoints of tracks with startpoints of other tracks. 
// Connections are only valid if startpoints are later in time than the endpoints
// From all points with the minimal time distance, the one with minimum spatial distance is chosen.
// Input:
//    endpoints - Source point cloud
//    startpoints - Target point cloud
//    max_distance_gap - Matches with a higher distance are neglected.
//    int max_frame_gap - Maximum difference in time for the two points to be accepted
// Output:
//    out_match_IDs - Contains the matching from endpoints to startpoints: out_match_IDs[endpoint] == startpoint. If no match exists the value is -1.
void radius_search_gapClosing (PointCloud<double>& endpoints,
                               PointCloud<double>& startpoints,
                               double max_distance_gap,
                               int max_frame_gap,
                               vector<int>& out_match_IDs // Output: For every point in frame 1, ID of nearest neighbor in frame 2
                              );

// Starting at a point with local index "point" in frame "frame" this walks forward in time 
// through the frame-to-frame data gathering the global indices of all points 
// of the track into -> vector<int>& track_IDs.
// Must be called on the first point of a track to gather the whole track.
// The flag processed[frame][point] is set to true for all points visited.
// Input:
//    frame - the frame to look at
//    point - local index of the starting search point in the given frame
//    match_IDs[frame][point] - Index of the points neighbor within frame frame+1 (-1 for no neighbor)
//    global_index_offset[frame] - How many points are in the data before [frame]
// Input & Output:
//    processed[frame][point] - indicates wether a point was processed
//    track_IDs[trackPoint] - vector saving global indices of all points belonging to this track
void build_track(int frame,
                 int point, 
                 vector< vector<int> >& match_IDs,
                 vector<int> global_index_offset,
                 vector< vector<bool> >& processed,
                 vector<int>& track_IDs);
    
private:
    PointCloud<double>& input_pts;
    int nr_frames;
    vector<int> frame_pt_count;
    vector<int> global_index_offset;
    vector< vector<int> > match_IDs;
    vector< vector<int> > all_tracks;
    
    typedef KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, PointCloud<double> > ,
            PointCloud<double>,
            3 /* dim */
            > my_kd_tree_t;
};








// Go through the data frame by frame and perform nearest neighbor search to the next frame.
    // Point pairs farther then "max_linking_dist" appart are neglected.
    void NearestNeighborTracker::frameToframeLinking(double max_linking_dist)
    {
        trackerInit(); // Initialise the tracker        
        
        cout << " Frame to frame linking.. (max_linking_dist=" << max_linking_dist << ")" << endl;
        match_IDs.clear();
		match_IDs.resize(nr_frames-1);

        int t1_pt_idx = 0;
        int t2_pt_idx = t1_pt_idx + frame_pt_count[0];
        for(int iF=0; iF<(nr_frames-1); ++iF) // For every frame except the last
		{		
				// Map points in frame to pointcloud
			PointCloud<double> pts_frame_t1( &input_pts(0, t1_pt_idx),  input_pts.nRows, frame_pt_count[iF]);
			PointCloud<double> pts_frame_t2( &input_pts(0, t2_pt_idx),  input_pts.nRows, frame_pt_count[iF+1] );    

			if( input_pts(0, t2_pt_idx) - input_pts(0, t1_pt_idx) == 1 ) // Only link data one frame apart
			{
				// Nearest neighbor search
				nn_search(pts_frame_t1, pts_frame_t2, max_linking_dist, match_IDs[iF]);
			}
			else
			{
				match_IDs[iF].resize( pts_frame_t1.nr_points, -1); // set not match for all points
			}


            // Shift startpoints by one frame
            t1_pt_idx += frame_pt_count[iF]; 
            t2_pt_idx += frame_pt_count[iF+1];
        }
    }
    
    // Must be called after frameToframeLinking. Builds the "all_tracks" data from the frame2frame
    // data, which contains for each track the global indices of all particle positions.
    // tracks below "min_track_length" are not saved.
    void NearestNeighborTracker::buildTracks(int min_track_length)
    {
        all_tracks.clear();
        
        cout << " Building tracks from frame-to-frame data.. (min_track_length=" << min_track_length << ")" << endl;
        vector< vector<bool> > processed(nr_frames);
        for(int iF = 0; iF<nr_frames; ++iF)
        {
            processed[iF].resize(frame_pt_count[iF], false); // Initialize as unprocessed
        }
        
        int overall_track_cnt = 0;
        int overall_valid_track_cnt = 0;
        for(int iF = 0; iF<nr_frames; ++iF) // For every frame
        {
            for(int iP=0; iP < frame_pt_count[iF]; ++iP) // For every point within the frame
            {
                // Skip already processed points
                if(processed[iF][iP] == true)
                    continue;
                
                vector<int> track_IDs;
                build_track(iF, iP, match_IDs, global_index_offset, processed, track_IDs);
                
                overall_track_cnt++;
                if(track_IDs.size() >= min_track_length)
                {
                    overall_valid_track_cnt++;
                    all_tracks.push_back(std::move(track_IDs));
                }
            }
        }
        
        cout << " -> Keeping " << overall_valid_track_cnt << " (of " << overall_track_cnt << ") tracks above minimum track length." << endl;
        
    }

    
    // Performs gap closing. Make sure frameToframeLinking and buildTracks is executed first!
    void NearestNeighborTracker::gapClosing(double max_distance_gap, int max_frame_gap)
    {
		cout << " Gap closing.. (max_distance_gap=" << max_distance_gap << ", max_frame_gap=" << max_frame_gap << ")" << endl; 

        // Build Point clouds holding the startpoints and endpoints of tracks
        int nr_tracks = all_tracks.size();
        int nr_dimensions = 4; // Here we copy only frame,x,y,z
        PointCloud<double> startpoints(nr_dimensions, nr_tracks);
        PointCloud<double> endpoints(nr_dimensions, nr_tracks);
        
        for(int iTrack=0; iTrack<nr_tracks; ++iTrack)
        {
            int startpointID = all_tracks[iTrack].front();
            int endpointID   = all_tracks[iTrack].back();
            for(int iParam = 0; iParam<nr_dimensions; ++iParam)
            {
                startpoints(iParam,iTrack) = input_pts(iParam,startpointID);
                endpoints(iParam,iTrack)   = input_pts(iParam,endpointID);
            }
        }
		
		// Perform gap closing
		vector<int> match_IDs;
		vector<bool> toRemove(nr_tracks, false);
		
		radius_search_gapClosing (endpoints, startpoints, max_distance_gap, max_frame_gap, match_IDs);
		
	    // Concatenate matching tracks
		int nr_gaps = 0;
		vector<bool> processed(nr_tracks, false);
		vector< vector<int> > tracks_to_concat(nr_tracks); // saves for every track the indices of tracks to concatenate

		for (int iTrack = 0; iTrack<match_IDs.size(); ++iTrack)
		{
			// Skip already tracks already connected
			if( processed[iTrack] == true)
				continue;
			processed[iTrack] = true;
		   
		   // Check if the current track A has a connected track M. If A->B, than we have to check if M has further connections B->C->.. and so on		   
  		   int match_track = match_IDs[iTrack];
		   while(match_track != -1)
		   {
			   tracks_to_concat[iTrack].push_back(match_track); // save ->B to be connected
			   toRemove[match_track] = true; // mark B's data to be removed after copy into this track
			   ++nr_gaps;

			   // If the matched track B was already processed, it was already gap-connected to all later tracks B->C->D
			   // In that case we need to copy B's connections then stop.
			   if(processed[match_track])
			   {
				   tracks_to_concat[iTrack].insert(tracks_to_concat[iTrack].end(), tracks_to_concat[match_track].begin(), tracks_to_concat[match_track].end());
				   tracks_to_concat[match_track].clear(); // Clear all of B's connections as we dont want to copy any data later
				   break;
			   }
			   processed[match_track] = true;

			   match_track = match_IDs[match_track]; // See if the end of this track was also connected. Exists: B->C ?		   
		   }
	  		  
		}

		// Here we copy for every track 'iTrack' the data of all track IDs saved in tracks_to_concat[iTrack]
		for (int iTrack = 0; iTrack<tracks_to_concat.size(); ++iTrack)
		{		
			// For efficiency reasons, we pre-reserve memory for the incoming new data before copying it
			int new_size = all_tracks[iTrack].size();
			for( int iAddTrack = 0; iAddTrack < tracks_to_concat[iTrack].size(); ++iAddTrack)
			{
			    new_size += all_tracks[tracks_to_concat[iTrack][iAddTrack]].size();
		    }
			all_tracks[iTrack].reserve(new_size);

		   // Add all data to the current track
		   for( int iAddTrack = 0; iAddTrack < tracks_to_concat[iTrack].size(); ++iAddTrack)
		   {
			  // Concatenate the vectors of the tracks
			  all_tracks[iTrack].insert(all_tracks[iTrack].end(), all_tracks[tracks_to_concat[iTrack][iAddTrack]].begin(), all_tracks[tracks_to_concat[iTrack][iAddTrack]].end());
		   }	
		}
		
		// Delete all tracks B whose data was transfered to other tracks (A->B)
		// Note: We have to go from end to start of the vector, because erasing an element invalidates all later indices
		for (int iTrack = nr_tracks-1; iTrack>=0; --iTrack)
		{
			// Remove connected tracks
			if(toRemove[iTrack])
				all_tracks.erase(all_tracks.begin() + iTrack);
		}

		cout << " -> Closed " << nr_gaps << " gaps. Remaining tracks: " << all_tracks.size() << endl;
    }
    
    
    int NearestNeighborTracker::findNrFrames()
    {
        //// Find the number of frames
        //int min_frame = input_pts(0,0);
        //int max_frame = min_frame;
        //for(unsigned int iP = 1; iP<input_pts.nr_points; ++iP)
        //{
        //    if (input_pts(0,iP) < min_frame)
        //        min_frame = input_pts(0,iP);
        //    if (input_pts(0,iP) > max_frame)
        //        max_frame = input_pts(0,iP);
        //}
        //return max_frame-min_frame+1;        
		
		int nr_full_frames = 1;
		int nr_empty_frames = 0;

		int curr_frame = input_pts(0,0);
		for(unsigned int iP = 1; iP<input_pts.nr_points; ++iP)
		{
			int next_frame = input_pts(0,iP);
			if(curr_frame != next_frame)
			{
				++nr_full_frames;
				nr_empty_frames += next_frame-curr_frame-1;
			}
			
			curr_frame = next_frame;
		}

		cout << " Processing " << nr_full_frames << " frames with " << input_pts.nr_points << " points. (" << nr_empty_frames << " frames empty). On average " << input_pts.nr_points/nr_full_frames << " points per frame." << endl;   

		return nr_full_frames;
    }
    
    
    // Initialize the tracker. 
    // a) This counts the number of points per frame -> frame_pt_count
    // b) Counts the cumulative number of points before each frame -> global_index_offset
    void NearestNeighborTracker::trackerInit()
    {
        frame_pt_count.resize(nr_frames);
        global_index_offset.resize(nr_frames);
        
        // --  Find #points per frame iF and #points up to (excluding) frame iF --
        unsigned int pt_idx = 0, check_idx = 0;
        int first_frame = input_pts(0, pt_idx);
		/*int last_frame  = first_frame + nr_frames -1;*/
		
		int curr_frame = first_frame;
        int pt_cnt = 0;        
        for( int iF=0; iF<nr_frames; iF++) // For every frame
        {
            pt_cnt = 0;
            // While point is from the same frame, count up
			
            while(input_pts(0, pt_idx) == curr_frame  &&  pt_idx <  input_pts.nr_points)
            {
                pt_idx++;
                pt_cnt++;
            };

		    frame_pt_count[iF] = pt_cnt;
            global_index_offset[iF] = pt_idx - pt_cnt;

			//// If there are frames without detection, like 1,1,2, (!) 4,5 ..
			//// curr_frame = iF+first_frame  before the loop starts.
			//int next_frame = input_pts(0, pt_idx);
   //         for(int jump=0; jump<(next_frame-curr_frame-1); ++iF, ++jump)
   //         {
			//	frame_pt_count[iF+1] = 0;
			//	global_index_offset[iF+1] = pt_idx- pt_cnt;
   //         };
            

//         cout << "Frame " << iF << "  with points " << pt_cnt << endl;
            curr_frame = input_pts(0, pt_idx);
        };
    }
    
    
// Nearest neighbor search connecting points in the two input clouds.
// Input:
//    pts_frame_1 - Source point cloud
//    pts_frame_2 - Target point cloud
//    max_linking_distance - Matches with a higher distance are neglected.
// Output:
//    out_match_IDs - Contains the matching from source to target: out_match_IDs[src_point] == tgt_match_point. If no match exists, the value is -1.
void NearestNeighborTracker::nn_search (PointCloud<double>& pts_frame_1,
                  PointCloud<double>& pts_frame_2,
                  double max_linking_dist,
                  vector<int>& out_match_IDs // Output: For every point in frame 1, ID of nearest neighbor in frame 2
                  )
{   
  double max_dist_sq = max_linking_dist*max_linking_dist;
    
  // Variables we need in the algorithm
  const unsigned int nr_pts_frame_1 = pts_frame_1.nr_points;
  const unsigned int nr_pts_frame_2 = pts_frame_2.nr_points;
  int k_neigh = 1;   // Number of neighbors to search for
 
  // Storage for current best match and its distance to the point
  vector<double> f2_distances(nr_pts_frame_2, -1);
  vector<int> f2_matches(nr_pts_frame_2, -1);
  
  // construct a kd-tree for frame t+1:
  my_kd_tree_t  tree(3 /*dim*/, pts_frame_2, KDTreeSingleIndexAdaptorParams(30 /* max leaf */) );
  tree.buildIndex();
    
    
    // Process all points in frame t, search nearest neighbor in frame t+1
  size_t nn_index = 0;
  double nn_dist_sq = 0;
  
  for (unsigned int i = 0; i < nr_pts_frame_1; ++i)
  {
      double query_pt[3] = {pts_frame_1(1, i), pts_frame_1(2, i), pts_frame_1(3, i)};
	      
      tree.knnSearch(&query_pt[0], k_neigh, &nn_index, &nn_dist_sq);
      
	  // If the nearest neighbor point was not assigned until now (distance<0) or the query point is a better match (closer)
	  // assign the query point (frame 1) as the nearest neighbor (of the point in frame 2)
      // also check if distance is below max_linking_distance
      if ( ((f2_distances[nn_index] < 0) || (nn_dist_sq < f2_distances[nn_index])) && (nn_dist_sq <=max_dist_sq))
      {
          f2_matches[nn_index] = i;
          f2_distances[nn_index] = nn_dist_sq;
      }
      
  }
    
  // Output vector
  out_match_IDs.clear(); 
  out_match_IDs.resize(nr_pts_frame_1, -1); // Every point in cloud1 gets a match or value -1
  
  for (unsigned int iPt = 0; iPt<nr_pts_frame_2; ++iPt)
  {   
   if (f2_matches[iPt] < 0)
        continue;
   
   out_match_IDs[f2_matches[iPt]] = iPt;
  };
}


// Radius search connecting endpoints of tracks with startpoints of other tracks. 
// Connections are made to the closest startpoint in time within the search radius.
// From all points with the minimal time distance, the one with minimum spatial distance is chosen.
// Input:
//    endpoints - Source point cloud
//    startpoints - Target point cloud
//    max_distance_gap - Matches with a higher distance are neglected.
//    int max_frame_gap - Maximum difference in time for the two points to be accepted
// Output:
//    out_match_IDs - Contains the matching from endpoints to startpoints: out_match_IDs[endpoint] == startpoint. If no match exists the value is -1.
void NearestNeighborTracker::radius_search_gapClosing (PointCloud<double>& endpoints,
                               PointCloud<double>& startpoints,
                               double max_distance_gap,
                               int max_frame_gap,
                               vector<int>& out_match_IDs // Output: For every point in frame 1, ID of nearest neighbor in frame 2
                              )
{
    double max_radius_sq = max_distance_gap*max_distance_gap;
    
    // Variables we need in the algorithm
    const unsigned int nr_endpoints = endpoints.nr_points;
    const unsigned int nr_startpoints = startpoints.nr_points;
    
    // Storage for current best match and its distance to the point
	vector<int> f2_framediff(nr_startpoints,-1);
    vector<int> f2_matches(nr_startpoints, -1);
    
    // construct a kd-tree for frame t+1:
    my_kd_tree_t  tree(3 /*dim*/, startpoints, KDTreeSingleIndexAdaptorParams(30 /* max leaf */) );
    tree.buildIndex();
    
	nanoflann::SearchParams params;
    params.sorted = false; // we don't need the results sorted by their distance, as we prioretize time over spatial distance
	
	// Process all endpoints search closest (in time!) startpoint that pops up later
    // From all points with the minimal time distance, the one with minimum spatial distance is chosen.
	for (unsigned int iPt = 0; iPt < nr_endpoints; ++iPt)
	{
		double query_pt[3] = {endpoints(1, iPt), endpoints(2, iPt), endpoints(3, iPt)};
		int endpoint_frame = endpoints(0, iPt);

		// Perform a radius search with the maximum allowed gap-bridging distance
		std::vector<std::pair<size_t,double> >   ret_matches;
		const size_t nMatches = tree.radiusSearch(&query_pt[0],max_radius_sq, ret_matches, params);

		int min_frame_diff = numeric_limits<int>::max();
        double min_distance = numeric_limits<double>::max();
		int closest_candidate_in_time  = -1;
		// For all points within the radius, find the closest in time
		for (size_t iMatch=0;iMatch<nMatches;iMatch++)
		{
			size_t index = ret_matches[iMatch].first;			
            double nn_distance = ret_matches[iMatch].second;
            
			// Skip those startpoints found that are earlier in time than the endpoint 
			// or where the difference in time is greater than the specified max_frame_gap
			int startpoint_frame = startpoints(0, index);
			int frame_diff = startpoint_frame-endpoint_frame;
			if( (frame_diff<= 0) || (frame_diff-1 > max_frame_gap) )
				continue;

			// 1) Check if this point is closer in time then the ones before
			if(frame_diff<min_frame_diff)
			{
				min_frame_diff = frame_diff;
				closest_candidate_in_time = index;
                min_distance = nn_distance;
			}
            // 2) If has the same time distance as the current best match,
            // check if it is spatially closer
            if( (frame_diff == min_frame_diff) && (nn_distance < min_distance))
            {
//                 min_frame_diff = frame_diff; // already equal
                closest_candidate_in_time = index;
                min_distance = nn_distance;  
            }
		}

		// If no candidate was found, check next point
		if(closest_candidate_in_time == -1)
			continue;

		// If the neighbor point was not assigned until now (frame_diff<0) or the query point is a better match (closer in time)
		// assign the query point (endpoint) as the nearest neighbor (of the startpoint)
		if ( ((f2_framediff[closest_candidate_in_time] < 0) || (min_frame_diff < f2_framediff[closest_candidate_in_time])))
		{
			f2_matches[closest_candidate_in_time] = iPt;
			f2_framediff[closest_candidate_in_time] = min_frame_diff;
		}

	}
    
    // Output vector
    out_match_IDs.clear();
    out_match_IDs.resize(nr_endpoints, -1); // Every endpoint gets a matching startpoint ID or value -1
    
    for (unsigned int iPt = 0; iPt<nr_startpoints; ++iPt)
    {
        if (f2_matches[iPt] < 0)
            continue;
        
        out_match_IDs[f2_matches[iPt]] = iPt;
    };
}


// Starting at a point with local index "point" in frame "frame" this walks forward in time 
// through the frame-to-frame data gathering the global indices of all points 
// of the track into -> vector<int>& track_IDs.
// Must be called on the first point of a track to gather the whole track.
// The flag processed[frame][point] is set to true for all points visited.
// Input:
//    frame - the frame to look at
//    point - local index of the starting search point in the given frame
//    match_IDs[frame][point] - Index of the points neighbor within frame frame+1 (-1 for no neighbor)
//    global_index_offset[frame] - How many points are in the data before [frame]
// Input & Output:
//    processed[frame][point] - indicates wether a point was processed
//    track_IDs[trackPoint] - vector saving global indices of all points belonging to this track
void NearestNeighborTracker::build_track(int frame, int point, vector< vector<int> >& match_IDs, vector<int> global_index_offset, vector< vector<bool> >& processed, vector<int>& track_IDs)
{   
    processed[frame][point] = true; // Mark point as processed
    track_IDs.push_back(point + global_index_offset[frame]); // Save its global index
    
    // Follow each match recursively and keep count using curr_length
    if(frame< match_IDs.size()) // Check if we reached the last frame (which does not possess tracking information)
    {
        int match_pt = match_IDs[frame][point];
        if(match_pt > -1) // if match exists
        {
            build_track(frame+1, match_pt, match_IDs, global_index_offset, processed, track_IDs);
        }
    }
}

