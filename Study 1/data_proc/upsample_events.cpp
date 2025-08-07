#include <Rcpp.h>
using namespace Rcpp;

// Function to convert onset and offset times into binary presence at frame level
// Arguments:
//   onset_time: NumericVector containing onset times in seconds
//   offset_time: NumericVector containing offset times in seconds
// Returns:
//   DataFrame: binary presence vector at frame and timestamp levels
// [[Rcpp::export]]
DataFrame convertToFrames(NumericVector onset_time, NumericVector offset_time) {
  // Determine the maximum offset time
  double max_offset = max(offset_time);
  
  // Calculate the total number of frames based on the maximum offset time and frame rate
  int num_frames = ceil(max_offset * 29.97);
  
  // Create an output vector to store the binary presence at frame level
  NumericVector presence(num_frames, 0.0);
  NumericVector timestamp(num_frames, 0.0);
  
  // Loop through the onset and offset times
  for (int i = 0; i < onset_time.length(); i++) {
    // Convert onset and offset times to frame indices
    int start_frame = ceil(onset_time[i] * 29.97);
    int end_frame = ceil(offset_time[i] * 29.97);
    
    // Set the presence to 1 for frames within the onset and offset range
    for (int frame = start_frame; frame <= end_frame; frame++) {
      presence[frame] = 1.0;
    }
  }
  
  // Calculate the timestamp for each frame
  for (int frame = 0; frame < num_frames; frame++) {
    timestamp[frame] = frame / 29.97;
  }
  
  // Create a DataFrame with presence and timestamp columns
  DataFrame output = DataFrame::create(Named("presence") = presence,
                                       Named("timestamp") = timestamp);
  
  // Return the binary presence vector at frame level
  return output;
}