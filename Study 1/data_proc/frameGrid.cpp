#include <Rcpp.h>
using namespace Rcpp;

// Function to generate a vector of frames between begin_frame and end_frame
// Arguments:
//   begin_frame: the starting frame value
//   end_frame: the ending frame value
// Returns:
//   NumericVector: a vector containing all the frames between begin_frame and end_frame (inclusive)
// [[Rcpp::export]]
NumericVector generateFrames(int begin_frame, int end_frame) {
  
  // Calculate the number of frames
  int num_frames = end_frame - begin_frame + 1;
  
  // Create an output vector of appropriate size
  NumericVector output(num_frames);
  
  // Fill the output vector with frame values
  for (int i = 0; i < num_frames; i++) {
    output[i] = begin_frame + i;
  }
  
  return output;
}