#include <Rcpp.h>
using namespace Rcpp;

// Function to convert onset and offset times into binary presence at frame level
// Arguments:
//   onset_time: NumericVector containing onset times in seconds
//   offset_time: NumericVector containing offset times in seconds
// Returns:
//   DataFrame: binary presence dataframe at frame level with timestamp column
DataFrame convertToFrames(NumericVector onset_time, NumericVector offset_time) {
  double max_offset = max(offset_time);
  int num_frames = ceil(max_offset * 29.97);
  
  NumericVector presence(num_frames, 0.0);
  NumericVector timestamp(num_frames, 0.0);
  
  for (int i = 0; i < onset_time.length(); i++) {
    int start_frame = ceil(onset_time[i] * 29.97);
    int end_frame = ceil(offset_time[i] * 29.97);
    
    for (int frame = start_frame; frame <= end_frame; frame++) {
      if (frame >= 0 && frame < num_frames) {
        presence[frame] = 1.0;
      }
    }
  }
  
  for (int frame = 0; frame < num_frames; frame++) {
    timestamp[frame] = frame / 29.97;
  }
  
  DataFrame output = DataFrame::create(Named("presence") = presence,
                                       Named("timestamp") = timestamp);
  
  return output;
}

// Function to subset dataframe, convert times to frames, and return the resulting dataframe
// Arguments:
//   data: input dataframe
//   category: category to filter on
// Returns:
//   DataFrame: resulting dataframe with converted frames and timestamps
DataFrame cat2frame(DataFrame data, String category) {
  
  NumericVector onset_col = data["onset"];
  NumericVector offset_col = data["offset"];
  CharacterVector category_col = data["category"];
  
  std::vector<double> subset_onset;
  std::vector<double> subset_offset;
  
  for (int i = 0; i < category_col.size(); i++) {
    if (category_col[i] == category) {
      subset_onset.push_back(onset_col[i]);
      subset_offset.push_back(offset_col[i]);
    }
  }
  
  NumericVector onset_subset = wrap(subset_onset);
  NumericVector offset_subset = wrap(subset_offset);
  
  DataFrame result = convertToFrames(onset_subset, offset_subset);
  
  return result;
}

// Function to generate temporal profile of changes in one category relative to another category
// Arguments:
//   data: input dataframe
//   alignment_frame: alignment frame type ("onset" or "offset")
//   alignment_category: alignment category string
//   window_start: number of frames before alignment time
//   window_end: number of frames after alignment time
//   target_var: target variable category string
// Returns:
//   DataFrame: alignment times, frame timestamps, frame index, and target presence
// [[Rcpp::export]]
DataFrame generateTemporalProfile(DataFrame data, String alignment_frame, String alignment_category,
                                  int window_start, int window_end, String target_var) {
  NumericVector alignment_times;
  
  if (alignment_frame == "onset") {
    NumericVector onset = data["onset"];
    CharacterVector category = data["category"];
    for (int i = 0; i < onset.size(); i++) {
      if (category[i] == alignment_category) {
        alignment_times.push_back(onset[i]);
      }
    }
  } else if (alignment_frame == "offset") {
    NumericVector offset = data["offset"];
    CharacterVector category = data["category"];
    for (int i = 0; i < offset.size(); i++) {
      if (category[i] == alignment_category) {
        alignment_times.push_back(offset[i]);
      }
    }
  } else {
    Rcpp::stop("Invalid alignment frame type");
  }
  
  NumericVector frame_timestamps;
  IntegerVector align_indices;
  IntegerVector frame_indices;
  
  DataFrame target_frame = cat2frame(data, target_var);
  
  for (int i = 0; i < alignment_times.size(); i++) {
    double align_time = alignment_times[i];
    int start_frame = ceil(align_time * 29.97) + window_start;
    int end_frame = ceil(align_time * 29.97) + window_end;
    
    for (int frame = start_frame; frame <= end_frame; frame++) {
      double timestamp = frame / 29.97;
      frame_timestamps.push_back(timestamp);
      align_indices.push_back(i + 1);
      frame_indices.push_back(frame - start_frame + window_start);
    }
  }
  
  NumericVector target_presence(frame_timestamps.size(), 0.0);
  
  NumericVector target_timestamp = target_frame["timestamp"];
  LogicalVector target_presence_col = target_frame["presence"];
  
  for (int i = 0; i < frame_timestamps.size(); i++) {
    double timestamp = frame_timestamps[i];
    int index = -1;
    
    for (int j = 0; j < target_timestamp.size(); j++) {
      if (std::abs(target_timestamp[j] - timestamp) < 1e-6) {
        index = j;
        break;
      }
    }
    
    if (index != -1) {
      target_presence[i] = target_presence_col[index];
    }
  }
  
  DataFrame result = DataFrame::create(Named("timestamp") = frame_timestamps,
                                       Named("align_index") = align_indices,
                                       Named("frame_index") = frame_indices,
                                       Named("target_presence") = target_presence);
  
  return result;
}