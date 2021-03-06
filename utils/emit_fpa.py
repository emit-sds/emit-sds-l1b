# David R Thompson
import numpy as np

# OSF seam channels
osf_seam_positions = ((187,189),)

# Number of basis vectors used to describe EMIT nonlinearity
linearity_nbasis = 2

# The columns on either side of the FPA are masked.
last_masked_col_left, first_masked_col_right = 9, 1272

# The first and last represent the extrema of the rows
# containing good data (inclusive, zero-indexed)
first_valid_row, last_valid_row =  6, 333

# Subframes have 328 rows
valid_rows = last_valid_row - first_valid_row + 1

# Rows that are not masked, columns not significantly
# Vignetted
first_illuminated_row, last_illuminated_row = 20, 320
first_illuminated_column, last_illuminated_column = 24, 1265

# The range of elements that is distributed
first_distributed_column, last_distributed_column = 24, 1265
first_distributed_row, last_distributed_row = 26, 313

# EMIT FPA size
native_rows, native_columns = 480, 1280

# Define masked rows and columns
masked_rows = np.concatenate((np.arange(first_valid_row, first_illuminated_row, dtype=int),
            np.arange(last_illuminated_row+1, last_valid_row+1, dtype=int)),
            axis=0)
masked_cols = np.concatenate((np.arange(0, last_masked_col_left+1, dtype=int),
            np.arange(first_masked_col_right, native_columns, dtype=int)),
            axis=0)

# These columns used for stray light checks
vignetted_cols = np.concatenate((np.arange(last_masked_col_left+1, first_illuminated_column, dtype=int),
                 np.arange(last_illuminated_column+1, first_masked_col_right, dtype=int)),axis=0)


# EMIT frames can be in native format or in subframe (328 row) format.
# This function extracts a subframe from a native format frame
def frame_extract(frame, clip_columns = False):
  if frame.shape[1] != native_columns:
     raise IndexError('All frames should have '+str(native_columns)+' columns')
  if frame.shape[0] != native_rows:
     raise IndexError('Native frames should have '+str(native_rows)+' rows')
  frame = frame[first_valid_row:(last_valid_row+1),:]
  if clip_columns:
      frame = frame[:,first_illuminated_column:(last_illuminated_column+1)]
  return frame
 

# EMIT frames can be in native format or in subframe (328 row) format.
# This function makes sure that all frames have native format by 
# embedding subframes inside some padding.
def frame_embed(frame):
  if frame.shape[1] != native_columns:
     raise IndexError('All frames should have '+str(native_columns)+' columns')
  if frame.shape[0] == native_rows:
     return frame
  if frame.shape[0] != valid_rows:
     raise IndexError('Invalid number of rows')
  embedded = np.zeros((native_rows, native_columns))
  embedded[first_valid_row, last_valid_row+1] = frame
  return embedded    
