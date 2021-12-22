# David R Thompson
import numpy as np

# The columns on either side of the FPA are masked.
last_masked_col_left, first_masked_col_right = 9, 1272

# The first and last represent the extrema of the rows
# containing good data (inclusive, zero-indexed)
first_valid_row, last_valid_row =  6, 333

# Subframes have 328 rows
valid_rows = last_valid_row - first_valid_row + 1

# Rows that are not masked 
first_illuminated_row, last_illuminated_row = 20, 320

# Columns that are not significantly vignetted
first_illuminated_column, last_illuminated_column = 24, 1265

# EMIT FPA size
rows, columns = 480, 1280

# Define masked rows and columns
masked_rows = np.concatenate((np.arange(first_valid_row, first_illuminated_row, dtype=int),
            np.arange(last_illuminated_row+1, last_valid_row+1, dtype=int)),
            axis=0)
masked_cols = np.concatenate((np.arange(0, last_masked_col_left+1, dtype=int),
            np.arange(first_masked_col_right, columns, dtype=int)),
            axis=0)

# These columns used for stray light checks
vignetted_cols = np.concatenate((np.arange(last_masked_col_left+1, first_illuminated_column, dtype=int),
                 np.arange(last_illuminated_column+1, first_masked_col_right, dtype=int)),axis=0)


# EMIT frames can be in native format or in subframe (328 row) format.
# This function extracts a subframe from a native format frame
def frame_extract(frame):
  if frame.shape[1] != columns:
     raise IndexError('All frames should have '+str(columns)+' columns')
  if frame.shape[0] != rows:
     raise IndexError('Native-format frames should have '+str(rows)+' rows')
  return frame[first_valid_row:(last_valid_row+1),:]
 

# EMIT frames can be in native format or in subframe (328 row) format.
# This function makes sure that all frames have native format by 
# embedding subframes inside some padding.
def frame_embed(frame):
  if frame.shape[1] != columns:
     raise IndexError('All frames should have '+str(columns)+' columns')
  if frame.shape[0] == rows:
     return frame
  if frame.shape[0] != valid_rows:
     raise IndexError('Invalid number of rows')
  embedded = np.zeros((rows, columns))
  embedded[first_valid_row, last_valid_row+1] = frame
  return embedded    
