# David R Thompson
import numpy as np
import os, os.path, sys
import json


class FPA:

  def __init__(self, filepath=None):

      if filepath is None:
          raise ValueError('Configuration file required')
          #basedir = os.path.abspath(os.path.split(__file__)[0])+'/../'
          #filepath = basedir+'config/emit.json'

      with open(filepath,'r') as fin:
          config_dict = json.load(fin)
          for key, val in config_dict.items():
              setattr(self,key,val) 

      
      # Adjust local filepaths where needed
      my_directory, my_executable = os.path.split(os.path.abspath(__file__))
      for fi in dir(self):
          if '_file' in fi:
              path = getattr(self,fi)
              if path[0] != '/':
                  path = os.path.join(my_directory, path)
                  setattr(self,fi,path)


      # Define masked rows 
      if self.first_illuminated_row>=0 and self.last_illuminated_row>=0:
          self.masked_rows = np.concatenate((
                np.arange(0, self.first_illuminated_row, dtype=int),
                np.arange(self.last_illuminated_row+1, self.native_rows, dtype=int)),
               axis=0)
      else:
          self.masked_rows = []

      # Define masked columns
      if self.first_masked_col_right >=0 and self.last_masked_col_left>=0:
          self.masked_cols = np.concatenate((
                np.arange(0, self.last_masked_col_left+1, dtype=int),
                np.arange(self.first_masked_col_right, self.native_columns, dtype=int)),
               axis=0)
      else:
          self.masked_cols = []

      # These columns used for stray light checks
      self.vignetted_cols = np.concatenate((np.arange(self.last_masked_col_left+1, 
              self.first_illuminated_column, dtype=int),
             np.arange(self.last_illuminated_column+1, 
               self.first_masked_col_right, dtype=int)),axis=0)


# EMIT frames can be in native format or in subframe (328 row) format.
# This function extracts a subframe from a native format frame
def frame_extract(frame, fpa, clip_columns = False):
  if frame.shape[1] != fpa.native_columns:
     raise IndexError('All frames should have '+str(fpa.native_columns)+' columns')
  if frame.shape[0] != fpa.native_rows:
     raise IndexError('Native frames should have '+str(fpa.native_rows)+' rows')
  frame = frame[fpa.first_valid_row:(fpa.last_valid_row+1),:]
  if clip_columns:
      frame = frame[:,fpa.first_illuminated_column:(fpa.last_illuminated_column+1)]
  return frame
 

# EMIT frames can be in native format or in subframe (328 row) format.
# This function makes sure that all frames have native format by 
# embedding subframes inside some padding.
def frame_embed(frame, fpa):
    if frame.shape[1] != fpa.valid_columns:
       raise IndexError('All frames should have '+str(fpa.valid_columns)+' columns')
    if frame.shape[0] == fpa.valid_rows:
       raise IndexError('Invalid number of rows: %i'%frame.shape[0])
    embedded = np.zeros((fpa.native_rows, fpa.native_columns))
    embedded[fpa.first_valid_row, fpa.last_valid_row+1] = frame
    return embedded    
