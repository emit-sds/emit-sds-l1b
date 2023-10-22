#! /usr/bin/env python
#
#  Copyright 2020 California Institute of Technology
#
# EMIT Radiometric Calibration code
# Author: Philip G Brodrick, philip.brodrick@jpl.nasa.gov

import os, sys
import numpy as np
import logging

OBC_DARK1     = 2
OBC_SCIENCE   = 3
OBC_DARK2     = 4
OBC_BRIGHTMED = 5
OBC_BRIGHTHI  = 6
OBC_LASER     = 7
NOTFOUND = -1


# def bytes_from_file(filename, chunksize=8192):
#     with open(filename, "rb") as f:
#         while True:
#             chunk = f.read(chunksize)
#             if chunk:
#                 for b in chunk:
#                     yield b
#             else:
#                 break

# # example:
# for b in bytes_from_file('filename'):
#     do_stuff_with(b)



# with open("myfile", "rb") as f:
#     byte = f.read(1)
#     while byte != b"":
#         # Do stuff with byte.
#         byte = f.read(1)


def read_frames(raw_file, num_frames_to_read, num_channels, num_samples, start_line=0):
    NC = num_channels
    NS = num_samples
    NFRAME = NS * NC
    NRAW   = NFRAME + NS # (1 extra metadata channel)
    RAWBYTES = NRAW * 2
    FRAMEBYTES = NFRAME * 2

    num_read = 0

    obcv_prev = NOTFOUND # to keep track of state transitions
    
    img_size = os.path.getsize(raw_file)
    if img_size % RAWBYTES != 0:
        logging.warning(f'file "{raw_file}" contains truncated frames')
    
    nl_max = int(float(img_size)/RAWBYTES)
    if start_line < 0: # indexing from EOF
        start_line = nl_max+start_line

    frames = np.zeros((num_frames_to_read, NC, NS),dtype=np.float32)
    frame_obcv = np.zeros(num_frames_to_read,dtype=np.uint16)
    frame_meta = np.zeros([num_frames_to_read,2],dtype=np.int64)

#     frames = np.zeros((num_frames_to_read, NS, NC),dtype=np.int16)
#     frame_obcv = np.ones(num_frames_to_read,dtype=np.int32)*NOTFOUND
#     frame_meta = np.zeros([num_frames_to_read,2],dtype=np.int64)
    

    if start_line > nl_max:
        logging.debug('start_line > max frames in file, nothing to read')
        return None, None, num_read, None

    with open(raw_file, 'rb') as f:
        if start_line > 0:
            f.seek(start_line*RAWBYTES,os.SEEK_CUR)
            nl_max = nl_max-start_line
            
        while num_read < nl_max and num_read < num_frames_to_read:
            # read frame header
            buf = np.fromfile(f, count=NS, dtype='<u2') #Little-endian uint-16
            if len(buf) > 0:
                obcv  = buf[321] >> 8                
                logging.debug(f'num_read: {num_read}, start_line+num_read: {start_line+num_read}, obcv: {obcv}, obcv_prev: {obcv_prev}')

                if obcv == OBC_DARK2 and obcv_prev == OBC_SCIENCE:
                    tn_idx = start_line+num_read-1
                    logging.debug(f'OBC_SCIENCE->OBC_DARK2 transition detected at index {tn_idx}')
                    break
                obcv_prev = obcv

                #gps_state = buf[340]
                #print(GPS_STATUS_MSG[gps_state])
                #if gps_state==HxBAD:
                #    warn('Bad GPS msg found')
                clock = format_clock_words(buf[0],buf[1])
                count,countN = extract_fc(buf[160]),extract_fc(buf[319])
                logging.debug(f'clock: {clock}, count: {count}, countN: {countN}')
                    
                frame_meta[num_read,:]  = [clock, count]
                frame_obcv[num_read] = obcv

                frames[num_read,...] = np.fromfile(f, count=NC*NS, dtype='<u2').reshape((NC,NS)).astype(float)

                num_read = num_read+1

            else:
                break

        frame_meta = frame_meta[:num_read,:]
        frame_obcv = frame_obcv[:num_read]
        
    if np.all(frame_obcv == 0):
        logging.debug(f'no OBC_SCIENCE frames read in range [{start_line}, {start_line+num_read}]')

    # Return in BIL format        
    return frames, frame_meta, num_read, frame_obcv

# def read_frames(raw_file, num_frames_to_read, num_channels, num_samples, start_line=0):
#     NC = num_channels
#     NS = num_samples
#     NFRAME = NS * NC
#     NRAW   = NFRAME + NS # (1 extra metadata channel)
#     RAWBYTES = NRAW * 2
#     FRAMEBYTES = NFRAME * 2

#     num_read = 0

#     obcv_prev = NOTFOUND # to keep track of state transitions
    
#     img_size = os.path.getsize(raw_file)
#     if img_size % RAWBYTES != 0:
#         logging.warning(f'file "{raw_file}" contains truncated frames')
    
#     nl_max = int(float(img_size)/RAWBYTES)
#     if start_line < 0: # indexing from EOF
#         start_line = nl_max+start_line

#     frames = np.zeros((num_frames_to_read, NS, NC),dtype=np.int16)
#     frame_obcv = np.ones(num_frames_to_read,dtype=np.int32)*NOTFOUND
#     frame_meta = np.zeros([num_frames_to_read,2],dtype=np.int64)

#     if start_line > nl_max:
#         logging.debug('start_line > max frames in file, nothing to read')
#         return None, None, num_read, None

#     with open(raw_file, 'rb') as f:
#         if start_line > 0:
#             f.seek(start_line*RAWBYTES,os.SEEK_CUR)
#             nl_max = nl_max-start_line
            
#         while num_read < nl_max and num_read < num_frames_to_read:
#             # read frame header
#             buf = np.fromfile(f, count=NS, dtype='<u2')
#             if len(buf) > 0:
#                 obcv  = buf[321] >> 8                
#                 logging.debug(f'num_read: {num_read}, start_line+num_read: {start_line+num_read}, obcv: {obcv}, obcv_prev: {obcv_prev}')

#                 if obcv == OBC_DARK2 and obcv_prev == OBC_SCIENCE:
#                     tn_idx = start_line+num_read-1
#                     logging.debug(f'OBC_SCIENCE->OBC_DARK2 transition detected at index {tn_idx}')
#                     break
#                 obcv_prev = obcv

#                 #gps_state = buf[340]
#                 #print(GPS_STATUS_MSG[gps_state])
#                 #if gps_state==HxBAD:
#                 #    warn('Bad GPS msg found')
#                 clock = format_clock_words(buf[0],buf[1])
#                 count,countN = extract_fc(buf[160]),extract_fc(buf[319])
#                 logging.debug(f'clock: {clock}, count: {count}, countN: {countN}')
                    
#                 frame_meta[num_read,:]  = [clock, count]
#                 frame_obcv[num_read] = obcv

#                 frames[num_read,...] = np.fromfile(f, count=NC*NS, dtype='<u2').reshape((NC,NS)).T

#                 num_read = num_read+1

#             else:
#                 break

#         frame_meta = frame_meta[:num_read,:]
#         frame_obcv = frame_obcv[:num_read]
        
#     if np.all(frame_obcv == NOTFOUND):
#         logging.debug(f'no OBC_SCIENCE frames read in range [{start_line}, {start_line+num_read}]')

#     # Return in BIL format
#     frames = np.transpose(frames,(0,2,1))
        
#     return frames, frame_meta, num_read, frame_obcv


def read_frames_metadata(raw_file, num_frames_to_read, num_channels, num_samples, start_line=0):
    NC = num_channels
    NS = num_samples
    NFRAME = NS * NC
    NRAW   = NFRAME + NS # (1 extra metadata channel)
    RAWBYTES = NRAW * 2
    FRAMEBYTES = NFRAME * 2

    num_read = 0

    obcv_prev = NOTFOUND # to keep track of state transitions
    
    img_size = os.path.getsize(raw_file)
    if img_size % RAWBYTES != 0:
        logging.warning(f'file "{raw_file}" contains truncated frames')
    
    nl_max = int(float(img_size)/RAWBYTES)
    if start_line < 0: # indexing from EOF
        start_line = nl_max+start_line

    if start_line > nl_max:
        logging.debug('start_line > max frames in file, nothing to read')
        return None, num_read, None

    frame_obcv = np.zeros(num_frames_to_read,dtype=np.int32)
    frame_meta = np.zeros([num_frames_to_read,2],dtype=np.int64)
    with open(raw_file, 'rb') as f:
        if start_line > 0:
            f.seek(start_line*RAWBYTES,os.SEEK_CUR)
            nl_max = nl_max-start_line
            
        while num_read < nl_max and num_read < num_frames_to_read:
            # read frame header
            buf = np.fromfile(f, count=NS, dtype='<u2')
            if len(buf) > 0:
                obcv  = buf[321] >> 8                
                logging.debug(f'num_read: {num_read}, start_line+num_read: {start_line+num_read}, obcv: {obcv}, obcv_prev: {obcv_prev}')

                if obcv == OBC_DARK2 and obcv_prev == OBC_SCIENCE:
                    tn_idx = start_line+num_read-1
                    logging.debug(f'OBC_SCIENCE->OBC_DARK2 transition detected at index {tn_idx}')
                    break
                obcv_prev = obcv

                #gps_state = buf[340]
                #print(GPS_STATUS_MSG[gps_state])
                #if gps_state==HxBAD:
                #    warn('Bad GPS msg found')
                clock = format_clock_words(buf[0],buf[1])
                count,countN = extract_fc(buf[160]),extract_fc(buf[319])
                logging.debug(f'clock: {clock}, count: {count}, countN: {countN}')
                    
                frame_meta[num_read,:]  = [clock, count]
                frame_obcv[num_read] = obcv
                num_read = num_read+1
                f.seek(FRAMEBYTES,os.SEEK_CUR)                
            else:
                break

        frame_meta = frame_meta[:num_read,:]
        frame_obcv = frame_obcv[:num_read]
        
    if np.all(frame_obcv == NOTFOUND):
        logging.debug(f'no OBC_SCIENCE frames read in range [{start_line}, {start_line+num_read}]')
        
    return frame_meta, num_read, frame_obcv


def format_clock_words(msw,lsw):
    return (np.int64(msw)<<16)+lsw

def extract_fc(word):
    # applies 14 bit mask to extract frame count from 16-bit pps/gps word
    return np.bitwise_and(np.uint32(2**14*-1), np.int32(word))

