############################################################################
#Project 2ACE
#PI: Prof. Lili Qiu @ UT Austin/MSRA
#Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
#Copyright @ The University of Texas at Austin, 2023
#Partially inherited from Teng Wei/Song Wang/Jingqi Huang/Xinyu Zhang @ UCSD
############################################################################

import subprocess as sp
#from scipy import io
import json
import os
from random import randint
import copy
import numpy as np
import scipy.io as sio

def set_beam_num(file_name, beam_num):
	try:
		sp.check_call(['./wil6210_brd_mod', '-set_beam_num=%d' % beam_num, \
			'-in_brd=%s' % file_name, '-out_brd=%s' % file_name])
	except Exception as e:
		print(e)

def set_beam(file_name, module_index, sector_type, sector_index, gain, phase, amplify):
	try:
		sp.check_call(['./wil6210_brd_mod', '-set_pattern', \
			'-in_brd=%s' % file_name, '-out_brd=%s' % file_name, '-module_index=%d' % module_index, '-sector_type=%d' % sector_type, '-sector_index=%d' % sector_index, \
			'-mag=%s' % gain, '-phase=%s' % phase, '-amp=%s' % amplify])
	except Exception as e:
		print(e)

def enable_modules(file_name, module_index, isTx):
	try:
		for m_idx in module_index:
			sp.check_call(['./wil6210_brd_mod','-in_brd=%s' % file_name, '-out_brd=%s' % file_name, '-enable_module', '-module_index=%d' % m_idx, '-sector_type=%d' % isTx])
	except Exception as e:
		print(e)

def disable_modules(file_name, module_index, isTx):
		try:
			for m_idx in module_index:
				sp.check_call(['./wil6210_brd_mod','-in_brd=%s' % file_name, '-out_brd=%s' % file_name, '-disable_module', '-module_index=%d' % m_idx, '-sector_type=%d' % isTx])
		except Exception as e:
			print(e)

#90 degree
grouping = [[[1,2,3,4],[5,6,7,8],[17,18,19,20],[21,22,23,24]],[[1,2],[3,4],[5,7],[6,8],[17,18],[19,20],[21,23],[22,24]]]
calibration_bit =  [0,2,3,0,0,3,0,3,0,0,0,0,0,0,0,0,1,0,0,3,0,3,0,0,0,0,0,0,0,0,0,0]
calibration_bit2 = [0,2,3,0,0,3,0,3,0,0,0,0,0,0,0,0,1,0,0,3,0,3,0,0,0,0,0,0,0,0,0,0]

mag = ['77777777000000007777777700000000']

amp = '66666666'
brd_name_in = 'wil6210_sparrow_plus_rx.brd'

cb_num = 160

separation = [32, 96, 160]
active_ant = list(np.array([1,2,3,4,5,6,7,8,17,18,19,20,21,22,23,24])-1)
phase = ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
sector_per_cb = 1
cb_all = []
cb_all_actual = []

for i in range(cb_num):

    beam_num = 1

    brd_name_out = './codebook_brd/multires_16ant_rx/wil6210_rx_cb' + str(i+1) + '.brd' 
    os.system('cp ' + brd_name_in + ' ' + brd_name_out)

    cb_cur = []
    cb_cur_actual = []
    for j in range(sector_per_cb):
        curphase = copy.deepcopy(phase)
        inferphase = copy.deepcopy(phase)
        if i<separation[0]:
           for k in range(4):
               cur_bit = randint(0,3)#'0'
               for _ in range(4):
                    randint(0,3)
               for l in range(4):
                   actuall_bit = cur_bit - calibration_bit[grouping[0][k][l]-1]
                   if actuall_bit < 0:
                        actuall_bit += 4
                   curphase[grouping[0][k][l]-1] = str(actuall_bit)
                   inferphase[grouping[0][k][l]-1] = str(cur_bit) 
           cb_cur.append(''.join(inferphase))
           cb_cur_actual.append(''.join(curphase))
        elif i<separation[1] and i >= separation[0]:
           for k in range(8):
               cur_bit = randint(0,3)
               for _ in range(4):
                    randint(0,3)
               for l in range(2):
                   actuall_bit = cur_bit - calibration_bit2[grouping[1][k][l]-1]
                   if actuall_bit < 0:
                        actuall_bit += 4
                   curphase[grouping[1][k][l]-1] = str(actuall_bit)
                   inferphase[grouping[1][k][l]-1] = str(cur_bit) 
           cb_cur.append(''.join(inferphase))
           cb_cur_actual.append(''.join(curphase))
        elif i<separation[2] and i >= separation[1]:
           for k in range(16):
               cur_bit = randint(0,3)
               for _ in range(4):
                    randint(0,3)
               actuall_bit = cur_bit
               if actuall_bit < 0:
                    actuall_bit += 4
               curphase[active_ant[k]] = str(actuall_bit)
               inferphase[active_ant[k]] = str(cur_bit) 
           cb_cur.append(''.join(inferphase))
           cb_cur_actual.append(''.join(curphase))

    #cur_mag = list(mag[0])
    #cur_mag[i] = '7'
    #cur_mag = ''.join(cur_mag)
    #cur_mag = [cur_mag]
    #print(cur_mag)

    #set_beam_num(brd_name, beam_num)

    # clear previous codebook entrys
    for sec_idx in range(0, sector_per_cb):
        set_beam(brd_name_out, 0, 0, sec_idx, '00000000000000000000000000000000', '00000000000000000000000000000000', '00000000')

    # set entry values
    for sec_idx in range(0, sector_per_cb):
        cur_sec = sec_idx
        print('Module:%d, Sector:%d' % (0, cur_sec))
        set_beam(brd_name_out, 0, 0, cur_sec, mag[sec_idx], cb_cur_actual[sec_idx], amp)

    cb_all.extend(cb_cur)
    cb_all_actual.extend(cb_cur_actual)

with open('./codebook_brd/multires_16ant_rx/multires_16_rx.txt', 'w') as f:
     for sector in cb_all:
         f.write("%s\n" % sector)

with open('./codebook_brd/multires_16ant_rx/multires_16_rx_actual.txt', 'w') as f:
     for sector in cb_all_actual:
         f.write("%s\n" % sector)