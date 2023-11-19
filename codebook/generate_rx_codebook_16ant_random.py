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


mag = ['77777777000000007777777700000000']

amp = '66666666'
brd_name_in = 'wil6210_sparrow_plus_rx.brd'

cb_num = 64
sector_per_cb = 1
cb_all = []

for i in range(cb_num):

    beam_num = 1

    brd_name_out = './codebook_brd/random_16ant_rx/wil6210_rx_16ant_cb' + str(i+1) + '.brd' 
    os.system('cp ' + brd_name_in + ' ' + brd_name_out)

    cb_cur = []
    for j in range(sector_per_cb):
        sector_cur = ''
        for k in range(32):
            if k == 0:
                sector_cur += '0'
                continue
            sector_cur += str(randint(0,3))
        cb_cur.append(sector_cur)

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
        set_beam(brd_name_out, 0, 0, cur_sec, mag[sec_idx], cb_cur[sec_idx], amp)

    cb_all.extend(cb_cur)

with open('./codebook_brd/random_16ant_rx/random_codebook_16ant_rx.txt', 'w') as f:
     for sector in cb_all:
         f.write("%s\n" % sector)