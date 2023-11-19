############################################################################
#Project 2ACE
#PI: Prof. Lili Qiu @ UT Austin/MSRA
#Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
#Copyright @ The University of Texas at Austin, 2023
#Partially inherited from Teng Wei/Song Wang/Jingqi Huang/Xinyu Zhang @ UCSD
############################################################################

from random import randint
import subprocess as sp
import json
import os

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

cb_all = []
mag = '77777777000000007777777700000000'

pha = ['00222220220202023311111311313131',
'00222220220202023311111311313131',
'00222220220202023311111311313131',
'00222220220202023311000310313030',
'00222220220202022211000210212020',
'00222220220202022200000200202020',
'00222220220202022200000200202020',
'00222220220202021133000130131010',
'00222220220202021122333123121313',
'00331110130301010022333023020303',
'00331110130301010011222012010202',
'00331110130301013300222302303232',
'00331110130301012233222232232222',
'00331110130301012222111221222121',
'00000000000000001111111111111111',
'00000000000000000000000000000000',
'00000000000000000000000000000000',
'00000000000000003333333333333333',
'00113330310103032222333223222323',
'00113330310103032211222212212222',
'00113330310103031100222102101212',
'00113330310103030033222032030202',
'00113330310103030022111021020101',
'00222220220202023322111321323131',
'00222220220202023311000310313030',
'00222220220202022200000200202020',
'00222220220202022200000200202020',
'00222220220202022233000230232020',
'00222220220202021133000130131010',
'00222220220202021133333133131313',
'00222220220202021133333133131313',
'00222220220202021133333133131313']

amp = '66666666'

for sec_idx in range(len(pha)):
    os.system('cp wil6210_sparrow_plus_tx.brd ./codebook_brd/sweeping_phi_16ant/wil6210_rx' + str(sec_idx+1) + '.brd')
    brd_name_out = './codebook_brd/sweeping_phi_16ant/wil6210_rx'+ str(sec_idx+1) +'.brd'
    set_beam(brd_name_out, 0, 0, 0, '00000000000000000000000000000000', '00000000000000000000000000000000', '00000000')

for sec_idx in range(len(pha)):
    brd_name_out = './codebook_brd/sweeping_phi_16ant/wil6210_rx'+ str(sec_idx+1) +'.brd'
    print('Module:%d, Sector:%d' % (0, sec_idx))
    set_beam(brd_name_out, 0, 0, 0, mag, pha[sec_idx], amp)

with open('./codebook_brd/sweeping_phi_16ant/codeboook_16ant_sweeping_phi.txt', 'w') as f:
     for sector in pha:
        f.write("%s\n" % sector)