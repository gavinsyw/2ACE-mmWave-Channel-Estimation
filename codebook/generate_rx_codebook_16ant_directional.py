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

pha = ['02323321301012102323120300020200',
'02323321301012102323110300020200',
'02322321301011102323110300020200',
'02322321301011101312110333320230',
'02322321301011131212110333320230',
'02322321301011131212010333313230',
'02322321301011130201010223313230',
'02322321301011130130003223203123',
'02022211311111133130303112202123',
'02022211211111133023303102132113',
'02031211211110132012332002121002',
'02031211211110122312232031010002',
'02031211211110121201221321300331',
'02131101221210120230121210333331',
'02101101121210120123120200232220',
'02100101121213113012010133122210',
'02100101121213112001010023111200',
'02200031131213112330003312000103',
'02200031031313111223303302330133',
'02213031031312110212332231223022',
'02213031031312100101232131212022',
'02213031031312103001221121102311',
'02313321001012102030121010131311',
'02313321301012102323121000030301',
'02322321301011101312110300020200',
'02322321301011131212010333320230',
'02322321301011131201010223313230',
'02322321301011130201010223213220',
'02322321301011130101003223203123',
'02322321301011130130003223202123',
'02322311301011130130003212202123',
'02322311311011130130003212202123']

amp = '66666666'

for sec_idx in range(len(pha)):
    os.system('cp wil6210_sparrow_plus_tx.brd ./codebook_brd/directional_16ant/wil6210_rx' + str(sec_idx+1) + '.brd')
    brd_name_out = './codebook_brd/directional_16ant/wil6210_rx'+ str(sec_idx+1) +'.brd'
    set_beam(brd_name_out, 0, 0, 0, '00000000000000000000000000000000', '00000000000000000000000000000000', '00000000')

for sec_idx in range(len(pha)):
    brd_name_out = './codebook_brd/directional_16ant/wil6210_rx'+ str(sec_idx+1) +'.brd'
    print('Module:%d, Sector:%d' % (0, sec_idx))
    set_beam(brd_name_out, 0, 0, 0, mag, pha[sec_idx], amp)

with open('./codebook_brd/directional_16ant/codeboook_16ant_directional.txt', 'w') as f:
     for sector in pha:
         f.write("%s\n" % sector)