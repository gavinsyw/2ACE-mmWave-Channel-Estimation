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
'03211321321300322002210331312200',
'02313312010310322031301001220113',
'02131132320012312013013032113210',
'03211321321301322103320002022310',
'00222220220202021133333133131313',
'00222220220202022200000200202020',
'03210321332011021321203231201300',
'02323212120031132021331012230230',
'02132132030130122020113003321031',
'03201021022312033221021133100132',
'00222220220202022200000200202020',
'00331110130301012233111231232121',
'00330110023333232122121210110300',
'03030101323221121010010133303222',
'03003001223321123332300012221112',
'00103030333032222221233302111112',
'00113330310103032211333213212323',
'00113330310103032211333213212323',
'00110330021111212322323230330100',
'01010303121223323030030311101222',
'01001003221123321112100032223332',
'00301010111012222223211102333332',
'00331110130301012233111231232121',
'00222220220202022200000200202020',
'01230123112033023123201213203100',
'02121232320013312023113032210210',
'02312312010310322020331001123013',
'01203023022132011223023311300312',
'00222220220202022200000200202020',
'00222220220202021133333133131313',
'01233123123100122002230113132200',
'02131132030130122013103003220331',
'02313312120032132031031012331230',
'01233123123103122301120002022130',
'00222220220202023311111311313131']

amp = '66666666'

for sec_idx in range(len(pha)):
    os.system('cp wil6210_sparrow_plus_tx.brd ./codebook_brd/sweeping_thetaNphi_16ant/wil6210_rx' + str(sec_idx+1) + '.brd')
    brd_name_out = './codebook_brd/sweeping_thetaNphi_16ant/wil6210_rx' + str(sec_idx+1) +'.brd'
    set_beam(brd_name_out, 0, 0, 0, '00000000000000000000000000000000', '00000000000000000000000000000000', '00000000')

for sec_idx in range(len(pha)):
    brd_name_out = './codebook_brd/sweeping_thetaNphi_16ant/wil6210_rx' + str(sec_idx+1) +'.brd'
    print('Module:%d, Sector:%d' % (0, sec_idx))
    set_beam(brd_name_out, 0, 0, 0, mag, pha[sec_idx], amp)

with open('./codebook_brd/sweeping_thetaNphi_16ant/codeboook_16ant_sweeping_thetaNphi.txt', 'w') as f:
     for sector in pha:
         f.write("%s\n" % sector)