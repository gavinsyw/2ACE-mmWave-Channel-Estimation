####################################################################
#Project 2ACE
#PI: Prof. Lili Qiu @ UT Austin/MSRA
#Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
#Copyright @ The University of Texas at Austin, 2023
####################################################################

import numpy as np
from itertools import combinations
import matlab.engine
import time
import subprocess
import socket
import json
import ast
from multiprocessing import Pool
import subprocess as sp
import os
import struct
import scipy.io as sio
import sys
import copy
from datetime import date
from codebook_library import codebook_generator, svd_beamformer, set_beam_num, set_beam, collect_ACO_tx, collect_ACO_rx, twos_comp, fetch_rss, rss2csi, get_ACO_codebook_bit,generate_tx_codebook

if __name__ == '__main__':

	ip = '0.0.0.0'
	ap_ipaddr = '192.168.0.10'
	port = 10002

	num_ant = str(16)
	active_ant = list(np.array([1,2,3,4,5,6,7,8,17,18,19,20,21,22,23,24])-1)
	
	env = "indoor"

	BF_on = True

	total_num_probe = int(num_ant)*int(num_ant)*4

	num_round = 64

	repeat = 40

	methods = ['Sector Level Sweeping','A2','A2 w/ Multiresolution','A2 nuclear','PLOMP','PLGAMP','PhaseLift', 'ACO']

	num_method = len(methods)

	num_round_multires = 160

	num_round_directional = 32

	num_round_sweeping_thetaNphi = 36

	num_round_sweeping_phi = 32

	rss = np.zeros((num_round, 62))

	rss_multires = np.zeros((num_round_multires, 62))

	rss_directional = np.zeros((num_round_directional, num_round_directional))

	rss_sweeping_thetaNphi = np.zeros((num_round_sweeping_thetaNphi, num_round_sweeping_thetaNphi))

	rss_sweeping_phi = np.zeros((num_round_sweeping_phi, num_round_sweeping_phi))

	M = [4, 36, 121, 225, 361, 529, 784, 1024]

	compare_bf_rss = np.zeros((repeat, len(M), num_method, 64))

	H_recover = np.zeros((repeat, len(M), num_method-4, int(num_ant)*int(num_ant)), dtype = np.complex128)

	file_out = '../result/training_rss_' + num_ant + 'x' + num_ant + '_' + str(date.today()) + '_' + env + '_test.mat'

	codebook_path = '../codebook/codebook_mat/random_probe_cb_' + num_ant + 'x' + num_ant + '.mat'

	directional_codebook_path = '../codebook/codebook_mat/directional_codebook_' + num_ant + 'x' + num_ant + '.mat'

	multires_codebook_path = '../codebook/codebook_mat/random_probe_cb_' + num_ant + 'x' + num_ant + '_multires.mat'

	bfresult_path = '../result/beamforming_' + num_ant + 'x' + num_ant + '_' + str(date.today()) + '_' + env + '_test.mat'

	print(file_out + '\n')
	print(codebook_path + '\n')
	print(directional_codebook_path + '\n')

	print('Matlab Engine Start')
	eng = matlab.engine.start_matlab()

	print('Please leave the environment! Start in 10s')
	time.sleep(10)

################################################################################################################################################################################################
	'''
	start thetaNphi sweeping probing
	'''
	print('Start sweeping thetaNphi probing')
	os.system('echo "123456" | sshpass -p "123456" scp ../codebook/codebook_brd/sweeping_thetaNphi_' + num_ant + 'ant/wil6210_tx.brd utwncg@' + ap_ipaddr + ':wil6210.brd')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo cp ~/wil6210.brd /lib/firmware/wil6210.brd"')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 down"')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 up"')
	for i in range(num_round_sweeping_thetaNphi):
		os.system('echo "123456" | sudo -S -k cp ../codebook/codebook_brd/sweeping_thetaNphi_' + num_ant + 'ant/wil6210_rx' + str(i+1) + '.brd /lib/firmware/wil6210.brd')
		os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
		os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
		time.sleep(1)

		dump_count = 0
		cur_rss = fetch_rss(ip, port, ap_ipaddr, dump_count)
		cur_rss = np.median(cur_rss, axis = 0)

		cur_rss[cur_rss>1000] = 0
		cur_rss = cur_rss*0.0652-74.3875
		print(cur_rss[2:38])
		rss_sweeping_thetaNphi[i,:] = cur_rss[2:38]

		print('\nFinished round ' + str(i+1) + '/'  + str(num_round_sweeping_thetaNphi) + '\n')

		#check phased array and baseband chipset temperature
		try:
			content = open('/sys/kernel/debug/ieee80211/phy0/wil6210/temp').readlines()
			mac_temp = float(content[0].split(' ')[-1].strip('\n'))
			radio_temp = float(content[1].split(' ')[-1].strip('\n'))
			if mac_temp > 70 or radio_temp > 62.5:
				os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
				time.sleep(20)
				os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
			else:
				print('Temperature test passes\n')
				print('Current Temperature is MAC: ' + str(mac_temp) + ' RADIO: ' + str(radio_temp) + '\n')
		except Exception as e:
			print(e)

	sio.savemat(file_out,{'rss':rss, 'rss_directional':rss_directional, 'rss_sweeping_phi':rss_sweeping_phi, 'rss_sweeping_thetaNphi':rss_sweeping_thetaNphi})

################################################################################################################################################################################################
	'''
	start phi sweeping probing
	'''
	print('Start sweeping phi probing')
	os.system('echo "123456" | sshpass -p "123456" scp ../codebook/codebook_brd/sweeping_phi_' + num_ant + 'ant/wil6210_tx.brd utwncg@' + ap_ipaddr + ':wil6210.brd')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo cp ~/wil6210.brd /lib/firmware/wil6210.brd"')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 down"')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 up"')
	for i in range(num_round_sweeping_phi):
		os.system('echo "123456" | sudo -S -k cp ../codebook/codebook_brd/sweeping_phi_' + num_ant + 'ant/wil6210_rx' + str(i+1) + '.brd /lib/firmware/wil6210.brd')
		os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
		os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
		time.sleep(1)

		dump_count = 0
		cur_rss = fetch_rss(ip, port, ap_ipaddr, dump_count)
		cur_rss = np.median(cur_rss, axis = 0)

		cur_rss[cur_rss>1000] = 0
		cur_rss = cur_rss*0.0652-74.3875
		print(cur_rss[2:34])
		rss_sweeping_phi[i,:] = cur_rss[2:34]

		print('\nFinished round ' + str(i+1) + '/'  + str(num_round_sweeping_phi) + '\n')

		#check phased array and baseband chipset temperature
		try:
			content = open('/sys/kernel/debug/ieee80211/phy0/wil6210/temp').readlines()
			mac_temp = float(content[0].split(' ')[-1].strip('\n'))
			radio_temp = float(content[1].split(' ')[-1].strip('\n'))
			if mac_temp > 70 or radio_temp > 62.5:
				os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
				time.sleep(20)
				os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
			else:
				print('Temperature test passes\n')
				print('Current Temperature is MAC: ' + str(mac_temp) + ' RADIO: ' + str(radio_temp) + '\n')
		except Exception as e:
			print(e)

	sio.savemat(file_out,{'rss':rss, 'rss_directional':rss_directional, 'rss_sweeping_phi':rss_sweeping_phi, 'rss_sweeping_thetaNphi':rss_sweeping_thetaNphi})

################################################################################################################################################################################################
	'''
	start directional probing
	'''
	print('Start directional probing for PLOMP and PLGAMP')
	os.system('echo "123456" | sshpass -p "123456" scp ../codebook/codebook_brd/directional_' + num_ant + 'ant/wil6210_tx.brd utwncg@' + ap_ipaddr + ':wil6210.brd')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo cp ~/wil6210.brd /lib/firmware/wil6210.brd"')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 down"')
	os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 up"')
	for i in range(num_round_directional):
		os.system('echo "123456" | sudo -S -k cp ../codebook/codebook_brd/directional_' + num_ant + 'ant/wil6210_rx' + str(i+1) + '.brd /lib/firmware/wil6210.brd')
		os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
		os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
		time.sleep(1)

		dump_count = 0
		cur_rss = fetch_rss(ip, port, ap_ipaddr, dump_count)
		cur_rss = np.median(cur_rss, axis = 0)

		cur_rss[cur_rss>1000] = 0
		cur_rss = cur_rss*0.0652-74.3875
		print(cur_rss[2:34])
		rss_directional[i,:] = cur_rss[2:34]

		print('\nFinished round ' + str(i+1) + '/'  + str(num_round_directional) + '\n')

		#check phased array and baseband chipset temperature
		try:
			content = open('/sys/kernel/debug/ieee80211/phy0/wil6210/temp').readlines()
			mac_temp = float(content[0].split(' ')[-1].strip('\n'))
			radio_temp = float(content[1].split(' ')[-1].strip('\n'))
			if mac_temp > 70 or radio_temp > 62.5:
				os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
				time.sleep(20)
				os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
			else:
				print('Temperature test passes\n')
				print('Current Temperature is MAC: ' + str(mac_temp) + ' RADIO: ' + str(radio_temp) + '\n')
		except Exception as e:
			print(e)

	sio.savemat(file_out,{'rss':rss, 'rss_directional':rss_directional, 'rss_sweeping_phi':rss_sweeping_phi, 'rss_sweeping_thetaNphi':rss_sweeping_thetaNphi})

	cb_directional = sio.loadmat(directional_codebook_path)['cb']
	cb_train_directional_cur_amp = matlab.double(np.abs(cb_directional).tolist())
	cb_train_directional_cur_angle = matlab.double(np.angle(cb_directional).tolist())
	rss_train_directional = matlab.double(rss_directional.tolist())

	'''
	Start Directional Estimation
	'''
	print('Start Directional Estimation')
	H_amp_cur, H_angle_cur = eng.channel_recovery_ADMM_v2_simulation_directional(int(num_ant), \
		int(num_ant), cb_train_directional_cur_amp, cb_train_directional_cur_angle, rss_train_directional, eng.double(3), nargout=2)
	H_amp_cur = np.array(H_amp_cur)
	H_angle_cur = np.array(H_angle_cur)  
	H_recover_directional = np.squeeze(H_amp_cur*np.exp(1j*H_angle_cur))

################################################################################################################################################################################################
	'''
		Start Random Probing
	'''
	print('Start first round random probing for phaselift')
	try:
		for i in range(num_round):
			os.system('echo "123456" | sudo -S -k cp ../codebook/codebook_brd/random_' + num_ant + 'ant_rx/wil6210_rx_' + num_ant + 'ant_cb' + str(i+1) + '.brd /lib/firmware/wil6210.brd')
			os.system('echo "123456" | sshpass -p "123456" scp ../codebook/codebook_brd/random_' + num_ant + 'ant_tx/wil6210_tx_cb' + str(i+1) + '.brd utwncg@' + ap_ipaddr + ':wil6210.brd')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo cp ~/wil6210.brd /lib/firmware/wil6210.brd"')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 down"')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 up"')
			os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
			os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
			time.sleep(1)
			
			dump_count = 0
			cur_rss = fetch_rss(ip, port, ap_ipaddr, dump_count)
			cur_rss = np.median(cur_rss, axis = 0)

			cur_rss[cur_rss>1000] = 0
			cur_rss = cur_rss*0.0652-74.3875
			print(cur_rss)
			rss[i,:] = cur_rss[2:]

			print('\nFinished round ' + str(i+1) + '/'  + str(num_round) + '\n')
			sio.savemat(file_out,{'rss':rss})

			#check phased array and baseband chipset temperature
			try:
				content = open('/sys/kernel/debug/ieee80211/phy0/wil6210/temp').readlines()
				mac_temp = float(content[0].split(' ')[-1].strip('\n'))
				radio_temp = float(content[1].split(' ')[-1].strip('\n'))
				if mac_temp > 70 or radio_temp > 62.5:
					os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
					time.sleep(20)
					os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
				else:
					print('Temperature test passes\n')
					print('Current Temperature is MAC: ' + str(mac_temp) + ' RADIO: ' + str(radio_temp) + '\n')
			except Exception as e:
				print(e)

		sio.savemat(file_out,{'rss':rss, 'rss_multires':rss_multires, 'rss_directional':rss_directional, 'rss_sweeping_phi':rss_sweeping_phi, 'rss_sweeping_thetaNphi':rss_sweeping_thetaNphi})

	except (Exception, np.AxisError) as e:
		print(e)

	#start estimation with random probing
	print('start sleeping!')
	#time.sleep(10)
	rss = np.reshape(np.transpose(rss),[-1,1])
	cb = sio.loadmat(codebook_path)['cb']

	for _ in range(5):
		rand_ind = list(np.random.permutation(num_round*62))
		rss = rss[rand_ind,:]
		cb = cb[rand_ind,:]

	rss_train = rss[:,:]
	cb_train = cb[:,:]

	rss_train_ori = copy.deepcopy(rss_train)
	rss_train = matlab.double(rss_train.tolist())
	cb_train_cur_amp = matlab.double(np.abs(cb_train).tolist())
	cb_train_cur_angle = matlab.double(np.angle(cb_train).tolist())

	'''
		Start phaselift estimation 
	'''
	print('Start phaselift Estimation')
	H_amp_cur, H_angle_cur = eng.channel_recovery_ADMM_v2_simulation_phaselift(int(num_ant), int(num_ant), cb_train_cur_amp, cb_train_cur_angle, rss_train, eng.double(3), nargout=2)
	H_amp_cur = np.array(H_amp_cur)
	H_angle_cur = np.array(H_angle_cur)            
	H_recover_phaselift =  np.squeeze(H_amp_cur*np.exp(1j*H_angle_cur))

################################################################################################################################################################################################
	'''
		start multires probing
	'''
	print('Start random probing for A2 and Multiresolution')
	try:
		for i in range(num_round_multires):
			os.system('echo "123456" | sudo -S -k cp ../codebook/codebook_brd/multires_' + num_ant + 'ant_rx/wil6210_rx_cb' + str(i+1) + '.brd /lib/firmware/wil6210.brd')
			os.system('echo "123456" | sshpass -p "123456" scp ../codebook/codebook_brd/multires_' + num_ant + 'ant_tx/wil6210_tx_cb' + str(i+1) + '.brd utwncg@' + ap_ipaddr + ':wil6210.brd')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo cp ~/wil6210.brd /lib/firmware/wil6210.brd"')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 down"')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 up"')
			os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
			os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
			time.sleep(1)
			
			dump_count = 0
			cur_rss = fetch_rss(ip, port, ap_ipaddr, dump_count)
			cur_rss = np.median(cur_rss, axis = 0)

			cur_rss[cur_rss>1000] = 0
			cur_rss = cur_rss*0.0652-74.3875
			print(cur_rss)
			rss_multires[i,:] = cur_rss[2:]

			print('\nFinished round ' + str(i+1) + '/'  + str(num_round_multires) + '\n')

			#check phased array and baseband chipset temperature
			try:
				content = open('/sys/kernel/debug/ieee80211/phy0/wil6210/temp').readlines()
				mac_temp = float(content[0].split(' ')[-1].strip('\n'))
				radio_temp = float(content[1].split(' ')[-1].strip('\n'))
				if mac_temp > 70 or radio_temp > 62.5:
					os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
					time.sleep(20)
					os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
				else:
					print('Temperature test passes\n')
					print('Current Temperature is MAC: ' + str(mac_temp) + ' RADIO: ' + str(radio_temp) + '\n')
			except Exception as e:
				print(e)

		sio.savemat(file_out,{'rss':rss, 'rss_directional':rss_directional, 'rss_multires':rss_multires, 'rss_sweeping_phi':rss_sweeping_phi, 'rss_sweeping_thetaNphi':rss_sweeping_thetaNphi})

	except (Exception, np.AxisError) as e:
		print(e)

	print('start sleeping!')
	#time.sleep(10)
	resolution_thresh = [1984, 3968, 3968]
	rss_train_multires = np.reshape(rss_multires,[-1,1])
	cb_multires = sio.loadmat(multires_codebook_path)['cb']

	# 90 degree
	calibration_bit = np.array([0,2,3,0,0,3,0,3,1,0,0,3,0,3,0,0])
	calibration_bit2 = np.array([0,2,3,0,0,3,0,3,1,0,0,3,0,3,0,0])
	calibration_bit3 = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])	

	for _ in range(500):
		rand_ind = list(np.random.permutation(resolution_thresh[0]))
		rss_train_multires[:resolution_thresh[0],:] = rss_train_multires[rand_ind,:]
		cb_multires[:resolution_thresh[0],:] = cb_multires[rand_ind,:]

		rand_ind = list(resolution_thresh[0]+np.random.permutation(resolution_thresh[1]))
		rss_train_multires[resolution_thresh[0]:(resolution_thresh[0]+resolution_thresh[1]),:] = rss_train_multires[rand_ind,:]
		cb_multires[resolution_thresh[0]:(resolution_thresh[0]+resolution_thresh[1]),:] = cb_multires[rand_ind,:]

		rand_ind = list(resolution_thresh[0]+resolution_thresh[1]+np.random.permutation(resolution_thresh[2]))
		rss_train_multires[(resolution_thresh[0]+resolution_thresh[1]):,:] = rss_train_multires[rand_ind,:]
		cb_multires[(resolution_thresh[0]+resolution_thresh[1]):,:] = cb_multires[rand_ind,:]

	rss_train = rss_train_multires[5952:,:]
	cb_train = cb_multires[5952:,:]

	rss_train_ori = copy.deepcopy(rss_train)
	rss_train = matlab.double(rss_train.tolist())
	cb_train_cur_amp = matlab.double(np.abs(cb_train).tolist())
	cb_train_cur_angle = matlab.double(np.angle(cb_train).tolist())

	cb_train_multires_cur_amp = matlab.double(np.abs(cb_multires).tolist())
	cb_train_multires_cur_angle = matlab.double(np.angle(cb_multires).tolist())
	rss_train_multires = matlab.double(rss_train_multires.tolist())

################################################################################################################################################################################################

	'''
	ACO
	'''
	print('collect ACO')
	# set rx side as quasi-omnidirectional
	os.system('echo "123456" | sudo -S -k cp ../codebook/wil6210_sparrow_plus_rx.brd /lib/firmware/wil6210.brd')
	os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
	os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')

	# get Tx side CSI
	h_tx = collect_ACO_tx([ip], ap_ipaddr, port)[:, active_ant]
	wt_aco = get_ACO_codebook_bit(h_tx)

	brd_name_tx = generate_tx_codebook(wt_aco, active_ant)

	#get Rx side CSI
	h_rx = collect_ACO_rx(ip, ap_ipaddr, port, brd_name_tx, num_ant = 16)
	wr_aco = get_ACO_codebook_bit(h_rx)

	print([wt_aco, wr_aco])

	H_ACO = np.concatenate([h_tx,h_rx], axis =0)

################################################################################################################################################################################################
	'''
		start beamforming
	'''
	print('Start A2/multiresolution estimation and beamforming')
	for r in range(repeat):
		H_amp_cur, H_angle_cur = eng.channel_recovery_ADMM_v2_simulation_A2only(int(num_ant), int(num_ant), cb_train_cur_amp, cb_train_cur_angle, rss_train, eng.double(r+1), nargout=2)
		H_amp_cur = np.array(H_amp_cur)
		H_angle_cur = np.array(H_angle_cur)            
		H_recover_random = np.squeeze(H_amp_cur*np.exp(1j*H_angle_cur))

		H_amp_cur, H_angle_cur = eng.channel_recovery_ADMM_v2_simulation_multiresolution(int(num_ant), int(num_ant), cb_train_multires_cur_amp, cb_train_multires_cur_angle, rss_train_multires, eng.double(r+1), nargout=2)
		H_amp_cur = np.array(H_amp_cur)
		H_angle_cur = np.array(H_angle_cur)            
		H_recover_multires = np.squeeze(H_amp_cur*np.exp(1j*H_angle_cur))

		H_amp_cur, H_angle_cur = eng.channel_recovery_ADMM_v2_simulation_A2nuclear(int(num_ant), int(num_ant), cb_train_cur_amp, cb_train_cur_angle, rss_train, eng.double(r+1), nargout=2)
		H_amp_cur = np.array(H_amp_cur)
		H_angle_cur = np.array(H_angle_cur)            
		H_recover_nuclear = np.squeeze(H_amp_cur*np.exp(1j*H_angle_cur))
		
		print(np.shape(H_ACO))
		print(np.shape(H_recover_random))
		print(np.shape(H_recover_multires))
		print(np.shape(H_recover_directional))
		print(np.shape(H_recover_phaselift))

		H_recover_temp = np.concatenate([np.expand_dims(H_recover_multires,1), np.expand_dims(H_recover_random,1), np.expand_dims(H_recover_nuclear,1), np.expand_dims(H_recover_phaselift,1)], axis = 1)

		H_recover[r, :, : ,:] = np.squeeze(H_recover_temp)

		for i in range(len(M)):
			if i in [0,1]:
				brd_name_rx, brd_name_tx = codebook_generator(np.squeeze(H_recover_temp[i,:,:]), np.squeeze(H_recover_directional[i,:,:]), wt_aco, wr_aco, active_ant, int(num_ant), int(num_ant), compensation = calibration_bit)
			elif i in [2,3]:
				brd_name_rx, brd_name_tx = codebook_generator(np.squeeze(H_recover_temp[i,:,:]), np.squeeze(H_recover_directional[i,:,:]), wt_aco, wr_aco, active_ant, int(num_ant), int(num_ant), compensation = calibration_bit2)
			else:
				brd_name_rx, brd_name_tx = codebook_generator(np.squeeze(H_recover_temp[i,:,:]), np.squeeze(H_recover_directional[i,:,:]), wt_aco, wr_aco, active_ant, int(num_ant), int(num_ant), compensation = calibration_bit3)

			os.system('echo "123456" | sshpass -p "123456" scp ' + brd_name_tx + ' utwncg@' + ap_ipaddr + ':wil6210.brd')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo cp ~/wil6210.brd /lib/firmware/wil6210.brd"')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 down"')
			os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 up"')
			time.sleep(1)

			count = 0
			for brd in brd_name_rx:
				os.system('echo "123456" | sudo -S -k cp ' + brd + ' /lib/firmware/wil6210.brd')
				os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
				os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
				time.sleep(0.1)

				dump_count = 0
				cur_rss = fetch_rss(ip, port, ap_ipaddr, dump_count)
				cur_rss = np.median(cur_rss, axis = 0)

				cur_rss[cur_rss>1000] = 0
				cur_rss = cur_rss*0.0652-74.3875
				print(cur_rss)
				compare_bf_rss[r, i, count, :] = cur_rss
				count += 1

	sio.savemat(bfresult_path,{'rss_bf':compare_bf_rss, 'rss_sweeping_phi':rss_sweeping_phi, 'rss_sweeping_thetaNphi':rss_sweeping_thetaNphi, 'rss_directional': rss_directional, 'H_recover':H_recover, 'H_recover_directional': H_recover_directional,'H_ACO': H_ACO})