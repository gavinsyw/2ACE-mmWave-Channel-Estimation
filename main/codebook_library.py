############################################################################
#Project 2ACE
#PI: Prof. Lili Qiu @ UT Austin/MSRA
#Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
#Copyright @ The University of Texas at Austin, 2023
#Partially inherited from Teng Wei/Song Wang/Jingqi Huang/Xinyu Zhang @ UCSD
############################################################################

import numpy as np
import scipy.io as sio
import copy
import json
import os
import socket
import ast
from multiprocessing import Pool
import subprocess as sp
import struct
import time

def set_beam_num(file_name, beam_num):
    try:
        sp.check_call(['../codebook/wil6210_brd_mod', '-set_beam_num=%d' % beam_num, \
            '-in_brd=%s' % file_name, '-out_brd=%s' % file_name])
    except Exception as e:
        print(e)

def set_beam(file_name, module_index, sector_type, sector_index, gain, phase, amplify):
    try:
        sp.check_call(['../codebook/wil6210_brd_mod', '-set_pattern', \
            '-in_brd=%s' % file_name, '-out_brd=%s' % file_name, '-module_index=%d' % module_index, '-sector_type=%d' % sector_type, '-sector_index=%d' % sector_index, \
            '-mag=%s' % gain, '-phase=%s' % phase, '-amp=%s' % amplify])
    except Exception as e:
        print(e)

def enable_modules(file_name, module_index, isTx):
    try:
        for m_idx in module_index:
            sp.check_call(['../codebook/wil6210_brd_mod','-in_brd=%s' % file_name, '-out_brd=%s' % file_name, '-enable_module', '-module_index=%d' % m_idx, '-sector_type=%d' % isTx])
    except Exception as e:
        print(e)

def disable_modules(file_name, module_index, isTx):
    try:
        for m_idx in module_index:
            sp.check_call(['../codebook/wil6210_brd_mod','-in_brd=%s' % file_name, '-out_brd=%s' % file_name, '-disable_module', '-module_index=%d' % m_idx, '-sector_type=%d' % isTx])
    except Exception as e:
        print(e)

def twos_comp(val, bits):
    """compute the 2's complement of int value val"""
    for i in range(len(val)):
        if (val[i] & (1 << (bits - 1))) != 0: # if sign bit is set e.g., 8bit: 128-255
            val[i] = val[i] - (1 << bits)       # compute negative value
    return val  

def svd_beamformer(H):
    ant_tx, ant_rx = np.shape(H)
    u, s, wr_opt = np.linalg.svd(H)
    u, s, wt_opt = np.linalg.svd(np.transpose(H))
    
    wr_quant_angle = -np.transpose(np.around(np.angle(wr_opt)/(np.pi/2))*(np.pi/2))
    wt_quant_angle = -np.transpose(np.around(np.angle(wt_opt)/(np.pi/2))*(np.pi/2))
    
    wr_quant = np.exp(1j*wr_quant_angle)
    wt_quant = np.exp(1j*wt_quant_angle)
    
    rss = []
    beam_idx = []
    for i in range(ant_tx):
        for j in range(ant_rx):
            signal = np.abs(np.matmul(wt_quant[:,i],np.matmul(H,wr_quant[:,j])))**2
            rss.append(10*np.log10(signal*1000))
            beam_idx.append([i,j])
    
    idx = np.argmax(np.array(rss))
    tx_idx = beam_idx[idx][0]
    rx_idx = beam_idx[idx][1]
            
    wr_quant_out = wr_quant[:,rx_idx]
    wr_quant_out = np.around((np.angle(wr_quant_out))/(np.pi/2))
    wr_quant_out[wr_quant_out<0] = wr_quant_out[wr_quant_out<0]+4
    wr_quant_out[wr_quant_out==4] = 0
    wr_out = ''
    for bit in wr_quant_out:
        wr_out = wr_out + str(int(bit))
        
    wt_quant_out = wt_quant[:,tx_idx] 
    wt_quant_out = np.around((np.angle(wt_quant_out))/(np.pi/2))
    wt_quant_out[wt_quant_out<0] = wt_quant_out[wt_quant_out<0]+4
    wt_quant_out[wt_quant_out==4] = 0
    wt_out = ''
    for bit in wt_quant_out:
        wt_out = wt_out + str(int(bit))
    return wr_out, wt_out

def svd_beamformer_compensation(H, offset):
    ant_tx, ant_rx = np.shape(H)
    u, s, wr_opt = np.linalg.svd(H)
    u, s, wt_opt = np.linalg.svd(np.transpose(H))
    
    wr_quant_angle = -np.transpose(np.around(np.angle(wr_opt)/(np.pi/2))*(np.pi/2))
    wt_quant_angle = -np.transpose(np.around(np.angle(wt_opt)/(np.pi/2))*(np.pi/2))
    
    wr_quant = np.exp(1j*wr_quant_angle)
    wt_quant = np.exp(1j*wt_quant_angle)
    
    rss = []
    beam_idx = []
    for i in range(ant_tx):
        for j in range(ant_rx):
            signal = np.abs(np.matmul(wt_quant[:,i],np.matmul(H,wr_quant[:,j])))**2
            rss.append(10*np.log10(signal*1000))
            beam_idx.append([i,j])
    
    idx = np.argmax(np.array(rss))
    tx_idx = beam_idx[idx][0]
    rx_idx = beam_idx[idx][1]
            
    wr_quant_out = wr_quant[:,rx_idx]*np.exp(-1j*offset)
    wr_quant_out = np.around((np.angle(wr_quant_out))/(np.pi/2))
    wr_quant_out[wr_quant_out<0] = wr_quant_out[wr_quant_out<0]+4
    wr_quant_out[wr_quant_out==4] = 0
    wr_out = ''
    for bit in wr_quant_out:
        wr_out = wr_out + str(int(bit))
        
    wt_quant_out = wt_quant[:,tx_idx]*np.exp(-1j*offset)
    wt_quant_out = np.around((np.angle(wt_quant_out))/(np.pi/2))
    wt_quant_out[wt_quant_out<0] = wt_quant_out[wt_quant_out<0]+4
    wt_quant_out[wt_quant_out==4] = 0
    wt_out = ''
    for bit in wt_quant_out:
        wt_out = wt_out + str(int(bit))
    return wr_out, wt_out

def generate_tx_codebook(pha, active_ant):

    gp_amp = '66666666'

    amp = ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
    phase = ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
    for idx, ant in enumerate(active_ant):
        phase[ant] = pha[idx]
        amp[ant] = '7'

    mag = [''.join(amp) for _ in range(64)]
    pha = [''.join(phase) for _ in range(64)]
    #write codebook
    brd_name_tx = '../codebook/codebook_brd/optimized_codebook/wil6210_optimized_tx.brd'
    beam_num_tx = len(pha)

    # clear previous tx codebook entrys
    for sec_idx in range(beam_num_tx):
        set_beam(brd_name_tx, 0, 1, sec_idx, '00000000000000000000000000000000', '00000000000000000000000000000000', '00000000')

    # clear previous tx codebook entrys
    for sec_idx in range(beam_num_tx):
        print('Module:%d, Sector:%d' % (0, sec_idx))
        set_beam(brd_name_tx, 0, 1, sec_idx, '40000000000000000000000000000000', '00000000000000000000000000000000', gp_amp)

    return brd_name_tx

def collect_ACO_rx(ip, ap_ipaddr, port, brd_name_tx, num_ant = 16):
    os.system('echo "123456" | sshpass -p "123456" scp ' + brd_name_tx + ' utwncg@' + ap_ipaddr + ':wil6210.brd')
    os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo cp ~/wil6210.brd /lib/firmware/wil6210.brd"')
    os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 down"')
    os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 up"')
    time.sleep(1)

    rss = np.zeros((num_ant*4, 1))
    for i in range((num_ant-1)*4):
        print('Now probing ' + str(i+1) + 'th codebook for rx side ACO')
        os.system('echo "123456" | sudo -S -k cp ../codebook/codebook_brd/ACO_rx/wil6210_rx' + str(i+1) + '.brd /lib/firmware/wil6210.brd')
        os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
        os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
        time.sleep(0.5)

        dump_count = 0
        cur_rss = fetch_rss(ip, port, ap_ipaddr, dump_count)
        cur_rss = np.median(cur_rss, axis = 0)

        cur_rss[cur_rss>1000] = 0
        cur_rss = cur_rss*0.0652-74.3875
        print(cur_rss)
        rss[i+4,:] = cur_rss[15]

    rss = np.expand_dims(np.squeeze(np.mean(rss, axis = 1)), axis = 0)
    csi = rss2csi(rss)
    return np.expand_dims(csi, axis = 0)

def codebook_generator(H_est, H_directional, wt_aco, wr_aco, active_ant, num_tx_ant, num_rx_ant,compensation = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])):

    wr = []
    wt = []
    for i in range(len(H_est)):
        H = np.reshape(H_est[i,:],[num_tx_ant,num_rx_ant])
        if i == 0:
            wr_cur, wt_cur = svd_beamformer_compensation(H, compensation*(np.pi/2))
        else:
            wr_cur, wt_cur = svd_beamformer(H)
        wr.append(wr_cur)
        wt.append(wt_cur)

    for i in range(len(H_directional)):
        H = np.reshape(H_directional[i,:],[num_tx_ant,num_rx_ant])
        wr_cur, wt_cur = svd_beamformer(H)
        wr.append(wr_cur)
        wt.append(wt_cur)

    wt.append(wt_aco)
    wr.append(wr_aco)

    gp_amp = '66666666'    
    amp = ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
    phase = ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
    amp_tx = ['14444444000000001444444400000000',
              '44144144000000004444444100000000',
              '44444144000000004444414400000000', 
              '14444144000000004444414400000000',
              '44444444000000001444444400000000',
              '55505055000000000555505500000000',
              '44141444000000001444444400000000',
              '05555015000000000555505500000000',
              '44441444000000001414444400000000',
              '44444141000000004444444400000000',
              '44144444000000004414444400000000',
              '44441444000000004444444400000000',
              '44144444000000001444444100000000',
              '44444441000000004444414400000000',
              '44441444000000001444444400000000',
              '44144444000000001444444400000000',
              '41144444000000004444444400000000',
              '44444444000000004444444100000000',
              '44444444000000004444444400000000',
              '44441444000000001444444400000000',
              '44144444000000004444414400000000',
              '44144444000000001414444400000000',
              '44441444000000001414444400000000',
              '44444441000000004444414400000000',
              '44444441000000004444414400000000',
              '44441444000000004444444400000000',
              '44444441000000004444414100000000',
              '44444444000000004444414100000000',
              '44444444000000004444414400000000',
              '14444444000000001444414400000000',
              '14444144000000004444414400000000',
              '44141444000000004444414400000000',
              '44444444000000001414444400000000',
              '44444444000000004444114100000000',
              '44441444000000001444444400000000',
              '44441444000000001444444400000000',
              '44444444000000004444414100000000',
              '44444444000000004444414400000000',
              '44444444000000001444414100000000',
              '44414444000000001444444400000000',
              '44114444000000001444444400000000',
              '14444144000000004444414400000000',
              '41444444000000001444444400000000',
              '44144144000000004444444400000000',
              '06666066000000006660601100000000',
              '44444144000000004444444400000000',
              '55115551000000005515555500000000',
              '44441444000000001444444400000000',
              '01717770000000000777017700000000',
              '44414444000000001444444400000000']

    phase_tx = ['22000121221010200031012121231323',
                '01122202113230233310230310220323',
                '22333121103233031002022210230112',
                '30012121113201202300223120221033',
                '32003222103332030320322100112031',
                '23311313123312010122131313110033',
                '02100030133320222001003033112320',
                '12233303221012313310330320221033',
                '32103320033313210213232102330133',
                '22333222133221322002033310223002',
                '02130101012102300122021223111003',
                '00221330332210311031033101220203',
                '32103223103321320223322123001310',
                '20021020123311300011001303003312',
                '03201101233320220123121100223321',
                '12233010213300100032012120231323',
                '10212111301003001301232332011213',
                '01232233123300231132123103003202',
                '00222333000222022132133222331020',
                '21032323311131330223232231223322',
                '01023131331030133012310212000322',
                '13101212223220121300232311303331',
                '20031213211130323001121110101200',
                '23010313023311303300020132333312',
                '30121021023300330011001333333202',
                '31032320201133220112132210111200',
                '23011010010200230011001332332201',
                '22300010231122013310301311221023',
                '03111333321110202032033211230213',
                '32100333203332332001033300123031',
                '23311020113201201233102013110033',
                '02123202331020130133031322110332',
                '21032112322213103101110031223032',
                '01222233013303131132123133032101',
                '20021213201023223011121110111200',
                '31033323312102000223232201333032',
                '23011020010200230022002303003212',
                '12303010221011303310301310221323',
                '32100333133320222001003333112320',
                '13211212233320220112121230222311',
                '20322313231020221223132330222311',
                '12300313302133120021031332333211',
                '32000333233332333113100010230002',
                '01122202113301203311330320231033',
                '32301012131010012223121111332211',
                '00122232000233132203133203103212',
                '30110113031112323320001301022233',
                '31033323312102000223232201333032',
                '22033132101002021213010122000303',
                '13211212233320220112121230222311']

    phase_rx = ['00000000000000000000000000000000']

    # Change Tx Amplitude and Phase:
    amp_tx_cur = copy.deepcopy(amp)
    if num_tx_ant == 8:
        #idx = [4, 6, 8, 15, 20, 22, 25, 31]
        idx = [0, 1, 4, 6, 16, 17, 20, 22]

        for ind in idx:
            amp_tx_cur[ind] = '7'

        for i in range(len(amp_tx)):
            cur_amp = list(amp_tx[i])
            for ind in range(32):
                if ind not in idx:
                    cur_amp[ind] = '0'
            amp_tx[i] = ''.join(cur_amp)

        for i in range(len(wt)):
            phase_tx_cur = copy.deepcopy(phase)
            for j in range(len(idx)):
                phase_tx_cur[idx[j]] = wt[i][j]
            phase_tx_cur = ''.join(phase_tx_cur)
            phase_tx.append(phase_tx_cur)
            phase_tx.append(phase_tx_cur)
            amp_tx.append(''.join(amp_tx_cur))
            amp_tx.append(''.join(amp_tx_cur))

    elif num_tx_ant == 16:
        idx = active_ant
        idx_directional = active_ant

        for ind in idx:
            amp_tx_cur[ind] = '7'

        for i in range(len(amp_tx)):
            cur_amp = list(amp_tx[i])
            for ind in range(32):
                if ind not in idx and ind != 0:
                    cur_amp[ind] = '0'
            amp_tx[i] = ''.join(cur_amp)

        for i in range(len(wt)):
            phase_tx_cur = copy.deepcopy(phase)
            if i < len(wt)-2:
                for j in range(len(idx)):
                    phase_tx_cur[idx[j]] = wt[i][j]
            else:
                for j in range(len(idx_directional)):
                    phase_tx_cur[idx_directional[j]] = wt[i][j]
            phase_tx_cur = ''.join(phase_tx_cur)
            phase_tx.append(phase_tx_cur)
            phase_tx.append(phase_tx_cur)
            amp_tx.append(''.join(amp_tx_cur))
            amp_tx.append(''.join(amp_tx_cur))

    elif num_tx_ant == 32:
        for i in range(len(wt)):
            amp_tx.append('77777777777777777777777777777777')
            amp_tx.append('77777777777777777777777777777777')
        phase_tx.extend(wt)
    else:
        print('Number of Tx antenna must be 8, 16 or 32')

    # Change Rx Amplitude and Phase
    amp_rx = copy.deepcopy(amp)
    if num_rx_ant == 8:
        #idx = [4, 6, 8, 15, 20, 22, 25, 31]
        idx = [0, 1, 4, 6, 16, 17, 20, 22]
        for ind in idx:
            amp_rx[ind] = '7'
        amp_rx = ''.join(amp_rx)
        for i in range(len(wr)):
            phase_rx_cur = copy.deepcopy(phase)
            for j in range(len(idx)):
                phase_rx_cur[idx[j]] = wr[i][j]
            phase_rx_cur = ''.join(phase_rx_cur)
            phase_rx.append(phase_rx_cur)
    elif num_rx_ant == 16:
        idx = active_ant
        idx_directional = active_ant

        for ind in idx:
            amp_rx[ind] = '7'
        amp_rx = ''.join(amp_rx)
        for i in range(len(wr)):
            phase_rx_cur = copy.deepcopy(phase)
            if i < len(wr)-2:
                for j in range(len(idx)):
                    phase_rx_cur[idx[j]] = wr[i][j]
            else:
                for j in range(len(idx_directional)):
                    phase_rx_cur[idx_directional[j]] = wr[i][j]
            phase_rx_cur = ''.join(phase_rx_cur)
            phase_rx.append(phase_rx_cur)
    elif num_rx_ant == 32:
        amp_rx = '77777777777777777777777777777777'
        phase_rx.extend(wr)
    else:
        print('Number of Tx antenna must be 8, 16 or 32')

    for _ in range(64-len(phase_tx)):
        amp_tx.append('77777777777777777777777777777777')
        phase_tx.append(''.join(phase))

    #write codebook
    brd_name_tx = '../codebook/codebook_brd/optimized_codebook/wil6210_optimized_tx.brd'
    brd_name_rx = ['../codebook/codebook_brd/optimized_codebook/wil6210_optimized_rx'+str(i+1)+'.brd' for i in range(len(phase_rx))]
    beam_num_tx = len(phase_tx)
    beam_num_rx = len(phase_rx)

    # clear previous tx codebook entrys
    for sec_idx in range(beam_num_tx):
        set_beam(brd_name_tx, 0, 1, sec_idx, '00000000000000000000000000000000', '00000000000000000000000000000000', '00000000')

    # clear previous tx codebook entrys
    for sec_idx in range(beam_num_tx):
        print('Module:%d, Sector:%d' % (0, sec_idx))
        set_beam(brd_name_tx, 0, 1, sec_idx, amp_tx[sec_idx], phase_tx[sec_idx], gp_amp)

    # clear previous rx codebook entrys
    for sec_idx in range(beam_num_rx):
        set_beam(brd_name_rx[sec_idx], 0, 0, 0, '00000000000000000000000000000000', '00000000000000000000000000000000', '00000000')

    # set rx entry values
    for sec_idx in range(beam_num_rx):
        print('Module:%d, Sector:%d' % (0, sec_idx))
        if sec_idx == 0:
            set_beam(brd_name_rx[sec_idx], 0, 0, 0,  '60000000000000000000000000000000', phase_rx[sec_idx], gp_amp)
        else:
            set_beam(brd_name_rx[sec_idx], 0, 0, 0, amp_rx, phase_rx[sec_idx], gp_amp)

    return brd_name_rx, brd_name_tx

def fetch_rss(server_address, port, ap_ipaddr, dump_count, num_measurement = 10):
    try:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        try:
            sock.connect((server_address, port))
        except Exception as e:
            print(e)
            input("Connection Failed ! Restart the dump server and press Enter to continue...")
            sock.connect((server_address, port))
        except ConnectionResetError as e:
            print(e)
            input("Connection Failed ! Restart the dump server and press Enter to continue...")
            sock.connect((server_address, port))

        buf_len = 10240

        data = {}
        data['cmd'] = 'per_beam_snr'
        data['args'] = {}

        m = json.dumps(data).encode('utf-8')

        measure_interval = 0.1

        info_file_input=[]
                
        for i in range(num_measurement):
            sock.sendall(m)

            d = b''
            while True:
                d_tmp = sock.recv(buf_len)
                d = d + d_tmp
                if not d_tmp or str(chr(d_tmp[-1])) == '}' or str(chr(d_tmp[-1])) == ']' or str(chr(d_tmp[-1])) == ')':
                    break

            snr_list = np.array(ast.literal_eval(d.decode("utf-8")))
            if np.any(snr_list) and len(snr_list[snr_list==0])<5:
                info_file_input.append(list(snr_list))
            time.sleep(measure_interval)
               
        info_file_input = np.array(info_file_input)
        info_file_input[info_file_input>1000] = twos_comp(info_file_input[info_file_input>1000],32)

        sock.close()

        if len(info_file_input) < 3:
            print('This dump bloody failed!!!!! Try bloody again!!!!!')
            dump_count += 1
            if dump_count > 10:
                print("Serious Dumping Error! Restart kernel and firmware controller!")
                time.sleep(30)
            os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 down"')
            os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sudo ifconfig wlp3s0 up"')
            os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
            os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
            time.sleep(1)
            info_file_input = fetch_rss(server_address, port, ap_ipaddr, dump_count)

        return info_file_input

    except Exception as e:
        print(e)

def rss2csi(rss):
    rss = np.reshape(rss,[16,4])
    csi_phase_temp = np.fft.fft(rss,4)
    csi_phase = np.angle(csi_phase_temp[:,1])
    gamma = csi_phase_temp[:,0]
    delta = np.abs(csi_phase_temp[:,1])
    csi_amp = 0.5*(np.sqrt(gamma + 2*delta) - np.sqrt(gamma - 2*delta))
    #csi = csi_amp*np.exp(1j*csi_phase);
    return np.abs(csi_amp)*np.exp(1j*csi_phase)

def collect_ACO_tx(ip, ap_ipaddr, port):
    os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sh ~/Documents/ACO/load_aco_codebook_upper.sh"')

    time.sleep(0.5)

    rss_input_args = []
    dump_count = 0
    for i in range(len(ip)):
        rss_input_args.append((ip[i], port,ap_ipaddr, dump_count))

    with Pool(len(ip)) as p:
        rss1 = p.starmap(fetch_rss,rss_input_args)

    rss1 = np.array(np.mean(np.squeeze(rss1), axis = 0))

    os.system('sshpass -p "123456" ssh utwncg@' + ap_ipaddr + ' "echo "123456" | sh ~/Documents/ACO/load_aco_codebook_lower.sh"')

    time.sleep(0.5)

    with Pool(len(ip)) as p:
        rss2 = p.starmap(fetch_rss,rss_input_args)
    
    rss2 = np.array(np.mean(np.squeeze(rss2), axis = 0))

    if len(ip) == 1:
        rss1 = np.expand_dims(rss1, axis = 0)
        rss2 = np.expand_dims(rss2, axis = 0)

    rss2[:,0] = rss1[:,3]

    rss1[rss1>1000] = 0
    rss1 = rss1*0.0652-74.3875
    rss1 = 10**(rss1/10)/1000

    csi1 = []
    for i in range(len(ip)):
        csi_temp = rss2csi(rss1[i,:])
        if len(ip) == 1:
            csi_temp = np.expand_dims(csi_temp, axis = 0)
        csi_temp[:,0] = np.sqrt(rss1[:,1])
        csi1.append(np.squeeze(csi_temp))
    csi1 = np.array(csi1)

    rss2[rss2>1000] = 0
    rss2 = rss2*0.0652-74.3875
    rss2 = 10**(rss2/10)/1000

    csi2 = []
    for i in range(len(ip)):
        csi_temp = rss2csi(rss2[i,:])
        csi2.append(np.squeeze(csi_temp))
    csi2 = np.array(csi2)

    csi = np.concatenate([csi1,csi2], axis=1)
    return csi

def get_ACO_codebook_bit(h):
    w = np.around((np.angle(h.conj()))/(np.pi/2))
    w[w<0] = w[w<0]+4
    w[w==4] = 0
    w_out = ''
    for bit in w[0]:
        w_out = w_out + str(int(bit))
    return w_out

def collect_ACO_full(num_ant, ip, port):
    csi = np.zeros((num_ant,num_ant),dtype=np.complex_)
    for i in range(num_ant):
        os.system('echo "123456" | sudo -S -k cp ./random_codebook_airfide/rx_codebook_1ant/wil6210_rx_1ant_cb' + str(i+1) + '.brd /lib/firmware/wil6210.brd')
        print('wil6210_rx_1ant_cb' + str(i+1))
        os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 down')
        os.system('echo "123456" | sudo -S -k ifconfig wlp3s0 up')
        time.sleep(1)
        csi[:,i] = np.squeeze(np.transpose(collect_ACO(ip,ap_ipaddr, port)))

    return csi