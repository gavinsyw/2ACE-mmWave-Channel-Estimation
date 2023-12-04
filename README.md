# MobiHoc23 2ACE Source Code
Source Code of "2ACE: Spectral Profile-driven Multi-resolutional Compressive Sensing for mmWave Channel Estimation" (ACM MobiHoc ’23)

#### PI: [Dr. Lili Qiu](https://www.cs.utexas.edu/~lili/) @ The University of Texas at Austin & Microsoft Research Asia Shanghai <br>
#### Developer: [Yiwen Song](https://gavinsyw.github.io/) @ Carnegie Mellon University, [Changhan Ge]() @ The University of Texas at Austin, Dr. Yin Zhang

Abstract: Channel estimation is critical to millimeter-wave capability. Unlike sub-6 GHz WiFi, commercial-off-the-shelf 60 GHz WiFi devices adopt single RF-chain and can only report the combined received signal strength (RSS) instead of the antenna-wise channel state information (CSI). Therefore, recovering the CSI using a limited number of RSS measurements is important but faces the following challenges: (i) solving a non-convex objective is hard and computationally heavy, (ii) the estimation error is high with insufficient RSS measurements, and (iii) channel fluctuates dynamically. To jointly tackle them, we propose 2ACE, an Accelerated and Accurate Channel Estimation approach using spectral profile-driven multiresolutional compressive sensing. Our thorough experiments show that 2ACE yields 2-8 dB reduction in CSI estimation error, 1-5 dB improvement in beamforming performance, and $5^{\circ} -10^{\circ}$ reduction in angle-of-departure estimation error over the existing schemes. <br>

Read our paper [here](https://doi.org/10.1145/3565287.3610252)

# Warning & Disclaimer
Hardware modifications in this design will void your devices' warranty. The code requires root access to the operating system and will hijack the network interface, which may cause data leakage and irreversible damage to your OS and hardware. Users try our code at their own risk and responsibility. Meanwhile, we are not responsible for problems caused by third party software involved in our code.

# Simulation

### Generate mmWave CSI trace using Wireless Insite
* Wireless Insite version tested: v3.4.4.14
* A valid Wireless Insite license (with X3D feature) is required. Please contact [Remcom](https://www.remcom.com/wireless-insite-em-propagation-software) for details.

# Testbed 
### Equipments
* Devices: Two Acer TravelMate P658-m serving as AP (TX) and STA (RX), respectively.
* Network Card: [NGFF595A-L](https://lian-mueller.com/qualcomm-atheros-communications-distributor/ngff595a-l.html) Qualcomm Atheros QCA6320-based 802.11ad WiGig Baseband M.2 Module (Sparrrow).
  + Qualcomm Atheros QCA6335 based [NGFF695A-L](https://lian-mueller.com/qualcomm-atheros-communications-distributor/ngff695a-l.html) (Sparrow Plus) in [M3](http://m3.ucsd.edu/sdr/) should also work but we never tried.
* Phased Array: Qualcomm Atheros QCA6310-based 802.11ad WiGig RF Antenna Module (32-antenna URA, the same one as the single antenna array module used in [M3](http://m3.ucsd.edu/sdr/) and [MultiLoc](https://dl.acm.org/doi/pdf/10.1145/3498361.3538920)).
  + Acer TravelMate P658-m orginally comes with [NGFF595A-L-Ant](https://lian-mueller.com/qualcomm-atheros-communications-distributor/ngff595a-l-ant.html) (32-antenna irregular array, the same one used in [ACO](https://dl.acm.org/doi/abs/10.1145/3241539.3241576)). We replaced the phased array with a 32-ant URA for a better evaulation. 
  + 2ACE should be able to work with any phased array, including NGFF595A-L-Ant, since it does not consider any antenna structure. However, basedline approaches [PLOMP/PLGMP](https://github.com/yzhang417/SANBA-mmWave-SDR)/azimuth sweeping/azimuth+elevation sweeping may not work with NGFF595A-L-Ant since they require antenna array structure information to construct directional beams.

### Dependencies
* OS tested: Ubuntu 16.04/18.04/20.04 (Ubuntu 20.04 will need to install python2 to support the legacy features).
* ```OpenSSH``` server/client and ```sshpass``` required on both machines.
* Python 3.8+ environment required, need to install numpy and scipy packages.
* Software required: Matlab (tested with 2021a, 2021b, and 2022a)
  + A valid matlab license is required.
  + require signal processing and phased array toolboxes
  + require matlab engine for python, detailed installation instruction can be found [here](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)
* A custom patched wil6210 firmware (driver of 802.11ad card) for ```wil6210.fw``` version 6.2.0.40 and a phased-array controller ```wil6210_server-2.2.0``` server are required, which are the same ones used in [X-array](https://dl.acm.org/doi/10.1145/3372224.3380882). These two modules togehter enable dumping received signal strength (RSS) from network interface. Meanwhile, a codebook module ```wil6210_brd_mod``` is needed to generate beamforming codebook. We are not the author of these third-party software and do not have the permission to release them. User have to request them from X-array.
  + Note that the "wil6210-server" used in [M3](https://github.com/RenjieZhao/M-Cube-Hostcmds/tree/main) does not work for our code since it does not support RSS dumping. The server should include ```per_beam_snr``` module.
  + Note that the patched firmware "wil6210.fw" used in [ACO](https://dl.acm.org/doi/abs/10.1145/3241539.3241576) does not work for our code.
* Your AP and STA should be connected via either ethernet port or sub-6 WiFi link, and the corresponding network interfaces should have IPv4 addresses under the same subnet, say ```AP: 192.168.0.10/24, STA: 192.168.0.66/24```.
* Do not connect your STA and AP via 60GHz WiGig link. During the excution of code, the WiGig interface will be hijacked and no TCP link would be able to setup between STA and AP. 

### Detailed steps
#### Make sure you already read throught the warning/disclaimer, and all the detailed steps before your experiment.
#### Preliminary (Optional but recommended)
* Hardware phase offset calibration: If you are not using the same phased array as ours, you will have to calibrate the hardware phase offset on each antenna following the approach in [M3](http://m3.ucsd.edu/sdr/) or [PLOMP/PLGAMP](https://dl.acm.org/doi/abs/10.1145/3323679.3326532). After the calibration, replace the ```./codebook/hardware_phaseoffset.mat``` file with your own measurement.
* RSSI to RSS mapping: Commercial WiGig card normaly only report an integer received signal strength indicator (RSSI) instead of absolute received signal strength (RSS) in dBm unit. However, 2ACE and all baselines will need RSS to recover CSI matrix. If you are not using the same baseband chipset NGFF595A-L as ours, you will have to calibrate the function for mapping RSSI to RSS.

#### Run with our codebook
* Make sure you already calibrate your experiment settings following the preliminary.
* Install all the dependencies on your AP and STA machine.
* Clone this repo to your STA machine by
  ```
  git clone https://github.com/ChanghanGe/MobiHoc23_2ACE_Artifact.git
  cd MobiHoc23_2ACE_Artifact
  ```
* You have to figure out the name of your WiGig interface first. You can check all your network interface by hit ```ifconfig -a```. In our experiment, the name of the WiGig interface is ```wlp3s0```. It varies by machine. If your WiGig interface is not ```wlp3s0```, you have to change it in ```./firmware/ap_mmwave.conf```,  ```./main/main.py``` and ```./firmware/load_csi_firmware.sh```
* One of your machine have to be on AP mode. Normally wireless interfaces of a linux machine are set as STA mode by default. To make it as an AP, run
  ```
  sudo cp ./firmware/ap_mmwave.conf /etc/wpa_supplicant/
  ```
  Then reboot your AP. It should then work at channel 2 of 802.11ad (60.48GHz). This code works for linux kernel version 4.15.0-142.
  If the rebooting does not make your laptop an AP, please ask ChatGPT or Stackoverflow.
* Put the ```wil6210_brd_mod``` under the ```./codebook/``` folder, put the patched ```wil6210.fw``` under the ```./firmware/patched_fw/``` folder and put the patched ```wil6210_server-2.2.0``` under the ```./firmware/``` folder. Then run the following commands.
  ```
  chmod a+x ./codebook/wil6210_brd_mod
  chmod a+x ./firmware/patched_fw/wil6210.fw
  chmod a+x ./firmware/wil6210_server-2.2.0
  ```
* Copy ```./codebook/codebook_brd/ACO``` folder to your AP machine and put it under ```~/Documents/```
* Open ```main.py``` in ```./main/``` folder. Change the ```ap_ipaddr, username, password``` variable to your own AP's ip address, username and password, respectively.
* Make sure the 10002 port of STA is idle. Otherwise, please change ```port``` variable in ```main.py``` and ```./firmware/load_csi_firmware.sh```. 
* Open a separate terminal and go to ```./frimware/```. Then ```sh load_csi_firmware.sh```. The server will then start to listen to command at the designated port.
* Run ```python3 main.py``` in ```./main/```.
* Wait for a while and the result will be stored in ```./result/```

#### Run with your own codebook
* This requires the codebook module ```wil6210_brd_mod``` mentioned above.
* After you get this module, play with the code inside ```./codebook/```.

# Q&A
* Why does it take too long to run the code? Why does solving the CSI matrix seem to take forever?
  + 2ACE solve the CSI matrix very fast, normally within a second for 16x16 CSI matrices. However, the baseline approach PhaseLift normally get stuck at somewhere during convex optimization. Meanwhile, PLOMP and PLGAMP also take long since their first step is PhaseLift. The larger the CSI matrix, the longer the algorithms take to solve it. PhaseLift/PLGAMP/PLOMP will take hours to solve a 32x32 CSI matrix with 4096 training probes.
  + For a fast run, you may disable the code related to PLOMP/PLGAMP/PhaseLift in simulation main.m, and comment the corresponding blocks in ```main.py``` for testbed.
* Why does the RSS reading suddenly stop? What should I do in this situation?
  + RSS dump failure has multiple root causes. Error message "This dump bloody failed!!!!! Try bloody again!!!!!" is a very common one. A simple reason might be that the network interface kernel is not ready while the RSS reading command is triggered.
  + Error message "Serious Dumping Error! Restart kernel and firmware controller!" means you are in big trouble.
    - If ```SIOCSIFFLAGS: Timer expired``` pops up at the same time, the error is highly caused by a corrupted codebook file ```.brd``` on either AP or STA side. 
    - If ```NOT_SUPPORTED_EVENTID detected!!``` or ```event 0x1900 no response.``` are displayed, the error is highly like caused by a wrong version of firmware. Our testbed code only works with firmware version 6.2.0.40 with a custom patch.
    - In some situation, like overheating, the NIC will shut off the WiGig interface automatically for self-protection. 
  + Error message "Connection Failed ! Restart the dump server and press Enter to continue..." means the controller server ```wil6210_server-2.2.0``` breaks unexpectedly or you forgot to start it. Resart the controller server will help.
  + The authors were as desperate as you when they saw these errors. You can feel it.
* Why does RSS value fluctuates a lot, even if with the same codebook？
  + mmWave link over commercial WiGig device is vulnerable to noise in the environment. As we observed, the RSS of the same codebook obtained under the same experiment setting can differ by as high as 2dB. The noise comes from multiple sources and fluctuates with temperature and humidity from the experiment environment and also the temperature on board (phased array and basedband chipset will overheat and the measurement quality goes down). We already made some efforts on this issue. There is a part of the code monitoring the temperature of phased array and basedband chipset by reading sensor data from linux kernel ```/sys/kernel/debug/ieee80211/phy0/wil6210/temp```. If the onboard temperature is too high (mac_temp > $70^{\circ}$ C or radio_temp > $62.5^{\circ}$ C), code will sleep for 20s to let the chipset cool down. Meanwhile, when we were doing hardware modifications, the antenna arrays are mounted on a large aluminium board with thermal paste in between. 
  + One way to mitigate this issue is to dump the RSS many times and take the average (mean after removing outlier). You may increase the number of measurements in ```fetch_rss``` function in ```./main/codebook_library.py```.
  + Migrate our design to a fancy high-end SDR.
  + You may read more about hardware overheating issue [here](https://ieeexplore.ieee.org/document/9259381).
* Can we try different number of antennas and larger CSI matrix, say 32x32?
  + Yes, you can try whatever number of antennas in simulation. However, we do not recommend to try 32 antenna setting on testbed. 
    - As we mentioned above, the larger the matrix, the longer the algorithms take to solve it.
    - Open all 32 antennas on phased array, especially on the RX side, will heat up the phased array and baseband chipset very quickly. The URA we used on RX side can only sustain for a few minutes with all 32 antennas open. Open all 32 antenna on NGFF595A-L-Ant will kill the WiGig kernel immediately and may burn out the phased array.
* Why suddenly nothing works and the OS become slow and less responsive?
  + Reboot your machine. If problems still exist, unfortunately your OS environment might be corrupted and you may need to reinstall the OS.
  + There is a possibility that the baseband chipset or phased array burn out. Good luck.
* The codes may exhibit a poor structure with limited comments. Please accept the apology from the authors in advance.

# Cite 2ACE
#### Bibtex
  ```
  @inproceedings{10.1145/3565287.3610252,
  author = {Song, Yiwen and Ge, Changhan and Qiu, Lili and Zhang, Yin},
  title = {2ACE: Spectral Profile-driven Multi-resolutional Compressive Sensing for mmWave Channel Estimation},
  year = {2023},
  isbn = {9781450391658},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3565287.3610252},
  html={https://dl.acm.org/doi/10.1145/3565287.3610252},
  doi = {10.1145/3565287.3610252},
  booktitle = {Proceedings of the Twenty-Fourth International Symposium on Theory, Algorithmic Foundations, and Protocol Design for Mobile Networks and Mobile Computing (ACM MobiHoc '23)},
  location = {Washington DC, USA},
  series = {MobiHoc '23},
  }
  ```

#### ACM Reference Format
Yiwen Song, Changhan Ge, Lili Qiu, Yin Zhang. 2023. 2ACE: Spectral Profile-driven Multi-resolutional Compressive Sensing for mmWave Channel Estimation. In Twenty-fourth International Symposium on Theory, Algorithmic Foundations, and Protocol Design for Mobile Networks and Mobile Computing (ACM MobiHoc '23), October 23–26, 2023, Washington, DC, USA. ACM, New York, NY, USA, 10 pages. https://doi.org/10.1145/3565287.3610252

# Acknowledgement
This work is supported in part by NSF Grant [CNS-2008824](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2008824&HistoricalAwards=false) and [CNS-2107037](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2107037&HistoricalAwards=false). We appreciate the insightful feedback from ACM MobiHoc 2023 anonymous reviewers.

#### The code in this repository is partly inherited from the following sources:
* Yi Zhang, Kartik Patel, Sanjay Shakkottai, and Robert W. Heath Jr. 2019. Side-information-aided Non-coherent Beam Alignment Design for Millimeter Wave Systems. In MobiHoc '19: The Twentieth ACM International Symposium on Mobile Ad Hoc Networking and Computing (ACM MobiHoc '19), July 02-05, 2019, Catania, Italy. ACM, New York, NY, USA, 10 pages. https://github.com/yzhang417/SANBA-mmWave-SDR

* Renjie Zhao, Timothy Woodford, Teng Wei, Kun Qian, and Xinyu Zhang. 2020. M-Cube: a millimeter-wave massive MIMO software radio. In Proceedings of the 26th Annual International Conference on Mobile Computing and Networking (ACM MobiCom '20). Association for Computing Machinery, New York, NY, USA, Article 15, 1–14. DOI:https://doi.org/10.1145/3372224.3380892

* Song Wang, Jingqi Huang, Xinyu Zhang, Hyoil Kim, and Sujit Dey. 2020. X-Array: approximating omnidirectional millimeter-wave coverage using an array of phased arrays. In Proceedings of the 26th Annual International Conference on Mobile Computing and Networking (ACM MobiCom '20). Association for Computing Machinery, New York, NY, USA, Article 5, 1–14. https://doi.org/10.1145/3372224.3380882

* Kun QIAN, Xinyu ZHANG, Zheng YANG, Yunhao LIU, AoD-based localization with cots millimeter-wave devices,  SCIENTIA SINICA Informationis, Volume 51, Issue 1, 2021, Pages 122-, ISSN 1674-7267, https://doi.org/10.1360/SSI-2019-0135.
(http://www.sciengine.com/doi/10.1360/SSI-2019-0135)

We greatly appreciate the sources above but we are not responsible for any potential problem caused by them.

# Contact
If you have any question regarding 2ACE, please contact Yiwen Song at yiwensong@cmu.edu.
