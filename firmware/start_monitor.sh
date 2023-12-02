############################################################################
#Project 2ACE
#PI: Prof. Lili Qiu @ UT Austin/MSRA
#Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
#Copyright @ The University of Texas at Austin, 2023
#Partially inherited from Teng Ghufran Baig @ UT Austin
############################################################################

sudo modprobe -vvv wil6210
sudo rfkill unblock all
sudo ifconfig wlp3s0 down
sudo iw wlp3s0 set type monitor
sudo iw wlp3s0 set monitor control
sudo ifconfig wlp3s0 up
sudo iw wlp3s0 set power_save off
sudo iw wlp3s0 set freq 60480
sudo ifconfig wlp3s0 down
sudo ifconfig wlp3s0 up
