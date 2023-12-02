############################################################################
#Project 2ACE
#PI: Prof. Lili Qiu @ UT Austin/MSRA
#Participants: Yiwen Song @ CMU, Changhan Ge @ UT Austin, Yin Zhang
#Copyright @ The University of Texas at Austin, 2023
#Partially inherited from Teng Ghufran Baig @ UT Austin
############################################################################

sudo fuser -k 10002/tcp
sudo cp ./patched_fw/wil6210.fw /lib/firmware/wil6210.fw
echo "123456" | sudo -S -k sh ./start_monitor.sh
sudo python2 ./wil6210_server-2.2.0 10002
