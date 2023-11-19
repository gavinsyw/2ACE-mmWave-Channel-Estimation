echo "123456" | sudo -S -k cp ~/Documents/ACO/wil6210_ACO_16ant_upper.brd /lib/firmware/wil6210.brd
echo "123456" | sudo -S -k ifconfig wlp3s0 down
echo "123456" | sudo -S -k ifconfig wlp3s0 up
