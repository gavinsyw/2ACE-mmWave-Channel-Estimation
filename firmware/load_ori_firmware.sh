sudo fuser -k 10002/tcp
sudo cp ./ori_fw/wil6210.fw /lib/firmware/wil6210.fw
echo "123456" | sudo -S -k sh ./start_monitor.sh
#sudo python2 ./wil6210_server-2.2.0 10002
