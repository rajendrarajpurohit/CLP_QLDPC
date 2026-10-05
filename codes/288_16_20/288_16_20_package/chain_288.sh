#!/bin/bash
cd /home/user1/Rajendra/Nova/frontier
while pgrep -f "certify_nolimit.py D2_10" > /dev/null; do sleep 60; done
echo "216 process gone at $(date), launching 288" >> chain_288.log
ID=W10_m56-100_1_n288k16d20w9
for r in "1 4" "5 8" "9 12" "13 16"; do
  set -- $r
  nohup python3 /home/user1/Rajendra/CLP/keep/sectors.py $ID 288 20 --side Z --first $1 --last $2 --dir . --threads 6 --ram 20 > ${ID}_Z$1-$2.log 2>&1 &
  nohup python3 /home/user1/Rajendra/CLP/keep/sectors.py $ID 288 20 --side X --first $1 --last $2 --dir . --threads 6 --ram 20 > ${ID}_X$1-$2.log 2>&1 &
done
sleep 30
echo "288 chunks running: $(pgrep -fc "sectors.py $ID")" >> chain_288.log
while pgrep -f "sectors.py $ID" > /dev/null; do sleep 120; done
echo "288 sectors finished at $(date)" >> chain_288.log
python3 /home/user1/Rajendra/CLP/keep/witness.py $ID 288 20 --dir . > ${ID}_witness.out 2>&1
python3 /home/user1/Rajendra/CLP/keep/verify.py  $ID 288 20 --dir . > ${ID}_verify.txt  2>&1
echo "288 witness and verify done at $(date)" >> chain_288.log
