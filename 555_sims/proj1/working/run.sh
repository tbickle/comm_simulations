g++ -c ./source/comm.cpp -std=c++11
g++ -c ./source/base.cpp -std=c++11
g++ -c ./source/fsm.cpp -std=c++11
g++ -c ./source/viterbi.cpp -std=c++11
g++ -c ./source/constellations.cpp -std=c++11
g++ -o sim comm.o base.o fsm.o viterbi.o constellations.o
rm comm.o base.o fsm.o viterbi.o constellations.o
./sim 0 2

if false
then

ts=$(date +%T)
echo $(date) > ./output/$ts
# PAM tests
echo PAM Tests:
echo 2
./sim 0 2 >> ./output/$ts
echo 4
./sim 0 4 >> ./output/$ts
echo 8
./sim 0 8 >> ./output/$ts
echo 16
./sim 0 16 >> ./output/$ts
# PSK tests
echo PSK Tests:
echo 2
./sim 1 2 >> ./output/$ts
echo 4
./sim 1 4 >> ./output/$ts
echo 8
./sim 1 8 >> ./output/$ts
echo 16
./sim 1 16 >> ./output/$ts
# QAM tests
echo QAM Tests:
echo 2
./sim 2 2 >> ./output/$ts
echo 4
./sim 2 4 >> ./output/$ts
echo 8
./sim 2 8 >> ./output/$ts
echo 16
./sim 2 16 >> ./output/$ts

fi

