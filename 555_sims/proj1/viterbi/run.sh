g++ -c -g ./source/base.cpp -std=c++11
g++ -c -g ./source/fsm.cpp -std=c++11
g++ -c -g ./source/viterbi.cpp -std=c++11
g++ -c -g ./source/cc_test.cpp -std=c++11
g++ -o VA base.o fsm.o viterbi.o cc_test.o
rm base.o fsm.o viterbi.o cc_test.o
echo $(date) > ./notes/out
./VA >> ./notes/out
./VA

