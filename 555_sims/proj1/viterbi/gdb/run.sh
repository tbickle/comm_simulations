g++ -g -c base.cpp -std=c++11
g++ -g -c fsm.cpp -std=c++11
g++ -g -c viterbi.cpp -std=c++11
g++ -g -c cc_test.cpp -std=c++11
g++ -o VA base.o fsm.o viterbi.o cc_test.o
#rm base.o fsm.o viterbi.o cc_test.o
#echo $(date) > ./notes/out
#./VA >> ./notes/out
#./VA

