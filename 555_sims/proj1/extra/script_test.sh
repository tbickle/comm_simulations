g++ -o test constellations.cpp -std=c++11
# PAM tests
./test 0 2 > output
./test 0 4 >> output

