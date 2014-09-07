g++ -o test constellations.cpp -std=c++11
ts=$(date +%T)

echo $(date) > ./output/$ts

# PAM tests
echo PAM Tests:
echo 2
./test 0 2 >> ./output/$ts
echo 4
./test 0 4 >> ./output/$ts

if false
then

echo 8
./test 0 8 >> ./output/$ts
echo 16
./test 0 16 >> ./output/$ts

# PSK tests
echo PSK Tests:
echo 2
./test 1 2 >> ./output/$ts
echo 4
./test 1 4 >> ./output/$ts
echo 8
./test 1 8 >> ./output/$ts
echo 16
./test 1 16 >> ./output/$ts

# QAM tests
echo QAM Tests:
echo 2
./test 2 2 >> ./output/$ts
echo 4
./test 2 4 >> ./output/$ts
echo 8
./test 2 8 >> ./output/$ts
echo 16
./test 2 16 >> ./output/$ts

fi

