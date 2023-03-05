#COMMAND="./Rss_nn_test $1 runtime-config private-$1.pem"
COMMAND="./fpSum $1 runtime-config private-$1.pem 1"

SLEEPTIME=(0 8 4 0)
size=(4 6 8 10 12 14 16 18 20)

for j in ${size[@]}
do
    for k in {1..2}
    do
        echo "COMMAND: $COMMAND 60 $j 52 11 32"
        eval "$COMMAND 60 $j 52 11 32"
        sleep ${SLEEPTIME[$1]}
    done
done
