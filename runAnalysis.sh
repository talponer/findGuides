#!/bin/bash

for NAME in $(cut -f1 targetGenes.txt); do
##for NAME in HARS TAF1; do
    echo $NAME
    ID=$(grep "^$NAME" targetGenes.txt | cut -f2)
    MUT=$(grep "^$NAME" targetGenes.txt | cut -f3)
    REG=$(grep "^$NAME" targetGenes.txt | cut -f4)
    if [ "$MUT" == "-" ]; then
	MUTATION=""
    else
	MUTATION="-m $MUT"
    fi
    LOGFILE=$NAME\_putativeGuides.log
    OUTFILE=$NAME\_putativeGuides.txt
    if [ ! -f $OUTFILE ]; then
	./bin/findGuides.py -a -t $MUTATION -r $REG -i $ID 2> $LOGFILE | sort -k3,3nr -k2,2n > $OUTFILE
    fi
    #sleep 1m
done

# get a negative control sequence
for NAME in $(cut -f1 targetGenes.txt); do
#for NAME in AARS KARS HARS TAF1; do
    echo $NAME
    ID=$(grep "^$NAME" targetGenes.txt | cut -f2)
    MUT=$(grep "^$NAME" targetGenes.txt | cut -f3)
    REG=$(grep "^$NAME" targetGenes.txt | cut -f4)
    if [ "$MUT" == "-" ]; then
	MUTATION=""
    else
	MUTATION="-m $MUT"
    fi
    LOGFILE=$NAME\_putativeControlGuides.log
    OUTFILE=$NAME\_putativeControlGuides.txt
    ./bin/findGuides.py -c -t -i $ID 2> $LOGFILE | sort -k3,3nr -k2,2n > $OUTFILE 
    ##sleep 1m
done

