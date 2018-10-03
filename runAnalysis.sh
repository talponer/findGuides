#!/bin/bash

for NAME in $(cut -f1 targetGenes.txt); do
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
    python ./bin/findGuides.py $MUTATION -r $REG -i $ID 2> $LOGFILE | sort -k3,3nr -k2,2n > $OUTFILE
done
