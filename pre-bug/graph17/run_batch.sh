#!/bin/bash
COUNTER=0
while [  $COUNTER -lt 20 ];
do
	./main
	let "COUNTER = COUNTER + 1"
done   
