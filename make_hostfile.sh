#!/bin/bash

if (( $# != 2 )); then
	echo "Usage: ./make_hostfile ROOM SLOTS"
	exit 1
fi

if [[ $1 == "403" ]]; then
    nmap -sP 132.227.113.160/27 | grep ppti | cut -d ' ' -f 5 | cut -d '.' -f 1 | grep -v gw | sed "s/\$/ slots=$2/" > hostfile
elif [[ $1 == "401" ]]; then
    nmap -sP 132.227.113.128/27 | grep ppti | cut -d ' ' -f 5 | cut -d '.' -f 1 | grep -v gw | sed "s/\$/ slots=$2/" > hostfile
else
    echo "localhost slots=4" > hostfile
fi
