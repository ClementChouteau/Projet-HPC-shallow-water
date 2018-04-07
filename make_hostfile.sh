#!/bin/bash

if [[ $1 == "403" ]]; then
    nmap -sP 132.227.113.160/27 | grep ppti | cut -d ' ' -f 5 | cut -d '.' -f 1 | grep -v gw | sed 's/$/ slots=4/' > hostfile
else
    echo "localhost slots=4" > hostfile
fi
