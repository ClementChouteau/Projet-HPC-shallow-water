if [[ $1 == "403" ]]; then
    nmap -sP 132.227.114.160/27 | grep ppti | cut -d ' ' -f 5 | cut -d '.' -f 1 | sed 's/$/ slots=1/' > hostfile
else
    echo "localhost slots=4" > hostfile
fi