nmap -sP 132.227.114.225/27 | grep ppti | cut -d ' ' -f 5 | cut -d '.' -f 1 | grep -vE "302-06|-gw" | sed 's/$/ slots=1/' > hostfile
