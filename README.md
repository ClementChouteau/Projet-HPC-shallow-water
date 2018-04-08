# [HPC] Projet : shallow water model
### Clément Chouteau & Romain Mekarni

## Compiler

- `make` : lance la compilation des différentes versions. Crée un fichier `hostfile`.

## Configurer

- `make hostfile` crée un fichier `hostfile` adapté pour les salles 401 (par défaut) et 403
- `./make_hostfile.sh ROOM SLOTS` exemple : `./make_hostfile.sh 0 4` retourne `localhost slots=4`

## Tester

- `make times` : exécute les 4 versions parallèles sur un test (test2 par défaut) et retourne un fichier `times.csv`. 

Dans `parallel` il y a les tests suivants :
- `test2` : test de référence performances en parallèle `-x 8192 -y 8192 -t 20`
- `test3` : export de référence parallèle `-x 512 -y 512 -t 40 --export --export-path ~/`
- `big_test` : `-x 32768 -y 32768 -t 40` à paralléliser sur beaucoup de noeuds (> 16)
- `small_test` : `-x 2048 -y 2048 -t 20` à paralléliser sur une machine locale avec plus de 2 coeurs
- `testSave` : tester un export custom. `make io` crée et ouvre la dernière image png de la simulation.

Les tests non de référence peuvent être customisés :
- `make NODES=4 TEST=testSave EXPORT_STEP=5 times`
- `make NODES=64 TEST=big_test times`
- `cd parallel; make NODES=4 MODE="--block --async" SIZE_TEST=2048 T_TEST=40 testSave`

## Nettoyer

- `make clean` ou `make mrproper`