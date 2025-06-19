#!/bin/sh

autossh -M 0 -N -i {SSH_TUNNEL_KEY} -o ServerAliveInterval=60 -o ServerAliveCountMax=3 -L 3306:{SSH_TUNNEL_DESTINATION_HOST_PORT} {SSH_TUNNEL_BASTION_USER_HOST} &

python /gwas-indexing/main.py
