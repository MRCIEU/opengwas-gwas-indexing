#!/bin/sh

autossh -M 0 -N -i "${SSH_TUNNEL_KEY}" \
    -o StrictHostKeyChecking=no -o ServerAliveInterval=60 -o ServerAliveCountMax=3 \
    -L 127.0.0.1:3306:"${SSH_TUNNEL_DESTINATION_HOST_PORT}" "${SSH_TUNNEL_BASTION_USER_HOST}" > /tmp/autossh.log 2>&1 &

python /gwas-indexing/main.py
