#!/bin/bash
#$ -N Lmon3
#$ -q free64
#$ -ckpt blcr

\time -o timeout /data/users/swflynn/DGB/cluster/src/dgb_cluster < input > out
