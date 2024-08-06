#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

WORK_PATH=${1}
SUBTASK_MEM=${2}        # memory required per job
SUBTASK_NPROCS=${3}     # number of cores required per job
IFS="," read -ra SUBTASK_PATHS <<< "${4}"   # list of run directories
IFS="," read -ra SUBTASK_LOGS <<< "${5}"    # list of worker counts
IFS="," read -ra NODES <<< "${6}"           # list of nodes for running
ACTIVATION_HOOK=${7}    # activation hook

echo "Working directory: ${WORK_PATH}"
echo "Subtask memory: ${SUBTASK_MEM}"
echo "Subtask nprocs: ${SUBTASK_NPROCS}"
echo "Subtask paths: ${SUBTASK_PATHS[@]}"
echo "Subtask logs: ${SUBTASK_LOGS[@]}"
echo "Nodes: ${NODES[@]}"
printf "Activation hook: ---\n${ACTIVATION_HOOK}\n---\n"

# Determine how many workers to put on each node, based on subtask specs
echo "Determining node capacities based on subtask specs..."
SSHLOGINS=()
SSHLOGIN=""
for node in "${NODES[@]}"; do
    node_mem_kb=$(ssh ${node} "grep MemTotal /proc/meminfo" | awk '{print $2}')
    node_mem=$((node_mem_kb / 1000000))
    node_nprocs=$(ssh ${node} "nproc --all")
    node_cap1=$((node_mem / SUBTASK_MEM))
    node_cap2=$((node_nprocs / SUBTASK_NPROCS))
    node_nwork=$((node_cap1 < node_cap2 ? node_cap1 : node_cap2))
    echo "Node ${node}: Memory=${node_mem} | Nprocs=${node_nprocs} | NWorkers=${node_nwork}"
    SSHLOGINS+=("${node_nwork}/${node}")
done
SSHLOGIN=$(IFS=,; echo "${SSHLOGINS[*]}")

LOG_FILE="{1}/{2}"
IS_RUNNING_FILE="${LOG_FILE}_IS_RUNNING"
RUN_COMMAND="automech run -p {1} &> ${LOG_FILE}"
CHECK_COMMAND="automech check-log -p ${LOG_FILE}"

parallel --sshlogin ${SSHLOGIN} "
    cd ${WORK_PATH};
    touch ${IS_RUNNING_FILE};
    eval ${ACTIVATION_HOOK@Q} &&
    printf \"Host: \$(hostname)\n| Working directory: \${PWD}\n| Command: ${RUN_COMMAND}\n\" &&
    ${RUN_COMMAND};
    rm ${IS_RUNNING_FILE};
    ${CHECK_COMMAND}
" ::: ${SUBTASK_PATHS[*]} :::+ ${SUBTASK_LOGS[*]}
