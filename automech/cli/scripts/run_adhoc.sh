#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

WORK_PATH=${1}
SUBTASK_MEM=${2}        # memory required per job
SUBTASK_NPROCS=${3}     # number of cores required per job
IFS="," read -ra SUBTASK_PATHS <<< "${4}"   # list of run directories
IFS="," read -ra NODES <<< "${5}"           # list of nodes for running
ACTIVATION_HOOK=${6}    # activation hook

echo "NODES: ${NODES[*]}"
echo "SUBTASK_PATHS: ${SUBTASK_PATHS[*]}"

# Determine how many workers to put on each node, based on job memory and nprocs
echo "Determining node capacities based on job memory ${SUBTASK_MEM} and nprocs ${SUBTASK_NPROCS}..."
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
echo "Running with --sshlogin ${SSHLOGIN}"

echo "WORK_PATH: ${WORK_PATH}"
echo "ACTIVATION HOOK:"
echo "${ACTIVATION_HOOK}"

run_command="automech run -p {1} &> {1}/out.log"
check_command="automech check-log -p {1}"
parallel --sshlogin ${SSHLOGIN} "
    cd ${WORK_PATH} &&
    eval ${ACTIVATION_HOOK@Q} &&
    echo \"Host: \$(hostname) | Working directory: \${PWD} | Command: ${run_command}\" &&
    ${run_command} &&
    ${check_command}
" ::: ${SUBTASK_PATHS[*]}
