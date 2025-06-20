#!/usr/bin/env bash

# Provides a less verbose way of executing subtask runner from a node, which would
# otherwise look like this:
#
#   pixi run node csed-0022 sub.log "automech subtasks run csed-00{20..22}"
#
# All options and flags are passed along to the subtask runner.

SUBTASK_LOG="sub.log"
ARG_OFFSET=1
while getopts ":o:" opt; do
  case $opt in
    o)
      SUBTASK_LOG=$OPTARG
      ARG_OFFSET=$((ARG_OFFSET+2))
      ;;
  esac
done

cd ${INIT_CWD:-$(pwd)}

ARGS=$(printf '%q ' "${@:ARG_OFFSET}")
COMMAND="automech subtasks run ${ARGS}"

pixi run node ${@: -1} ${SUBTASK_LOG} "${COMMAND}"
