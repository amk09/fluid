#!/bin/bash

PARTITION="gpu"
GRES="gpu:1"
TIME="02:00:00"
MEM="4G"
CPUS="1"
CONSTRAINT=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --time)
      TIME="$2"
      shift 2
      ;;
    --mem)
      MEM="$2"
      shift 2
      ;;
    --cpus)
      CPUS="$2"
      shift 2
      ;;
    --profile)
      CONSTRAINT="--constraint=gtx_2080_ti|titan_rtx|rtx_a6000|l40"
      shift
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Launch an interactive job with specified options
srun \
  --partition=gpus \
  --gres=gpu:1 \
  --time="$TIME" \
  --mem="$MEM" \
  --cpus-per-task="$CPUS" \
  $CONSTRAINT \
  --pty \
  bash