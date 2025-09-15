#!/bin/bash
#SBATCH --job-name=main
#SBATCH --output=load.log
#SBATCH --time=01:00:00
#SBATCH --account=PAS2967  # Specify PI's account here
#SBATCH --ntasks=4                   # Number of parallel tasks (processes)
#SBATCH --cpus-per-task=48            # CPUs allocated per task (adjust as needed)
#SBATCH --constraint=48core



# Your job's commands go here
source ~/.bashrc
conda activate fragment


# Get list of allocated nodes
nodes=$(scontrol show hostnames $SLURM_JOB_NODELIST)
head_node=$(echo $nodes | awk '{print $1}')
head_ip=$(srun -N1 -n1 -w $head_node hostname --ip-address)

echo "Head node: $head_node with IP: $head_ip"

# Start Ray head on first node
srun -N1 -n1 -w $head_node ray start --head \
    --node-ip-address=$head_ip \
    --port=6379 \
    --num-cpus=40 &

sleep 10  # give head some time to start

# Start Ray workers on the other nodes
for node in $nodes; do
    if [ "$node" != "$head_node" ]; then
        srun -N1 -n1 -w $node ray start \
            --address $head_ip:6379 \
            --num-cpus=40 &
    fi
done

wait



python load2.py