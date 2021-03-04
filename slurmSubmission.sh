if [ ! -d slurmOut/ ];
then
				mkdir slurmOut
fi

snakemake --use-envmodules --snakefile Snakefile --cluster-config slurmConfig.json -R all --latency-wait 60 --cluster "sbatch -J {rule} -o slurmOut/slurm-%j.out -e slurmOut/slurm-%j.err -N1 -n {cluster.threads} --time {cluster.time} --mem={cluster.mem} -A {cluster.account}" --jobs 500
