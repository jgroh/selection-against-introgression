import numpy as np
reps = np.arange(20)
gens = ["gen002","gen005", "gen010", "gen050", "gen100", "gen500", "gen750", "gen1000"]



rule all:
	input:
		"results/multi-scale/ancestry_master.txt"

rule simulation:
	input:
	output:
		expand("results/multi-scale/replicate{{rep}}_{gen}.trees", rep=reps, gen=gens)
	shell:
		"slim -d seed={wildcards.rep} code/multi-scale.slim"

rule ancestry_calc:
	input:
		"results/multi-scale/{rep_gen}.trees"
	output:
		"results/multi-scale/{rep_gen}_ancestry.txt"
	shell:
		"python code/ancestry.py {input} > {output}"

rule aggregate:
	input:
		expand("results/multi-scale/replicate{rep}_{gen}_ancestry.txt", rep=reps, gen=gens)
	output:
		"results/multi-scale/ancestry_master.txt"
	shell:
		"cat {input} > {output}"
