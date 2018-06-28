#!/bin/bash
. ~kjyi/src/parse
PARSE $@ << EOF
#FLAG	    	    VARNAME		DEFAULT		COMMENT
-j|--cores|--jobs	thread		1			qsub_ppn
-s|--snakefile		snakefile	./Snakefile	
EOF
# snakemake wrapper code
SNAKEMAKE=/home/users/cjyoon/anaconda3/bin/snakemake
case $1 in
	setup)
		ln -s $(readlink -f $0) ~/bin/snakemake
		exit 0
		;;
	create_templates)
		setup_dir=$(dirname $(readlink -f $0))
		cp $setup_dir/Snakefile \
		   $setup_dir/sampleConfig.yaml \
		   $setup_dir/sampleConfig.yaml \
		   $setup_dir/sampleConfig.yaml \

		rm snakemake.sh
		exit 0
		;;
	draw)
		echo draw_command her
        $SNAKEMAKE --forceall --dag | dot -Tpng > dag.png
		exit 0
		;;
	dry)
		echo drycommand here
        $SNAKEMAKE -np ${@:2}
		exit 0
		;;
	qsub)
		cat << EOF | qsub -e /dev/null -o /dev/null
#!/bin/bash
#PBS -N $(basename $snakefile)
#PBS -l nodes=1:ppn=$((thread + 1))
cd \$PBS_O_WORKDIR
$SNAKEMAKE ${@:2} &> snakemake_$(basename snakefile).out
EOF
		exit 0
		;;
esac
$SNAKEMAKE $@
