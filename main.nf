#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Telosealer pipeline to correctly finish T2T contigs from a flye assembly path for species with short telomeres

println "Project : $workflow.projectDir"


// Workflow processes

process getMainEdges {
    label "telosealer"
    cpus 1

    input:
        file gfa
    
    output:
        path "telolist", emit: main_edges
    
    """
    awk '/^S/ && length(\$3) > ${params.min_main_edge_length}{print \$2}' $gfa > telolist
    """
}

process get_teloreads {
    label "telosealer"
    input:
        file fastq
    output:
        path "allTelomeres.csv", emit: telomereCsv
    
    """
    seqkit fq2fa $fastq | NCRF ${params.telomotif} --minlength=${params.minMotifLength} --stats=events > telomeres.ncrf

    cat telomeres.ncrf\
    | grep -B 2 --no-group-sep "^${params.telomotif}+" \
    | grep -A 1 --no-group-sep -P "mRatio=9[0-9]|mRatio=100" | grep -v "^#" | sed 's/ \\([0-9]*\\)-\\([0-9]*\\) / \\1 \\2 /' \
    | awk '\$2 - \$5 < 100{print \$1","\$1"_fTelo"}'\
    > allTelomeres.csv

    cat telomeres.ncrf\
    | grep -B 2 --no-group-sep "^${params.telomotif}-"\
    | grep -A 1 --no-group-sep -P "mRatio=9[0-9]|mRatio=100" | grep -v "^#" | sed 's/ \\([0-9]*\\)-\\([0-9]*\\) / \\1 \\2 /' \
    | awk '\$4 < 100{print \$1","\$1"_rTelo"}'\
    >> allTelomeres.csv
    """
}

process map_to_graph {
    label "telosealer"
    cpus 16
    input:
        file gfa
        file fastq
        file telomereCsv
    output:
        path "*gaf", emit: gaf
    """
    cat $fastq | sed -e 's/ .*//' > reads.fastq
    GraphAligner ${params.GraphAligner_options} -f reads.fastq -g $gfa -a readsVgraph.gaf -t 16
    cat $telomereCsv | while read csvrow ; do 
        rname=\$(echo \$csvrow | cut -d ',' -f 1) 
        nrname=\$(echo \$csvrow | cut -d ',' -f 2)
        sed -i -e "s/\$rname\t/\$nrname\t/" readsVgraph.gaf
        done
    """
}

process seal {
    label 'telosealer'
    input:
        file gaf
        file fastq
        file gfa
        each edge
    output:
        path "*.fasta", emit: contigFasta
    """
echo "${edge}" | head -n 1 > temp.txt
bashedge=\$(head -n 1 temp.txt)
echo \$bashedge
grep -P "fTelo.*\\t<\${bashedge}[<>]" $gaf | awk '\$3 < 100 &&  \$4 > 1000' | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq | seqkit seq -rp > templ.fq
grep -P "fTelo.*\\t<\${bashedge}\\t" $gaf | awk '\$3 < 100 &&  \$4 > 1000 && \$8 < 10000' | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq| seqkit seq -rp >> templ.fq

grep -P "rTelo.*[0-9]>\${bashedge}" $gaf | awk '\$2 - \$4 < 100 && \$4 - \$3 > 1000' | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq >> templ.fq
grep -P "rTelo.*\\t>\${bashedge}" $gaf | awk '\$2 - \$4 < 100 &&  \$4 - \$3 > 1000 && \$8 < 10000' | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq >> templ.fq

if [[ \$(cat templ.fq | wc -l) -gt 0 ]] ; then
    spoa -r 0 templ.fq | sed -e 's/>/>left/' > \${bashedge}_left.fa
else
    touch \${bashedge}_left.fa
fi

grep -P "rTelo.*\\t<\${bashedge}\\t" $gaf | awk '\$2 - \$4 < 100 &&  \$4 - \$3 > 1000 && \$7 - \$9 < 10000' | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq | seqkit seq -rp  > tempr.fq
grep -P "rTelo.*[0-9]<\${bashedge}" $gaf | awk '\$2 - \$4 < 100 && \$4 - \$3 > 1000' | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq | seqkit seq -rp  > tempr.fq

grep -P "fTelo.*\\t>\${bashedge}[<>]" $gaf | awk '\$3 < 100 &&  \$4 > 1000' | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq >> tempr.fq
grep -P "fTelo.*\\t>\${bashedge}\\t" $gaf | awk '\$3 < 100 &&  \$4 > 1000 && \$7 - \$9 < 10000' | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq >> tempr.fq

if [[ \$(cat tempr.fq | wc -l) -gt 0 ]] ; then
    spoa -r 0 tempr.fq | sed -e 's/>/>right/' > \${bashedge}_right.fa
else
    touch \${bashedge}_right.fa
fi

grep -P "S\\t\${bashedge}\\t" $gfa | awk '{print ">" \$2"\\n" \$3}' > middle.fa
minimap2 -cx asm5 middle.fa \${bashedge}_left.fa \${bashedge}_right.fa | sort -k10,10nr | awk '!a[\$1]++' > endsVmiddle.paf
if [[ \$(cat endsVmiddle.paf | wc -l ) -gt 0 ]] ; then
    python $workflow.projectDir/scripts/merge_edges.py --middle middle.fa --left \${bashedge}_left.fa --right \${bashedge}_right.fa --paf endsVmiddle.paf > consensus.fasta
else
    cp middle.fa consensus.fasta
fi
    """
}

workflow {
    main:
        gfa = file(params.gfa)
        fastq = file(params.fastq)
        main_edges = getMainEdges(gfa)
        main_edge_array = main_edges
            .splitText()
        teloreads = get_teloreads(fastq)
        gaf = map_to_graph(gfa,fastq,teloreads)
        sealedContigs = seal(gaf,fastq,gfa,main_edge_array)
}