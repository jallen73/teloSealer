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
    | grep -B 2 "^TTAGGG+" \
    | grep -A 1 "mRatio=9[0-9]\|mRatio=100" | grep -v "^#" | sed 's/ \([0-9]*\)-\([0-9]*\) / \1 \2 /' \

    """
}

process map_to_graph {
    label "telosealer"
    input:
        file gfa
        file fastq
        file telomereCsv
    output:
        path "*gaf", emit: gaf
    """
    cat $fastq | sed -e 's/ .*//' > reads.fastq
    GraphAligner ${params.GraphAligner_options} -f reads.fastq -g gfa -a readsVgraph.gaf
    cat $telomereCsv | while read csvrow ; do 
        rname=\$(echo \$csvrow | cut -c ',' -f 1) 
        nrname=\$(echo \$csvrow | cut -c ',' -f 2)
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
    grep -P "fTelo.*\t>${edge}" $gaf | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq > templ.fq
    grep -P "rTelo.*[0-9]<${edge}" $gaf | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq | seqkit seq -rc >> templ.fq
    spoa -r 0 templ.fq > ${edge}_left.fa

    grep -P "rTelo.*\t<${edge}" $gaf | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq > tempr.fq
    grep -P "fTelo.*[0-9]>${edge}" $gaf | cut -d '_' -f 1 | fgrep -f - -A 3 --no-group-sep $fastq | seqkit seq -rc >> tempr.fq
    spoa -r 0 tempr.fq > ${edge}_right.fa
    
    awk '/^S/ && \$2 == ${edge}{print ">" \$2 "\n" \$3}' > middle.fa 
    minimap2 -cx asm5 middle.fa ${edge}_left.fa ${edge}_right.fa | sort -k10,10nr | awk '!a[\$1]++' > endsVmiddle.paf
    python $workflow.projectDir/nfdir//scripts/merge_edges.py --edges middle.fa ${edge}_left.fa ${edge}_right.fa --paf endsVmiddle.paf
    """
}

workflow {
    main:
        gfa = file(params.gfa)
        fastq = file(params.fastq)
        main_edges = getMainEdges(gfa)
        main_edge_array = Channel
    .fromPath('/some/path/*.txt')
    .splitText()
        teloreads = get_teloreads(fastq)
}