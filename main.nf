#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process findTeloReads {
    label "TeloID"
    cpus 1
   
    input: 
      file fastq

    output: 
      path "telomeric_read_names.list", emit: telo_read_names
      
    """
      catfishq $fastq \
      | seqkit fq2fa \
      | NCRF $params.telomere_sequence --stats=events --minlength=45 \
      | sed -e 's/[A-Za-z]*=//g' -e 's/ \\([0-9]*\\)-\\([0-9]*\\) / \\1 \\2 /' \
      | tr -d "%" \
      | awk '/^#/{p=0}/^#/ && \$4>90{p=1} !/^TTAGGG/ && !/^#/ && p>0 && \$2>4000 && ((/TAACCC/&& \$4 < 30)||(/TTAGGG/&&\$2-\$5 < 30)) {print \$1}' \
      > telomeric_read_names.list 
    """
}

// Telosealer pipeline to correctly finish T2T contigs from a flye assembly path for species with short telomeres

// Workflow processes

process getMainEdges {
    label "telosealer"
    cpus 1

    input:
        file gfa
    
    output:
        array, emit: main_edges
    
    """
    awk '/^S/ && length(\$3) > ${params.min_main_edge_length}{print \$2}' $gfa
    """
}

process get_teloreads {
    label "telosealer"
    input:
        file fastq
    output:
        path "allTelomeres.csv", emit: telomereCsv
    
    """
    seqkit fq2fa $fastq | NCRF ${params.telomotif} --minlength=${params.minlength} --stats=events > telomeres.ncrf
    cat telomeres.ncrf | //awk madness//  > allTelomeres.csv
    cat telomeres.ncrf | //awk madness// >> allTelomeres.csv
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
    python $//nfdir//scripts/merge_edges.py --edges middle.fa ${edge}_left.fa ${edge}_right.fa --paf endsVmiddle.paf
    """
}

workflow {
    main:
        gfa = file($params.gfa)
        main_edges = getMainEdges(gfa)
}