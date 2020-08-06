//set params

params.dev = false
ID = params.ID_file
Overlaps = params.overlap_file
reads_dir = params.reads_dir
out_dir = params.out_dir

/*
make samplefile for singletons
*/

process singleton{

    publishDir "$out_dir/sampleFiles"

    output:
    file 'singletons.tsv' into singleton_file

    """
    #!/usr/bin/env Rscript

    library("tidyverse", quietly = T)

    singleton.patterns = list.files("${reads_dir}") %>%
        str_extract(".*-[:alpha:](?=_S)") %>%
        unique()

    write(singleton.patterns,"")
    singleton.table = tibble(
        target = "pfama1",
        sample = singleton.patterns,
        rep1 = singleton.patterns
        )

    write.table(singleton.table,
    "singletons.tsv",
    quote = F,
    sep = "\\t",
    row.names = F,
    col.names = F)
    """

    }

/*
make sample file for replicates.
*/

process replicates{

    publishDir "$out_dir/sampleFiles"

    output:
    file 'replicates.tsv' into replicates_file

    """
    #!/usr/bin/env Rscript

    library("tidyverse", quietly = T)

    replicate.patterns = list.files("${reads_dir}") %>%
        str_extract(".*(?=-[:alpha:]_S)") %>%
        unique()

    replicate.table = tibble(
        target = "pfama1",
        sample = replicate.patterns) %>%
        rowwise() %>%
        mutate(
        rep1 = list.files("${reads_dir}", pattern = sample) %>%
            str_extract(".*[:alpha:](?=_S)") %>% unique() %>% .[1],
        rep2 = list.files("${reads_dir}", pattern = sample) %>%
            str_extract(".*[:alpha:](?=_S)") %>% unique() %>% .[2]
        )

    write(show(replicate.table),"")
    write.table(replicate.table,
    "replicates.tsv",
    quote = F,
    sep = "\\t",
    row.names = F,
    col.names = F)
    """

    }



process collapse_singleton_runs{

    publishDir "$out_dir/Seekdeep"

    input:
    file sngt_file from singleton_file


    output:
    //stdout echoer
    file 'singeltons/*' into singleton_out

    when:
    !params.dev

    """
    #!/usr/bin/env bash

    SeekDeep setupTarAmpAnalysis \
        --outDir singletons \
        --inputDir ${reads_dir} \
        --idFile ${ID} \
        --overlapStatusFnp ${Overlaps} \
        --samples ${sngt_file} \

    cd singletons
    ./runAnalysis.sh 4
    """
    }



process collapse_replicate_runs{

    publishDir "$out_dir/Seekdeep"

    input:
    file rep_file from replicates_file


    output:
    file 'replicates/*'

    when:
    !params.dev

    """
    #!/usr/bin/env bash

    SeekDeep setupTarAmpAnalysis \
        --outDir replicates \
        --inputDir ${reads_dir} \
        --idFile ${ID} \
        --overlapStatusFnp ${Overlaps} \
        --samples ${rep_file} \

    cd replicates
    ./runAnalysis.sh 4
    """
    }
