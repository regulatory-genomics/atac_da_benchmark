nextflow.enable.dsl=2

params.outdir = "results"
include { da_peakvi as peakvi } from '../software/peakvi.nf'
include { da_scada as scada } from '../software/scaDA.nf'
include { da_pacs as pacs } from '../software/pacs.nf'
include { da_snapatac2 as snapatac2 } from '../software/snapatac2.nf'
include { da_archr as archr } from '../software/archr.nf'
include { da_deseq2 as deseq2 } from '../software/chraccr.nf'
include { json; genBenchId } from '../../common/utils.gvy'
include { pr_curve; output_metrics } from '../../common/metrics.nf'
workflow bench {
    take: datasets
    main:
    da =snapatac2(datasets) | concat(
        peakvi(datasets),
        archr(datasets),
        pacs(datasets),
        deseq2(datasets)
      //  scada(datasets)
	)
	| map({ [json(it[0]).data_name, it] })
        | combine(
                datasets | map({ [json(it[0]).data_name, it[1]] }),
                by: 0
        )
        | map ({ [genBenchId(it[1][0]), it[1][1], it[2]] })
    da
        |pr_curve |toSortedList|output_metrics
}
