nextflow.enable.dsl=2

params.outdir = 'results'

include { bench_atac_simulated  } from './ATAC/main.nf'
include { parseYamlString } from './common/utils.gvy'

datasets = Channel.from(
    parseYamlString(file(params.input_data).text)
)


workflow {
    if (params.assays == null || params.assays.contains('ATAC')) {
	bench_atac_simulated(datasets)
    }
    
}
