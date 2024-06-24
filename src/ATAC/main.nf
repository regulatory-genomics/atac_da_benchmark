nextflow.enable.dsl=2

params.outdir = 'results'

include { download_dataset } from '../common/download3'
include { bench as DA } from './benchmark/DA'
include { pr_curve } from '../common/metrics.nf'
workflow bench_atac {
    take: metadata

    main:
        data = metadata
            | filter {it.assay == "ATAC" }
            | download_dataset 

        if (params.components == null || params.components.contains('DA')) {
            DA(data)
        }


}


workflow bench_atac_simulated {
    take: metadata

    main:
	data = metadata
	    | filter {it.assay == "simulating" }
	    | download_dataset
        if (params.components == null || params.components.contains('DA')) {
	res = DA(data)
//	ress= res| map({ [json(it[0]).data_name, it] })

  //      ress.view { println("${it}")} 
      //  pr_curve(res)
        //   input_data.view { println("${it}")}
           
 }
}
