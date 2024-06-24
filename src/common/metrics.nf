nextflow.enable.dsl = 2

include { json } from './utils.gvy'

params.resultDir = 'results'

process pr_curve {
    container 'jiaqifan/scvi-tools:v25'
    publishDir "${params.resultDir}", mode: 'copy'
    tag "${json(metadata).data_name}"
    input:
    tuple val(metadata), path("${json(metadata).data_name}_${json(metadata).method}_da_region.csv")
    output:
    //tuple path("${json(metadata).data_name}_${json(metadata).method}_Precision_Recall_Curve.pdf"), path("${json(metadata).data_name}_${json(metadata).method}_metrics.tsv"), stdout
    stdout
    """
    #!/usr/local/bin/python
    from sklearn.metrics import precision_recall_curve, auc
    from sklearn.metrics import precision_score, recall_score
    import pandas as pd
    import matplotlib.pyplot as plt
    import json
    import sys
    metadata = json.loads('${metadata}')
    metrics = {}
    original_stdout = sys.stdout
    sys.stdout = sys.stderr
    y_true = pd.read_csv("/workspace/"+metadata['data_name']+"_groundTruth.csv")
    y_pred = pd.read_csv(metadata['data_name']+"_"+metadata['method']+"_da_region.csv")
    results = []
    threshold = [0.01,0.05,0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8,0.9,1]
    def calculate_precision_recall(pred_genes, true_genes):
          # Convert lists to sets if not already sets
          pred_genes = set(pred_genes)
          true_genes = set(true_genes)

          # Calculate the intersection for true positive
          true_positive = pred_genes.intersection(true_genes)

          # Calculate the union for union_genes
          union_genes = pred_genes.union(true_genes)

          # Calculate precision and recall
          precision = len(true_positive) / len(pred_genes) if pred_genes else 0
          recall = len(true_positive) / len(true_genes) if true_genes else 0

          return precision, recall, union_genes
    if 'pval' in y_pred.columns:
        y_pred['padj'] = y_pred['pval']

    if 'prob_da' in y_pred.columns:
        y_pred['padj'] = 1 - y_pred['prob_da']
    for p in threshold:
      y_pred_sig=y_pred[(y_pred['padj']<p)]
      cell_types_t = set(y_true['cell_type'])
      cell_types_p = set(y_pred_sig['cell_type'])
      cell_types = cell_types_t.intersection(cell_types_p)
      for cell_type in cell_types:
          true_genes = set(y_true[y_true['cell_type'] == cell_type]['index'])
          pred_genes = set(y_pred_sig[y_pred_sig['cell_type'] == cell_type]['index'])
          # true_labels = [1 if index in true_genes else 0 for index in pred_genes]
          # pred_labels = [1] * len(pred_genes)  # Assuming all predicted are positive
          precision, recall, thresholds = calculate_precision_recall(pred_genes,true_genes)
          results.append((cell_type, precision, recall))
    results_df = pd.DataFrame(results, columns=['cell_type','Precision', 'Recall'])
    auc_pr=[]
    for cell_type, group in results_df.groupby('cell_type'):
        sorted_group = group.sort_values(by='Recall')
        pr_auc = auc(sorted_group['Recall'], sorted_group['Precision'])
        #   print(f"PR AUC for {cell_type}: {pr_auc:.2f}")
        auc_pr.append((cell_type,pr_auc))
    auc_pr_df = pd.DataFrame(auc_pr, columns=['cell_type','pr_auc'])
    metadata["metrics"] = auc_pr
    plt.figure(figsize=(8, 6))
    plt.plot(results_df['Recall'],results_df['Precision'],marker='o')
    plt.title('Precision-Recall curve by Cell Type')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.savefig(metadata['data_name']+"_"+metadata['method']+'_Precision_Recall_Curve.pdf')
    plt.close()
    auc_pr_df.to_csv(metadata['data_name']+"_"+metadata['method']+"_metrics.tsv")
    sys.stdout = original_stdout
    print(json.dumps(metadata), end='')
    """
}
process output_metrics {
    container 'jiaqifan/scvi-tools:v15'
    publishDir "${params.resultDir}", mode: 'copy'

    input:
        val(result)
    output:
        path('metrics.tsv')

    """
    #!/usr/local/bin/python
    import json
    import pandas as pd

    dicts = []
    for item in ${result}:
        metadata1 = item
        metric = metadata1['metrics']
        metric = dict(metric)
        metric['dataset'] = metadata1['data_name']
        metric['algorithm'] = metadata1['method']
        metric['group'] = metadata1['group'] if 'group' in metadata1 else 'main'
        dicts.append(metric)
    pd.DataFrame(dicts).to_csv("metrics.tsv", sep="\t", index=False, header=True)
    """
}

