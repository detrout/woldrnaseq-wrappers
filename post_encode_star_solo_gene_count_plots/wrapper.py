import logging
import pandas
import sys
from encoded_client.encoded import ENCODED, DCCValidator, make_attachment

filename_to_output_type = {
    "Aligned.sortedByCoord.out.bam": "alignments",
    "Gene_Unique_filtered.tar.gz": "sparse gene count matrix of unique reads",
    "Gene_EM_filtered.tar.gz": "sparse gene count matrix of all reads",
    "Gene_Unique_raw.tar.gz": "unfiltered sparse gene count matrix of unique reads",
    "Gene_EM_raw.tar.gz": "unfiltered sparse gene count matrix of all reads",
    "GeneFull_Ex50pAS_Unique_filtered.tar.gz": "sparse gene count matrix of unique reads",
    "GeneFull_Ex50pAS_EM_filtered.tar.gz": "sparse gene count matrix of all reads",
    "GeneFull_Ex50pAS_Unique_raw.tar.gz": "unfiltered sparse gene count matrix of unique reads",
    "GeneFull_Ex50pAS_EM_raw.tar.gz": "unfiltered sparse gene count matrix of all reads",
    "SJ_Unique_raw.tar.gz": "unfiltered sparse splice junction count matrix of unique reads",
}

def prepare_sc_count_matrix_qc_metric(
        config, metric_of, pct_mt_plot, gene_by_count_plot, genes_by_count_plot):
    sc_count_metric = {
        "assay_term_name": "single-cell RNA sequencing assay",
        "total_counts_vs_pct_mitochondria": make_attachment(pct_mt_plot),
        "total_counts_vs_genes_by_count": make_attachment(gene_by_count_plot),
        "counts_violin_plot": make_attachment(genes_by_count_plot),

        # run parameters not from the log file
        "quality_metric_of": metric_of,
        "step_run": config["quantification_step_run"],
        "lab": config["lab"],
        "award": config["award"],
    }
    return sc_count_metric


logger = logging.getLogger("post_encode_star_solo_gene_count_qc_plots")
logger.setLevel(logging.INFO)
if len(snakemake.log) > 0:
    log = snakemake.log[0]
    logger.addHandler(logging.FileHandler(log))
logger.addHandler(logging.StreamHandler(sys.stderr))

archive = str(snakemake.input.archive)
pct_mt_plot = str(snakemake.input.pct_mt_plot)
genes_by_count_plot = str(snakemake.input.genes_by_count_plot)
counts_violin_plot = str(snakemake.input.counts_violin_plot)
posted = str(snakemake.input.posted)

output = str(snakemake.output[0])

dry_run = snakemake.params.get("dry_run", False)

server = ENCODED(snakemake.params.submit_host)
uploaded = pandas.read_csv(posted)

logger.info("Using {} as QC file".format(archive))
output_type = filename_to_output_type[archive.name]
accession = uploaded[uploaded["output_type"] == output_type]["accession"].to_list()
assert len(accession) > 0 and not pandas.isnull(accession)

qc = prepare_sc_count_matrix_qc_metric(
    snakemake.config,
    accession,
    input.pct_mt_plot,
    input.genes_by_count_plot,
    input.counts_violin_plot
)

try:
    validator = DCCValidator(server)
    validator.validate(qc, "scrna_seq_counts_summary_quality_metric")
    if not dry_run:
        results = server.post_json("/scrna_seq_counts_summary_quality_metric/", qc)
    else:
        results = qc

    with open(output, "wt") as outstream:
        outstream.write(str(results))
except Exception as e:
    logger.error(e)
    raise e
