import pandas
from encoded_client.encoded import ENCODED, DCCValidator, make_attachment
from woldrnaseq.models import (
    load_star_solo_quality_metric,
)


def prepare_star_solo_qc(config, metric_of, filename, umi_plot):
    summary = load_star_solo_quality_metric(filename)

    star_solo_metrics = {
        "assay_term_name": "single-cell RNA sequencing assay",
        "barcode_rank_plot": make_attachment(umi_plot),
        "estimated_number_of_cells": summary["Estimated Number of Cells"],
        "fraction_of_unique_reads_in_cells": summary["Fraction of Unique Reads in Cells"],
        "mean_UMI_per_cell": summary["Mean UMI per Cell"],
        "mean_reads_per_cell": summary["Mean Reads per Cell"],
        "median_UMI_per_cell": summary["Median UMI per Cell"],
        "median_reads_per_cell": summary["Median Reads per Cell"],
        "number_of_reads": summary["Number of Reads"],
        "q30_bases_in_CB_UMI": summary["Q30 Bases in CB+UMI"],
        "q30_bases_in_rna_read": summary["Q30 Bases in RNA read"],
        "reads_mapped_to_genome_unique": summary["Reads Mapped to Genome: Unique"],
        "reads_mapped_to_genome_unique_and_multiple": summary["Reads Mapped to Genome: Unique+Multiple"],
        "reads_with_valid_barcodes": summary["Reads With Valid Barcodes"],
        "sequencing_saturation": summary["Sequencing Saturation"],
        "umis_in_cells": summary["UMIs in Cells"],

        # run parameters not from the log file
        "quality_metric_of": metric_of,
        "step_run": config["quantification_step_run"],
        "lab": config["lab"],
        "award": config["award"],
    }

    if "Mean Gene per Cell" in summary:
        star_solo_metrics.update({
            "mode": "Gene",
            "reads_mapped_to_gene_unique_and_multiple_gene": summary["Reads Mapped to Gene: Unique+Multiple Gene"],
            "reads_mapped_to_gene_unique_gene": summary["Reads Mapped to Gene: Unique Gene"],
            "unique_reads_in_cells_mapped_to_gene": summary["Unique Reads in Cells Mapped to Gene"],
            "mean_gene_per_cell": summary["Mean Gene per Cell"],
            "median_gene_per_cell": summary["Median Gene per Cell"],
            "total_gene_detected": summary["Total Gene Detected"],
        })
    elif "Mean GeneFull_Ex50pAS per Cell" in summary:
        star_solo_metrics.update({
            "mode": "GeneFull_Ex50pAS",
            "reads_mapped_to_genefull_ex50pas_unique_and_multiple_gene_ex50pas": summary["Reads Mapped to GeneFull_Ex50pAS: Unique+Multiple GeneFull_Ex50pAS"],
            "reads_mapped_to_genefull_ex50pas_unique_genefull_ex50pas": summary["Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS"],
            "unique_reads_in_cells_mapped_to_genefull_ex50pas": summary["Unique Reads in Cells Mapped to GeneFull_Ex50pAS"],
            "mean_genefull_ex50pas_per_cell": summary["Mean GeneFull_Ex50pAS per Cell"],
            "median_genefull_ex50pas_per_cell": summary["Median GeneFull_Ex50pAS per Cell"],
            "total_genefull_ex50pas_detected": summary["Total GeneFull_Ex50pAS Detected"],
        })
    elif "Mean GeneFull_Ex50pAS per Cell" in summary:
        star_solo_metrics.update({
            "mode": "GeneFull",
            "reads_mapped_to_genefull_unique_and_multiple_genefull": summary["Reads Mapped to GeneFull: Unique+Multiple GeneFull"],
            "reads_mapped_to_genefull_unique_genefull": summary["Reads Mapped to GeneFull: Unique GeneFull"],
            "unique_reads_in_cells_mapped_to_genefull": summary["Unique Reads in Cells Mapped to GeneFull"],
            "mean_genefull_per_cell": summary["Mean GeneFull_Ex50pAS per Cell"],
            "median_genefull_per_cell": summary["Median GeneFull_Ex50pAS per Cell"],
            "total_genefull_detected": summary["Total GeneFull_Ex50pAS Detected"],
        })
    else:
        raise ValueError("Unknown mode? {}".format(sorted(summary.keys())))

    return star_solo_metrics


logger = logging.getLogger("post_star_qc_metrics")
logger.setLevel(logging.INFO)
if len(snakemake.log) > 0:
    log = snakemake.log[0]
    logger.addHandler(logging.FileHandler(log))
logger.addHandler(logging.StreamHandler(sys.stderr))

gene_summary = snakemake.input.gene_summary
umi_plot = snakemake.input.umi_plot
posted = str(snakemake.input.posted)

output = snakemake.output[0]

dry_run = snakemake.params.get("dry_run", False)

server = ENCODED(snakemake.params.submit_host)
uploaded = pandas.read_csv(posted)

accession = uploaded[uploaded["file_format"] == "bam"]["accession"].to_list()
assert not pandas.isnull(accession), "Please submit files. Accession column is empty"

qc = prepare_star_solo_qc(snakemake.config, accession, gene_summary, umi_plot)

try:
    validator = DCCValidator(server)
    validator.validate(qc, "star_solo_quality_metric")
    if not dry_run:
        results = server.post_json("/star_solo_quality_metric/", qc)
    else:
        results = qc
    with open(posted, "wt") as outstream:
        outstream.write(str(results))
except Exception as e:
    logger.error(e)
    raise e
