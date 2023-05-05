
import logging
import os
import pandas
from pathlib import Path
import sys
from encoded_client.encoded import ENCODED, HTTPError
from encoded_client.submission import process_endpoint_files

assert len(snakemake.input) == len(snakemake.output), "List of files to process needs to be same length"

logger = logging.getLogger("post_encode_star_solo_files")
logger.setLevel(logging.INFO)
if len(snakemake.log) > 0:
    log = snakemake.log[0]
    logger.addHandler(logging.FileHandler(log))
logger.addHandler(logging.StreamHandler(sys.stderr))

submit_host = snakemake.params.submit_host
server = ENCODED(submit_host)

metadata_filename = str(snakemake.input.metadata)
dry_run = snakemake.params.get("dry_run", False)

# Keep track of output files we think are supposed to be created
upload_names = ["bam", "gene_unique_raw", "gene_multi_raw", "sj_unique_raw"]
upload_paths = []
for output_name in upload_names:
    output_value = snakemake.output.get(output_name)
    if output_value is not None:
        upload_paths.append(Path(str(output_value)))

# posted_filename is required
posted_filename = str(snakemake.output.posted)
upload_paths.append(Path(posted_filename))

metadata = pandas.read_csv(
    metadata_filename,
    dtype={
        "uuid": str,
        "accession": str
    },
    index_col=None)

# catch previous submission and add their accession to the metadata sheet
for i, row in metadata.iterrows():
    if pandas.isnull(row["accession"]):
        try:
            submitted_file = server.get_json("md5:{}".format(row["md5sum"]))
            metadata.at[i, "accession"] = submitted_file["accession"]
            metadata.at[i, "uuid"] = submitted_file["uuid"]
            upload_file = Path("{}.{}.upload".format(
                row["submitted_file_name"].replace("/", "_"), submit_host))
            print("Creating already uploaded file marker", upload_file)
            if not upload_file.exists():
                upload_file.touch()
        except HTTPError as e:
            if e.response.status_code != 404:
                logger.warning("Unexpected status code {}".format(e))

uploaded = process_endpoint_files(server, "/files/", metadata, dry_run=dry_run)

logger.info("Processed {} files".format(len(uploaded)))
logger.info(metadata)
metadata.to_csv(posted_filename, index=False)

creation_errors = False
for filename in upload_paths:
    if not filename.exists():
        logger.info("Unable to find {}".format(filename))
        creation_errors = True

if creation_errors:
    for name in sorted(os.listdir()):
        logger.info("{} {}".format(os.getcwd(), name))

    curdir = Path(metadata_filename).parent
    for name in sorted(os.listdir()):
        logger.info("{} {}".format(curdir, name))
