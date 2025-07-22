import hail as hl
from gnomad.sample_qc.ancestry import pc_project
import argparse
import time

def project_individuals(project_mt, pca_loadings):
    """
    Project samples into predefined PCA space
    :param project_mt: matrix table of samples to project 
    :param pca_loadings: existing PCA space of reference samples 
    """
    ht_projections = pc_project(project_mt, pca_loadings)
    ht_projections = ht_projections.transmute(**{f'PC{i}': ht_projections.scores[i - 1] for i in range(1, 51)})

    return ht_projections

if __name__ == "__main__":
    overall_start = time.time()
    parser = argparse.ArgumentParser(description="Project data to HGDP+1kG PCA space using HAIL")
    parser.add_argument("--input", type=str,  help="Input genotype PLINK file to process")
    parser.add_argument("--reference", type=str,  help="path to the gnomad_loadings HT file, e.g. gnomad_loadings.ht")
    parser.add_argument("--output", type=str, help="Output file")

    args = parser.parse_args()

    QUERYDATA = args.input
    OUTPUT = args.output
    gnomad = args.reference


    # initialise the pyspark
    hl.init(spark_conf={'spark.driver.memory': '128g'}, default_reference='GRCh38')

    # load in PLINK data and convert to HAIL matrixtable (it doesnt have to be PLINK, feel free to edit, but be sure to convert to HAIL)

    genotype_hail = './genotype_hail.mt'
    hl.import_plink(bed=f'{QUERYDATA}.bed',
                    bim=f'{QUERYDATA}.bim',
                    fam=f'{QUERYDATA}.fam').write(genotype_hail, overwrite=True)

    genotype_mt = hl.read_matrix_table(genotype_hail)

    # load in the gnomad_loadings which is downloaded from the github
    gnomad_loadings_path = gnomad
    loadings_ht = hl.read_table(gnomad_loadings_path)
    print(loadings_ht.count())

    ### sanity check on overlapping variants
    genotype_common = genotype_mt.semi_join_rows(loadings_ht)
    overlap_count = genotype_common.count()[0]

    print("OVERLAPPING VARIANTS:\t", overlap_count)

    if overlap_count == 0:
        raise ValueError("No overlapping variants found between genotype and gnomAD loadings.")
    elif overlap_count < 5445:
        print("WARNING: Less than 50% overlapping variants found, make sure to align reference alleles and check the data quality.")
    elif overlap_count < 8167:
        print("WARNING: Less than 75% overlapping variants found, this may affect the analysis.")

    genotype_pca = project_individuals(genotype_mt, loadings_ht)

    genotype_pca = genotype_pca.to_pandas()
    genotype_pca = genotype_pca.rename(columns={"s":"ID"})
    genotype_pca.to_csv(f"{OUTPUT}_pca.tsv", sep="\t", index=False)

    # Calculate total elapsed time
    overall_end = time.time()
    total_seconds = overall_end - overall_start

    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60

    parts = []
    if hours > 0:
        parts.append(f"{hours} hour{'s' if hours > 1 else ''}")
    if minutes > 0:
        parts.append(f"{minutes} minute{'s' if minutes > 1 else ''}")
    if seconds > 0 or not parts:
        parts.append(f"{seconds:.4f} second{'s' if seconds > 1 else ''}")

    total_str = " ".join(parts)
    print(f"Total elapsed time: {total_str}", flush=True)