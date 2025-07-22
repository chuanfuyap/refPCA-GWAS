import pandas as pd
import numpy as np

import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
import argparse

import multiprocessing as mp
import logging
import time
import sys


def init_worker_log():
    """Set up logging for each worker process."""
    worker_name = mp.current_process().name
    print(f"Initializing worker: {worker_name}", flush=True)  # Debug print to confirm worker initialization
    logging.basicConfig(
        filename=f'worker_{worker_name}.log',
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logging.debug(f"Worker {worker_name} started.")

def readplinkraw(raw, variants):
    "function to read raw genotype count data produced by PLINK for chunk of variants out of the full list"
    cols_to_read = ["FID"]
    cols_to_read.extend(variants)
    raw_df = pd.read_csv(raw, usecols=cols_to_read, sep="\t")

    return raw_df

def read_loco_file(loco_path, chr_num):
    """
    Read a specific REGENIE Step 1 LOCO file and return a DataFrame.
    """
    loco_df = pd.read_csv(loco_path, sep="\s+",
                        skiprows=lambda x: x != 0 and x not in [chr_num])
    loco_df = loco_df.rename(columns={"FID_IID":"CHR"})
    loco_df.set_index("CHR", inplace=True)
    loco_df = loco_df.T.reset_index().rename(columns={"index":"FID"})
    # need 
    loco_df["FID"] = loco_df["FID"].apply(lambda x: x.split("_")[0])
    loco_df.columns = ["FID", f"CHR{chr_num}"]
    
    return loco_df

def process_fit_linear_model(args):
    """Wrapper for fit_model with logging."""
    try:
        logging.debug(f"Processing fit_model with variant: {args[3]}...")  # Log a summary of args
        return fit_linear_model(*args)
    except Exception as e:
        logging.error(f"Error in fit_model with args {args}: {e}", exc_info=True)
        raise

def process_fit_interaction_model(args):
    """Wrapper for fit_interaction_model with logging."""
    try:
        logging.debug(f"Processing fit_interaction_model with variant: {args[3]}...")  # Log a summary of args
        return fit_interaction_model(*args)
    except Exception as e:
        logging.error(f"Error in fit_interaction_model with args {args}: {e}", exc_info=True)
        raise

def fit_linear_model(abt, covar_names, phenotype, variant, chr_count, model_type="linear"):
    "fit multiple regression model for a single variant"
    abt = abt.rename(columns={variant: "SNP", })

    # Build formula string
    covar_str = " + ".join(covar_names)
    # Choose model
    if model_type.lower() == "logistic":
        # Build the formula string
        formula = f"{phenotype} ~ SNP" + (f" + {covar_str}" if covar_str else "")
        model = smf.glm(formula=formula, family=sm.families.Binomial(), data=abt, offset=abt[chr_count])
    elif model_type.lower() == "linear":
        formula = f"{phenotype} ~ SNP + {chr_count}" + (f" + {covar_str}" if covar_str else "")
        model = smf.ols(formula=formula, data=abt)


    # fit the model   
    result = model.fit()
    # extract results
    output = pd.DataFrame([{
        "VARIANT": variant,
        "MAF": 1 - abt["SNP"].mean() / 2, 
        "LOG10P": -np.log10(result.pvalues["SNP"]),
        "BETA": result.params["SNP"],
        "SE": result.bse["SNP"],
        "LogLikelihood": result.llf,
    }])

    return output

def fit_interaction_model(abt, covar_names, phenotype, variant, chr_count, model_type="linear"):
    "fit multiple regression model for a single variant with interaction terms with PC1/2/3"
    abt = abt.rename(columns={variant: "SNP", })

    # create interaction terms in the formula
    interaction_terms = "SNP:PC1 + SNP:PC2 + SNP:PC3"
    covar_str = " + ".join(covar_names)

    # Choose model type
    if model_type.lower() == "logistic":
        # Build the formula string
        formula = f"{phenotype} ~ SNP + {interaction_terms}" + (f" + {covar_str}" if covar_str else "")
        model = smf.glm(formula=formula, family=sm.families.Binomial(), data=abt, offset=abt[chr_count])
    elif model_type.lower() == "linear":
        # Build the formula string
        formula = f"{phenotype} ~ SNP + {interaction_terms} + {chr_count}" + (f" + {covar_str}" if covar_str else "")
        model = smf.ols(formula=formula, data=abt)


    result = model.fit()
    # extract results
    output = pd.DataFrame([{
        "VARIANT": variant,
        "MAF": 1 - abt["SNP"].mean() / 2,
        "LOG10P": -np.log10(result.pvalues["SNP"]),
        "SNPxPC1_LOG10P": -np.log10(result.pvalues["SNP:PC1"]),
        "SNPxPC2_LOG10P": -np.log10(result.pvalues["SNP:PC2"]),
        "SNPxPC3_LOG10P": -np.log10(result.pvalues["SNP:PC3"]),
        "BETA": result.params["SNP"],
        "SE": result.bse["SNP"],
        "SNPxPC1_BETA": result.params["SNP:PC1"],
        "SNPxPC1_SE": result.bse["SNP:PC1"],
        "SNPxPC2_BETA": result.params["SNP:PC2"],
        "SNPxPC2_SE": result.bse["SNP:PC2"],
        "SNPxPC3_BETA": result.params["SNP:PC3"],
        "SNPxPC3_SE": result.bse["SNP:PC3"],
        "LogLikelihood": result.llf,
    }])

    # extract covmat, flatten it and append to the output dataframe
    covmat = result.cov_params().loc[["SNP", "SNP:PC1", "SNP:PC2", "SNP:PC3"], ["SNP", "SNP:PC1", "SNP:PC2", "SNP:PC3"]]
    covmat = covmat.values.flatten()
    covmat = pd.DataFrame(covmat, index=[f"CovMat_{i}" for i in range(16)]).T
    output = pd.concat([output, covmat], axis=1)  
    
    return output

def prep_abt(genotype, covariates, loco_pgs, phenotype, chunk):
    "loads in data and covariates along with LOCO and phenotype data to build abt for stats model"
    covar_df = covariates.copy()
    loco_df = loco_pgs.copy()
    pheno_df = phenotype.copy()

    if chunk >1:
        variants_chunked = np.array_split(genotype.columns[1:], chunk)
        abt_list = []
        for variants_subset in variants_chunked:
            # merge
            genotype_df = genotype[["FID"] + list(variants_subset)]
            abt = genotype_df.merge(covar_df, on="FID")
            abt = abt.merge(loco_df, on="FID")
            abt = abt.merge(pheno_df, on="FID")
            abt = abt.set_index("FID")
            # get the covariates and pheno names for model
            print("ABT shape after merge for analysis:", abt.shape)
            abt_list.append(abt)
            
        return abt_list, variants_chunked
    
    else:
        genotype_df = genotype
        variants = genotype_df.columns[1:]  # Exclude the FID column

        # merge
        abt = genotype_df.merge(covar_df, on="FID")
        abt = abt.merge(loco_df, on="FID")
        abt = abt.merge(pheno_df, on="FID")
        abt = abt.set_index("FID")
        # get the covariates and pheno names for model

        print("ABT shape after merge for analysis:", abt.shape)

        return [abt], [variants]

def analyse_snps(raw_df, covar_df, loco_df, pheno_df, chr_count, covar_names, pheno_name, num_cores, sub_chunk_size, model_type="linear"):
    "loads in data and covariates along with LOCO and phenotype data and fits the model for every variant in snps"
    # fix the FID column to ensure they can merge
    abt_file_list, snps_batches = prep_abt(raw_df, covar_df, loco_df, pheno_df, chunk=sub_chunk_size)
    colnames = list(covar_names.copy())
    colnames.extend([pheno_name, chr_count])
    
    all_standard_outputs = []
    all_interaction_outputs = []

    for abt_file, snps in zip(abt_file_list, snps_batches):
        print("Starting Standard Model...")
        logging.info("Starting Standard Model analysis...")
        with mp.Pool(processes=num_cores) as pool:
            standard_args = [(abt_file[colnames+[variant]], covar_names, pheno_name, variant, chr_count, model_type) for variant in snps]
            try:
                batch_standard_out = pool.map(process_fit_linear_model, standard_args)
            except Exception as e:
                logging.error(f"Error during multiprocessing (Standard Model): {e}", exc_info=True)
                raise

        all_standard_outputs.extend(batch_standard_out)

        print("Starting Interaction Model...")
        logging.info("Starting Interaction Model analysis...")

        with mp.Pool(processes=num_cores) as pool:
            interaction_args = [(abt_file[colnames+[variant]], covar_names, pheno_name, variant, chr_count, model_type) for variant in snps]
            try:
                batch_interaction_out = pool.map(process_fit_interaction_model, interaction_args)
            except Exception as e:
                logging.error(f"Error during multiprocessing (Interaction Model): {e}", exc_info=True)
                raise
        all_interaction_outputs.extend(batch_interaction_out)

    # Concat results
    standard_out = pd.concat(all_standard_outputs, ignore_index=True)
    standard_out["CHR"] = chr_count[3:]

    interaction_out = pd.concat(all_interaction_outputs, ignore_index=True)
    interaction_out["CHR"] = chr_count[3:]

    interaction_out["LRT"] = 2*(interaction_out["LogLikelihood"] - standard_out["LogLikelihood"])
    interaction_out["LRT_P"] = stats.chi2.sf(interaction_out["LRT"], 3)

    return standard_out, interaction_out


if __name__ == "__main__":
    overall_start = time.time()

    parser = argparse.ArgumentParser(description="Run REGENIE Interaction Analysis (REGINTEST) analysis")
    parser.add_argument("--chr", type=int, required=True, help="Chromosome number to process")
    parser.add_argument("--num-cores", type=int, default=4, help="Number of cores for multiprocessing")
    parser.add_argument("--chunksize", type=int, default=1000, help="Number of variants to load per chunk")
    parser.add_argument("--subchunksize", type=int, default=1, help="Number of chunks to split variants for processing")

    # Optionally allow overriding file paths
    parser.add_argument("--covar-path", type=str, )
    parser.add_argument("--pheno-path", type=str, )
    parser.add_argument("--loco-path", type=str,)
    parser.add_argument("--rawfile-path", type=str, )
    parser.add_argument("--variant-list-path", type=str,)
    parser.add_argument("--output-path", type=str, )
    parser.add_argument("--model-type", type=str, choices=["linear", "logistic"], default="linear")

    args = parser.parse_args()

    CHR = args.chr # Chromosome number to process
    NUM_CORES = args.num_cores # Number of cores to use for multiprocessing
    CHUNKSIZE= args.chunksize  # Number of variants to load in each chunk 
    SUBCHUNKSIZE = args.subchunksize  # Number of chunks to split the variants into for processing
    
    ###########################################################################################################################################################
    # Configure logging at the global level
    logging.basicConfig(
        filename='script_chr{}.log'.format(CHR),  # Log to a single file
        level=logging.DEBUG,  # Log everything (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logging.debug("Logging is set up.")
    ###########################################################################################################################################################
    # sets up parameters for analysis 

    covar_path = args.covar_path
    pheno_path = args.pheno_path
    loco_path = args.loco_path
    rawfile_path = args.rawfile_path
    variant_list_path = args.variant_list_path
    output_path = args.output_path.format(chr=CHR)
    model_type = args.model_type

    print("MULTIPROCESSING ENABLED WITH {} CORES".format(NUM_CORES), flush=True)
    print("--------------------------------------------------------------------------", flush=True)
    print("ANALYSING SNPs IN CHR{}\n".format(CHR), flush=True)
    print("MODEL TYPE: {}".format(model_type), flush=True)
    
    print("LOADING COVARIATES AND PHENOTYPES AND LOCO-PGS", flush=True)
    print("--------------------------------------------------------------------------", flush=True)
    print(f"COVARIATE FILE: {covar_path}", flush=True)
    print(f"LOCO-PGS FILE: {loco_path}", flush=True)
    print(f"PHENOTYPE FILE: {pheno_path}", flush=True)

    start = time.time()
    loco_file = read_loco_file(loco_path, CHR)
    print("LOCO file loaded with shape:", loco_file.shape, flush=True)
    pheno_file = pd.read_csv(pheno_path, sep="\t").drop(columns=[ "IID"])
    print("Phenotype file loaded with shape:", pheno_file.shape, flush=True)
    covar_file = pd.read_csv(covar_path, sep="\t").drop(columns=[ "IID"])    
    print("Covariate file loaded with shape:", covar_file.shape, flush=True)
    variant_list = list(pd.read_csv(variant_list_path, header=None, sep="\t")[0].values)

    # do some preprecessing to ensure the FID columns are strings since UKBB sample IDs are numbers
    covar_file["FID"] = covar_file["FID"].astype("str")
    loco_file["FID"] = loco_file["FID"].astype("str")
    pheno_file["FID"] = pheno_file["FID"].astype("str")

    covariates = covar_file.columns[1:]
    pheno_name = pheno_file.columns[1]

    end = time.time()

    print(f"Elapsed time for loading: {end - start:.4f} seconds", flush=True)

    print("--------------------------------------------------------------------------", flush=True)
    print("LOADING SNPs AND ANALYSING", flush=True)
    print("SNPs list file: {}".format(variant_list_path), flush=True)
    print("SNPs raw file: {}".format(rawfile_path), flush=True)
    print("Total number of variants to process: {}".format(len(variant_list)), flush=True)
    print("Processing SNPs in chunks of size: {}".format(CHUNKSIZE), flush=True)
    print("Number of chunks to process: {}".format(len(variant_list) // CHUNKSIZE + (1 if len(variant_list) % CHUNKSIZE > 0 else 0)), flush=True)
    print("--------------------------------------------------------------------------", flush=True)
    all_standard_outputs = []
    all_interaction_outputs = []

    # Process the variants in chunks
    for ix, i in enumerate(range(0, len(variant_list), CHUNKSIZE)):
        chunk = variant_list[i:i + CHUNKSIZE]
        start = time.time()

        # Process this chunk
        print(f"Processing chunk {ix} from {i} to {i+len(chunk)-1}", flush=True)
        allele_count = readplinkraw(rawfile_path, variants=chunk)
        allele_count['FID'] = allele_count['FID'].astype("str")

        standard_out, interaction_out = analyse_snps(allele_count, covar_file, loco_file, pheno_file, "CHR{}".format(CHR), covariates, pheno_name, NUM_CORES, SUBCHUNKSIZE, model_type)

        all_standard_outputs.append(standard_out)
        all_interaction_outputs.append(interaction_out)

        end = time.time()
        print(f"Elapsed time for chunk {ix}: {end - start:.4f} seconds", flush=True)

    print("SAVING ANALYSIS OUTPUT FOR CHR{}".format(CHR), flush=True)
    pd.concat(all_standard_outputs).to_csv(f"{output_path}_chr{CHR}_standard.txt", sep="\t", index=False)
    pd.concat(all_interaction_outputs).to_csv(f"{output_path}_chr{CHR}_interaction.txt", sep="\t", index=False)

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