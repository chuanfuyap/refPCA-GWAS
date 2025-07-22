import numpy as np
import pandas as pd
from numpy.linalg import inv
from scipy.stats import chi2
from tqdm import tqdm
import argparse
import time
import sys

def stack_covariance_matrices(*matrices):
    """
    Stack multiple covariance matrices into a block-diagonal matrix.

    Parameters
    ----------
    matrices : tuple of np.ndarray
        Each element is a 2D numpy array (square matrix),
        representing a covariance matrix from one cohort/study.

    Returns
    -------
    block_diag : np.ndarray
        A larger 2D numpy array whose diagonal blocks are the input matrices.
    """
    total_size = sum(mat.shape[0] for mat in matrices) # determine how big the block diagonal matrix should be
    block_diag = np.zeros((total_size, total_size)) # initialize block diagonal matrix with zeroes
    current_index = 0

    # loop through and fill in the block diagonal matrix by indexing into the correct position
    for mat in matrices:
        size = mat.shape[0]
        block_diag[current_index:current_index + size, current_index:current_index + size] = mat
        current_index += size
    return block_diag

def create_stacked_identity(slopes):
    """
    Create a stacked identity matrix corresponding to the slopes layout.

    Parameters
    ----------
    slopes : np.ndarray
        A 2D numpy array of shape (K, P):
        - K = number of cohorts/studies.
        - P = number of parameters (e.g., 4: BETA, SNPxPC1, SNPxPC2, SNPxPC3).

    Returns
    -------
    stacked_identity : np.ndarray
        A 2D numpy array of shape (K*P, P), which vertically stacks
        an identity matrix of size (P, P) K times.
    """
    K, P = slopes.shape # K is the number of cohorts/studies, P is the number of parameters
    stacked_identity = np.vstack([np.eye(P) for _ in range(K)])
    return stacked_identity

def synregslope(slopes, covariance_matrices):
    """
    Generalized Least Squares (GLS) estimator for regression slopes.

    Parameters
    ----------
    slopes : np.ndarray
        A 2D numpy array of shape (K, P) containing slope estimates from each
        cohort/study for a single variant. For example, P = 4 for:
        [BETA, SNPxPC1_BETA, SNPxPC2_BETA, SNPxPC3_BETA].
    covariance_matrices : list of np.ndarray
        Each element is a 2D numpy array of shape (P, P) representing the
        covariance of the slope estimates in the corresponding cohort.

    Returns
    -------
    beta_hat : np.ndarray
        A 1D numpy array of length P containing the meta-analysed slope estimates.
    beta_cov : np.ndarray
        A 2D numpy array of shape (P, P) representing the covariance matrix
        of these meta-analysed slope estimates.
    """
    b_flat = slopes.flatten() # flatten the slopes into a 1D array (vector of Betas)
    Sigma = stack_covariance_matrices(*covariance_matrices) # stack the covariance matrices into a block diagonal matrix
    W = create_stacked_identity(slopes) # create a stacked identity matrix
    Sigma_inv = inv(Sigma) 
    W_T_Sigma_inv = W.T @ Sigma_inv # matrix multiplication, W transpose times Sigma inverse, part of the GLS formula
    beta_hat = inv(W_T_Sigma_inv @ W) @ W_T_Sigma_inv @ b_flat
    beta_cov = inv(W_T_Sigma_inv @ W)
    return beta_hat, beta_cov

def chi_square_value(beta_hat, beta_cov):
    """
    Compute chi-square statistic for the composite hypothesis beta = 0.

    Parameters
    ----------
    beta_hat : np.ndarray
        A 1D numpy array of length P with the estimated slopes.
    beta_cov : np.ndarray
        The corresponding 2D covariance matrix of beta_hat.

    Returns
    -------
    chi_squared_stat : float
        The chi-square test statistic for the hypothesis H0: beta = 0.
    """
    if not isinstance(beta_hat, np.ndarray):
        beta_hat = np.array(beta_hat)
    chi_squared_stat = beta_hat.T @ inv(beta_cov) @ beta_hat

    return chi_squared_stat


def meta_analyse_variants(dataframes):
    """
    Perform meta-analysis of multiple dataframes, each containing:
      - A 'MarkerName' column (unique variant identifier).
      - 4 slope columns: 'BETA', 'SNPxPC1_BETA', 'SNPxPC2_BETA', 'SNPxPC3_BETA'.
      - 16 CovMat_* columns: 'CovMat_0' ... 'CovMat_15', forming a 4x4 covariance.

    Parameters
    ----------
    dataframes : list of pd.DataFrame
        A list of pandas DataFrames. Each DataFrame must have the above columns
        and can contain multiple variants (rows).

    Returns
    -------
    meta_df : pd.DataFrame
        A DataFrame with meta-analysed results. For each variant:
          - 'MarkerName'
          - 'meta_BETA', 'meta_SNPxPC1_BETA', 'meta_SNPxPC2_BETA', 'meta_SNPxPC3_BETA'
          - 'meta_CovMat_0' ... 'meta_CovMat_15'
          - 'chi_square'
    """
    # 1. Identify the set of variants present in ALL dataframes
    marker_sets = [set(df['MarkerName']) for df in dataframes]
    union_variants = set.union(*marker_sets)

    results = []

    # 2. We wrap the for-loop in tqdm to show a progress bar
    for variant in tqdm(union_variants, desc="Meta-Analysing Variants",file=sys.stdout, disable=False):
        # Extract the slopes and covariances for this variant from each DataFrame
        slopes_list = []
        cov_matrices = []

        for df in dataframes:
            rows = df.loc[df['MarkerName'] == variant]
            if len(rows) == 0:
                continue  # this DataFrame does not have the variant, skip to the next one
            row = rows.iloc[0]
            slopes_list.append([
                row['BETA'],
                row['SNPxPC1_BETA'],
                row['SNPxPC2_BETA'],
                row['SNPxPC3_BETA']
            ])

            cov_entries = [row[f'CovMat_{i}'] for i in range(16)]
            cov_matrix = np.array(cov_entries).reshape((4, 4))
            cov_matrices.append(cov_matrix)

        # **Skip** if fewer than 2 DataFrames contributed data
        if len(slopes_list) < 2:
            continue
        
        n_studies = len(slopes_list)

        # Convert the collected slopes into a numpy array of shape (K, 4), it gets flattened later in the `synregslope` function
        slopes_array = np.array(slopes_list)

        # 3. Perform synregslope to get combined Beta estimates & covariance
        beta_hat, beta_cov = synregslope(slopes_array, cov_matrices)

        # 4. Compute the chi-square statistic of just the interaction terms
        chi_sq = chi_square_value(beta_hat, beta_cov) # this is computed for all the coefficients, checking if SNP/SNPxPC1/SNPxPC2/SNPxPC3 are all zero
        het_chi_sq = chi_square_value(beta_hat[1:], beta_cov[[1, 2, 3]][:, [1, 2, 3]]) # this is computed for the interaction terms only, checking if SNPxPC1/SNPxPC2/SNPxPC3 are all zero

        # Flatten beta_hat
        meta_BETA = beta_hat[0]
        meta_SNPxPC1_BETA = beta_hat[1]
        meta_SNPxPC2_BETA = beta_hat[2]
        meta_SNPxPC3_BETA = beta_hat[3]

        # Flatten the 4x4 covariance matrix into 16 columns
        beta_cov_flat = beta_cov.flatten()  # shape: (16,)
        cov_dict = {f'meta_CovMat_{i}': beta_cov_flat[i] for i in range(16)}

        # Gather everything into a dictionary for appending to results
        results_dict = {
            'MarkerName': variant,
            'meta_BETA': meta_BETA,
            'meta_SNPxPC1_BETA': meta_SNPxPC1_BETA,
            'meta_SNPxPC2_BETA': meta_SNPxPC2_BETA,
            'meta_SNPxPC3_BETA': meta_SNPxPC3_BETA,
            'chi_square': chi_sq,
            'het_chi_square': het_chi_sq,
            'n_studies': n_studies
        }
        results_dict.update(cov_dict)

        results.append(results_dict)

    # 5. Convert the list of result dicts to a DataFrame
    meta_df = pd.DataFrame(results)
    return meta_df
    
if __name__ == "__main__":
    overall_start = time.time()

    parser = argparse.ArgumentParser(description="Run Multivariate Meta-Analysis of REGINTEST output")
    parser.add_argument("--input-filelist", type=str, required=True, help="Txt file containing the list of summary stats files, one per line")
    parser.add_argument("--output-file", type=str, default="meta_results.tsv", help="Output file for meta-analysis results")

    args = parser.parse_args()
    # Read the input file list


    # Load all the dataframes from the file list

    with open(args.input_filelist, 'r') as f:
        file_list = f.read().splitlines()
    dataframes = []
    for file in file_list:
        print("Processing file:", file, flush=True)
        df = pd.read_csv(file, sep="\t")  # Assuming tab-separated values
        dataframes.append(df)
        print(f"Loaded {len(df)} rows from {file}", flush=True)

    # Perform the meta-analysis
    print("Starting meta-analysis...", flush=True)
    meta_results = meta_analyse_variants(dataframes)
    # Print the number of variants and studies
    print(f"Number of variants meta-analysed: {len(meta_results)}", flush=True)
    # Heterogeneity statistics
    meta_results['het_p'] = chi2.sf(meta_results["het_chi_square"], 3)
    significant_het = meta_results[meta_results['het_p']<0.05].shape[0]
    print(f"Number of variants with significant heterogeneity (p < 0.05): {significant_het}", flush=True)

    # Save the results to a file
    output_filename = args.output_file
    meta_results.to_csv(output_filename, sep="\t", index=False)

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