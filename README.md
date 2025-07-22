# refPCA-GWAS
Scripts and necessary variant list to run GWAS without labelling samples based on ancestry or ethnicity. 

Overall process:
1. run extract_fix_ref.sh (edit the contents to the right path)
    - this extracts variants from a PLINK file and fixes the reference allele to match the HGDP+1kG reference panel (which I have been calling gnomad for convenience)
2. run pca_projection.sh (edit the contents to the right path)
    - this calls project_gnomad.py
3. now you can run REGENIE with the additional 20 PCA components as covariates in both STEP 1 and STEP 2. 

OPTIONAL to run REGINTEST:
- this test for interaction of a SNP with the first 3 PCA components, which tests for heterogeneity of effect across ancestry groups.
- run regintest.sh (edit the contents to the right path)
    - this calls regintest.py

Requirements to run PCA projection:
- PLINK (1 or 2 should work)
- HAIL (necessary for PCA projection)

Requirements to run REGINTEST:
- statsmodels=0.14.4
- pandas=2.3.0
- scipy=1.11.3"
- numpy=1.26.1"
- patsy=1.0.1
