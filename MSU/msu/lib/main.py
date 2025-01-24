import pandas as pd
import os
from suda_new import find_unique_genes, find_minimal_unique_subsets, calculate_fk_msu_fm_suda, calculate_all_subsets


def main():
    # Input and output file paths
    input_file = "C:/Users/HP/PycharmProjects Dataset/sudaaaAlgo/MSU/msu/lib/datasets/cleaned_binary_gene_fpkm_matrix_dataset.csv"
    output_file = "C:/Users/HP/PycharmProjects Dataset/sudaaaAlgo/MSU/msu/lib/fk_msu_results.csv"

    # Check if the input file exists
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"The input file '{input_file}' does not exist. Please verify the path.")

    # Load dataset
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        raise ValueError(f"Error loading input file: {e}")

    # Rename the first column to "Cancer" for clarity
    if 'Unnamed: 0' in df.columns:
        df.rename(columns={"Unnamed: 0": "Cancer"}, inplace=True)
    elif 'Cancer' not in df.columns:
        raise ValueError("The dataset does not contain the expected first column for 'Cancer'.")

    # Step 1: Identify unique genes
    print("Step 1: Identifying unique genes for each cancer type...")
    unique_genes = find_unique_genes(df)
    print("Unique Genes for each cancer type:\n", unique_genes)  # Debugging: Check unique genes

    # Step 2: Calculate all possible subsets (from `suda_new.py`)
    print("Step 2: Calculating all possible subsets...")
    all_subsets = {}
    try:
        all_subsets = calculate_all_subsets(unique_genes, sample_size=1000)
        print("All Subsets for each cancer type:")
        for cancer, subsets in all_subsets.items():
            print(f"\nCancer Type: {cancer}")
            for subset in subsets:
                print(subset)
    except MemoryError:
        print("MemoryError occurred while calculating all subsets. Skipping this step.")

    # Step 3: Calculate minimal unique subsets
    print("Step 3: Calculating minimal unique subsets...")
    minimal_subsets = find_minimal_unique_subsets(unique_genes)
    print("Minimal Unique Subsets for each cancer type:\n", minimal_subsets)

    # Step 4: Calculate FK, fM, SUDA, and MSU values
    print("Step 4: Calculating FK, fM, SUDA, and MSU values...")
    fk_values, msu_values, fm_values, suda_values = calculate_fk_msu_fm_suda(unique_genes, minimal_subsets)

    # Step 5: Calculate total SUDA scores
    total_suda = sum(suda_values.values())

    # Step 6: Define file-level risk score (DIS)
    dis_score = 100  # Example value; replace with actual DIS if available

    # Step 7: Compile results into a DataFrame
    print("\nCompiling results into a DataFrame...")
    results = []
    for cancer in unique_genes.keys():
        unique_genes_str = ', '.join(unique_genes.get(cancer, [])) if unique_genes.get(cancer) else 'None'
        all_subsets_str = ', '.join(
            ['{' + ', '.join(subset) + '}' for subset in all_subsets.get(cancer, [])]) if all_subsets.get(
            cancer) else '{}'
        minimal_subsets_str = ', '.join(
            ['{' + ', '.join(subset) + '}' for subset in minimal_subsets.get(cancer, [])]) if minimal_subsets.get(
            cancer) else '{}'

        fk = fk_values.get(cancer, 0)
        msu = msu_values.get(cancer, 0)
        fm = fm_values.get(cancer, 0)
        suda = suda_values.get(cancer, 0)



        # Calculate DIS-SUDA
        dis_suda = (dis_score / total_suda) * suda if total_suda > 0 else 0

        results.append(
            [cancer, unique_genes_str, all_subsets_str, minimal_subsets_str, fk, msu, fm, suda, dis_suda])

    # Define DataFrame columns
    results_df = pd.DataFrame(results,
                              columns=['Cancer', 'Unique Genes', 'All Subsets', 'Minimal Unique Subsets', 'FK', 'MSU',
                                       'fM',
                                       'SUDA', 'DIS-SUDA'])

    # Step 8: Save results to a CSV file
    try:
        results_df.to_csv(output_file, index=False)
        print(f"\nResults successfully saved to '{output_file}'.")
    except Exception as e:
        raise IOError(f"Error saving results to file: {e}")

    # Print the results for the first 10 cancer types only
    print("\nResults for the first 10 cancer types:")
    print(results_df.head(10))


if __name__ == "__main__":
    main()