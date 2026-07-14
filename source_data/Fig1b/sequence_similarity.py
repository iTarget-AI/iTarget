import numpy as np
import pandas as pd
import csv

# ============================
# Part 0: Generate FASTA from CSV
# ============================

def csv_to_fasta(csv_file, database_fasta, input_fasta, activity_col=3):
    """
    Read a CSV file with columns: uniprot_id, sequence, ..., activity (1 for database, 0 for input)
    and write two FASTA files.
    """
    with open(csv_file, 'r') as csvfile, open(database_fasta, 'w') as db_f, open(input_fasta, 'w') as in_f:
        reader = csv.reader(csvfile)
        header = next(reader)  # skip header
        # Assume columns: 0=ID, 1=sequence, 2=..., 3=activity
        for row in reader:
            uniprot_id = row[0]
            sequence = row[1]
            activity = row[3]  # '1' or '0'
            if activity == '1':
                db_f.write(f'>{uniprot_id}\n{sequence}\n')
            elif activity == '0':
                in_f.write(f'>{uniprot_id}\n{sequence}\n')
            else:
                # ignore or raise
                pass

# ============================
# Part 1: Build similarity matrix (original first script)
# ============================

def read_fasta_ids(path):
    """Extract protein IDs from a FASTA file (first word after '>')."""
    ids = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def build_matrix(result_csv, input_fasta, db_fasta, output_csv):
    """
    Build similarity matrix from DIAMOND CSV output (header, comma-separated).
    Columns: qseqid, sseqid, pident, ...
    Values: pident, 0 if no hit.
    """
    db_ids = read_fasta_ids(db_fasta)
    input_ids = read_fasta_ids(input_fasta)

    if not input_ids or not db_ids:
        print("ERROR: FASTA files are empty. Check CSV→FASTA step.")
        return

    col_idx = {pid: i for i, pid in enumerate(db_ids)}
    row_idx = {pid: i for i, pid in enumerate(input_ids)}

    n_rows, n_cols = len(input_ids), len(db_ids)
    mat = np.zeros((n_rows, n_cols), dtype=np.float32)

    matched = 0
    missing_query = set()
    missing_subject = set()

    with open(result_csv, "r") as f:
        # skip header
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 3:
                continue
            q, s, pident_str = parts[0].strip(), parts[1].strip(), parts[2].strip()
            try:
                pident = float(pident_str)
            except ValueError:
                continue

            ri = row_idx.get(q)
            ci = col_idx.get(s)
            if ri is None:
                missing_query.add(q)
                continue
            if ci is None:
                missing_subject.add(s)
                continue

            mat[ri, ci] = pident
            matched += 1

    print(f"Matrix built: {matched} hits inserted, "
        f"{len(missing_query)} query IDs absent, "
        f"{len(missing_subject)} subject IDs absent")

    with open(output_csv, "w", newline="") as f:
        f.write("," + ",".join(db_ids) + "\n")
        for i in range(n_rows):
            f.write(input_ids[i] + "," + ",".join(f"{v:.1f}" for v in mat[i]) + "\n")

    print(f"Matrix saved to {output_csv}")


def compute_max_per_row(input_csv, output_csv):
    """
    Read a CSV, keep the first column as ID, compute row-wise maximum
    over all remaining columns, and save as a two-column CSV ('id', 'max_value').
    """
    df = pd.read_csv(input_csv)
    id_col = df.iloc[:, 0]          # first column as ID
    data_cols = df.iloc[:, 1:]      # remaining columns for max
    max_vals = data_cols.max(axis=1)

    result = pd.DataFrame({
        'id': id_col,
        'max_value': max_vals
    })
    result.to_csv(output_csv, index=False)


def count_ranges(input_csv, output_csv):
    """
    Read a CSV with 'max_value' column and count occurrences in predefined intervals.
    Intervals: <=20, 20~30, ..., 80~90, >90   (right-closed, so <=20 includes 0 and negatives)
    Additionally, add a separate row '<20' (strictly less than 20) as a sub‑category.
    """
    df = pd.read_csv(input_csv)
    vals = df['max_value']

    # Define bins (left‑open, right‑closed) with -inf to include all values <=20
    bins = [-float('inf'), 20, 30, 40, 50, 60, 70, 80, 90, float('inf')]
    labels = ['<=20', '20~30', '30~40', '40~50', '50~60',
            '60~70', '70~80', '80~90', '>90']

    # right=True gives (a,b] so (-inf,20] includes 0 and 20
    df['range'] = pd.cut(vals, bins=bins, labels=labels, right=True)

    # Count and reorder from high to low
    counts = df['range'].value_counts()
    order = ['>90', '80~90', '70~80', '60~70', '50~60',
            '40~50', '30~40', '20~30', '<=20']
    counts = counts.reindex(order)

    counts_df = counts.reset_index()
    counts_df.columns = ['range', 'count']


    counts_df.to_csv(output_csv, index=False)

    # Verification print (English)
    print(f"Verification: count of '<=20' is {counts.loc['<=20']} ")


# ============================
# Main execution
# ============================

if __name__ == '__main__':
    # Define input file names (all in current directory)
    csv_file = "UniProt-Reviewed-Activity.csv"
    diamond_result = "Human_DIAMOND_Result_1e-5.csv"  # adjust as needed
    # Output names
    db_fasta = "database.fasta"
    in_fasta = "input.fasta"
    matrix_csv = "Similarity_Matrix.csv"
    max_csv = "Max_Identity_per_Protein.csv"
    range_csv = "Max_Identity_Distribution.csv"

    # Step 1: Generate FASTA files
    csv_to_fasta(csv_file, db_fasta, in_fasta)
    print(f"Generated {db_fasta} and {in_fasta} from {csv_file}")

    # Step 2: Build similarity matrix
    build_matrix(diamond_result, in_fasta, db_fasta, matrix_csv)

    # Step 3: Compute max per row
    compute_max_per_row(matrix_csv, max_csv)

    # Step 4: Count ranges
    count_ranges(max_csv, range_csv)

    print("All tasks completed.")