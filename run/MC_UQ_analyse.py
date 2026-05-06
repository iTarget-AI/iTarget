import pandas as pd
import numpy as np
import glob
from scipy.special import softmax
import re
import os

def load_mc_predictions(folder_path):
    """
    Load MC Dropout predictions from CSV files
    
    Args:
        folder_path: Path to folder containing prediction CSV files
        
    Returns:
        logits: Array of shape (n_samples, n_dropouts, 2)
        true_labels: Array of true labels
    """
    files = glob.glob(f"{folder_path}/*.csv")
    all_logits, true_labels = [], None
    
    for f in files:
        df = pd.read_csv(f)
        if true_labels is None:
            true_labels = df['y_true'].values
        
        # Parse logits from string format
        logits_batch = []
        for pred_str in df['y_pred']:
            numbers = re.findall(r"[-+]?\d*\.\d+|\d+", str(pred_str))
            logit0 = float(numbers[0])
            logit1 = float(numbers[1])
            logits_batch.append([logit0, logit1])
        all_logits.append(logits_batch)
    
    return np.array(all_logits).transpose(1, 0, 2), true_labels

def calculate_mean_entropy_uncertainty(logits):
    """
    Calculate uncertainty based on mean entropy of MC Dropout predictions
    
    Args:
        logits: Array of shape (n_samples, n_dropouts, 2)
        
    Returns:
        mean_probs: Mean predicted probabilities for positive class
        mean_entropy: Mean entropy across dropout samples (uncertainty measure)
    """
    # Convert logits to probabilities using softmax
    probs = softmax(logits, axis=2)
    pos_probs = probs[:, :, 1]  # Probability of positive class
    
    # Mean prediction across dropout samples
    mean_probs = np.mean(pos_probs, axis=1)
    
    # Calculate entropy for each dropout sample
    # H(p) = -p * log(p) - (1-p) * log(1-p)
    entropy_per_dropout = -np.sum(probs * np.log(probs + 1e-15), axis=2)
    
    # Average entropy across dropout samples (epistemic + aleatoric uncertainty)
    mean_entropy = np.mean(entropy_per_dropout, axis=1)
    
    return {
        'mean_probs': mean_probs,
        'mean_entropy': mean_entropy
    }

def save_predictions(true_labels, mean_probs, mean_entropy, output_path):
    """
    Save predictions and uncertainty estimates to CSV
    
    Args:
        true_labels: Array of true labels
        mean_probs: Mean predicted probabilities
        mean_entropy: Mean entropy uncertainty
        output_path: Path to save output CSV file
    """
    predictions_df = pd.DataFrame({
        'true_label': true_labels,
        'predicted_label': (mean_probs > 0.5).astype(int),
        'mean_positive_prob': mean_probs,
        'mean_entropy_uncertainty': mean_entropy
    })
    predictions_df.to_csv(output_path, index=False)
    print(f"Predictions saved to: {output_path}")

def mc_uq_calculate(folder_path, output_file):
    """
    Main function: load data and compute mean entropy uncertainty
    
    Args:
        folder_path: Path to folder containing prediction CSV files
        output_file: Path to save output CSV file
    """
    print(f"Loading data from: {folder_path}")
    logits, true_labels = load_mc_predictions(folder_path)
    print(f"Loaded: {logits.shape[1]} MC samples, {len(true_labels)} instances")
    
    print("Computing mean entropy uncertainty...")
    results = calculate_mean_entropy_uncertainty(logits)
    # {
    # 'mean_probs': mean_probs, # （n_samples, )
    # 'mean_entropy': mean_entropy # （n_samples, )
    # }    

    print(f"Saving results to: {output_file}")
    save_predictions(true_labels, results['mean_probs'], 
                    results['mean_entropy'], output_file)
    
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder_path", type=str)
    parser.add_argument("--output_file", type=str)
    params = parser.parse_args()
    print(vars(params))
    # Configuration
    folder_path = params.folder_path  # Update with your data path
    output_file = params.output_file  # Update with output path
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Run analysis
    results = mc_uq_calculate(folder_path, output_file)
