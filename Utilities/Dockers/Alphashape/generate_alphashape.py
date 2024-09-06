import argparse
import pickle
import pandas as pd
import alphashape

# Parse command-line arguments
parser = argparse.ArgumentParser(description = "Generate and save alphashape with the given training data and the alpha parameter.")
parser.add_argument("-t", "--training-data", required = True)
parser.add_argument("-a", "--alpha", default = 8.0)
args = parser.parse_args()

# Read training data
# Expected input format: tsv with columns PC1 and PC2
df_train = pd.read_csv(args.training_data, sep = '\t', header = 0)
training_points = list(zip(df_train.PC1, df_train.PC2))

# Generate alphashape
alpha_shape = alphashape.alphashape(training_points, float(args.alpha))

# Pickle alphashape and save to file
with open("alphashape.pickle", "wb") as outfile:
	pickle.dump(alpha_shape, outfile)
