import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix

# List of .bam files corresponding to different cell lineages
bam_files = ["file1.bam", "file2.bam", "file3.bam"]

# Extract data from .bam files into a Pandas DataFrame
data = []
cell_lineages = []
for bam_file in bam_files:
    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam.fetch():
        read_id = read.query_name
        cpg_methylation_status = None
        chromosome = read.reference_name
        read_pair = read.is_read1
        cytosine_conversions = 0
        cytosine_retentions = 0
        for tag, value in read.get_tags():
            if tag == "ZW":
                cpg_methylation_status = value
            if tag == "ZC":
                cytosine_conversions = value
            if tag == "ZR":
                cytosine_retentions = value
        if cpg_methylation_status:
            data.append([read_id, cpg_methylation_status, chromosome, "read1" if read_pair else "read2", cytosine_conversions, cytosine_retentions])
            cell_lineages.append(bam_file)

data = pd.DataFrame(data, columns=["read_id", "methylation_status", "chromosome", "read_pair", "conversions", "retentions"])
cell_lineages = pd.Series(cell_lineages)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(data, cell_lineages, test_size=0.2, random_state=0)

# Fit a GradientBoostingClassifier to the training data
gbc = GradientBoostingClassifier()
gbc.fit(X_train, y_train)

# Predict cell lineage for the test data
y_pred = gbc.predict(X_test)

# Plot the confusion matrix
cm = confusion_matrix(y_test, y_pred)
plt.imshow(cm, interpolation="nearest", cmap=plt.cm.Blues)
plt.colorbar()
tick_marks = np.arange(len(bam_files))
plt.xticks(tick_marks, bam_files, rotation=45)
plt.yticks(tick_marks, bam_files)
plt.xlabel("Predicted cell lineage")
plt.ylabel("True cell lineage")
plt.title("Confusion matrix")
plt.tight_layout()
plt.show()

