import csv
import sys

import matplotlib.pyplot as plt


filename = sys.argv[1] if len(sys.argv) > 1 else "uniform_samples.csv"

with open(filename, newline="", encoding="utf-8") as sample_file:
    samples = csv.reader(sample_file)
    next(samples, None)
    values = [float(row[0]) for row in samples if row]

plt.hist(values, bins=20, range=(0, 1), density=True, edgecolor="black")
plt.axhline(1, color="red", linestyle="--", label="Expected density")
plt.xlabel("Value")
plt.ylabel("Probability density")
plt.title(f"Uniform sample ({len(values)} values)")
plt.legend()
plt.tight_layout()
plt.show()
